"""Common evaluation and frozen-field transfer; no prediction is acceptance."""
from dataclasses import dataclass, asdict
from time import perf_counter, process_time
import hashlib
import json
import numpy as np

from .geometry import make_l_mesh, subset_mesh, bracket_pads, bracket_fixed_dofs, bracket_load_basis
from .fem import solve, stress_metrics
from .optimization import geometry_status
from .verification import bracket_verify
from .envelope import compliance_envelope


@dataclass(frozen=True)
class Candidate:
    name: str
    physical: np.ndarray
    representation: str  # density or binary

    def __post_init__(self):
        values=np.array(self.physical,dtype=float,copy=True)
        if values.ndim != 1 or not np.isfinite(values).all() or np.any((values<0)|(values>1)):
            raise ValueError('Physical density must be a finite [0,1] vector')
        if self.representation not in ('density','binary'):
            raise ValueError('Unknown representation')
        if self.representation=='binary' and not np.isin(values,[0,1]).all():
            raise ValueError('Binary candidate must contain only zero or one')
        values.setflags(write=False)
        object.__setattr__(self,'physical',values)

    @property
    def content_hash(self):
        return hashlib.sha256(self.physical.astype('<f8').tobytes()+self.representation.encode()).hexdigest()


def transfer_frozen(source_mesh, values, target_mesh):
    """Exact inheritance on origin-aligned square grids; no refiltering."""
    values=np.asarray(values)
    if values.shape != (len(source_mesh.cells),):
        raise ValueError('Source coefficient shape mismatch')
    for grid in (source_mesh,target_mesh):
        expected=(grid.cell_indices+.5)*grid.h
        if not np.allclose(grid.centroids,expected,rtol=0,atol=1e-10):
            raise ValueError('Transfer requires origin-aligned cell indices and geometry')
        corner=grid.cell_indices[:,None,:]*grid.h+grid.h*np.array([[0,0],[1,0],[1,1],[0,1]])
        if not np.allclose(grid.nodes[grid.cells],corner,rtol=0,atol=1e-10):
            raise ValueError('Transfer requires regular square cells')
    ratio=source_mesh.h/target_mesh.h
    if ratio<1 or not np.isclose(ratio,round(ratio),rtol=0,atol=1e-10):
        raise ValueError('Target must be an integer subdivision of the source')
    lookup={tuple(ij):i for i,ij in enumerate(source_mesh.cell_indices)}
    parent=np.floor(target_mesh.centroids/source_mesh.h).astype(int)
    try: indices=np.array([lookup[tuple(ij)] for ij in parent])
    except KeyError as exc:raise ValueError('Target extends outside frozen source') from exc
    return values[indices].copy()


def evaluate_candidate(candidate, source_mesh, material, benchmark, h=2., independent=False,
                       include_stress=True):
    started,cpu=perf_counter(),process_time()
    if candidate.physical.shape != (len(source_mesh.cells),):
        raise ValueError('Candidate/source shape mismatch')
    mesh=make_l_mesh(h)
    rho=transfer_frozen(source_mesh,candidate.physical,mesh)
    pads=bracket_pads(mesh)
    volume=float(material.thickness*(mesh.areas@rho))
    total=float(material.thickness*mesh.areas.sum())
    volume_ok=volume <= benchmark.volume_fraction*total+1e-6*total
    pads_ok=bool(np.all(rho[pads]==1))
    problem_signature=hashlib.sha256(json.dumps({'material':asdict(material),'benchmark':asdict(benchmark)},sort_keys=True).encode()).hexdigest()
    source_hash=hashlib.sha256(source_mesh.nodes.tobytes()+source_mesh.cells.tobytes()).hexdigest()
    record={'candidate':candidate.name,'hash':candidate.content_hash,
            'problem_signature':problem_signature,'source_geometry_hash':source_hash,
            'representation':candidate.representation,'mesh_h_mm':h,
            'source_h_mm':source_mesh.h,'volume_mm3':volume,
            'mass_kg':volume*material.density,'volume_fraction':volume/total,
            'volume_budget_fraction':benchmark.volume_fraction,
            'volume_feasible':bool(volume_ok),'pads_preserved':pads_ok,
            'geometry_status':None,'status':'pending','solver_passed':False,
            'optimization_feasible':False,'independent_verified':False,
            'physical_validation':'not_assessed','manufacturability':'not_assessed',
            'response_feasibility':'not_defined_no_response_limit',
            'stress_status':'diagnostic_only','factorizations':0,'rhs_solves':0}
    if candidate.representation=='binary':
        geo=geometry_status(mesh,rho)
        record['geometry_status']=geo
        if not geo['valid']:
            record.update(status='geometry_failed',wall_seconds=perf_counter()-started,
                          cpu_seconds=process_time()-cpu)
            return record,{}
        active=rho.astype(bool)
        mesh=subset_mesh(mesh,active)
        rho=rho[active]
        moduli=np.full(len(mesh.cells),material.E)
    else:
        moduli=material.E*(benchmark.emin_ratio+(1-benchmark.emin_ratio)*rho**benchmark.p)
    record['q4_elements']=len(mesh.cells)
    record['q4_ndof']=mesh.ndof
    forces=benchmark.load*bracket_load_basis(mesh)
    solved=solve(mesh,moduli,forces,bracket_fixed_dofs(mesh),material.nu,material.thickness)
    nominal=solved.u[:,1].reshape(-1,2)
    magnitudes=np.linalg.norm(nominal,axis=1)
    material_nodes=np.unique(mesh.cells[rho>=.9])
    pad_nodes=np.unique(mesh.cells[bracket_pads(mesh)])
    gram=forces.T@solved.u
    envelope=compliance_envelope(gram,*np.deg2rad([-20,20]))
    free=~bracket_pads(mesh)
    grey=float(np.sum(mesh.areas[free]*4*rho[free]*(1-rho[free]))/max(mesh.areas[free].sum(),1))
    record.update(status='evaluated',solver_passed=bool(solved.diagnostics['passed']),
                  optimization_feasible=bool(volume_ok and pads_ok and solved.diagnostics['passed']),
                  compliance_n_mm=float(solved.compliance[1]),
                  load_weighted_displacement_mm=float(solved.compliance[1]/benchmark.load),
                  max_displacement_all_nodes_mm=float(magnitudes.max()),
                  max_material_displacement_mm=float(magnitudes[material_nodes].max()) if len(material_nodes) else None,
                  max_pad_displacement_mm=float(magnitudes[pad_nodes].max()),
                  small_deformation_review_trigger=bool(magnitudes[pad_nodes].max()/120>.01),
                  worst_compliance_n_mm=envelope.maximum,
                  worst_angle_degrees=float(np.rad2deg(envelope.angle_radians)),
                  load_gram=gram.tolist(),grey_indicator=grey,
                  primary_analysis_seconds=solved.wall_seconds,
                  factorizations=solved.factorizations,rhs_solves=solved.rhs_solves,
                  diagnostics=solved.diagnostics)
    arrays={'nodes':mesh.nodes,'cells':mesh.cells,'rho':rho,'u_basis':solved.u,
            'forces_basis':forces}
    if include_stress:
        stress=stress_metrics(mesh,moduli,solved.u[:,1],material.nu,material.thickness,densities=rho)
        record['stress']={k:v for k,v in stress.items() if not isinstance(v,np.ndarray)}
        arrays.update({k:v for k,v in stress.items() if isinstance(v,np.ndarray)})
    if independent:
        checked=bracket_verify(mesh,moduli,material,load=benchmark.load)
        ev=compliance_envelope(checked.gram,*np.deg2rad([-20,20]))
        d=checked.diagnostics
        energy=max(d['energy_error_relative'])
        force=float(np.linalg.norm(np.asarray(d['resultant_balance']),axis=0).max())/benchmark.load
        moment=float(np.max(np.abs(d['moment_balance'])))/(benchmark.load*120)
        passed=bool(max(d['free_residual_relative'])<=1e-8 and energy<=1e-8
                    and force<=1e-7 and moment<=1e-7 and d['symmetry_error']<=1e-12)
        record['independent']={'element':'P2_triangle','compliance_n_mm':float(checked.compliance[1]),
                               'max_pad_displacement_mm':checked.max_pad_displacement,
                               'worst_compliance_n_mm':ev.maximum,
                               'worst_angle_degrees':float(np.rad2deg(ev.angle_radians)),
                               'q4_relative_compliance_difference':abs(float(solved.compliance[1]-checked.compliance[1]))/abs(float(checked.compliance[1])),
                               'q4_relative_pad_difference':abs(record['max_pad_displacement_mm']-checked.max_pad_displacement)/checked.max_pad_displacement,
                               'diagnostics':d,'solver_passed':passed,'wall_seconds':checked.wall_seconds,
                               'triangles':len(checked.mesh.cells),'ndof':checked.mesh.ndof}
        record['factorizations']+=1
        record['rhs_solves']+=2
        arrays.update(p2_nodes=checked.mesh.nodes,p2_cells=checked.mesh.cells,p2_u_basis=checked.u,
                      p2_forces_basis=checked.forces)
    record.update(wall_seconds=perf_counter()-started,cpu_seconds=process_time()-cpu)
    return record,arrays


def verification_summary(records):
    """Accept only the complete h2,h1,h.5 fixed-field sequence plus P2 checks."""
    if not records:
        return {'passed':False,'status':'resolution_unverified','reason':'No records'}
    if any(r.get('status')=='geometry_failed' for r in records):
        return {'passed':False,'status':'geometry_failed','reason':'Invalid binary geometry'}
    identity=('candidate','hash','representation','source_h_mm','problem_signature','source_geometry_hash')
    if any(tuple(r.get(k) for k in identity)!=tuple(records[0].get(k) for k in identity) for r in records):
        return {'passed':False,'status':'unverified','reason':'Mixed frozen designs or problem definitions'}
    if len({r.get('mesh_h_mm') for r in records})!=len(records):
        return {'passed':False,'status':'unverified','reason':'Duplicate mesh levels'}
    if any(not r.get('pads_preserved') for r in records):
        return {'passed':False,'status':'unverified','reason':'Preserved pads not established'}
    by_h={r['mesh_h_mm']:r for r in records}
    if any(h not in by_h for h in (2.,1.,.5)):
        return {'passed':False,'status':'resolution_unverified','reason':'Required mesh sequence incomplete'}
    if any(not r.get('solver_passed') or not r.get('independent',{}).get('solver_passed') for r in records):
        return {'passed':False,'status':'unverified','reason':'Solver/geometry/independent audit failed'}
    coarse,middle,fine=(by_h[h] for h in (2.,1.,.5))
    if any(not np.isfinite(r[key]) or r[key]<=0 or not np.isfinite(r['independent'][key]) or r['independent'][key]<=0
           for r in records for key in ('compliance_n_mm','max_pad_displacement_mm')):
        return {'passed':False,'status':'unverified','reason':'Invalid numerical metrics'}
    def change(a,b,key):return abs(a[key]-b[key])/max(abs(b[key]),1e-12)
    c=change(middle,fine,'compliance_n_mm');d=change(middle,fine,'max_pad_displacement_mm')
    pc=change(middle['independent'],fine['independent'],'compliance_n_mm')
    pd=change(middle['independent'],fine['independent'],'max_pad_displacement_mm')
    agree_c=fine['independent']['q4_relative_compliance_difference']
    agree_d=fine['independent']['q4_relative_pad_difference']
    previous_c=change(coarse,middle,'compliance_n_mm')
    previous_p=change(coarse['independent'],middle['independent'],'compliance_n_mm')
    trend_ok=bool(c<=previous_c+1e-10 and pc<=previous_p+1e-10)
    passed=bool(c<.01 and d<.02 and pc<.01 and pd<.02 and agree_c<.01 and agree_d<.02 and trend_ok)
    return {'passed':passed,'status':'numerically_verified' if passed else 'resolution_unverified',
            'q4_compliance_h1_to_h0_5':c,'q4_pad_h1_to_h0_5':d,
            'p2_compliance_h1_to_h0_5':pc,'p2_pad_h1_to_h0_5':pd,
            'finest_compliance_agreement':agree_c,'finest_pad_agreement':agree_d,
            'q4_compliance_h2_to_h1':previous_c,'p2_compliance_h2_to_h1':previous_p,
            'compliance_refinement_trend_improves':trend_ok,
            'physical_validation':'not_assessed'}
