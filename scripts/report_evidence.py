"""Render archived measurements; never rerun or alter solver evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from prometheus.geometry import make_l_mesh


def separated_best(order,metric):
    """2% reporting margin reflects two 1% resolution gates, not an error bound."""
    if not order:return None
    if len(order)>1 and abs(metric(order[1])-metric(order[0]))/max(abs(metric(order[1])),abs(metric(order[0])))<=.02:
        return None
    return order[0]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',nargs='+',required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    out=Path(args.output)
    out.mkdir(parents=True,exist_ok=True)
    runs=[(Path(r),json.loads((Path(r)/'summary.json').read_text())) for r in args.runs]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    root,main_run=runs[0]
    mesh=make_l_mesh(2.)
    fig,axs=plt.subplots(2,3,figsize=(12,6.7),layout='constrained')
    for ax,(name,records) in zip(axs.flat,main_run['evaluations'].items()):
        rho=np.load(root/'candidates'/name/'physical_source.npy')
        grid=np.full((40,60),np.nan)
        grid[mesh.cell_indices[:,1],mesh.cell_indices[:,0]]=rho
        ax.imshow(grid,origin='lower',extent=(0,120,0,80),vmin=0,vmax=1,cmap='Greys',interpolation='nearest')
        first=records[0]
        ax.set(title=f"{name} | volume {first.get('volume_fraction',float('nan')):.3f}\n{main_run['verification'][name]['status']}",xlabel='x [mm]',ylabel='y [mm]')
        ax.set_aspect('equal')
    fig.suptitle('PROMETHEUS-01: frozen h=2 mm candidate fields; numerical models only',fontsize=14)
    fig.savefig(out/'candidate_fields.png',dpi=180);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for run,data in runs:
        hist=json.loads((run/'optimization_1'/'history.json').read_text())
        label=f"volume {data['evaluations']['B3'][0]['volume_budget_fraction']:.2f}"
        for ax,x,xlabel in [(axs[0],[r['evaluation'] for r in hist],'Factorizations / evaluated designs'),
                           (axs[1],[r.get('elapsed_wall_seconds',r['cumulative_solve_seconds']) for r in hist],'Elapsed optimization wall time [s]')]:
            ax.plot(x,[r['compliance'] for r in hist],label=label)
            ax.set(xlabel=xlabel,ylabel='Q4 nominal compliance [N mm]',yscale='log')
            ax.grid(alpha=.2);ax.legend()
    fig.suptitle('SIMP histories: same algorithm and 300-design cap')
    fig.savefig(out/'optimization_histories.png',dpi=180);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    for (name,rr),color in zip(main_run['evaluations'].items(),colors):
        rr=[r for r in rr if 'independent' in r]
        if not rr:continue
        for ax,key,ylabel in zip(axs,['compliance_n_mm','max_pad_displacement_mm'],['Compliance [N mm]','Maximum pad displacement [mm]']):
            ax.plot([r['mesh_h_mm'] for r in rr],[r[key] for r in rr],'-o',label=f'{name} Q4',color=color)
            ax.plot([r['mesh_h_mm'] for r in rr],[r['independent'][key] for r in rr],'--x',label=f'{name} P2',color=color)
            ax.set(xlabel='Analysis h [mm]; frozen original field',ylabel=ylabel,yscale='log')
            ax.set_xticks([.5,1,2]);ax.set_xlim(2.1,.4);ax.grid(alpha=.2)
    axs[1].legend(ncol=2,fontsize=8,loc='center left',bbox_to_anchor=(.01,.63))
    fig.suptitle('Fixed-design refinement: compliance and displacement are distinct metrics')
    fig.savefig(out/'refinement.png',dpi=180);plt.close(fig)
    comparisons=[]
    for path,data in runs:
        ev=data['evaluations'];vf=ev['B3'][0]['volume_budget_fraction']
        binary=['B1','B3_binary']
        coarse=[n for n in binary if ev[n][0].get('optimization_feasible')]
        verified=[n for n in binary if data['verification'][n]['passed'] and ev[n][-1].get('optimization_feasible')]
        coarse=sorted(coarse,key=lambda n:ev[n][0]['compliance_n_mm'])
        verified=sorted(verified,key=lambda n:ev[n][-1]['independent']['compliance_n_mm'])
        item={'run':str(path),'volume_fraction':vf,'equal_mass_binary_candidates':binary,
              'coarse_valid_order':coarse,'independently_verified_order':verified,
              'selected_coarse':separated_best(coarse,lambda n:ev[n][0]['compliance_n_mm']),
              'selected_verified':separated_best(verified,lambda n:ev[n][-1]['independent']['compliance_n_mm']),
              'tie_reporting_relative_margin':.02,'tie_rule_note':'Postprocessing reporting convention, not a certified error bound; no winner for <=2% separation',
              'candidate_count':len(ev),'mass_compliant_model_feasible_count':sum(r[0].get('optimization_feasible',False) for r in ev.values()),
              'valid_binary_geometry_count':sum(bool(r[0].get('geometry_status',{}).get('valid')) for r in ev.values() if r[0].get('representation')=='binary'),
              'binary_geometry_count':sum(r[0].get('representation')=='binary' for r in ev.values()),
              'mass_compliant_independently_verified_count':sum(data['verification'][n]['passed'] and r[-1].get('optimization_feasible',False) for n,r in ev.items()),
              'response_false_feasible_rate':None,'response_false_feasible_reason':'No response limit or cheap feasibility classifier defined',
              'partial_audit_budgets':[]}
        for policy,order in [('fixed_ID',binary),('coarse_objective',coarse+[n for n in binary if n not in coarse])]:
            for budget in [1,2]:
                attempted=order[:budget]
                accepted=[n for n in attempted if n in verified]
                accepted=sorted(accepted,key=lambda n:ev[n][-1]['independent']['compliance_n_mm'])
                selected=separated_best(accepted,lambda n:ev[n][-1]['independent']['compliance_n_mm'])
                item['partial_audit_budgets'].append({'policy':policy,'attempt_budget':budget,'attempted':attempted,'accepted':accepted,'selected':selected,
                    'archived_audit_wall_seconds':sum(r.get('wall_seconds',0) for n in attempted for r in ev[n]),
                    'factorizations':sum(r.get('factorizations',0) for n in attempted for r in ev[n]),
                    'rhs_solves':sum(r.get('rhs_solves',0) for n in attempted for r in ev[n]),
                    'timing_note':'Conditional audit-only retrospective replay, not total algorithm cost or a separately timed policy. Common generation and screening are excluded from these increments. Full attempt includes mandatory terminal verification; counts include rerun coarse solves with no reuse credit.',
                    'common_prerequisites':{'generation_factorizations':data['optimization'][0]['evaluations'],
                        'generation_wall_seconds':data['optimization'][0]['wall_seconds'],
                        'screening_factorizations':sum(bool(ev[n][0].get('solver_passed')) for n in binary),
                        'screening_rhs_solves':2*sum(bool(ev[n][0].get('solver_passed')) for n in binary),
                        'screening_primary_analysis_seconds':sum(ev[n][0].get('primary_analysis_seconds',0) for n in binary),
                        'note':'Same generation and geometry/Q4 screening of both candidates charged to both policies; invalid geometry requires no Q4 solve. Screening wall overhead was not separately instrumented, so no end-to-end policy speedup is claimed.'}})
        comparisons.append(item)
    (out/'comparisons.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    # Illustrative fixed cases are evaluations of the archived quadratic,
    # not substitutes for the exact stationary-angle envelope.
    angle_rows=[]
    for path,data in runs:
        for name,records in data['evaluations'].items():
            for row in records:
                for family,metrics in [('Q4',row),('P2',row.get('independent',{}))]:
                    if 'load_gram' not in metrics:continue
                    q=np.asarray(metrics['load_gram'])
                    values={}
                    for degrees in (-20,0,20):
                        theta=np.deg2rad(degrees);a=np.array([np.sin(theta),np.cos(theta)])
                        values[str(degrees)]=float(a@q@a)
                    angle_rows.append({'run':str(path),'candidate':name,'h_mm':row['mesh_h_mm'],
                                       'family':family,'compliance_n_mm_by_angle_degrees':values,
                                       'status':'discrete_model_diagnostic'})
    (out/'fixed_angle_diagnostics.json').write_text(json.dumps(angle_rows,indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'run_summary_sha256':{str(path):hashlib.sha256((path/'summary.json').read_bytes()).hexdigest() for path,_ in runs},
        'tie_rule':'Explicit reporting-only 2% relative separation margin, not a rigorous uncertainty interval; introduced during postprocessing and not used to tune designs.'},indent=2)+'\n')
    print(json.dumps({'figures':str(out),'comparisons':[{k:c[k] for k in ['volume_fraction','selected_coarse','selected_verified']} for c in comparisons]},indent=2))


if __name__=='__main__':main()
