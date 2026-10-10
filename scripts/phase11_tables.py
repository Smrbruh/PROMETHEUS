"""Assemble traceable Phase 1.1 tables and a convergence diagnostic figure."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time

ROOT=Path(__file__).resolve().parents[1]


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path,value):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def csv_save(path,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('x',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def main(output,refinement):
    start=time.perf_counter();cpu=time.process_time();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    ref=Path(refinement).resolve();catalog=read(ROOT/'results/phase11/geometry/catalog.json')['entries']
    entries={e['catalog_key']:e for e in catalog};labels=read(ref/'labels.json');responses=read(ref/'responses.json')
    costs=read(ref/'costs.json');imports=read(ref/'imported_evidence.json');rows=[]
    for result in labels:
        key=result['catalog_key'];angle=result['angle_degrees'];entry=entries[key];g=entry['geometry'];label=result['label']
        history=[r for r in responses if r['catalog_key']==key and r['angle_degrees']==angle]
        valid=[r for r in history if r['status']=='evaluated' and r['p2_compliance_n_mm'] is not None]
        finest=min(valid,key=lambda r:r['mesh_h_mm']) if valid else {}
        q=label.get('q4_refinement',{});p=label.get('p2_refinement',{})
        estimate=label.get('reference_estimate') or {}
        rows.append({'catalog_key':key,'candidate_name':entry['candidate_name'],'volume_setting':entry['volume_setting'],
            'angle_degrees':angle,'volume_fraction':g['volume_fraction'],'components':g['edge_component_count'],
            'local_point_contacts':g['local_point_contact_count'],'minimum_axis_run_mm':g['minimum_axis_run_mm'],
            'minimum_ligament_mm':None,'mass_minus_budget_kg':g['mass_minus_budget_kg'],
            'historical_nominal_passed':(entry['historical_label'] or {}).get('passed'),
            'numerical_status':label['numerical_status'],'geometry_status':label['geometry_status'],
            'constraint_status':label['constraint_status'],'classification':label['classification'],
            'eligible_reference':label['eligible_reference'],'reference_mesh_h_mm':estimate.get('mesh_h_mm'),
            'reference_compliance_n_mm':estimate.get('compliance_n_mm'),
            'reference_pad_mm':estimate.get('max_pad_displacement_mm'),
            'last_successful_pair_h_mm':finest.get('mesh_h_mm'),
            'last_observed_p2_compliance_n_mm':finest.get('p2_compliance_n_mm'),
            'last_observed_p2_pad_mm':finest.get('p2_pad_mm'),
            'last_observed_q4_compliance_change':finest.get('q4_compliance_change'),
            'last_observed_p2_compliance_change':finest.get('p2_compliance_change'),
            'last_observed_q4_pad_change':finest.get('q4_pad_change'),
            'last_observed_p2_pad_change':finest.get('p2_pad_change'),
            'last_observed_same_h_compliance_gap':finest.get('same_h_compliance_gap'),
            'last_observed_same_h_pad_gap':finest.get('same_h_pad_gap'),
            'terminal_h_mm':min(r['mesh_h_mm'] for r in history),
            'q4_last_compliance_change':q.get('last_compliance_change'),
            'p2_last_compliance_change':p.get('last_compliance_change'),
            'q4_last_pad_change':q.get('last_pad_change'),'p2_last_pad_change':p.get('last_pad_change'),
            'finest_compliance_gap':label.get('finest_compliance_agreement'),'finest_pad_gap':label.get('finest_pad_agreement'),
            'reference_reason':'; '.join(label.get('reasons',[])),
            'historical_evaluation_wall_seconds':sum(i.get('historical_wall_seconds') or 0 for i in imports if i['catalog_key']==key),
            'new_attempt_wall_seconds':sum(c['wall_seconds'] for c in costs['attempts'] if c['catalog_key']==key),
            'cost_scope':'Basis cost shared across three angles; do not sum repeated candidate costs across angle rows',
            'array_path':entry['array_path'],'source_input_sha256':entry['source_input_sha256'],'geometry_sha256':entry['geometry_sha256']})
    save(out/'candidate_results.json',rows);csv_save(out/'candidate_results.csv',rows)
    schedule=[]
    for key,entry in entries.items():
        history={r['mesh_h_mm']:r for r in responses if r['catalog_key']==key and r['angle_degrees']==0}
        terminal=history[min(history)]
        for h in (2.,1.,.5,.25,.125,.0625):
            if h in history:status=history[h]['status'];performed=True
            elif not entry['geometry']['geometry_valid']:status='not_attempted_geometry_rejected';performed=False
            elif terminal['status']!='evaluated':status='not_attempted_after_terminal_failure';performed=False
            else:status='not_required_all_declared_angles_passed';performed=False
            schedule.append({'catalog_key':key,'mesh_h_mm':h,'status':status,'attempt_record_exists':performed,
                'note':'Conditional protocol: unattempted levels have no numerical response or convergence claim.'})
    save(out/'refinement_schedule.json',schedule);csv_save(out/'refinement_schedule.csv',schedule)
    thresholds=[]
    for e in read(ROOT/'results/phase11/geometry/thresholds.json')['entries']:
        g=e['diagnostics']
        thresholds.append({'candidate_id':e['candidate_id'],'volume_setting':e['volume_setting'],
            'threshold':e['rule']['threshold'],'volume_fraction':g['volume_fraction'],
            'mass_minus_budget_kg':g['mass_minus_budget_kg'],'mass_change_from_density_kg':g['mass_change_from_density_kg'],
            'components':g['edge_component_count'],'local_point_contacts':g['local_point_contact_count'],
            'minimum_axis_run_mm':g['minimum_axis_run_mm'],'minimum_ligament_mm':None,
            'geometry_valid':g['geometry_valid'],'within_budget':g['budget_feasible'],
            'numerical_status':'not_evaluated','eligible_reference':False,'array_path':e['array_path'],
            'array_file_sha256':e['array_file_sha256'],'generation_wall_seconds':e['generation_wall_seconds']})
    save(out/'threshold_results.json',thresholds);csv_save(out/'threshold_results.csv',thresholds)
    ranking=[]
    for volume in (.35,.40,.45):
        for angle in (-20.,0.,20.):
            subset=[r for r in rows if r['volume_setting']==volume and r['angle_degrees']==angle and r['constraint_status']=='within_budget']
            observed=sorted([r for r in subset if r['last_observed_p2_compliance_n_mm'] is not None],key=lambda r:r['last_observed_p2_compliance_n_mm'])
            eligible=sorted([r for r in subset if r['eligible_reference']],key=lambda r:r['reference_compliance_n_mm'])
            ties=[]
            if eligible:
                best=eligible[0]['reference_compliance_n_mm']
                ties=[r['catalog_key'] for r in eligible if (r['reference_compliance_n_mm']-best)/best<=.02]
            ranking.append({'volume_setting':volume,'angle_degrees':angle,
                'observed_unqualified_order':[r['catalog_key'] for r in observed],
                'observed_meshes_mm':[r['last_successful_pair_h_mm'] for r in observed],
                'eligible_order':[r['catalog_key'] for r in eligible],'within_2_percent_of_best':ties,
                'note':'Observed order may mix mesh spacings and is not an eligibility or certified continuum-ranking claim; tie margin is reporting only.'})
    save(out/'rankings.json',ranking)
    # The old evaluator counts only returned P2 results. Its residual exception
    # occurs after LU.solve, so charge this known completed work separately.
    ledger=[]
    for c in costs['attempts']:
        row={k:v for k,v in c.items() if k not in ('preflight',)}
        path=ROOT/c['record_path'] if c.get('record_path') else None
        record=read(path) if path and path.exists() else {}
        exc=record.get('exception',{})
        row['category']='new_refinement';row['known_extra_p2_factorizations']=0;row['known_extra_p2_rhs_solves']=0
        if exc.get('type')=='ArithmeticError' and exc.get('message')=='P2 free-DOF equilibrium residual exceeds 1e-8':
            row['known_extra_p2_factorizations']=1;row['known_extra_p2_rhs_solves']=2
            row['accounting_note']='Residual gate follows a returned P2 LU solve of both basis RHS; original evaluator counters stop at Q4. Added stage charges do not imply P2 acceptance.'
        elif exc.get('type')=='MemoryError' and 'factor.U.diagonal()' in exc.get('traceback','') and 'verification.py' in exc.get('traceback',''):
            row['known_extra_p2_factorizations']=1
            row['accounting_note']='Traceback locates memory failure while exposing P2 factor.U for pivot inspection, after splu returned and before factor.solve. Charge one P2 factorization and zero P2 RHS; no pivot or P2 response is inferred.'
        row['charged_factorizations']=None if c['factorizations'] is None else c['factorizations']+row['known_extra_p2_factorizations']
        row['charged_rhs_solves']=None if c['rhs_solves'] is None else c['rhs_solves']+row['known_extra_p2_rhs_solves']
        row['solve_count_scope']='completed executed stages' if record.get('status')=='evaluated' or row['known_extra_p2_factorizations'] else 'confirmed returned stages; other incurred work may be unknown'
        if c['status']=='resource_preflight_deferred':row['solve_count_scope']='preflight: exactly zero solves'
        if exc and not row['known_extra_p2_factorizations']:
            row['confirmed_factorization_lower_bound']=row['charged_factorizations']
            row['confirmed_rhs_lower_bound']=row['charged_rhs_solves']
            row['charged_factorizations']=None;row['charged_rhs_solves']=None
        row['partial_nominal_q4_compliance_n_mm']=record.get('compliance_n_mm') if record.get('status')!='evaluated' else None
        row['partial_nominal_q4_pad_mm']=record.get('max_pad_displacement_mm') if record.get('status')!='evaluated' else None
        ledger.append(row)
    for c in read(ROOT/'results/phase11/reproduction.json')['comparisons']:
        cost=c['replay_cost'];ledger.append({'category':'baseline_replay','catalog_key':c['replay'],
            'wall_seconds':cost['wall_seconds'],'worker_cpu_seconds':cost['cpu_seconds'],
            'charged_factorizations':cost['completed_factorizations'],'charged_rhs_solves':cost['completed_rhs_solves'],
            'peak_working_set_bytes':cost['peak_memory_bytes'],'status':'reproduced'})
    generation=read(ROOT/'results/phase11/geometry/generation_cost.json')
    for path,name in [(ROOT/'results/phase11/kinematics/complete/cost.json','kinematics_complete')]:
        c=read(path);ledger.append({'category':name,'wall_seconds':c['wall_seconds'],'worker_cpu_seconds':c['cpu_seconds'],
            'charged_factorizations':0,'charged_rhs_solves':0,'status':'completed'})
    ledger.append({'category':'geometry_generation','wall_seconds':generation['generation_wall_seconds'],
        'worker_cpu_seconds':generation['generation_cpu_seconds'],'charged_factorizations':0,'charged_rhs_solves':0,'status':'completed'})
    ledger.append({'category':'angle_projection','wall_seconds':costs['projection_wall_seconds'],
        'charged_factorizations':0,'charged_rhs_solves':0,'status':'completed'})
    validation=read(ROOT/'results/phase11/validation_replay.json')['measured_cost']
    ledger.append({'category':'validation_replay','wall_seconds':validation.get('wall_seconds'),
        'charged_factorizations':validation['factorizations'],'charged_rhs_solves':validation['rhs_solves'],'status':'completed'})
    aborted=ROOT/'results/phase11/kinematics/failure.json'
    if aborted.exists():
        c=read(aborted);ledger.append({'category':'kinematics_csv_failure','wall_seconds':c.get('wall_seconds'),
            'worker_cpu_seconds':c.get('cpu_seconds'),'charged_factorizations':0,'charged_rhs_solves':0,'status':'failed_retained',
            'accounting_note':'See preserved failure record for exact timing scope; no mechanics solves.'})
    save(out/'evaluation_costs.json',ledger);csv_save(out/'evaluation_costs.csv',ledger)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    keys=['v035_B1','v040_B3_binary','v040_B3_threshold','v045_B3_binary','v045_B3_threshold','v035_repair']
    fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    for ax,key in zip(axes.flat,keys):
        rr=sorted([r for r in responses if r['catalog_key']==key and r['angle_degrees']==0],key=lambda r:-r['mesh_h_mm'])
        for prefix,style in [('q4','o-'),('p2','s-')]:
            points=[r for r in rr if r[prefix+'_compliance_change'] is not None]
            ax.loglog([r['mesh_h_mm'] for r in points],[100*r[prefix+'_compliance_change'] for r in points],style,label=prefix.upper())
        ax.axhline(1,color='black',ls='--',lw=.8,label='1% gate');ax.invert_xaxis();ax.grid(True,which='both',alpha=.25)
        ax.set_title(key);ax.set_xlabel('Mesh spacing h [mm]');ax.set_ylabel('Successive compliance change [%]')
        failed=[r for r in rr if r['status']!='evaluated']
        if failed:ax.text(.03,.06,f"Terminal h={failed[-1]['mesh_h_mm']:g}:\n{failed[-1]['status']}",transform=ax.transAxes,fontsize=8)
    axes.flat[0].legend(fontsize=8);fig.suptitle('Frozen binary geometry: nominal compliance refinement\nSuccessive changes are diagnostics, not exact error bounds')
    fig.savefig(out/'compliance_refinement.png',dpi=160);fig.savefig(out/'compliance_refinement.svg');plt.close(fig)
    save(out/'table_build.json',{'wall_seconds':time.perf_counter()-start,'cpu_seconds':time.process_time()-cpu,
        'factorizations':0,'rhs_solves':0,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                  [ROOT/'results/phase11/geometry/catalog.json',ref/'labels.json',ref/'responses.json',ref/'costs.json']}})
    print(json.dumps({'candidate_rows':len(rows),'threshold_rows':len(thresholds),'cost_rows':len(ledger),'output':str(out)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='results/phase11/tables')
    p.add_argument('--refinement',default='results/phase11/refinement');args=p.parse_args();main(args.output,args.refinement)
