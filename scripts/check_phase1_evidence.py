"""Check archived evidence integrity without rerunning numerical solves."""
from pathlib import Path
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from prometheus.evaluation import verification_summary


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8'))


def main():
    paths=['phase1_v040_verified','phase1_v035','phase1_v045','phase1_h1_reoptimization']
    checked=[]
    for name in paths:
        folder=ROOT/'results'/name
        manifest=read(folder/'manifest.json');summary=read(folder/'summary.json')
        for rel,expected in manifest['source_sha256'].items():
            saved=folder/'source_snapshot'/rel
            if not saved.exists():saved=folder/'source_snapshot'/Path(rel).name
            assert digest(saved)==expected,(name,rel,'source hash')
        for rel,expected in manifest['specification_sha256'].items():
            assert digest(folder/'specification_snapshot'/rel)==expected,(name,rel,'spec hash')
        optimizations=summary['optimization']
        if isinstance(optimizations,dict):optimizations=[optimizations]
        for number,opt in enumerate(optimizations,1):
            sub=folder/f'optimization_{number}'
            hist=read(sub/'history.json')
            assert len(hist)==opt['evaluations']
            assert hist[-1]['compliance']==opt['final_compliance_n_mm']
            assert max(r['volume_violation'] for r in hist)<=.032
            for kind in ('raw','physical'):
                values=np.load(sub/f'{kind}.npy')
                assert hashlib.sha256(values.astype('<f8').tobytes()).hexdigest()==opt[f'{kind}_sha256']
        groups=summary['evaluations']
        if isinstance(groups,list):groups={groups[0]['candidate']:groups}
        for candidate,rows in groups.items():
            verdict=summary['verification'] if name.endswith('reoptimization') else summary['verification'][candidate]
            assert verification_summary(rows)==verdict,(name,candidate,'verification')
            for row in rows:
                saved=read(folder/'candidates'/candidate/f"h{row['mesh_h_mm']:g}.json")
                assert saved==row,(name,candidate,'individual record mismatch')
                assert row['independent_verified']==verdict['passed']
                field=folder/'candidates'/candidate/f"h{row['mesh_h_mm']:g}_fields.npz"
                if field.exists():
                    with np.load(field) as arrays:
                        volume=float(arrays['rho'].sum()*row['mesh_h_mm']**2*manifest['material']['thickness'])
                        assert abs(volume-row['volume_mm3'])<1e-8
                        assert arrays['u_basis'].shape==arrays['forces_basis'].shape
        cost=summary['cost']
        assert cost['evaluation_factorizations']==sum(r.get('factorizations',0) for rows in groups.values() for r in rows)
        assert cost['evaluation_rhs_solves']==sum(r.get('rhs_solves',0) for rows in groups.values() for r in rows)
        checked.append({'run':name,'summary_sha256':digest(folder/'summary.json'),'candidate_groups':len(groups),
                        'source_and_specification_hashes_match':True,'arrays_counts_and_statuses_match':True})
    tests=ET.parse(ROOT/'results'/'tests_final.xml').getroot()
    suites=list(tests.iter('testsuite'))
    assert all(int(s.get('failures','0'))==0 and int(s.get('errors','0'))==0 for s in suites)
    docs=[ROOT/'README.md']+[ROOT/'research'/f'{name}.md' for name in
        ['IMPLEMENTATION_REPORT','NUMERICAL_VALIDATION','EXPERIMENT_PROTOCOL','BASELINE_RESULTS','KNOWN_LIMITATIONS','PHASE1_DECISION']]
    for path in docs:
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            assert (path.parent/target.split('#')[0]).exists(),(path.name,target,'missing local link')
    report={'checks_passed':True,'scope':'Artifact integrity, not a new mechanics solve or physical validation',
            'runs':checked,'tests':sum(int(s.get('tests','0')) for s in suites),
            'test_record_sha256':digest(ROOT/'results'/'tests_final.xml'),
            'report_links_checked':len(docs),'script_sha256':digest(Path(__file__))}
    (ROOT/'results'/'evidence_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
