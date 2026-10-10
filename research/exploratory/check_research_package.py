"""Structural document QA only; not scientific or online source validation."""
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = '''RESEARCH_CHARTER LITERATURE_REVIEW COMPETITOR_ANALYSIS
MATHEMATICAL_FOUNDATION RESEARCH_HYPOTHESES BASELINE_EXPERIMENT ARCHITECTURE
VALIDATION_PROTOCOL IMPLEMENTATION_ROADMAP REFERENCES DECISION_MEMO
RESEARCH_STATUS'''.split()


def main():
    issues, ids, sources = [], set(), []
    for name in REQUIRED:
        if not (ROOT/(name+'.md')).is_file():
            issues.append('Missing required file: '+name)
    for path in sorted((ROOT/'sources').glob('*.json')):
        data = json.loads(path.read_text(encoding='utf-8-sig'))
        if path.name == 'bibliography_index.json':
            continue
        for record in data:
            rid = record['id']
            if rid in ids:
                issues.append('Duplicate source identifier: '+rid)
            ids.add(rid)
            for key in ['title','authors','year','venue','url','access',
                        'contribution','supported_claim','limitations']:
                if not record.get(key):
                    issues.append(f'{rid}: empty metadata {key}')
            if not record['url'].startswith('https://'):
                issues.append(rid+': non-HTTPS canonical source')
            sources.append(record)
    links = 0
    for path in sorted(ROOT.rglob('*.md')):
        body = path.read_text(encoding='utf-8-sig')
        if '\ufffd' in body:
            issues.append(str(path.relative_to(ROOT))+': replacement character')
        for rid in re.findall(r'\b[LMRT][0-9]{2}\b',body):
            if rid not in ids:
                issues.append(str(path.relative_to(ROOT))+': unknown source '+rid)
        for target in re.findall(r'\]\(([^\s)]+)\)',body):
            if '://' in target or target.startswith('#'):
                continue
            target = target.split('#')[0]
            if target:
                links += 1
                if not (path.parent/target).exists():
                    issues.append(f'{path.relative_to(ROOT)}: broken local link {target}')
    result = {'scope':'structural document QA, not scientific validation',
              'required_documents':len(REQUIRED),
              'markdown_files':len(list(ROOT.rglob('*.md'))),
              'source_records':len(sources), 'local_links_checked':links,
              'issues':sorted(set(issues)), 'passed':not issues}
    (ROOT/'exploratory'/'package_check_results.json').write_text(
        json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if issues:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
