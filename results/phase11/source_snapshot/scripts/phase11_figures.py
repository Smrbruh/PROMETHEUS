"""Render the immutable response table with explicit, legible mesh ticks."""
from pathlib import Path
import argparse
import hashlib
import json
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,FixedLocator,FixedFormatter


def main(output):
    start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    source=root/'results/phase11/refinement/responses.json';rows=json.loads(source.read_text())
    keys=['v035_B1','v040_B3_binary','v040_B3_threshold','v045_B3_binary','v045_B3_threshold','v035_repair']
    fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    for ax,key in zip(axes.flat,keys):
        data=sorted([r for r in rows if r['catalog_key']==key and r['angle_degrees']==0],key=lambda r:-r['mesh_h_mm'])
        for prefix,style in [('q4','o-'),('p2','s-')]:
            points=[r for r in data if r[prefix+'_compliance_change'] is not None]
            ax.plot([r['mesh_h_mm'] for r in points],[100*r[prefix+'_compliance_change'] for r in points],style,label=prefix.upper())
        ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xlim(1.15,.105)
        ax.xaxis.set_major_locator(FixedLocator([1,.5,.25,.125]));ax.xaxis.set_major_formatter(FixedFormatter(['1','0.5','0.25','0.125']))
        ax.xaxis.set_minor_formatter(NullFormatter());ax.tick_params(labelsize=9)
        ax.axhline(1,color='black',ls='--',lw=.8,label='1% gate');ax.grid(True,which='both',alpha=.22)
        ax.set_title(key,fontsize=11);ax.set_xlabel('Mesh spacing h [mm]');ax.set_ylabel('Successive compliance change [%]',fontsize=9)
        failed=[r for r in data if r['status']!='evaluated']
        if failed:
            r=failed[-1];label='memory stop' if r['status']=='memory_limit' else 'P2 verification failed'
            ax.text(.03,.04,f"h={r['mesh_h_mm']:g}: {label}",transform=ax.transAxes,fontsize=8,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
    axes.flat[0].legend(fontsize=8);fig.suptitle('Frozen binary geometry: nominal compliance refinement\nSuccessive changes are diagnostics, not exact error bounds',fontsize=13)
    fig.savefig(out/'compliance_refinement.png',dpi=160);fig.savefig(out/'compliance_refinement.svg');plt.close(fig)
    with (out/'manifest.json').open('x',encoding='utf-8') as f:
        json.dump({'source':str(source.relative_to(root)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                   'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'wall_seconds':time.perf_counter()-start,'factorizations':0,'rhs_solves':0,
                   'note':'A legibility revision; original table-export figures remain preserved.'},f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='results/phase11/figures')
    main(p.parse_args().output)
