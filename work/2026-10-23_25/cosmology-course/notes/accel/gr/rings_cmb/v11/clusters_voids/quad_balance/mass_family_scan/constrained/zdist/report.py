import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from fitz import L,Y,S
B='../';R=json.load(open('fitz_results.json'));IR=json.load(open(B+'interp_results.json'));FF={u['name']:u for u in json.load(open(B+'fullfree_results.json'))['universes']}
P=np.loadtxt('../../../../reproduce_cmb/planckTT.txt');PL={int(r[0]):r[1:] for r in P}
def pchi(Dl,a=2,b=30):
    return float(sum(((Dl[l]-PL[l][0])/(PL[l][2] if Dl[l]>PL[l][0] else PL[l][1]))**2 for l in range(a,b+1)))
rows=[]
for k,v in R.items():
    src=IR.get(k) or IR.get('path10_RED' if k=='grid_C3_0.1353_f20000' else k) or FF.get(k[3:])
    f=FF.get(k[3:]) or next((u for u in FF.values() if u['name']==k),None)
    Dl=np.array(src['Dl']);old=src.get('chi2') if 'chi2' in src else src['mig']['chi2']
    old360=f['mig360']['chi2'] if f else None
    rows.append(dict(name=k,ratio=v['ratio'],exp=v['exp']['chi2'],expR2=v['exp']['R2'],lC=v['exp']['lC'],bC=v['exp']['bC'],amp=v['exp']['amp'],
      full=v['full']['chi2'],fullR2=v['full']['R2'],th=v['full']['theta'],z0exp=v['z0']['exp']['chi2'],z0full=v['z0']['full']['chi2'],old=old,old360=old360,
      cmb=pchi(Dl),D=Dl[2:7].tolist()))
    fig,ax=plt.subplots(figsize=(7,4));ax.errorbar(L,Y,S,fmt='o',ms=3,color='k',label='Migkas 2020 Fig.6')
    ax.plot(L,v['exp']['A'],'b-',label=f"expansion only, real z: chi2 {v['exp']['chi2']:.2f}")
    ax.plot(L,v['full']['A'],'r-',label=f"expansion + bulk-flow terms: chi2 {v['full']['chi2']:.2f}")
    ax.set_xlabel('Galactic longitude (deg)');ax.set_ylabel('A/A_all');ax.set_title(k+f"  (z->0 ratio {v['ratio']:.4f})",fontsize=9);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(f'figs/{k}_AAall.png',dpi=80);plt.close(fig)
json.dump(rows,open('report_rows.json','w'))
def fm(x):return '-' if x is None else (f'{x:.3g}' if abs(x)>=1e4 or abs(x)<0.01 else f'{x:.2f}')
for key in ['full','exp']:
    print('\n## ranked by',key);print('|rank|universe|ratio z->0|exp chi2|exp R2|CoE l,b (exp)|full chi2|full R2|full lC,bC,bf,vp,lpp,sb|binned-z0 exp/full|old z->0 chi2 (l0 -60..100 / 0..360)|CMB chi2 2-30|D2..D6|')
    for i,r in enumerate(sorted(rows,key=lambda r:r[key])):
        t=r['th'];print(f"|{i+1}|{r['name']}|{r['ratio']:.4f}|{r['exp']:.2f}|{r['expR2']:.3f}|{r['lC']:.1f},{r['bC']:.1f}|{r['full']:.2f}|{r['fullR2']:.3f}|{t[0]:.1f},{t[1]:.1f},{t[2]:.4f},{t[3]:.0f},{t[4]:.1f},{t[5]:.2f}|{r['z0exp']:.2f}/{r['z0full']:.2f}|{fm(r['old'])} / {fm(r['old360'])}|{fm(r['cmb'])}|{', '.join(fm(x) for x in r['D'])}|")
