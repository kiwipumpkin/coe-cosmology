import json,glob,os,numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
rows={r['name']:r for r in json.load(open('../zdist/report_rows.json'))}
R=[]
for f in sorted(glob.glob('res/*.json')):
    k=os.path.basename(f)[:-5];o=json.load(open(f));r=rows.get(k,{})
    R.append(dict(name=k,G=o['G'],F=o['F'],ratio=o['ratio'],old=r.get('old'),zdist=r.get('full'),zdist_exp=r.get('exp'),cmb=r.get('cmb'),D=(r.get('D') or [None]*5)[:5]))
json.dump([{a:({c:d for c,d in b.items() if c!='A'} if a in 'GF' else b) for a,b in x.items()} for x in R],open('ring_table.json','w'),indent=0)
f=lambda v,n=3:'—' if v is None else (f'{v:.{n}g}')
for key in 'FG':
    print('\n##',key);S=sorted(R,key=lambda x:x[key]['chi2'])
    print('rank|universe|chi2|R2|l0|z|R(DA,Gpc)|bf|vp|lpp|sb|z0 chi2|zdist chi2|CMB chi2|D2|D3|D4|D5|D6')
    for i,x in enumerate(S):
        g=x[key];print(i+1,x['name'],f(g['chi2'],4),f(g['R2']),f(g['l0'],4),f(g['z']),f(g['DA_Gpc']),*( [f(g['bf']),f(g['vp'],4),f(g['lpp'],4),f(g['sb'])] if key=='F' else ['','','','']),f(x['old'],4),f(x['zdist'] if key=='F' else x['zdist_exp'],4),f(x['cmb']),*[f(d) for d in x['D']],sep='|')
M=__import__('csv').DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv'));M=list(M)
L=np.array([float(r['l_deg']) for r in M]);Y=np.array([float(r['A_over_Aall']) for r in M]);lo=[float(r['A_err_low']) for r in M];hi=[float(r['A_err_high']) for r in M]
for x in sorted(R,key=lambda x:x['F']['chi2'])[:4]+[y for y in R if y['name'] in('path00_REF','grid_C3_0.1353_f20000')]:
    o=json.load(open('res/%s.json'%x['name']));fig,ax=plt.subplots(figsize=(8,4.5))
    ax.errorbar(L,Y,yerr=[lo,hi],fmt='o',ms=3,c='k',label='Migkas 2020 Fig.6')
    ax.plot(L,o['F']['A'],'r-',label=f"ring fit + Migkas terms χ²={o['F']['chi2']:.2f}, z={o['F']['z']:.3g}, ℓ0={o['F']['l0']:.1f}")
    ax.plot(L,o['G']['A'],'b--',label=f"ring geometry only χ²={o['G']['chi2']:.2f}, z={o['G']['z']:.3g}, ℓ0={o['G']['l0']:.1f}")
    ax.set_xlabel('galactic longitude l (deg)');ax.set_ylabel('A/A_all');ax.set_title(x['name']+f"  CMB χ²(ℓ2–30)={f(x['cmb'])}");ax.legend(fontsize=8);fig.tight_layout();fig.savefig('ring_%s.png'%x['name'],dpi=110);plt.close()
