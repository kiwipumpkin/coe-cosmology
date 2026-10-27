# Seed robustness of the web compensation (Route A eps=0.01, cap 1e22). Same metric definitions as dTT_rms_vs_angle_web.py.
# Includes the original seed 2000 run (p5_web_RA0.01_M1e22) plus run_seeds.sh seeds. Writes seed_metrics.json and seed_rms.png.
import json,glob,re,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase5/';NS=64;L=3*NS-1;ll=np.arange(L+1);RHOC=2.775e11*(68.64/100)**2
def stats(m):
    a=hp.map2alm(m,lmax=L,iter=3);cl=hp.alm2cl(a);nd=m-hp.alm2map(hp.almxfl(a,(ll==1).astype(float)),NS)-m.mean()
    return float(nd.std()),float(np.sqrt(5*cl[2]/(4*np.pi)))
fs=[O+'p5_web_RA0.01_M1e22.json']+sorted(glob.glob(O+'p5_web_RA0.01_M1e22_s*.json'))
rows={}
for f in fs:
    j=json.load(open(f));z=np.load(f[:-5]+'.npz');Rs=z['Resid'];assert np.all(np.isfinite(Rs))
    r,l2=stats(Rs);R=j['R']*1e3;V=4/3*np.pi*((j['RMAXF']*R)**3-(0.075*R)**3)/j['NB'];dl=z['M']/(RHOC*V)
    rows[j['seed']]=dict(rms_nodip=r,l2_amp=l2,comp_max=j['cancel_err_max'],comp_rms=j['cancel_err_rms'],frac_at_cap=j['frac_at_bound'],
        max_abs_delta=float(np.abs(dl).max()),n_delta_lt_m1=int((dl<-1).sum()))
k=list(next(iter(rows.values())).keys());summ={q:dict(mean=float(np.mean([v[q] for v in rows.values()])),std=float(np.std([v[q] for v in rows.values()],ddof=1)) if len(rows)>1 else 0.,
    min=float(min(v[q] for v in rows.values())),max=float(max(v[q] for v in rows.values()))) for q in k}
json.dump(dict(per_seed=rows,summary=summ,n=len(rows)),open(O+'web/seed_metrics.json','w'),indent=1)
for s,v in rows.items(): print(s,' '.join(f'{q}={v[q]:.3g}' for q in k))
for q in k: print(q,summ[q])
S=list(rows);y=[rows[s]['rms_nodip'] for s in S]
fig,ax=plt.subplots(figsize=(8,5));ax.semilogy(range(len(S)),y,'o',ms=9);ax.set_xticks(range(len(S)));ax.set_xticklabels(S)
ax.axhline(1e-3,color='gray',ls='-.',label='1e-3 target');ax.axhline(1e-5,color='k',ls=':',label='Planck ~1e-5')
ax.axhline(np.mean(y),color='C1',ls='--',label=f'mean {np.mean(y):.2e}');ax.set_xlabel('random seed (body placement)');ax.set_ylabel('dT/T rms, dipole removed')
ax.set_title('Web compensation vs seed: Route A eps=0.01, cap 1e22');ax.legend();fig.tight_layout();fig.savefig(O+'web/seed_rms.png',dpi=130)
