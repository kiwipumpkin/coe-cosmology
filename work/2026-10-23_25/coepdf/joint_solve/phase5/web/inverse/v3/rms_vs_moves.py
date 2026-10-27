# v3 residual vs number of greedy changes (each = one body moved and/or re-massed). Residual = (bg + web sky) - Planck target (SMICA l=2..130),
# dipole removed, from per-cycle snapshots snaps/s<seed>_m<nnnn>.npz written by inverse3.py (move 0 = pure uniform random web).
# Reference lines: v2 seed-1 final (../inv_s1.json rms_nodip_end), earlier shell-concentrated web, Planck target rms (l<=40 v2 target and l<=130 v3 target).
# Also v3 residual restricted to l=2..40 (comparable band to v2). Usage: python3 rms_vs_moves.py [dir/]. Fails on NaN.
import sys,glob,json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
W='/workspace/coepdf/joint_solve/phase5/web/inverse/v3/'+(sys.argv[1] if len(sys.argv)>1 else '');P5='/workspace/coepdf/joint_solve/phase5/'
def fin(x,n):
    if not np.all(np.isfinite(x)):raise SystemExit('NONFINITE '+n)
    return x
T=fin(np.load(P5+'web/inverse/v3/target_Tobs_l130.npy'),'T');f40=np.zeros(192);f40[2:41]=1
def b40(m):return float(hp.alm2map(hp.almxfl(hp.map2alm(m,lmax=191,iter=3),f40),64,lmax=191).std())
v2=json.load(open(P5+'web/inverse/inv_s1.json'))['rms_nodip_end'];old=np.load(P5+'p5_web_RA0.01_M1e22.npz');old_nd=float(hp.remove_dipole(old['Resid']-old['Tobs']).std())
fig,ax=plt.subplots(figsize=(9.5,6));out={}
seeds=sorted({f.split('/')[-1].split('_')[0] for f in glob.glob(W+'snaps/s*_m*.npz')})
for s in seeds:
    fs=sorted(glob.glob(W+f'snaps/{s}_m*.npz'));nm=[];rd=[];r40=[]
    for f in fs:z=np.load(f);m=fin(z['Rmap'],f);nm.append(int(z['nmv']));rd.append(float(hp.remove_dipole(m-T).std()));r40.append(b40(m-T))
    ax.plot(nm,rd,'o-',ms=3,label=f'v3 {s}: residual, dipole removed (final {rd[-1]:.2e} at {nm[-1]})');ax.plot(nm,r40,'--',lw=1,label=f'v3 {s}: residual in l=2-40 only (final {r40[-1]:.2e})')
    out[s]=dict(n_moves=nm,resid_nodip=rd,resid_l2_40=r40)
ax.axhline(v2,color='m',ls=':',label=f'v2 seed 1 final (positions only, l<=60): {v2:.2e}')
ax.axhline(old_nd,color='g',ls=':',label=f'earlier shell-concentrated free-mass web: {old_nd:.2e}')
ax.axhline(1.64e-5,color='k',ls='-.',lw=0.8,label='rms of v2 Planck target (l<=40): 1.64e-5');ax.axhline(float(T.std()),color='k',ls='--',lw=0.8,label=f'rms of v3 Planck target (l<=130): {T.std():.2e}')
ax.set_yscale('log');ax.set_xlabel('number of greedy changes (LM refinements add none)');ax.set_ylabel('rms dT/T')
ax.set_title('v3 inverse web (position + mass free, fit to l<=130): residual vs changes');ax.legend(fontsize=7);fig.tight_layout();fig.savefig(W+'rms_vs_moves.png',dpi=130)
json.dump(dict(seeds=out,v2_s1=v2,old_web=old_nd,target_rms_l130=float(T.std())),open(W+'rms_vs_moves.json','w'),indent=1);print({s:(v['n_moves'][-1],v['resid_nodip'][-1],v['resid_l2_40'][-1]) for s,v in out.items()})
