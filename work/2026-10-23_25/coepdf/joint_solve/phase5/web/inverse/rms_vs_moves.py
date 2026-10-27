# (a) Residual vs number of moves, seed 1. Residual = (bg + web sky) - Planck low-l target Tobs (l=2..40 tapered, no monopole/dipole).
# Dipole-removed rms from the deterministic replay snapshots (replay/snap_s1_*.npz, written by inverse2_replay.py, which reproduces
# inv_s1.log cycle-by-cycle); move 0 = pure uniform random web. Also plots the log's cancel_rms (std of Resid-Tobs, dipole INCLUDED).
import glob,re,ast,json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
W='/workspace/coepdf/joint_solve/phase5/web/inverse/'
def fin(x,n):
    if not np.all(np.isfinite(x)):raise SystemExit('NONFINITE '+n)
    return x
z0=np.load(W+'replay/snap_s1_0000.npz');T=z0['Tobs'];B=z0['B'];nm,rd,rs=[],[],[]
for f in sorted(glob.glob(W+'replay/snap_s1_*.npz')):
    z=np.load(f);m=fin(z['Rmap'],f);nm.append(int(f[-8:-4]));rd.append(float(hp.remove_dipole(m-T).std()));rs.append(float(hp.remove_dipole(m).std()))
# NOTE: replay files are overwritten per n_moves; repeated LM-only cycles at the same n_moves keep the last one.
cy=[ast.literal_eval(l[6:]) for l in open(W+'inv_s1.log') if l.startswith('cycle ')]
fz=np.load(W+'inv_s1.npz');Rf=fin(fz['Resid'],'Resid');fin_nd=float(hp.remove_dipole(Rf-T).std());nf=json.load(open(W+'inv_s1.json'))['n_moves']
old=np.load('/workspace/coepdf/joint_solve/phase5/p5_web_RA0.01_M1e22.npz');old_nd=float(hp.remove_dipole(old['Resid']-old['Tobs']).std())
fig,ax=plt.subplots(figsize=(9,5.5))
ax.plot(nm,rd,'o-',ms=3,label='residual (sky − Planck target), dipole removed [replay]')
ax.plot([0]+[c['n_moves'] for c in cy],[float((z0['Rmap']-T).std())]+[c['cancel_rms'] for c in cy],'s--',ms=2,alpha=0.6,label='log cancel_rms (dipole included)')
ax.plot(nf,fin_nd,'r*',ms=12,label=f'final (n_moves={nf}): {fin_nd:.2e}')
ax.axhline(old_nd,color='g',ls=':',label=f'earlier shell-concentrated free-mass web: {old_nd:.2e}')
ax.axhline(float(hp.remove_dipole(T).std()),color='k',ls='-.',lw=0.8,label='rms of Planck target itself')
ax.set_yscale('log');ax.set_xlabel('number of moves (greedy moves; LM refinements do not add moves)');ax.set_ylabel('rms dT/T')
ax.set_title('Inverse problem seed 1: residual vs moves (move 0 = uniform random web, 3000 bodies)');ax.legend(fontsize=8);fig.tight_layout()
fig.savefig(W+'rms_vs_moves.png',dpi=130)
json.dump(dict(n_moves=nm,rms_nodip_resid=rd,rms_nodip_sky=rs,final_n=nf,final_rms_nodip_resid=fin_nd,old_web_rms_nodip_resid=old_nd),open(W+'rms_vs_moves.json','w'),indent=1)
print(nm[0],rd[0],nf,fin_nd,old_nd)
