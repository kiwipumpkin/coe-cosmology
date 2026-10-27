# (b) rms dT/T vs angle at stages 0 (pure random web), 100, 400 moves and final, seed 1. Same definitions as web/dTT_rms_vs_angle_web.py:
# left = rms in 15-deg psi bins (psi = healpix colatitude, 180 = CoE), dipole removed; right = sigma(theta)=sqrt(2[C(0)-C(theta)]), l>=2.
# Top row: residual (sky - Planck low-l target); bottom row: the sky itself (bg + web). Stage maps from replay/ snapshots (inverse2_replay.py).
import json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.special import eval_legendre
W='/workspace/coepdf/joint_solve/phase5/web/inverse/';NS=64;L=3*NS-1;ll=np.arange(L+1);th=np.linspace(0,180,361)
P=np.array([eval_legendre(l,np.cos(np.radians(th))) for l in ll]);psi=np.degrees(hp.pix2ang(NS,np.arange(hp.nside2npix(NS)))[0]);bins=np.arange(0,181,15)
def fin(x,n):
    if not np.all(np.isfinite(x)):raise SystemExit('NONFINITE '+n)
    return x
z0=np.load(W+'replay/snap_s1_0000.npz');T=z0['Tobs'];fz=np.load(W+'inv_s1.npz');nf=json.load(open(W+'inv_s1.json'))['n_moves']
st=[('0 (random)',z0['Rmap']),('100',np.load(W+'replay/snap_s1_0100.npz')['Rmap']),('400',np.load(W+'replay/snap_s1_0400.npz')['Rmap']),(f'final ({nf})',fz['Resid'])]
old=np.load('/workspace/coepdf/joint_solve/phase5/p5_web_RA0.01_M1e22.npz')
fig,ax=plt.subplots(2,2,figsize=(14,10));out={}
def draw(a0,a1,m,lab,c,ls='-'):
    nd=hp.remove_dipole(fin(m,lab));r=[np.sqrt(np.mean(nd[(psi>=b)&(psi<b+15)]**2)) for b in bins[:-1]];a0.step(bins,np.r_[r,r[-1]],ls,where='post',color=c,label=lab)
    cl=hp.anafast(nd,lmax=L);k=(2*ll+1)/(4*np.pi)*cl*(ll>=2);Ct=k@P;a1.plot(th[1:],np.sqrt(np.clip(2*(Ct[0]-Ct[1:]),0,None)),ls,color=c,label=lab);return float(nd.std())
for (lab,m),c in zip(st,plt.cm.viridis(np.linspace(0,0.9,4))):
    out[lab]=dict(resid=draw(ax[0,0],ax[0,1],m-T,f'{lab} moves',c),sky=draw(ax[1,0],ax[1,1],m,f'{lab} moves',c))
out['old_shell_web_resid']=draw(ax[0,0],ax[0,1],old['Resid']-old['Tobs'],'earlier shell-concentrated web','r',':')
out['planck_target']=draw(ax[1,0],ax[1,1],T,'Planck low-l target','k','--')
for i in range(2):
    for j in range(2):ax[i,j].set_yscale('log');ax[i,j].axhline(1e-5,color='k',ls=':',lw=0.7);ax[i,j].set_xlim(0,180);ax[i,j].legend(fontsize=7)
    ax[i,0].set_xlabel('ψ (deg; 180 = CoE)');ax[i,1].set_xlabel('θ (deg)');ax[i,0].set_ylabel('dT/T')
ax[0,0].set_title('residual (sky − Planck target): rms in 15° ψ bins, dipole removed');ax[0,1].set_title('residual σ(θ), l≥2')
ax[1,0].set_title('sky (bg + web): rms in 15° ψ bins, dipole removed');ax[1,1].set_title('sky σ(θ), l≥2')
fig.suptitle('Inverse problem seed 1: rms vs angle at 0 / 100 / 400 / final moves');fig.tight_layout();fig.savefig(W+'rms_vs_angle_stages.png',dpi=120)
json.dump(out,open(W+'rms_vs_angle_stages.json','w'),indent=1);print(json.dumps(out))
