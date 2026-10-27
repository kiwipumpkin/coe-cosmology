# Phase5 web test: dT/T rms vs angle BEFORE and AFTER the stochastic-web compensation bodies (p5_comp5.py output Resid = bg + bodies).
# Same definitions as routes/dTT_rms_vs_angle.py: panel 1 rms in 15-deg psi bins (dipole removed); panel 2 sigma(theta)=sqrt(2[C(0)-C(theta)]), l>=2.
# Also writes web_metrics.json: rms total / dipole-removed / l=2 amplitude of Resid, of Resid-Tobs (deviation from the Planck low-l target),
# comp max/rms, fraction of bodies at the MMAX cap, and body density contrast delta = M/(rho_crit(H0=68.64) * V_cell), V_cell = shell volume/NB.
import json,glob,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.special import eval_legendre
O='/workspace/coepdf/joint_solve/phase5/';NS=64;L=3*NS-1;ll=np.arange(L+1);th=np.linspace(0,180,361)
P=np.array([eval_legendre(l,np.cos(np.radians(th))) for l in ll]);psi=np.degrees(hp.pix2ang(NS,np.arange(hp.nside2npix(NS)))[0])
RHOC=2.775e11*(68.64/100)**2  # Msun/Mpc^3
def stats(m):
    a=hp.map2alm(m,lmax=L,iter=3);cl=hp.alm2cl(a);nd=m-hp.alm2map(hp.almxfl(a,(ll==1).astype(float)),NS)-m.mean()
    return a,cl,nd,dict(rms=float(m.std()),rms_nodip=float(nd.std()),l2_amp=float(np.sqrt(5*cl[2]/(4*np.pi))))
fig,ax=plt.subplots(1,2,figsize=(14,5.5));out={};bins=np.arange(0,181,15)
def draw(m,lab,c,ls):
    a,cl,nd,s=stats(m);r=[np.sqrt(np.mean(nd[(psi>=b)&(psi<b+15)]**2)) for b in bins[:-1]]
    ax[0].step(bins,np.r_[r,r[-1]],ls,where='post',color=c,label=lab)
    k=(2*ll+1)/(4*np.pi)*cl*(ll>=2);Ct=k@P;ax[1].plot(th[1:],np.sqrt(np.clip(2*(Ct[0]-Ct[1:]),0,None)),ls,color=c,label=lab);return s
cols=iter(plt.cm.tab10.colors)
for E in ['0.01','0.1']:
    c=next(cols);B=np.load(O+f'bg_RA_{E}_St_sm_ns64.npy');out[f'bg_RA{E}']=draw(B,f'A eps={E} bg',c,'-')
    for M,ls in (('1e22','--'),('1e24','-.'),('1e27',':')):
        f=O+f'p5_web_RA{E}_M{M}';
        try: z=np.load(f+'.npz');j=json.load(open(f+'.json'))
        except FileNotFoundError: continue
        Rs=z['Resid'];assert np.all(np.isfinite(Rs)),'nonfinite resid'
        s=draw(Rs,f'A eps={E} +web cap {M}',c,ls);d=stats(Rs-z['Tobs'])[3]
        R=j['R']*1e3;V=4/3*np.pi*((j['RMAXF']*R)**3-(0.075*R)**3)/j['NB'];dl=z['M']/(RHOC*V)
        s.update(vs_planck_target=d,comp_max=j['cancel_err_max'],comp_rms=j['cancel_err_rms'],frac_at_cap=j['frac_at_bound'],
          max_abs_delta=float(np.abs(dl).max()),min_delta=float(dl.min()),max_delta=float(dl.max()),n_void_delta_lt_m1=int((dl<-1).sum()),
          M_pos=j['M_pos'],M_neg=j['M_neg'],net_delta_over_shell=float(z['M'].sum()/(RHOC*V*j['NB'])));out[f'web_RA{E}_M{M}']=s
for i,t in enumerate(['rms dT/T in 15° ψ bins, dipole removed','σ(θ), l≥2']):
    ax[i].set_yscale('log');ax[i].axhline(1e-5,color='k',ls=':');ax[i].axhline(1e-3,color='gray',ls='-.');ax[i].set_xlim(0,180);ax[i].set_title(t)
ax[0].set_xlabel('ψ from anti-CoE (deg; 180=CoE)');ax[1].set_xlabel('θ (deg)');ax[0].set_ylabel('dT/T');ax[0].legend(fontsize=7)
fig.suptitle('Route A backgrounds before/after random cluster–void web compensation (lines: 1e-3 target, 1e-5 Planck)');fig.tight_layout()
fig.savefig(O+'web/dTT_rms_vs_angle_web.png',dpi=130);json.dump(out,open(O+'web/web_metrics.json','w'),indent=1);print(json.dumps(out,indent=1))
