# dT/T vs angle and rms vs angular scale for phase-5 backgrounds (all SURF=t).
# Definition: dT/T(n) = (1+z_mean)/(1+z(n)) - 1, with the mean taken over the full sky
# (solid-angle weighted), i.e. monopole removed. This equals the bg_*_sm_ns64.npy maps
# built by p2_bg5.py (sm = exp(-spline ln(1+z)) / mean - 1). psi = angle from the
# anti-CoE axis (psi=0 anti-CoE, psi=180 CoE), axially symmetric maps (Nside 64).
# Panel 1: dT/T(psi) from the ray traces (37 rays, 5 deg), monopole removed.
# Panel 2: rms in 15-deg psi bins of the dipole-removed map (solid-angle weighted).
# Panel 3: rms temperature difference vs separation theta,
#   sigma(theta) = sqrt(<[dT/T(n1)-dT/T(n2)]^2>) = sqrt(2[C(0)-C(theta)]),
#   C(theta) = sum_l (2l+1)/(4pi) C_l P_l(cos theta), C_l from the map (lmax=191);
#   with dipole (l>=1) and without dipole (l>=2).
import json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.special import eval_legendre
O='/workspace/coepdf/joint_solve/phase5/'
C=[('RA_0.1_St','Route A eps=0.1'),('RA_0.01_St','Route A eps=0.01'),('RA_0.001_St','Route A eps=1e-3'),
   ('RB_St','Route B (CatWISE slope)'),('K-1_St','Migkas K-1 (orig.)'),('VALpolyA_St','old poly (ref)')]
NS=64;pix_th,_=hp.pix2ang(NS,np.arange(hp.nside2npix(NS)));psi_pix=np.degrees(pix_th)
th=np.linspace(0,180,361);L=3*NS-1;ll=np.arange(L+1)
P=np.array([eval_legendre(l,np.cos(np.radians(th))) for l in ll])
fig,ax=plt.subplots(1,3,figsize=(18,5.5));out={}
for (T,lab),c in zip(C,plt.cm.tab10.colors):
    sm=np.load(O+f'bg_{T}_sm_ns64.npy');o=np.load(O+f'bg_{T}_trace.npz')['o'];o=o[np.argsort(o[:,0])]
    w=np.sin(np.radians(o[:,0]));inv=1/o[:,1];m=np.trapezoid(inv*w,o[:,0])/np.trapezoid(w,o[:,0])
    ax[0].plot(o[:,0],np.abs(inv/m-1),'o-',ms=3,color=c,label=lab)
    a=hp.map2alm(sm,lmax=L,iter=3);cl=hp.alm2cl(a)
    d1=hp.alm2map(hp.almxfl(a,(ll==1).astype(float)),NS);nd=sm-d1
    bins=np.arange(0,181,15);r=[np.sqrt(np.mean(nd[(psi_pix>=b)&(psi_pix<b+15)]**2)) for b in bins[:-1]]
    ax[1].step(bins,np.r_[r,r[-1]],where='post',color=c,label=lab)
    for lmin,ls in ((1,'-'),(2,'--')):
        k=(2*ll+1)/(4*np.pi)*cl*(ll>=lmin);Ct=k@P;s=np.sqrt(np.clip(2*(Ct[0]-Ct),0,None))
        ax[2].plot(th[1:],s[1:],ls,color=c,label=lab if lmin==1 else None)
        if lmin==2:s2=s
        else:s1=s
    out[T]=dict(rms=float(sm.std()),rms_nodip=float(nd.std()),onepz=[float(o[:,1].min()),float(o[:,1].max())],
               sig_dip_90=float(s1[180]),sig_nodip_90=float(s2[180]),sig_nodip_10=float(s2[20]),sig_nodip_max=float(s2.max()))
for i,t in enumerate(['|dT/T|(psi), monopole removed','rms dT/T in 15° psi bins, dipole removed','rms difference vs separation θ (solid: with dipole, dashed: no dipole)']):
    ax[i].set_yscale('log');ax[i].axhline(1e-5,color='k',ls=':');ax[i].axhline(1e-3,color='gray',ls='-.')
    ax[i].text(2,1.2e-5,'Planck 1e-5',fontsize=8);ax[i].text(2,1.2e-3,'comp. target 1e-3',fontsize=8)
    ax[i].set_xlim(0,180);ax[i].set_title(t,fontsize=10);ax[i].set_ylim(1e-6,3)
ax[0].set_xlabel('ψ from anti-CoE (deg; 180 = CoE)');ax[1].set_xlabel('ψ from anti-CoE (deg)');ax[2].set_xlabel('θ separation (deg)')
ax[0].set_ylabel('dT/T');ax[0].legend(fontsize=8)
fig.suptitle('Phase 5 backgrounds, SURF=t: dT/T = (1+z_mean)/(1+z) - 1');fig.tight_layout()
fig.savefig(O+'routes/dTT_rms_vs_angle.png',dpi=130);json.dump(out,open(O+'routes/dTT_rms_vs_angle.json','w'),indent=1)
print(json.dumps(out,indent=1))
