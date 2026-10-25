import sys,pickle; sys.path.insert(0,'..'); import numpy as np, healpy as hp, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S=pickle.load(open('low_S.pkl','rb')); T0=2.7255
u=lambda l,b: np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(15.166,30.0); dvec=u(264.021,48.253); beta=3362.08e-6/T0
P=dvec-(dvec@A)*A; P/=np.linalg.norm(P)
th=np.radians(np.arange(361)); n=np.outer(np.cos(th),A)+np.outer(np.sin(th),P)
I=hp.read_map('/workspace/planck/smica.fits',field=0,nest=True); M=hp.read_map('/workspace/planck/mask.fits',field=0,nest=True); ns=hp.get_nside(I)
cmb=np.empty(361); mf=np.empty(361)
for i,v in enumerate(n):
    pix=hp.query_disc(ns,v,np.radians(1.0),nest=True); cmb[i]=I[pix].mean(); mf[i]=1-M[pix].mean()
b=np.degrees(np.arcsin(n[:,2])); bad=(mf>0.5)|(np.abs(b)<10)
dip=beta*(n@dvec); full=cmb/T0+dip; full_c=full-full[:-1].mean()
rms=lambda x,m=None: np.sqrt(np.mean((x[:-1] if m is None else x[:-1][~m[:-1]])**2))
print("ring mean of full obs %.3e; masked/flagged directions %d of 360; toy angles flagged:"%(full[:-1].mean(),bad[:-1].sum()))
idx=np.where(bad[:-1])[0]; print([ (a,b_) for a,b_ in zip(idx[np.r_[True,np.diff(idx)>1]], idx[np.r_[np.diff(idx)>1,True]])])
print("obs RMS: dipole-only %.4e (unmasked %.4e); full %.4e (unmasked %.4e); CMB-only part %.4e (unmasked %.4e)"%(rms(dip),rms(dip,bad),rms(full_c),rms(full_c,bad),rms(cmb/T0),rms(cmb/T0,bad)))
for v in sorted(S): print("%g resid RMS all %.4e unmasked %.4e | core-minus-baseline all %.4e unmasked %.4e"%(v,rms(S[v][0]-full_c),rms(S[v][0]-full_c,bad),rms(S[v][2]-full_c),rms(S[v][2]-full_c,bad)))
P4=[0.0,0.00123,0.01,0.1]; cols=plt.cm.viridis(np.linspace(0,0.85,4)); x=np.arange(361)
def shade(a):
    for s,e in [(s,e) for s,e in zip(idx[np.r_[True,np.diff(idx)>1]], idx[np.r_[np.diff(idx)>1,True]])]: a.axvspan(s,e+1,color='0.85',lw=0)
fig,ax=plt.subplots(3,1,figsize=(9,10),sharex=True)
for c,v in zip(cols,P4):
    ls='--' if v==0.00123 else '-'; nm=f"{v:g}c"+(" (observed speed)" if v==0.00123 else "")
    ax[0].plot(x,full_c,color=c,ls=ls,lw=1.0,label=f"{nm}: obs RMS {rms(full_c):.2e}, calc RMS {rms(S[v][0]):.4f}, resid RMS {rms(S[v][0]-full_c):.4f} (unmasked {rms(S[v][0]-full_c,bad):.4f})")
    ax[1].plot(x,S[v][0],color=c,ls=ls,lw=1.6); ax[2].plot(x,S[v][0]-full_c,color=c,ls=ls,lw=1.2)
for a,t in zip(ax,['observed CMB ΔT/T, full sky (Planck 2018 SMICA + dipole, 1° disc mean)','calculated ΔT/T (LTB, with CoE mass)','residual = calculated − observed']):
    shade(a); a.set_ylabel('ΔT/T'); a.set_title(t,fontsize=9); a.axhline(0,color='k',lw=.4)
ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.12),fontsize=7.5,frameon=False)
ax[2].set_xlim(0,360); ax[2].set_xticks(range(0,361,45)); ax[2].set_xlabel('LTB a(r,t): space expanding (3c recession) + motion through space, CoE at 180 deg')
fig.suptitle('Great circle through CoE axis and CMB dipole apex; grey = common mask or |b|<10°',fontsize=9,y=0.995)
fig.tight_layout(); fig.savefig('../ltb_residual_stacked.png',dpi=150)
fig,ax=plt.subplots(figsize=(9,4.5)); shade(ax)
ax.plot(x,dip,'C3',lw=1.8,label=f"dipole only: RMS {rms(dip):.3e}"); ax.plot(x,full_c,'k',lw=1,label=f"full sky (SMICA + dipole): RMS {rms(full_c):.3e} (unmasked {rms(full_c,bad):.3e})")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_ylabel('observed ΔT/T'); ax.set_xlabel('toy angle on the great circle (0 = anti-CoE, 180 = CoE)')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=1,fontsize=8,frameon=False); fig.tight_layout(); fig.savefig('../ltb_observed_dipole_vs_full.png',dpi=150)
