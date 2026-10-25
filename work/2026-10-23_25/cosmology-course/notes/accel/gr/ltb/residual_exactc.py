import sys,pickle; sys.path.insert(0,'..'); import numpy as np, healpy as hp, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from grcore import stats, norm
S=pickle.load(open('low_S.pkl','rb')); X=pickle.load(open('exactc.pkl','rb')); T0=2.7255
u=lambda l,b: np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(15.166,30.0); dvec=u(264.021,48.253); beta=3362.08e-6/T0
P=dvec-(dvec@A)*A; P/=np.linalg.norm(P); th=np.radians(np.arange(361)); n=np.outer(np.cos(th),A)+np.outer(np.sin(th),P)
I=hp.read_map('/workspace/planck/smica.fits',field=0,nest=True); M=hp.read_map('/workspace/planck/mask.fits',field=0,nest=True); ns=hp.get_nside(I)
cmb=np.empty(361); mf=np.empty(361)
for i,v in enumerate(n): pix=hp.query_disc(ns,v,np.radians(1.0),nest=True); cmb[i]=I[pix].mean(); mf[i]=1-M[pix].mean()
bad=(mf>0.5)|(np.abs(np.degrees(np.arcsin(n[:,2])))<10); full=cmb/T0+beta*(n@dvec); full-=full[:-1].mean()
rms=lambda x,m=None: np.sqrt(np.mean((x[:-1] if m is None else x[:-1][~m[:-1]])**2))
V=[0.0,0.00123,0.01,0.1,0.5]; new={}
for v in V:
    a=norm(X[(0.0946,v,{0.0:1.4,0.00123:1.4,0.01:1.4,0.1:1.6,0.5:2.4}[v])]); b=norm(X[(0.0,v,{0.0:1.4,0.00123:1.4,0.01:1.4,0.1:1.6,0.5:2.4}[v])]); new[v]=(a,b,a-b)
for v in V:
    for nm,Z in (('about c (0.9999c)',S[v]),('exactly c',new[v])):
        a,b,c=Z; s=stats(a)
        print("%-8g %-18s core dip %.5f quad %.5f cold %.5f rms %.5f | base dip %.2e rms %.2e | added dip %.5f quad %.5f rms %.5f | resid rms all %.5f unmasked %.5f"%(v,nm,s[0],s[1],s[3],rms(a),stats(b)[0],rms(b),stats(c)[0],stats(c)[1],rms(c),rms(a-full),rms(a-full,bad)))
    print("   max |exact - about| with core: %.2e"%np.abs(new[v][0]-S[v][0]).max())
idx=np.where(bad[:-1])[0]; segs=list(zip(idx[np.r_[True,np.diff(idx)>1]], idx[np.r_[np.diff(idx)>1,True]]))
P4=[0.0,0.00123,0.01,0.1]; cols=plt.cm.viridis(np.linspace(0,0.85,4)); x=np.arange(361)
fig,ax=plt.subplots(2,1,figsize=(9,8),sharex=True)
for c,v in zip(cols,P4):
    a=new[v][0]; ls='--' if v==0.00123 else '-'; nm=f"{v:g}c"+(" (measured CMB motion)" if v==0.00123 else "")
    ax[0].plot(x,a,color=c,ls=ls,lw=1.6,label=f"{nm}: calculated RMS {rms(a):.4f}, residual RMS {rms(a-full):.4f} (unmasked {rms(a-full,bad):.4f})")
    ax[1].plot(x,a-full,color=c,ls=ls,lw=1.2)
for a_,t in zip(ax,['calculated ΔT/T, shell edge receding at exactly c (with CoE mass)','residual = calculated − observed (Planck 2018 SMICA + dipole)']):
    for s,e in segs: a_.axvspan(s,e+1,color='0.85',lw=0)
    a_.set_ylabel('ΔT/T'); a_.set_title(t,fontsize=9); a_.axhline(0,color='k',lw=.4)
ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.1),fontsize=7.5,frameon=False)
ax[1].set_xlim(0,360); ax[1].set_xticks(range(0,361,45)); ax[1].set_xlabel('LTB a(r,t): space expanding (3c recession) + motion through space, CoE at 180 deg')
fig.suptitle('grey = Planck common mask or |b| < 10°',fontsize=8,y=0.995); fig.tight_layout(); fig.savefig('../ltb_exactc_stacked.png',dpi=150)
