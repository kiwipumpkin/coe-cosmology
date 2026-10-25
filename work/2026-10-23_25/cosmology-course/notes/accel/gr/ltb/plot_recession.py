import sys,pickle; sys.path.insert(0,'..'); import numpy as np, healpy as hp, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from grcore import stats, norm
T0=2.7255; u=lambda l,b: np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(15.166,30.0); dvec=u(264.021,48.253); beta=3362.08e-6/T0
P=dvec-(dvec@A)*A; P/=np.linalg.norm(P); th=np.radians(np.arange(361)); n=np.outer(np.cos(th),A)+np.outer(np.sin(th),P)
I=hp.read_map('/workspace/planck/smica.fits',field=0,nest=True); M=hp.read_map('/workspace/planck/mask.fits',field=0,nest=True); ns=hp.get_nside(I)
cmb=np.empty(361); mf=np.empty(361)
for i,v in enumerate(n): pix=hp.query_disc(ns,v,np.radians(1.0),nest=True); cmb[i]=I[pix].mean(); mf[i]=1-M[pix].mean()
bad=(mf>0.5)|(np.abs(np.degrees(np.arcsin(n[:,2])))<10); full=cmb/T0+beta*(n@dvec); full-=full[:-1].mean()
rms=lambda x,m=None: np.sqrt(np.nanmean((x[:-1] if m is None else x[:-1][~m[:-1]])**2))
R=pickle.load(open('recession.pkl','rb')); X3=pickle.load(open('exactc.pkl','rb')); R6=pickle.load(open('recession6_fixedmass.pkl','rb'))
DT3={0.0:1.4,0.00123:1.4,0.1:1.6}
cases={}
for v in (0.0,0.00123,0.1):
    cases[('1.5c','fraction 0.0946',v)]=(R[(1.5,0.0946,v)][1],R[(1.5,0.0,v)][1],R[(1.5,0.0,v)][0])
    cases[('3c','fraction 0.0946',v)]=(X3[(0.0946,v,DT3[v])],X3[(0.0,v,DT3[v])],DT3[v])
    k=(0.0946/4,v)
    if k in R6 and not isinstance(R6[k],str): cases[('6c','same core mass as 3c (fraction 0.0237)',v)]=(R6[k][1],R[(6.0,0.0,v)][1],R6[k][0])
for k,(a,b,dt) in cases.items():
    nn=int(np.isnan(a).sum())
    if nn: a=a.copy(); bb=np.isnan(a); a[bb]=np.interp(np.where(bb)[0],np.where(~bb)[0],a[~bb]); cases[k]=(a,cases[k][1],dt)
    a=norm(a); b=norm(b); s=stats(np.nan_to_num(a)); c=a-b
    print(k,'view %.1f'%dt,'missing',nn,'core dip %.5f quad %.5f cold %.5f rms %.5f | base dip %.2e | added dip %.5f quad %.5f rms %.5f | resid all %.5f unmasked %.5f'%(s[0],s[1],np.nanmin(a),rms(a),stats(b)[0],stats(np.nan_to_num(c))[0],stats(np.nan_to_num(c))[1],rms(c),rms(a-full),rms(a-full,bad)))
idx=np.where(bad[:-1])[0]; segs=list(zip(idx[np.r_[True,np.diff(idx)>1]], idx[np.r_[np.diff(idx)>1,True]]))
plt.rcParams.update({'font.size':13}); fig,ax=plt.subplots(2,1,figsize=(12,10),sharex=True); x=np.arange(361); cols=['#1f77b4','#2ca02c','#d62728']
for c,lab in zip(cols,['1.5c','3c','6c']):
    k=[k for k in cases if k[0]==lab and k[2]==0.00123]
    if not k: continue
    a=norm(cases[k[0]][0]); ax[0].plot(x,a,color=c,lw=2.4,label=f"space carrying us at {lab}: calculated RMS {rms(a):.4f}, leftover RMS {rms(a-full):.4f}")
    ax[1].plot(x,a-full,color=c,lw=2.0)
for a_,t in zip(ax,['calculated ΔT/T at the measured CMB motion (0.00123c), shell edge at exactly c','leftover = calculated − observed (Planck 2018 SMICA + dipole)']):
    for s,e in segs: a_.axvspan(s,e+1,color='0.88',lw=0)
    a_.set_ylabel('ΔT/T'); a_.set_title(t,fontsize=12); a_.axhline(0,color='k',lw=.5)
ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.08),fontsize=12,frameon=False)
ax[1].set_xlim(0,360); ax[1].set_xticks(range(0,361,45)); ax[1].set_xlabel('direction on the ring (deg): 0 = away from CoE, 180 = toward CoE; grey = masked sky')
fig.tight_layout(); fig.savefig('../ltb_recession_stacked.png',dpi=150)
