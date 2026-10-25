import sys,pickle; sys.path.insert(0,'..'); import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S=pickle.load(open('low_S.pkl','rb'))
u=lambda l,b: np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(15.166,30.0)            # anti-CoE (toy 0 deg); CoE = -A (toy 180 deg)
dvec=u(264.021,48.253); beta=3362.08e-6/2.7255
P=dvec-(dvec@A)*A; P/=np.linalg.norm(P)   # toy 90 deg: the side toward the dipole apex
th=np.radians(np.arange(361)); n=np.outer(np.cos(th),A)+np.outer(np.sin(th),P)
psi=np.degrees(np.arccos(dvec@A)); print("angle anti-CoE to dipole apex %.2f deg -> apex at toy angle %.2f"%(psi,psi))
obs=beta*(n@dvec)
lb=lambda x:(np.degrees(np.arctan2(x[1],x[0]))%360,np.degrees(np.arcsin(x[2])))
print("pole of great circle (l,b)=(%.2f,%.2f)"%lb(np.cross(A,P)))
rms=lambda x: np.sqrt(np.mean(x[:-1]**2))
print("observed RMS %.3e"%rms(obs))
for v in sorted(S): print("%g resid RMS core %.4e | core-minus-baseline %.4e"%(v,rms(S[v][0]-obs),rms(S[v][2]-obs)))
P4=[0.0,0.00123,0.01,0.1]; cols=plt.cm.viridis(np.linspace(0,0.85,4)); x=np.arange(361)
fig,ax=plt.subplots(3,1,figsize=(9,10),sharex=True)
for c,v in zip(cols,P4):
    ls='--' if v==0.00123 else '-'; nm=f"{v:g}c"+(" (observed)" if v==0.00123 else "")
    ax[0].plot(x,obs,color=c,ls=ls,lw=1.6,label=f"{nm}: obs RMS {rms(obs):.2e}, calc RMS {rms(S[v][0]):.4f}, resid RMS {rms(S[v][0]-obs):.4f}")
    ax[1].plot(x,S[v][0],color=c,ls=ls,lw=1.6); ax[2].plot(x,S[v][0]-obs,color=c,ls=ls,lw=1.6)
for a,t in zip(ax,['observed CMB ΔT/T (Planck 2018 dipole only)','calculated ΔT/T (LTB, with CoE mass)','residual = calculated − observed']):
    a.set_ylabel('ΔT/T'); a.set_title(t,fontsize=9); a.axhline(0,color='k',lw=.4); a.axvline(psi,color='r',lw=.6)
ax[0].legend(loc='lower center',bbox_to_anchor=(0.5,1.12),ncol=1,fontsize=8,frameon=False)
ax[2].set_xlim(0,360); ax[2].set_xticks(range(0,361,45))
ax[2].set_xlabel('LTB a(r,t): space expanding (3c recession) + motion through space, CoE at 180 deg')
fig.suptitle('Great circle through the CoE axis and the CMB dipole apex (red line: apex at 80.5 deg)',fontsize=9,y=0.995)
fig.tight_layout(); fig.savefig('../ltb_residual_stacked.png',dpi=150)
