import sys,pickle; sys.path.insert(0,'..'); import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from grcore import stats, norm
d={}; [d.update(pickle.load(open(f,'rb'))) for f in ['ltb_0.0_0.1_0.3_0.5_0.9_0.99.pkl','ltb_0.0001_0.0003_0.00123_0.003_0.01_0.03.pkl']]
V=sorted(k for k in d if d[k][1]); C='0.0946'
rms=lambda x: np.sqrt(np.mean(x[:-1]**2))
S={}
for v in V:
    a=norm(d[v][1][0.0946][0]); b=norm(d[v][1][0.0][0]); S[v]=(a,b,a-b)
a0=S[0.0][0]
print("v0 dt | core dip quad cold rms | base dip quad cold rms | core-base dip quad cold rms | change vs v0=0: rms, dipole")
for v in V:
    row=[v,d[v][0]]
    for x in S[v]: dp,q,_,mn=stats(x); row+=[dp,q,mn,rms(x)]
    dx=S[v][0]-a0; row+=[rms(dx),stats(dx)[0]]
    print(" ".join("%.3e"%r for r in row))
print("noise floor: no-core v0=0 rms %.2e max %.2e"%(rms(S[0.0][1]),abs(S[0.0][1]).max()))
# (a)
x=np.array([v for v in V if v>0]); fig,ax=plt.subplots(figsize=(8,5))
ax.plot(x,[rms(S[v][0]) for v in x],'o-',label='with CoE mass'); ax.plot(x,[rms(S[v][2]) for v in x],'s-',label='added by mass (core − no-core)')
ax.axhline(rms(S[0.0][0]),ls='--',color='gray',label='v0 = 0 (with CoE mass)'); ax.axvline(1.23357e-3,color='r',lw=1,label='observed 369.82 km/s (Planck 2018)')
ax.set_xscale('log'); ax.set_xlim(1e-4,0.9); ax.set_xlabel('peculiar speed through expanding space, v0/c'); ax.set_ylabel('RMS of ΔT/T around the ring')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=8,frameon=False); fig.tight_layout(); fig.savefig('../ltb_rms_vs_speed_log.png',dpi=150)
th=np.arange(361); P=[0.0,0.00123,0.01,0.1]; cols=plt.cm.viridis(np.linspace(0,0.85,4))
fig,ax=plt.subplots(figsize=(9,5.5))
for c,v in zip(cols,P):
    dp,q,_,_=stats(S[v][0]); ax.plot(th,S[v][0],color=c,lw=1.6,ls='--' if v==0.00123 else '-',label=f"{v:g}c{' (observed)' if v==0.00123 else ''}: dip {dp:.4f}, quad {q:.4f}, RMS {rms(S[v][0]):.4f}")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_ylabel('ΔT/T')
ax.set_xlabel('LTB a(r,t): space expanding (3c recession) + motion through space, CoE at 180 deg')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=8,frameon=False); fig.tight_layout(); fig.savefig('../ltb_levels_lowspeed.png',dpi=150)
pickle.dump(S,open('low_S.pkl','wb'))
