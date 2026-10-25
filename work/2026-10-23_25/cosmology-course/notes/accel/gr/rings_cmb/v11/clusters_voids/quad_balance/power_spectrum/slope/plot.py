import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
r=json.load(open('slope/slope_release_rows.json'));old=json.load(open('shape/joint_shape_rows.json'))
P=sorted(set(x['p'] for x in r));S=sorted(set(x['near_logslope'] for x in r))
fig,ax=plt.subplots(1,3,figsize=(18,5.5))
for i,p in enumerate(P):
  g=sorted([x for x in r if x['p']==p],key=lambda x:x['near_logslope']);s=[x['near_logslope'] for x in g]
  ax[0].plot(s,[x['a2'] for x in g],'o-',label=f'curv p={p}');ax[1].semilogy(s,[x['D2'] for x in g],'o-',label=f'curv p={p}')
for a in ax[:2]: a.axvline(-0.04103732750534994,c='r',ls='--',label='CatWISE 3% lock');a.set_xlabel('near-observer d ln ρ / d ln r');a.legend(fontsize=8)
ax[0].axhline(0,c='k',lw=.5);ax[0].set_ylabel('a₂ (dT/T)');ax[1].axhline(226.1,c='gold',lw=2);ax[1].set_ylabel('D₂ [μK²]');ax[1].set_title('Planck D₂=226 (gold)')
cm=plt.cm.coolwarm
for j,s in enumerate(S):
  g=[x for x in r if x['near_logslope']==s];a2=np.array([x['a2'] for x in g]);a3=np.array([x['a3'] for x in g]);b,c=np.polyfit(a2,a3,1);xx=np.linspace(-0.05,0.065,2)
  ax[2].plot(a2,a3,'o',c=cm(j/(len(S)-1)));ax[2].plot(xx,c+b*xx,'-',c=cm(j/(len(S)-1)),lw=3 if abs(s+0.041)<1e-3 else 1,label=f'slope {s:+.3f}: a₃={c:+.2e}{b:+.3f}a₂')
ax[2].plot([o['a2'] for o in old],[o['a3'] for o in old],'k.',ms=4,label='earlier locked runs (curv/cub/soft)')
ax[2].plot(0,0,'k*',ms=15);ax[2].axhline(0,c='k',lw=.5);ax[2].axvline(0,c='k',lw=.5);ax[2].set_xlabel('a₂');ax[2].set_ylabel('a₃');ax[2].legend(fontsize=6.5);ax[2].set_title('a₃ vs a₂ lines: only slope≈0 passes through the origin')
plt.suptitle('Releasing the CatWISE near-slope: ln ρ = −s·ln x − p(x−1)², M=30, D=42');plt.savefig('slope_release_quad.png',dpi=110,bbox_inches='tight')
