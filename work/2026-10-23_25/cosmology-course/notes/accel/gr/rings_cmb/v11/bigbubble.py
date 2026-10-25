import numpy as np,v11,base as V,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
M=0.008532635720463873;R=0.5;F=v11.field(M,R);wo,t0,T,W,S=F
ps=np.arange(0,181,3.0);o,_=v11.trace(F,R,ps)
re=np.array([V.Rfun(te,w,M) for te,w in zip(o[:,2],o[:,4])]);rt=np.array([V.Rfun(t0,w,M) for w in o[:,4]]);phi=o[:,5]
print('max |r_emit-0.5| =',np.abs(re-R).max())
fig,ax=plt.subplots(figsize=(12,9));th=np.linspace(0,2*np.pi,300)
ax.fill(R*np.cos(th),R*np.sin(th),fill=False,hatch='///',color='orange');ax.plot(R*np.cos(th),R*np.sin(th),'k',lw=1.5,label='recombination surface R_rec = 0.5 D')
for s in (1,-1):
    P=s*phi
    for i in range(0,len(ps),3): ax.plot([0,rt[i]*np.cos(P[i])],[0,rt[i]*np.sin(P[i])],'k-',lw=.4,alpha=.5)
    ax.plot(re*np.cos(P),re*np.sin(P),'b.',ms=5,label='emission points (on R_rec)' if s>0 else None)
    ax.plot(rt*np.cos(P),rt*np.sin(P),'r-',lw=2,label='big bubble today: same matter moved out along its CoE radial line' if s>0 else None)
xo=V.Rfun(t0,wo,M);ax.plot(xo,0,'go',ms=9,label='us today (1.0 D)');ax.plot(0,0,'k+',ms=15,mew=2,label='CoE')
ax.plot([],[],'k-',lw=.4,label='radial lines from CoE')
ax.set_aspect('equal');ax.grid(alpha=.3);ax.legend(loc='upper left',fontsize=9)
ax.set_title('Big bubble check (v11, with core mass, R_rec=0.5 D, today t=0.924 D/c), rays every 3 deg; distances in D')
fig.savefig('bigbubble_fixed.png',dpi=110,bbox_inches='tight')
for i in range(0,len(ps),10): print(f"sky {ps[i]:5.0f}  CoE angle {np.degrees(phi[i]):6.1f}  emit r {re[i]:.3f} at t {o[i,2]:.3f}  shell today r {rt[i]:.3f}  dist from us {np.hypot(rt[i]*np.cos(phi[i])-xo,rt[i]*np.sin(phi[i])):.3f}")
