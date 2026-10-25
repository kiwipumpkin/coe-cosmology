import numpy as np,v11,base as V,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
M=0.008532635720463873;R=0.5;F=v11.field(M,R);wo,t0,T,W,S=F
X=lambda w:V.Rfun(t0,w,M)
tx=brentq(lambda t:V.Rfun(t,wo,M)-R,3e-3,t0)
times=[0.30,0.36,0.42,tx,0.569,0.687,t0];rho=0.25;rows=[]
for t in times:
    wr=brentq(lambda w:V.Rfun(t,w,M)-R,W[0]+1e-6,40);x=X(wo)-X(wr)   # our distance beyond the surface (today-distance units)
    small=big=0.
    if t>tx+1e-6:
        o,P=v11.trace(F,R,[0,90,180],t,keep=True);bd=[];sd=[]
        for (ps,tt,ww,pp),oo in zip(P,o):
            q=lambda w,p:np.hypot(X(w)*np.cos(p)-X(wo),X(w)*np.sin(p))
            bd.append(q(oo[4],oo[5]))
            j=np.argmin(abs(tt-tx)) if tt[-1]<tx else -1;sd.append(q(ww[j],pp[j]))
        big=float(np.mean([bd[0],bd[1],bd[1],bd[2]]));small=float(np.mean([sd[0],sd[1],sd[1],sd[2]]))
    rows.append((t,x,small,big))
fig,ax=plt.subplots(figsize=(8,20));
for i,(t,x,sm,bg) in enumerate(rows):
    y=-i*1.25;c=plt.Circle((x,y),rho,fill=False,lw=1.5);ax.add_patch(c)
    ph=np.linspace(0,2*np.pi,300);cx,cy=x+rho*np.cos(ph),y+rho*np.sin(ph);m=cx<0
    if m.any():
        from matplotlib.patches import Polygon
        xs=np.clip(cx,None,0);ax.add_patch(Polygon(np.c_[xs,cy],closed=True,fill=False,hatch='///',lw=0))
    ax.plot(x,y,'ko',ms=7);ax.annotate('',xy=(x+rho+0.25,y),xytext=(x,y),arrowprops=dict(arrowstyle='->',lw=1.2))
    if sm>0: ax.add_patch(plt.Circle((x,y),sm,fill=False,ls='--',color='g',lw=1.5))
    if bg>0: ax.add_patch(plt.Circle((x,y),bg,fill=False,color='r',lw=1.5))
    lab=f't = {t:.3f} D/c'+('  we cross' if abs(t-tx)<1e-6 else '')+(f'\nsee {sm:.2f} D, CMB at {bg:.2f} D' if bg>0 else '')
    ax.text(1.25,y,lab,va='center',fontsize=10)
ax.axvline(0,color='k',lw=2);ax.text(-0.9,0.6,'PLASMA',fontsize=13);ax.text(0.3,0.6,'LIGHT',fontsize=13)
ax.text(-1.0,-len(rows)*1.25+0.3,'hatched = still plasma\nblack circle = fixed patch of space around us (0.25 D)\ngreen dashed = how far we can see since crossing\nred = CMB distance (recombined edge), mean of 4 directions',fontsize=9)
ax.set_xlim(-1.1,2.3);ax.set_ylim(-len(rows)*1.25+0.1,1.0);ax.set_aspect('equal');ax.axis('off')
ax.set_title('Steps of our crossing (v11, with core mass, R_rec = 0.5 D)\ndistances in today\'s units of D; CoE to the left')
fig.savefig('steps_handdrawn_style.png',dpi=110,bbox_inches='tight')
for r in rows: print([round(v,3) for v in r])
