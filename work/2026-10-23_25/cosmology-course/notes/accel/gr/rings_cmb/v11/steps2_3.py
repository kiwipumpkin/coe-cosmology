import numpy as np,v11,base as V,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
M=0.008532635720463873;R=0.5;F=v11.field(M,R);wo,t0,T,W,S=F
tx=brentq(lambda t:V.Rfun(t,wo,M)-R,3e-3,t0);tr=0.30;ro=V.Rfun(tr,wo,M)
ph=np.linspace(0,2*np.pi,240);rings=[]
for rho in (0.05,0.10,0.15):
    px,py=ro+rho*np.cos(ph),rho*np.sin(ph);rr=np.hypot(px,py);an=np.arctan2(py,px)
    ws=np.array([brentq(lambda w:V.Rfun(tr,w,M)-r,W[0]+1e-6,40) for r in rr]);rings.append((ws,an))
times=[tr,tx,0.687,t0];cols=['0.4','C0','C2','C3']
fig,ax=plt.subplots(figsize=(15,8))
th=np.linspace(-np.pi/2,np.pi/2,200);ax.fill(np.r_[0,R*np.cos(th),0],np.r_[0,R*np.sin(th),0],fill=False,hatch='///',color='orange',lw=0)
ax.plot(R*np.cos(th),R*np.sin(th),'k-',lw=2);ax.text(0.05,0.05,'PLASMA',fontsize=12);ax.text(R+0.02,-0.55,'recombination\nsurface (0.5 D)',fontsize=9)
for a in (0,0.25,0.5):
    for s in (1,-1): ax.plot([0,2.6*np.cos(s*a)],[0,2.6*np.sin(s*a)],'k--',lw=.6,alpha=.5)
caps=[]
for t,c in zip(times,cols):
    for ws,an in rings:
        r=np.array([V.Rfun(t,w,M) for w in ws]);x,y=r*np.cos(an),r*np.sin(an)
        ax.plot(x,y,color=c,lw=1.6)
        ins=r<R
        if ins.any(): ax.plot(x[ins],y[ins],color=c,lw=4,alpha=.35)
    xo=V.Rfun(t,wo,M);ax.plot(xo,0,'o',color=c,ms=6);ax.text(xo,0.2+0.02*times.index(t),f't={t:.3f}',color=c,ha='center')
    if t>tx:
        o,P=v11.trace(F,R,np.arange(0,360,15.0),t,keep=True);bx,by,sx,sy=[],[],[],[]
        for (ps,tt,ww,pp),oo in zip(P,o):
            rb=V.Rfun(t,oo[4],M);bx.append(rb*np.cos(oo[5]));by.append(rb*np.sin(oo[5]))
            if tt[-1]<tx: j=np.argmin(abs(tt-tx))
            else: j=-1
            rs=V.Rfun(t,ww[j],M);sx.append(rs*np.cos(pp[j]));sy.append(rs*np.sin(pp[j]))
        ax.add_patch(plt.Circle((xo,0),t-tx,fill=False,ls='--',color=c,lw=1.3));ax.plot(np.r_[bx,bx[:1]],np.r_[by,by[:1]],':',color=c,lw=2)
        d=np.hypot(np.array(bx)-xo,by);caps.append((round(t,3),round(t-tx,3),{int(a):round(float(v),3) for a,v in zip(o[:,0],d) if int(a) in (0,45,90,135,180)}))
ax.plot(0,0,'k+',ms=16,mew=2);ax.text(0.02,-0.08,'CoE',fontsize=10)
ax.annotate('',xy=(2.3,-0.9),xytext=(1.6,-0.9),arrowprops=dict(arrowstyle='->',lw=1.5));ax.text(1.6,-0.85,'outward from CoE',fontsize=9)
from matplotlib.lines import Line2D
h=[Line2D([],[],color=c,lw=1.6,label=f'patch of space around us at t={t:.3f}') for t,c in zip(times,cols)]+[Line2D([],[],color='k',lw=4,alpha=.35,label='part of patch still plasma'),Line2D([],[],color='k',ls='--',label='small bubble: how far we see = circle c(t - t_x)'),Line2D([],[],color='k',ls=':',lw=2,label='big bubble: CMB edge (where seen matter is now)')]
ax.legend(handles=h,loc='upper left',fontsize=9);ax.set_aspect('equal');ax.set_xlim(-0.1,2.7);ax.set_ylim(-1.2,1.2)
ax.set_title(f'Our patch of space stretching outward from the CoE (v11, with core mass, R_rec=0.5 D, we cross at t={tx:.3f}; distances in D at each time)')
fig.savefig('steps_handdrawn_style_v2_v3.png',dpi=110,bbox_inches='tight');print(tx,ro,caps)
