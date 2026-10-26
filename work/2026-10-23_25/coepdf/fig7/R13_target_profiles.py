# rev7 remake of rev6 Fig 13 (target angular profiles matched by masses). Original script not found; rebuilt from the caption's method.
# Dimensionless in R (lengths / R_ls), so it does not depend on D. Potential on the sphere |n|=1:
#   Phi(theta) = M0/|n - 1.55 z| + sum_j m_j / |n - p_j|, M0 fixed (=1), p_j = 40 off-sphere sites (radii 0.5-0.9 and 1.15-1.4, polar angles 0-180 deg
#   in the x-z plane, axisymmetrised by averaging over azimuth), m_j signed (clump>0, void<0) fitted by ridge least squares so that
#   Phi(theta)-mean matches a*T(theta) for T = cos(theta/2), cos(1.5 theta), cos(2.5 theta), cos(2 theta).
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
th=np.radians(np.linspace(0,180,181));ph=np.linspace(0,2*np.pi,72,endpoint=False)
def ring(r,a):  # azimuth-averaged 1/|n-p| for p at radius r, polar angle a (an axisymmetric ring of mass)
    n=np.stack([np.sin(th)[:,None]*np.cos(ph),np.sin(th)[:,None]*np.sin(ph),np.cos(th)[:,None]+0*ph],-1);p=np.array([r*np.sin(a),0,r*np.cos(a)])
    p2=np.stack([r*np.sin(a)*np.cos(ph),r*np.sin(a)*np.sin(ph),r*np.cos(a)+0*ph],-1)
    return np.mean(1/np.linalg.norm(n[:,None,:,:]-p2[None,:,None,:],axis=-1),axis=(1,2))
sites=[(r,a) for r in (0.5,0.7,0.9,1.15,1.4) for a in np.radians(np.linspace(0,180,8))]
G=np.stack([ring(r,a) for r,a in sites],1);G-=G.mean(0);P0=ring(1.55,0.0);P0-=P0.mean()
T={'half wave cos(θ/2)':np.cos(th/2),'1½ waves cos(1.5θ)':np.cos(1.5*th),'2½ waves cos(2.5θ)':np.cos(2.5*th),'two waves cos(2θ)':np.cos(2*th)}
fig,ax=plt.subplots(2,4,figsize=(22,9));out={}
for i,(k,t) in enumerate(T.items()):
    t=t-t.mean();a=np.ptp(P0)/np.ptp(t);rhs=a*t-P0;lam=1e-6*np.trace(G.T@G)
    m=np.linalg.solve(G.T@G+lam*np.eye(G.shape[1]),G.T@rhs);fit=P0+G@m;r2=1-np.sum((fit-a*t)**2)/np.sum((a*t)**2)
    out[k]=dict(r2=float(r2),n_clump=int((m>0).sum()),n_void=int((m<0).sum()),sum_m=float(m.sum()))
    A=ax[0,i];A.plot(np.degrees(th),a*t,'k',lw=2.5,label='target');A.plot(np.degrees(th),fit,'r--',lw=2,label=f'fit, R²={r2:.4f}');A.plot(np.degrees(th),P0,'0.6',label='fixed mass alone')
    A.set_title(k);A.set_xlabel('θ from mass axis [deg]');A.legend(fontsize=8)
    B=ax[1,i];xs=np.array([r*np.sin(a_) for r,a_ in sites]);zs=np.array([r*np.cos(a_) for r,a_ in sites])
    B.scatter(xs,zs,s=20+300*np.abs(m)/np.abs(m).max(),c=np.where(m>0,'crimson','royalblue'),alpha=.7);B.plot(0,1.55,'k*',ms=14)
    B.add_patch(plt.Circle((0,0),1,fill=False,ls='--'));B.set_aspect('equal');B.set_xlim(-.2,1.6);B.set_ylim(-1.6,1.7);B.set_xlabel('ρ / R');B.set_ylabel('z / R')
    B.set_title(f'{out[k]["n_clump"]} clump rings (red), {out[k]["n_void"]} void rings (blue); ★ fixed mass at 1.55R',fontsize=8)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=90);json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1);print(out)
