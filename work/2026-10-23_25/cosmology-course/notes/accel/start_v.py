exec(open("expand.py").read().split("C=[];to=2.0")[0])
from scipy.integrate import solve_ivp
res={};to=3.0
def run(mu,uexp,v0):
    def tr(p0,u0):
        def f(t,y):
            p=y[2:]; v=p/np.sqrt(1+p@p); F=-mu*(y[:2]-Mp)/np.linalg.norm(y[:2]-Mp)**3
            return [*v,*F]
        g=1/np.sqrt(1-np.dot(u0,u0)); u0=[g*u0[0],g*u0[1]]
        return solve_ivp(f,(0,to),[*p0,*u0],dense_output=True,rtol=1e-11,atol=1e-13).sol
    yo=tr([0,0],[v0,0])(to); yo=np.r_[yo[:2],yo[2:]/np.sqrt(1+yo[2:]@yo[2:])]; cache={}
    def emit(s,e):
        key=(round(e[0],12),round(e[1],12))
        if key not in cache:
            G0=1/np.sqrt(1-v0*v0); up=uexp*e; den=1+v0*up[0]
            cache[key]=tr(list(e),[(up[0]+v0)/den, up[1]/(G0*den)])
        y=cache[key](s); vv=y[2:]/np.sqrt(1+y[2:]@y[2:])
        return y[:2],vv
    return T_of(yo[:2],yo[2:],emit,to,mu,0.0)

for v0 in (0.1,0.5):
  for mu in (0.0,0.1):
    for u in (0.0,0.9999):
      r=run(mu,u,v0); res[(v0,mu,u)]=r; print(v0,mu,u,"dip %.4f quad %.4f max %.4f min %.4f"%stats(r),flush=True)
  fig,ax=plt.subplots(figsize=(13,7.5)); x=np.degrees(th)
  for u,st,l in ((0.0,"-","not expanding"),(0.9999,"--","expanding at 0.9999c (≈c)")):
    r=res[(v0,0.1,u)]; d,q,_,_=stats(r); ax.plot(x,r,st,lw=3,label=f"{l}: dipole {d:.4f}, quadrupole {q:.4f}")
  ax.axhline(0,color="gray"); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
  ax.set_xlabel(f"direction (deg); GM/(Rc²) = 0.1, start {v0}c, slowed by mass at 180°",fontsize=15)
  ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.15),fontsize=14); plt.tight_layout(); plt.savefig(f"decel_v{v0}.png",dpi=120); plt.close()
