exec(open("expand.py").read().split("C=[];to=2.0")[0])
from scipy.integrate import solve_ivp
res={}
def run(mu,uexp,v0,to):
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

import sys
for v0 in (0.3,0.9,0.99,0.999):
  to=3.0/(1-v0)
  for mu in (0.0,0.1):
    for u in (0.0,0.9999):
      try:
        r=run(mu,u,v0,to); res[(v0,mu,u)]=r; print(v0,to,mu,u,"dip %.4f quad %.4f max %.4f min %.4f"%stats(r),flush=True)
      except Exception as ex: print(v0,mu,u,"FAIL",ex,flush=True)
import pickle; pickle.dump(res,open("startc.pkl","wb"))
