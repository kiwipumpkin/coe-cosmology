# exact no-mass (Minkowski/Milne) version of v8 for any redshift; observer at areal R=1 with w=uD
import numpy as np,sys
from scipy.optimize import brentq
uD=10**0.021927885054513485; t0=1/uD; g=np.sqrt(1+uD**2); to=g*t0; xo=np.array([uD*t0,0.])
def run(fr,q,psis):
    ts=fr*t0; out=[]
    for ps in psis:
        # direction in observer rest frame -> lab via aberration (observer moves along +x with v)
        v=uD/g; c=np.cos(np.radians(ps)); s=np.sin(np.radians(ps))
        cl=(c-v)/(1-v*c); sl=np.sqrt(max(0,1-cl*cl))*np.sign(s if s!=0 else 1)   # source direction in lab
        n=np.array([cl,sl])
        f=lambda L:(lambda t,x:np.log(np.sqrt(max(t*t-x@x,1e-300)))-np.log(ts)+q*np.log(np.hypot(*x)+1e-300))(to-L,xo+L*n)
        Ls=np.linspace(0,to,20001);fv=np.array([f(L) for L in Ls]);i=np.where(fv[:-1]*fv[1:]<0)[0]
        if not len(i): out.append((ps,)+(np.nan,)*4);continue
        L=brentq(f,Ls[i[0]],Ls[i[0]+1]);t=to-L;x=xo+L*n;tau=np.sqrt(t*t-x@x)
        k=np.r_[1.,-n];ue=np.r_[t,x]/tau;uo=np.r_[to,xo]/t0;mink=lambda a,b:a[0]*b[0]-a[1:]@b[1:]
        out.append((ps,mink(k,ue)/mink(k,uo),tau,t0-tau,np.hypot(*x)))
    return np.array(out)
if __name__=='__main__':
    fr,q=float(sys.argv[1]),float(sys.argv[2]); r=run(fr,q,np.arange(0,181,1.0)); np.save(f'milne_f{fr}_q{q}.npy',r)
    print(r[::30].round(5))
