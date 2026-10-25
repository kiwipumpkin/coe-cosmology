import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp; from scipy.optimize import brentq
th=np.radians(np.arange(0,361,1.0)); Mp=np.array([-2.0,0.0]); k=0.5
def T_of(xo,vo,emit,to,mu,lo):
    out=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)])
        te=brentq(lambda s: np.linalg.norm(xo-emit(s,e)[0])-(to-s),lo,to-1e-6)
        xe,ve=emit(te,e); n=(xo-xe)/np.linalg.norm(xo-xe)
        ge,go=1/np.sqrt(1-ve@ve),1/np.sqrt(1-vo@vo)
        g=np.sqrt((1-2*mu/np.linalg.norm(xe-Mp))/(1-2*mu/np.linalg.norm(xo-Mp))) if mu else 1.0
        out.append(go*(1-vo@n)/(ge*(1-ve@n))*g)
    d=np.array(out); return d/d.mean()-1
def stats(r): m=np.fft.rfft(r[:-1])/360; return 2*abs(m[1]),2*abs(m[2]),r.max(),r.min()
def plot(C,xl,fn):
    plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(figsize=(13,8))
    for lab,r in C: ax.plot(np.degrees(th),r,lw=3,label=lab)
    ax.axhline(0,color="gray"); ax.set_xlabel(xl); ax.set_ylabel("ΔT/T"); ax.set_xlim(0,360); ax.grid(alpha=.3)
    ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.13),ncol=4,fontsize=15); plt.tight_layout(); plt.savefig(fn,dpi=120); plt.close()
X=lambda t,a:(np.sqrt(1+(a*t)**2)-1)/a
V=lambda t,a:a*t/np.sqrt(1+(a*t)**2)
C=[];to=2.0
for a in (0.05,0.1,0.2,0.5):
    def emit(s,e,a=a):
        vb=V(s,a); G0=1/np.sqrt(1-vb*vb); up=k*vb*e; den=1+vb*up[0]
        return np.array([X(s,a),0])+(1+k*X(s,a))*e, np.array([(up[0]+vb)/den, up[1]/(G0*den)])
    r=T_of(np.array([X(to,a),0]),np.array([V(to,a),0]),emit,to,0,0.0); C.append((f"a = {a} c²/L",r)); print("acc",a,"dip %.4f quad %.4f max %.4f min %.4f"%stats(r),flush=True)
plot(C,"direction (deg); accelerating toward 0°, shell expanding at half the speed","expand_accel.png")
v=0.2;C=[]
for mu in (0.05,0.1,0.2,0.3):
    G=1/np.sqrt(1-v*v); u=k*v
    def emit(s,e):
        # exact Lorentz boost of the rest-frame expanding shell x'=(1+u tau)e at time tau
        ex=e[0]; tau=(s/G - v*ex)/(1+v*u*ex)   # solve t = G(tau + v x'_x) for tau
        xp=(1+u*tau)*e; x=np.array([G*(xp[0]+v*tau),xp[1]])
        up=u*e; den=1+v*up[0]; w=np.array([(up[0]+v)/den, up[1]/(G*den)])
        return x,w
    r=T_of(np.zeros(2),np.array([v,0]),emit,0.0,mu,-5.0); C.append((f"GM/(Rc²) = {mu}",r)); print("const",mu,"dip %.4f quad %.4f max %.4f min %.4f"%stats(r),flush=True)
plot(C,"direction (deg); 0.2c toward 0° (with thrust), expanding at 0.1c, mass at 180°","expand_const.png")
v0=0.3;to=3.0;C=[]
for mu in (0.0,0.02,0.05,0.1):
    def tr(p0,u0,mu=mu):
        f=lambda t,y:[y[2],y[3],*(-mu*(y[:2]-Mp)/np.linalg.norm(y[:2]-Mp)**3)]
        return solve_ivp(f,(0,to),[*p0,*u0],dense_output=True,rtol=1e-10,atol=1e-12).sol
    yo=tr([0,0],[v0,0])(to); cache={}
    def emit(s,e):
        key=(round(e[0],12),round(e[1],12))
        if key not in cache:
            G0=1/np.sqrt(1-v0*v0); up=k*v0*e; den=1+v0*up[0]
            cache[key]=tr(list(e),[(up[0]+v0)/den, up[1]/(G0*den)])
        y=cache[key](s); return y[:2],y[2:]
    r=T_of(yo[:2],yo[2:],emit,to,mu,0.0); C.append((f"GM/(Rc²) = {mu}",r)); print("decel",mu,"dip %.4f quad %.4f max %.4f min %.4f"%stats(r),flush=True)
plot(C,"direction (deg); start 0.3c toward 0°, expanding at 0.15c, slowed by mass at 180°","expand_decel.png")
