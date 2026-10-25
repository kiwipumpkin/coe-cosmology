# SR toy: emitting sphere, observer inside. Units: c=1, R(t_obs)=1. Observer at sphere centre at t_obs=0.
# Sight direction n' in observer frame (angle th from +x = direction of motion). Photon propagates along -n.
import numpy as np,json
from scipy.optimize import brentq
def addv(u,v):  # v measured in frame moving at u (3-vectors) -> lab
    u=np.asarray(u,float);v=np.asarray(v,float);uu=u@u
    if uu==0:return v
    g=1/np.sqrt(1-uu);vpar=(v@u)/uu*u;vperp=v-vpar
    return (u+vpar+vperp/g)/(1+u@v)
def ring(vc,xc,R,Rdot,vo,N=360):
    """vc(t),xc(t): centre velocity/position (lab); R(t),Rdot(t) in centre frame (Galilean radius in lab, toy);
    vo: observer lab velocity at t=0 (observer at origin). Returns T_obs/T_emit vs observer-frame sight angle."""
    th=np.radians(np.arange(N)*360/N);b=np.array([vo,0,0]);go=1/np.sqrt(1-vo*vo);out=[]
    for t in th:
        kp=-np.array([np.cos(t),np.sin(t),0])                                  # photon dir, observer frame
        kx=(kp[0]+vo)/(1+vo*kp[0]);ky=kp[1]/(go*(1+vo*kp[0]));k=np.array([kx,ky,0]);k/=np.linalg.norm(k)   # lab
        n=-k;f=lambda s:np.linalg.norm(n*s-xc(-s))-R(-s)
        s=brentq(f,1e-9,50);te=-s;xe=n*s;m=(xe-xc(te));m/=np.linalg.norm(m)
        ve=addv(vc(te),Rdot(te)*m);ge=1/np.sqrt(1-ve@ve)
        out.append(go*(1-b@k)/(ge*(1-ve@k)))   # nu_o/nu_e = (k.u_o)/(k.u_e), lab frame
    return np.array(out)

def stats(x):
    d=x/x.mean()-1;m=np.fft.rfft(d)/len(d)
    return dict(rms=float(np.sqrt(np.mean(d**2))),range=float(d.max()-d.min()),dipole=float(2*abs(m[1])),quad=float(2*abs(m[2])),
                dipole_sign=float(np.sign(m[1].real)),monopole_T_over_Temit=float(x.mean())),d
VEXP=0.5   # expanding: R = 1 + VEXP*t (R=VEXP*(t+2) i.e. sphere born at t=-2; R(0)=1)
def build(mode,v,a,expanding):
    Rf=(lambda t:1+VEXP*t) if expanding else (lambda t:1.0+0*t);Rd=(lambda t:VEXP) if expanding else (lambda t:0.0)
    z=lambda t:np.zeros(3)
    if mode==1:return dict(vc=z,xc=z,R=Rf,Rdot=Rd,vo=0.0)
    if mode==2:return dict(vc=z,xc=z,R=Rf,Rdot=Rd,vo=v)
    if mode==3:return dict(vc=lambda t:np.array([v,0,0]),xc=lambda t:np.array([v*t,0,0]),R=Rf,Rdot=Rd,vo=v)
    sg=1 if mode==4 else -1   # 4: acceleration along motion; 5: opposite
    return dict(vc=lambda t:np.array([v+sg*a*t,0,0]),xc=lambda t:np.array([v*t+0.5*sg*a*t*t,0,0]),R=Rf,Rdot=Rd,vo=v)
