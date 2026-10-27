# Phase5 routes (NEW Oct 26): cheap per-shell LTB evaluation of the density at recombination, same physics as dist2.build
# (common bang t=0, growing mode only, label x = today's radius in units of D, M(<x) from rho_t0(x)=shape(x) truncated at Rb).
# Returns rho(t_rec,x)/rho(t_rec,1) where t_rec = time our shell had a=1/1100. Mtot=frac*u^2/2, frac solved so H t0 = age 13.787 Gyr.
import sys,os,numpy as np
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B);import base as V
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
RV=14.30048287260589;H0=68.64;AGE=13.787;HT=AGE/(977.8/H0)
def tage(R,M,E):  # closed-form time to reach R from bang (outgoing branch)
    if M<=0:return R/np.sqrt(2*E)
    if abs(E)<1e-14*M/R:return np.sqrt(2*R**3/(9*M))
    a=M/(2*abs(E));b=a/np.sqrt(2*abs(E))
    if E>0:
        e=np.arccosh(1+R/a);return b*(np.sinh(e)-e)
    c=1-R/a
    if c<-1:return np.inf
    e=np.arccos(c);return b*(e-np.sin(e))
def our(u):  # frac for age: our shell has M=Mtot, E=u^2/2-Mtot, H=u at R=1
    f=lambda fr:tage(1.,fr*u*u/2,u*u/2-fr*u*u/2)*u-HT
    return brentq(f,1e-8,0.999)
def drec(D,LD,eta,xs,shape=None):
    u=H0*D*1000/299792.458;Rb=1.5*(1+RV/D)
    if shape is None:shape=lambda x:np.exp(-(np.abs(x)/LD)**eta)
    X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=shape(np.maximum(X,1e-12))*(X<=Rb)
    with np.errstate(all='ignore'):m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0)
    from scipy.interpolate import CubicSpline;cs=CubicSpline(np.log(X[1:]),np.log(np.maximum(m[1:],1e-300)));n1=np.exp(cs(0.0));  # FIX: cubic log-log M(x) (linear interp gave ~3e-3 spurious delta)
    fr=our(u);Mt=fr*u*u/2;t0=tage(1.,Mt,u*u/2-Mt)
    Mx=lambda x:Mt*np.exp(cs(np.log(x)))/n1
    def R(x,t):
        M=Mx(x)
        if M<=0:return x*t/t0
        lo=-M/x*(1-1e-12);hi=50*(x/t0)**2+10*M/x+10
        E=brentq(lambda e:tage(x,M,e)-t0,lo,hi,xtol=1e-15,rtol=1e-13)
        return V.Rfun(t,np.sign(E)*np.sqrt(2*abs(E)),M)
    tr=brentq(lambda t:R(1.,t)-1/1100.,1e-12,t0)
    def rho(x):
        h=1e-4*x;Rr=R(x,tr);Rp=(R(x+h,tr)-R(x-h,tr))/(2*h)
        return shape(x)*(x<=Rb)*x*x/(Rr*Rr*Rp)
    r1=rho(1.);return np.array([rho(x)/r1 for x in xs]),fr,tr,t0
if __name__=='__main__':
    import time;t=time.time()
    xs=[0.2,0.5,0.7,1.3,1.6,2.0,2.5]
    print(drec(27.737393,0.42296181187611975,3.068065,xs)[0],time.time()-t)
