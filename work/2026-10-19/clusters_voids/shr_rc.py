# Recombination when each shell's own density drops to a critical value. Density from exact LTB: rho = M'(x)/(4 pi R^2 R').
import sys,os,json,numpy as np
sys.path.insert(0,'../dipole_amp');os.environ.setdefault('TMIN','1e-8')
import dist2 as D
from scipy.optimize import brentq
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
mass,ztarget,k,out=sys.argv[1],sys.argv[2],float(sys.argv[3]),sys.argv[4]
LD=18159.848409348455/39.10676196793444;al=0.30000009809958916;DR=42/14.30048287260589;Rb=1.5*(1+1/D.DR) if hasattr(D,'DR') else 1.5*(1+1/DR)
mt=D.mfrac_table(LD,al,Rb,'soft');frac=4.362106762731466 if mass=='b' else 0.0
F=D.build(frac*0.5*D.u**2,mt);S,t0=F['S'],F['t0']
shape=lambda x:(1+(x/LD)**al)**-1.0     # plasma density today, CatWISE-locked soft law (relative units)
def rho(t,x):   # comoving mass conservation: rho(t,x) R^2 R' = rho_today(x) x^2  (exact for non-crossing shells)
    R=S(t,x)[0,0];Rp=S(t,x,dy=1)[0,0];return shape(x)*x**-k*x*x/(R*R*Rp)
# critical density: our shell (x=1) recombines when its radius is 0.5 D (matches old setup) or 1/1100 D (1+z ~ 1100)
Rus=0.5 if ztarget=='old' else 1/1100.
tus=brentq(lambda t:S(t,1.0)[0,0]-Rus,F['T'][0]*1.001,t0);rc=rho(tus,1.0)*float(os.environ.get('RCF','1'))
XX=np.linspace(0.6,1.6,201);te=np.array([brentq(lambda t:np.log(rho(t,x)/rc),F['T'][0]*1.0001,t0,xtol=1e-16,rtol=1e-14) for x in XX])
Rr=np.array([S(t,x)[0,0] for t,x in zip(te,XX)]);tf=interp1d(XX,te,fill_value='extrapolate')
def trace(psis):
    def AB(t,w):
        Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
        return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
    def rhs(l,y):
        t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B,Bt,Bw=AB(t,w)
        return [-pt,pw/A,pp/B,0.5*At/A**2*pw**2+0.5*Bt/B**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B**2*pp**2,0]
    hit=lambda l,y:y[0]-float(tf(y[1]));hit.terminal=True
    A0=AB(t0,1.0)[1];B0=AB(t0,1.0)[4];o=[]
    for ps in psis:
        pr=np.radians(ps);y0=[t0,1.0,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
        s=solve_ivp(rhs,(0,-200),y0,events=[hit],rtol=1e-10,atol=1e-13)
        y=s.y_events[0][0];o.append((ps,-y[3],y[0],t0-y[0],y[1],y[2],S(y[0],y[1])[0,0]))
    return np.array(o)
o=trace(np.arange(0,180.001,5))
np.savez(out,o=o,XX=XX,te=te,Rr=Rr,rc=rc,tus=tus,t0=t0,k=k,mass=mass)
print(out,'t0',t0,'te range',o[0,2],o[-1,2],'1+z',o[0,1],o[-1,1],flush=True)
