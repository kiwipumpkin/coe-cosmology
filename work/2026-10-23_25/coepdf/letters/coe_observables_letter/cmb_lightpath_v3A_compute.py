# v3A: CMB dT/T from the rev7 LTB famP background, using its own null geodesics and per-shell recombination exactly as built in
#   /workspace/coepdf/run_out/run9/realweb/los_flow/ltb.py (copied verbatim below up to the geodesic; nothing saved over its npz).
#   Recombination: shell x recombines when rho(t,x) = rho(t_us,1), t_us defined by S(t,1)=1/1100 (as built; T∝1/a at our shell).
#   Redshift: geodesic y=[t,w,phi,p_t,p_w,p_phi], dt/dlambda=-p_t; photon energy seen by comoving dust ∝ -p_t, so 1+z = p_t(src)/p_t(obs).
#   Observer comoving with the dust at x=1 (LTB radial coordinate in units of D); T_rec fixed => T_obs(n) ∝ 1/(1+z(n)).
#   This contains kinematic+gravitational+ISW/RS in one GR redshift (no split). Writes cmb_lightpath_v3A.npz (profile vs psi).
import sys,os,numpy as np,json
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');cwd=os.getcwd();os.chdir(B+'dipole_amp')
import dist2 as D
os.chdir(cwd)
from scipy.integrate import solve_ivp,cumulative_trapezoid
from scipy.optimize import brentq
from scipy.interpolate import interp1d
DR=42/14.30048287260589;Rb=1.5*(1+1/DR);CC=[-0.041,-0.035,0.1353,0.0];frac=20000.
shape=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=shape(np.maximum(X,1e-12))*(X<=Rb)
m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
F=D.build(frac*0.5*D.u**2,mt);S,t0=F['S'],F['t0']
def rho(t,x):R=S(t,x)[0,0];Rp=S(t,x,dy=1)[0,0];return shape(x)*x*x/(R*R*Rp)
tus=brentq(lambda t:S(t,1.0)[0,0]-1/1100.,F['T'][0]*1.001,t0);rc=rho(tus,1.0)
XX=np.linspace(0.98,1.02,81);te=np.array([brentq(lambda t:np.log(rho(t,x)/rc),F['T'][0]*1.0001,t0,xtol=1e-16,rtol=1e-14) for x in XX]);tf=interp1d(XX,te,fill_value='extrapolate')
def AB(t,w):
    Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
    return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B_,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B_,0.5*At/A**2*pw**2+0.5*Bt/B_**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B_**2*pp**2,0]
hit=lambda l,y:y[0]-float(tf(y[1]));hit.terminal=True
A0=AB(t0,1.0)[1];B0=AB(t0,1.0)[4];PS=np.arange(0,180.001,1.0);out=[]
for ps in PS:
    pr=np.radians(ps);y0=[t0,1.0,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
    s=solve_ivp(rhs,(0,-200),y0,events=[hit],rtol=1e-10,atol=1e-13)
    ye=s.y_events[0][0] if len(s.t_events[0]) else np.full(6,np.nan)
    w,ph=ye[1],ye[2];out.append([ps,ye[0],w,ph,ye[3],float(np.sqrt(1+w*w-2*w*np.cos(ph)))])
    print(ps,out[-1],flush=True)
o=np.array(out);np.savez('cmb_lightpath_v3A.npz',psi_deg=o[:,0],t_src=o[:,1],x_src=o[:,2],phi_src=o[:,3],pt_src=o[:,4],onez=o[:,4]/(-1.0)*-1,t0=t0,tus=tus,
  note='psi measured at observer from the +w direction (away from CoE = anti-CoE) as in ltb.py y0; 1+z = p_t(src)/p_t(obs), p_t(obs)=-1')
print('t0',t0,'tus',tus)
