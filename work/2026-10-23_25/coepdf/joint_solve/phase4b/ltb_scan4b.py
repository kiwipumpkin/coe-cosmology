# Phase4b (NEW) copy: output to phase4b; DSHAPE env sets the D used in the density shape (default 34.96 as phase4).
# Phase4 copy: frac scan near age anchor (Ht0~0.968).
# Phase3 A: famP poly LTB background vs kinematic speed u (=H0 D/c at the observer shell) and binding frac (Mtot = frac*u^2/2).
# Reports, in units D=1,c=1: H_loc, t0, H t0, Omega_local, eta_rec (conformal distance to the 1+z=1100 surface on the observer's shell, radial).
import sys,os,json,time,numpy as np
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');os.chdir(B+'dipole_amp')
import dist2 as D
from scipy.integrate import cumulative_trapezoid,quad
from scipy.optimize import brentq
DR=float(os.environ.get('DSHAPE','34.96'))/14.30048287260589;Rb=1.5*(1+1/DR);CC=[-0.041,-0.035,0.1353,0.0]
shape=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=shape(np.maximum(X,1e-12))*(X<=Rb)
m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
res=[]
for u,frac in [tuple(map(float,a.split(','))) for a in sys.argv[1:]]:
    D.u=u;t=time.time()
    F=D.build(frac*0.5*u**2,mt);S,t0=F['S'],F['t0']
    H=S(t0,1.,dx=1)[0,0]/S(t0,1.)[0,0]
    Mp=F['Mtot']*np.gradient(np.interp(np.array([0.999,1.0,1.001]),*mt),0.001)[1];Om=2*Mp/(S(t0,1.)[0,0]**2*S(t0,1.,dy=1)[0,0])/3/H**2
    a0=S(t0,1.)[0,0];tr=brentq(lambda tt:S(tt,1.)[0,0]-a0/1100.,F['T'][0]*1.001,t0)
    eta=quad(lambda tt:a0/S(tt,1.)[0,0],tr,t0,limit=400)[0]   # comoving (today) radial light distance, homogeneous-shell approx
    r=dict(u=u,frac=frac,H_loc=H,t0=t0,Ht0=H*t0,Omega=Om,eta_rec_D=eta,sec=time.time()-t);res.append(r);print(r,flush=True)
json.dump(res,open('/workspace/coepdf/joint_solve/phase4b/ltb_scan_%s.json'%os.environ.get('TAG','x'),'w'),indent=1)
