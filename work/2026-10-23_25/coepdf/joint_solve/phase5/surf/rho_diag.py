# Diagnostic: density and expansion state of each shell at t_rec (time our shell had a=1/1100), Migkas K-1 vs poly validation. Same build as famP_mig.py.
import sys,os,numpy as np
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');os.chdir(B+'dipole_amp')
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
RV=14.30048287260589
for name,Dg,frac,mode,LDX,ETA in [('K-1',27.737393,0.01847547274832157,'mig',0.42296181187611975,3.068065),('VALpoly',27.3385,0.0185,'poly',0,0)]:
    os.environ['U']=str(68.64*Dg*1000/299792.458/1.00016);import importlib,dist2 as D;importlib.reload(D)
    DR=Dg/RV;Rb=1.5*(1+1/DR);CC=[-0.041,-0.035,0.1353,0.0]
    shape=(lambda x:np.exp(-(np.abs(x)/LDX)**ETA)) if mode=='mig' else (lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC))))
    X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=shape(np.maximum(X,1e-12))*(X<=Rb)
    m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);F=D.build(frac*0.5*D.u**2,(X,m/np.interp(1.0,X,m)));S,t0=F['S'],F['t0']
    tus=brentq(lambda t:S(t,1.0)[0,0]-1/1100.,F['T'][0]*1.001,t0)
    rho=lambda t,x:shape(x)*x*x/(S(t,x)[0,0]**2*S(t,x,dy=1)[0,0]);rc=rho(tus,1.)
    print(name,'t_rec',tus,'t0',t0,'Rb',Rb)
    print(' x     shape(x)/shape(1)  rho(t_rec,x)/rho(t_rec,1)  a(t_rec,x)/a(t0,x)*1100  H_r/H_t(t0)')
    for x in [0.2,0.3,0.4,0.5,0.7,1.0,1.3,1.6,2.0,2.5]:
        R=S(tus,x)[0,0];a=R/S(t0,x)[0,0]*1100;Hr=S(t0,x,dx=1,dy=1)[0,0]/S(t0,x,dy=1)[0,0];Ht=S(t0,x,dx=1)[0,0]/S(t0,x)[0,0]
        print(f' {x:4.1f} {shape(x)/shape(1.):11.3e} {rho(tus,x)/rc:11.3e} {a:9.4f} {Hr/Ht:8.4f}')
