# Phase5 (NEW): CMB background built from the SAME density profile the Migkas fit uses.
# rho(r) ~ exp(-(r/L)^eta) (L from CatWISE lock delta=0.031k at depth d), truncated at Rb=1.5(1+RV/D) in units of D (prof_core rb_of).
# Normalisation: LTB/famP mass Mtot = frac*u^2/2 (D=c=1), u = <H0> D/c with <H0>=68.64. In these units prof_core's mc (= 2 G M(<D)/(D^3 H^2))
# maps to frac exactly (GM/(D c^2) = mc u^2/2), so frac=mc is the 'same normalisation'. Age 13.787 Gyr is held FIXED by solving frac so that
# H t0 = 13.787/(977.8/68.64) in the LTB build; frac/mc is reported as the consistency check between the two dynamics codes.
# Then: R_ls = eta_rec*D (radial conformal distance to 1+z=1100 on our shell, as ltb5.py), famP_mig trace -> ns64 bg map (as p2_bg5.py).
# Usage: python3 bg_mig.py TAG D mc eta L_Gpc [poly]   (poly = validation mode: old polynomial shape, frac fixed 0.0185)
import sys,os,json,time,subprocess,numpy as np,healpy as hp
O='/workspace/coepdf/joint_solve/phase5/';B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/'
sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');cwd=os.getcwd();os.chdir(B+'dipole_amp');import dist2 as DD;os.chdir(cwd)
from scipy.integrate import cumulative_trapezoid,quad
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
TAG=sys.argv[1];Dg,mc,eta,L=map(float,sys.argv[2:6]);POLY=len(sys.argv)>6 and sys.argv[6]=='poly'
RV=14.30048287260589;H0=68.64;AGE=13.787;HT=AGE/(977.8/H0);DR=Dg/RV;Rb=1.5*(1+1/DR);LDX=L/Dg
CC=[-0.041,-0.035,0.1353,0.0]
if POLY:shape=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
else:shape=lambda x:np.exp(-(np.abs(x)/LDX)**eta)
X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=shape(np.maximum(X,1e-12))*(X<=Rb)
m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
u=H0*Dg*1000/299792.458;DD.u=u;ck=O+f'bg_{TAG}_ckpt.json';C=json.load(open(ck)) if os.path.exists(ck) else {}
def ltb(frac):
    key='%.10g'%frac
    if key in C:return C[key]
    t=time.time();F=DD.build(frac*0.5*u**2,mt);S,t0=F['S'],F['t0'];H=S(t0,1.,dx=1)[0,0]/S(t0,1.)[0,0]
    a0=S(t0,1.)[0,0];tr=brentq(lambda tt:S(tt,1.)[0,0]-a0/1100.,F['T'][0]*1.001,t0);er=quad(lambda tt:a0/S(tt,1.)[0,0],tr,t0,limit=400)[0]
    r=dict(frac=frac,H_loc=H,t0=t0,Ht0=H*t0,eta_rec_D=er,shell_crossing=bool(F['shell_crossing']),n_dropped_inner=int(F['n_dropped_inner']),sec=time.time()-t)
    C[key]=r;json.dump(C,open(ck,'w'),indent=1);print(r,flush=True);return r
if POLY:fr=0.0185
else:
    g=lambda lf:ltb(10**lf)['Ht0']-HT;a=np.log10(mc);fa=g(a);step=-0.3 if fa<0 else 0.3;b=a+step;fb=g(b);n=0
    while fa*fb>0 and n<12:a,fa=b,fb;b=b+step;fb=g(b);n+=1   # Ht0 decreases with frac
    fr=10**brentq(g,min(a,b),max(a,b),xtol=1e-4) if fa*fb<0 else None
if fr is None:json.dump(dict(tag=TAG,fail='no frac gives age 13.787'),open(O+f'bg_{TAG}.json','w'));sys.exit('no age root')
r=ltb(fr);Rls=r['eta_rec_D']*Dg
# famP_mig trace + map (as p2_bg5.py)
# SURF env: rho (famP default, density-set surface) or a (1+z_local=1100)
tr=O+f'bg_{TAG}_trace.npz';Cd='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained'
if not os.path.exists(tr):
    env=dict(os.environ,SURF=os.environ.get('SURF','rho'),U=str(u/1.00016),DRV=str(DR),C1='-0.041',C2='-0.035',C3='0.1353',C4='0',LDX=str(LDX),ETA=str(eta))
    subprocess.run(['python3',O+'famP_mig.py','poly' if POLY else 'mig','0',str(fr),tr],env=env,check=True,cwd=Cd)
o=np.load(tr)['o'];o=o[np.argsort(o[:,0])]
def lb(l,b):l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
V=np.array(hp.pix2vec(64,np.arange(hp.nside2npix(64)))).T;psi=np.degrees(np.arccos(np.clip(V@lb(0,15),-1,1)))
sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));sm=sm/sm.mean()-1;np.save(O+f'bg_{TAG}_sm_ns64.npy',sm)
out=dict(tag=TAG,surf=os.environ.get('SURF','rho'),D=Dg,mc=mc,eta=eta,L=L,LD=LDX,poly=POLY,frac=fr,frac_over_mc=fr/mc,Ht0=r['Ht0'],age_Gyr=r['Ht0']*977.8/H0,eta_rec_D=r['eta_rec_D'],R_ls=Rls,
 shell_crossing=r['shell_crossing'],sm_rms=float(sm.std()),onepz_range=[float(o[:,1].min()),float(o[:,1].max())])
json.dump(out,open(O+f'bg_{TAG}.json','w'),indent=1);print('BG',out,flush=True)
