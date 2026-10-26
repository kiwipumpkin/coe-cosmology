# Phase3 (NEW): H0 set by kinematic speed U=H0*D/c (famP U env), D fixed. From phase2: rebuild CoE smooth CMB background (sm) at chosen D, parameterised (DRV=D/RV) via the original famP.py poly trace,
# then map it exactly as thin4_run.py did (psi from anti-CoE axis (0,15), sm=exp(-spline ln(1+z)) normalised).  rs_cc omitted (rms 1.3e-8 in hl_R4).
import os,sys,subprocess,numpy as np,healpy as hp
from scipy.interpolate import CubicSpline
D=float(sys.argv[1]);RV=14.30048287260589;C='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained'
H0=float(sys.argv[2]);U=H0/(299792.458/(D*1000))/1.00016;O='/workspace/coepdf/joint_solve/phase3/';f=O+f'p2_trace_D{D:g}_H{H0:g}.npz'
if not os.path.exists(f):
    env=dict(os.environ,U=str(U),DRV=str(D/RV),C1='-0.041',C2='-0.035',C3='0.1353',C4='0')
    subprocess.run(['python3',C+'/famP.py','poly','0','20000',f],env=env,check=True,cwd=C)
o=np.load(f)['o'];o=o[np.argsort(o[:,0])]
def lb(l,b):l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
ns=64;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;psi=np.degrees(np.arccos(np.clip(V@lb(0,15),-1,1)))
sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));sm=sm/sm.mean()-1
np.save(O+f'p2_sm_D{D:g}_H{H0:g}_ns64.npy',sm)
cl=hp.anafast(sm,lmax=6);print('D',D,'H0',H0,'U',U,'rms',sm.std(),'per-l rms',np.sqrt(cl*(2*np.arange(7)+1)/4/np.pi),flush=True)
ref=np.load('/workspace/coepdf/run_out/run9/realweb/los_flow_R/hl_R4_ns64.npz')['sm'];print('corr with hl_R4 sm',np.corrcoef(sm,ref)[0,1],'ratio rms',sm.std()/ref.std())
