import json,numpy as np,healpy as hp,sys
from scipy.interpolate import CubicSpline
CON='.';B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
sys.path.insert(0,B+'/../../rotfit');import os;cw=os.getcwd();os.chdir(B);import fit;os.chdir(cw)
ns=32;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T
aG=fit.vec(0,15);psi=np.degrees(np.arccos(np.clip(V@aG,-1,1)))
T0=2.7255e6;l=np.arange(31)
def Dl(x):return l*(l+1)*hp.anafast(np.asarray(hp.remove_dipole(x)),lmax=30)*T0**2/(2*np.pi)
Pl=np.loadtxt(B+'/planckTT.txt')
def chi(D,lmax=30):
  c=0
  for r in Pl:
    L=int(r[0])
    if 2<=L<=lmax:e=r[3] if D[L]>r[1] else r[2];c+=((D[L]-r[1])/e)**2
  return float(c)
def smooth(f):
  o=np.load(f)['o'];o=o[np.argsort(o[:,0])];sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));return Dl(sm/sm.mean()-1)
pl={int(r[0]):r[1] for r in Pl};LL=np.arange(2,31);P=np.array([pl[L] for L in LL])
def st(D):return dict(chi2=chi(D),chi2_2_6=chi(D,6),D2=float(D[2]),r=float(np.corrcoef(D[2:31],P)[0,1]),Dl=D.tolist())
th=np.loadtxt('/tmp/th.txt');TH=np.interp(l,th[:,0],th[:,1])
R={'red_inv_t0.5':st(smooth('inversion/runs/red_inv_t0.5.npz')),'red_base':st(smooth('runs_poly/f20000_p_-0.041_-0.035_0.1353_0.npz')),'LCDM':st(TH)}
json.dump(R,open('red_inv_t0.5_spectrum.json','w'),indent=1)
for k,v in R.items():print(k,{q:v[q] for q in('chi2','chi2_2_6','D2','r')})
I=json.load(open('inversion/inv1_results.json'))['red_inv_t0.5'];am=np.array(I['a_motion']);ag=np.array(I['a_grav'])
f=lambda a:l*(l+1)*a**2/(2*l+1)/(2*np.pi)
for nm,A in(('inv_model_Ag1083',I['params']['Ag']),('inv_model_Ag313',I['Ag_free']),('inv_motion_only',0)):
  R[nm]=st(f(am[:31]+A*ag[:31]));R[nm]['Ag']=A;print(nm,{q:R[nm][q] for q in('chi2','chi2_2_6','D2','r')})
json.dump(R,open('red_inv_t0.5_spectrum.json','w'),indent=1)
