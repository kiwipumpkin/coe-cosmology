import json,glob,sys,numpy as np,healpy as hp
from scipy.interpolate import CubicSpline,PchipInterpolator
from numpy.polynomial.legendre import Legendre
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;Q='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/'
ns=64;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));l=np.arange(31);mu=np.linspace(-1,1,4001)
M={float(k):v for k,v in json.load(open(Q+'power_spectrum/slope_refit/migkas_vs_ratio_dense.json')).items()}
for f in glob.glob('mig_*.json'):M.update({float(k):v for k,v in json.load(open(f)).items()})
rr=np.array(sorted(M));fc=PchipInterpolator(rr,[M[r]['chi2'] for r in rr])
def an(fn):
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1));Dl=l*(l+1)*hp.anafast(sm,lmax=30)*T0**2/(2*np.pi)
  y=-cs(np.degrees(np.arccos(mu)));y=y-0.5*np.trapezoid(y,mu)
  a=[(2*L+1)/2*np.trapezoid(y*Legendre.basis(L)(mu),mu) for L in (1,2,3,4)]
  rat=float(d['ratio']);return dict(run=fn,C=d['C'].tolist() if 'C' in d.files else None,ratio=rat,a=a,Dl=Dl.tolist(),chi2_interp=float(fc(rat)) if rr[0]<=rat<=rr[-1] else None)
if __name__=='__main__':
  R=[an(f) for f in sorted(glob.glob('runs_poly/*.npz'))]
  for r in R:print(r['C'],'rat %.5f a2 %+.3e D2 %.3g D3 %.3g D4 %.3g chi2 %s'%(r['ratio'],r['a'][1],r['Dl'][2],r['Dl'][3],r['Dl'][4],r['chi2_interp']))
  json.dump(R,open('poly_scan.json','w'))
