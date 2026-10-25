import json,sys,glob,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import Legendre
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/'
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T
psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));LMAX=20;l=np.arange(LMAX+1)
mu=np.linspace(-1,1,4001);rows=[]
for fn in glob.glob(B+'quad_balance/mass_family_scan/runs_rec/*.npz'):
  d=np.load(fn);o=d['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]))
  sm=np.exp(-cs(psi));sm=sm/sm.mean()-1;Cl0=hp.anafast(sm,lmax=LMAX)
  smd=np.asarray(hp.remove_dipole(sm));Dl=l*(l+1)*hp.anafast(smd,lmax=LMAX)*T0**2/(2*np.pi)
  A=np.exp(-cs(np.degrees(np.arccos(mu))));dl=A/np.trapezoid(A,mu)*2-1
  a=[(2*L+1)/2*np.trapezoid(dl*Legendre.basis(L)(mu),mu) for L in (1,2,3)]
  rows.append(dict(f=float(d['fB']),D2=float(Dl[2]),D2_planck=float(Dl[2]/226.1),a1=float(a[0]),a2=float(a[1]),a3=float(a[2]),D1_uK2=float(2*Cl0[1]*T0**2/(2*np.pi)),D3=float(Dl[3]),
   dt_emit_range=float(o[:,2].max()-o[:,2].min()),z0=float(o[0,1]),z180=float(o[-1,1])))
rows.sort(key=lambda r:r['f']);[print({k:(round(v,6) if abs(v)<1 else float('%.4g'%v)) for k,v in r.items()}) for r in rows]
json.dump(rows,open('recomb_sequential_quad.json','w'),indent=1)
f=[r['f'] for r in rows];fig,ax=plt.subplots(1,2,figsize=(14,5.5))
ax[0].semilogy(f,[r['D2'] for r in rows],'ko-',label='model D_2 (dipole removed)');ax[0].axhline(226.1,c='gold',lw=2,label='Planck D_2 = 226 μK²')
ax[0].axvline(0,c='gray',ls=':');ax[0].axvline(1,c='gray',ls=':');ax[0].set_xlabel('sequentiality f  (0 = simultaneous, 1 = physical surface)');ax[0].set_ylabel('D_2 [μK²]');ax[0].legend()
ax[1].plot(f,[r['a2'] for r in rows],'C3o-',label='a2 (P2 projection of A/A_all)');ax[1].plot(f,[r['a1'] for r in rows],'C0s--',label='a1 (dipole)');ax[1].axhline(0,c='k',lw=.5);ax[1].set_xlabel('f');ax[1].set_ylabel('dT/T Legendre coefficient');ax[1].legend()
fig.suptitle('Soft p=2, M=30, D=42: quadrupole vs recombination sequentiality');plt.savefig('../recomb_sequential_quad.png',dpi=110,bbox_inches='tight')
