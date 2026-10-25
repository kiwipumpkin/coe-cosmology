import json,sys,glob,re,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/'
web=np.load(B+'reproduce_cmb/webrec_D42.npy');ns=hp.npix2nside(len(web));V=np.array(hp.pix2vec(ns,np.arange(len(web)))).T
web=np.asarray(hp.remove_dipole(web-web.mean()))
psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)))
LMAX=600;l=np.arange(LMAX+1);Dl=lambda m:l*(l+1)*hp.anafast(m,lmax=LMAX)*T0**2/(2*np.pi)
P=np.loadtxt(B+'reproduce_cmb/planckTT.txt');pl=P[:,0].astype(int);PD=P[:,1];sig=(P[:,2]+P[:,3])/2
Dw=Dl(web);sel=(pl>=3)&(pl<=300)
res=[];best=None
for f in sorted(glob.glob(B+'quad_balance/mass_family_scan/runs_md/*.npz')):
  D,M=map(float,re.findall(r'D([\d.]+)_M([\d.]+)\.npz',f)[0])
  o=np.load(f)['o'];o=o[np.argsort(o[:,0])]
  sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1))
  Dt=Dl(sm+web);Ds=Dl(sm)
  chi=float(np.sum(((Dt[pl[sel]]-PD[sel])/sig[sel])**2)/sel.sum())
  r=dict(D_Gpc=D,M=M,D2=float(Dt[2]),D2_planck=float(Dt[2]/PD[pl==2][0]),D2_smooth=float(Ds[2]),D3=float(Dt[3]),chi2nu_l3_300=chi,zspread=float(o[:,1].max()/o[:,1].min()-1))
  res.append(r);print(r)
  if best is None or r['D2']<best[0]['D2']: best=(r,Dt,Ds)
json.dump(dict(rows=res,best_by_D2=best[0],web_chi2nu=float(np.sum(((Dw[pl[sel]]-PD[sel])/sig[sel])**2)/sel.sum()),note='web map webrec_D42 held fixed'),open('massdist_fit.json','w'),indent=1)
r,Dt,Ds=best;fig,ax=plt.subplots(1,2,figsize=(18,7))
for a,xs in zip(ax,('log','linear')):
  a.plot(l[2:],Dt[2:],'k',label=f"best fit M={r['M']:g} M_crit, D={r['D_Gpc']:g} Gpc: smooth+web");a.plot(l[2:],Ds[2:],'C3--',label='smooth only');a.plot(l[2:],Dw[2:],'C0:',label='web only')
  a.errorbar(pl,PD,yerr=[P[:,2],P[:,3]],fmt='.',c='gold',ms=3,elinewidth=.4,label='Planck 2018 TT')
  a.set_yscale('log');a.set_xscale(xs);a.set_xlim(2,LMAX);a.set_ylim(1e-1,1e9);a.set_xlabel('multipole l');a.set_ylabel('D_l [μK²]');a.legend(fontsize=9)
fig.suptitle('Mass–distance best fit vs Planck 2018 TT');plt.savefig('../power_spectrum_massdist_fit.png',dpi=110,bbox_inches='tight')
