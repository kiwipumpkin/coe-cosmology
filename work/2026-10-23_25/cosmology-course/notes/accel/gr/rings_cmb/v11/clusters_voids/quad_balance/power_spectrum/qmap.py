import sys,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/'
web=np.load(B+'reproduce_cmb/webrec_D42.npy');ns=hp.npix2nside(len(web));V=np.array(hp.pix2vec(ns,np.arange(len(web)))).T
o=np.load(B+'quad_balance/mass_family_scan/runs/soft_2_30.npz')['o'];o=o[np.argsort(o[:,0])]
psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));sm=sm/sm.mean()-1
sm=np.asarray(hp.remove_dipole(sm))
# isolate l=2 only
alm=hp.map2alm(sm,lmax=64)
f2=np.zeros_like(alm);ls,ms=hp.Alm.getlm(64)
f2[ls==2]=alm[ls==2]
q=hp.alm2map(f2,ns)*T0  # l=2 only map in uK
# full smooth in uK for reference
smu=sm*T0
# plot mollview of l=2 map, CoE on left: rotate so CoE axis points left. CoE direction fit.vec(0,15)
hp.mollview(q,title='Quadrupole (l=2 only) component of smooth CoE sky  [uK]',unit='uK',cmap='RdBu_r')
hp.graticule()
plt.savefig('quad_l2_map.png',dpi=110,bbox_inches='tight');plt.close()
# amplitude vs galactic longitude: average |q| and mean q over longitude bins
lon,lat=hp.pix2ang(ns,np.arange(len(q)),lonlat=True)
lon=np.mod(lon,360)
bins=np.arange(0,361,10);idx=np.digitize(lon,bins)-1
meanq=np.array([q[idx==i].mean() for i in range(36)])
rmsq=np.array([np.sqrt((q[idx==i]**2).mean()) for i in range(36)])
ctr=bins[:-1]+5
fig,ax=plt.subplots(figsize=(11,5))
ax.plot(ctr,meanq,'o-',color='crimson',label='mean l=2 signal')
ax.plot(ctr,rmsq,'s--',color='navy',label='rms l=2 amplitude')
ax.axhline(0,color='gray',lw=.6);ax.set_xlabel('galactic longitude l (deg)');ax.set_ylabel('quadrupole contribution (uK)')
ax.set_title('Where the quadrupole acts vs longitude');ax.legend();ax.grid(alpha=.3)
plt.tight_layout();plt.savefig('quad_vs_longitude.png',dpi=120);plt.close()
# find extrema of l=2 map
imax=np.argmax(q);imin=np.argmin(q)
print('l=2 map max %.0f uK at (l,b)='%q[imax],np.round(hp.pix2ang(ns,imax,lonlat=True),1))
print('l=2 map min %.0f uK at (l,b)='%q[imin],np.round(hp.pix2ang(ns,imin,lonlat=True),1))
print('CoE axis dir vec',np.round(fit.vec(0,15),3),'(l,b)=',np.round(hp.vec2ang(fit.vec(0,15),lonlat=True),1))
print('l=2 rms over sky uK',np.std(q))
