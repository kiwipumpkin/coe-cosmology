import json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,sys
from scipy.interpolate import CubicSpline
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;P=np.loadtxt('../../reproduce_cmb/planckTT.txt');Pl={int(l):np.sqrt((2*l+1)*2*np.pi*d/(l*(l+1))/(4*np.pi))/27.255 for l,d in P[:12,:2]};Pl[1]=3362/np.sqrt(3)/27.255
# recombination-era imprints (from reproduce_cmb/web_rec: Bernardeau lin + /1107 growth)
dV=np.array([-2.8644e-4,-6.6249e-4]);dC=np.array([1.21747e-3,1.32462e-3]);tV=np.log1p(dV)/3;tC=np.log1p(dC)/3
print('void dT/T range',tV,'cluster dT/T range',tC)
ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;aG=fit.vec(0,15);x=V@aG
o=np.load('../mass_family_scan/runs/soft_2_30.npz')['o'];o=o[np.argsort(o[:,0])]
base=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(np.degrees(np.arccos(np.clip(x,-1,1)))));base=hp.remove_dipole(base/base.mean()-1)
def per_l(m):cl=hp.anafast(m,lmax=6);l=np.arange(7);return np.sqrt((2*l+1)*cl/(4*np.pi))/1e-5
def zonal_a(m,l):  # signed Legendre amplitude about CoE axis, per-l rms units 1e-5
    from numpy.polynomial.legendre import legval;c=np.zeros(l+1);c[l]=1;return (2*l+1)*np.mean(m*legval(x,c))/np.sqrt(2*l+1)/1e-5
b=per_l(base);print('base per-l',np.round(b,2),'signed a20',zonal_a(base,2))
# (A) upper bound: entire sky tiled, hottest clusters in caps |x|>x0, deepest voids elsewhere (no overlap, filling factor 1)
best=None
for x0 in np.linspace(0,1,201):
    corr=np.where(np.abs(x)>x0,tC.max(),tV.min());corr-=corr.mean();a2=zonal_a(corr,2)
    if best is None or a2>best[1]:best=(x0,a2,corr)
x0,a2,corrA=best;tot=hp.remove_dipole(base+corrA);pA=per_l(tot);pc=per_l(corrA)
print(f'BOUND x0={x0:.3f} (cap half-angle {np.degrees(np.arccos(x0)):.1f} deg) correction a20=+{a2:.2f}e-5 = {a2/246.6*100:.1f}% of needed 246.6')
# (B) realistic packed disks: clusters r=0.30 Gpc, voids r=2.14 Gpc on the 12.0 Gpc surface (max sizes in generator), hex-ish Fibonacci packing
def disks(ang,mask_fn,N):
    i=np.arange(N)+.5;z=1-2*i/N;ph=np.pi*(1+5**.5)*i;c=np.c_[np.sqrt(1-z*z)*np.cos(ph),np.sqrt(1-z*z)*np.sin(ph),z]
    # rotate so z->aG
    k=np.cross([0,0,1],aG);s=np.linalg.norm(k);from scipy.spatial.transform import Rotation as R;c=c@R.from_rotvec(k/s*np.arcsin(s)).as_matrix().T
    return c[mask_fn(c@aG)]
corrB=np.zeros(len(x));rc=np.arcsin(0.30/12.0);rv=np.arcsin(2.14/12.0)
NC=int(4/(rc**2)*0.9);NV=int(4/(rv**2)*0.9)
cc=disks(rc,lambda u:np.abs(u)>x0,NC);vv=disks(rv,lambda u:np.abs(u)<=x0-np.sin(rv)*0.5,NV)
for c in cc:corrB[hp.query_disc(ns,c,rc)]+=tC.max()
for c in vv:corrB[hp.query_disc(ns,c,rv)]+=tV.min()
corrB-=corrB.mean();totB=hp.remove_dipole(base+corrB);pB=per_l(totB)
print('packed: Nclusters',len(cc),'Nvoids',len(vv),'correction a20',zonal_a(corrB,2))
R=dict(void_dTT=tV.tolist(),cluster_dTT=tC.tolist(),base_per_l_1e5=b.tolist(),needed_a20=246.6,
  bound=dict(cap_x0=x0,cap_halfangle_deg=float(np.degrees(np.arccos(x0))),a20_supplied=a2,frac_of_needed=a2/246.6,final_per_l_1e5=pA.tolist(),corr_per_l_1e5=pc.tolist()),
  packed=dict(N_clusters=len(cc),N_voids=len(vv),a20_supplied=zonal_a(corrB,2),final_per_l_1e5=pB.tolist(),corr_per_l_1e5=per_l(corrB).tolist()),
  planck_per_l_1e5={l:Pl[l] for l in range(1,7)})
json.dump(R,open('vc_results.json','w'),indent=1,default=float)
print(f"{'l':>2} {'base':>9} {'bound':>9} {'packed':>9} {'packedcorr':>10} {'Planck':>7} | final(packed) uK, xPlanck; corr-only xPlanck")
for l in range(1,7):print(f"{l:2d} {b[l]:9.2f} {pA[l]:9.2f} {pB[l]:9.2f} {per_l(corrB)[l]:10.3f} {Pl[l]:7.3f} | {pB[l]*27.255:8.1f} {pB[l]/Pl[l]:8.1f}  {per_l(corrB)[l]/Pl[l]:6.2f}")
fig=plt.figure(figsize=(18,12));mk=[((264.02,48.25),'Planck hot','r'),((84.02,-48.25),'Planck cold','b'),((180,-15),'CoE','m'),((0,15),'anti-CoE','g')]
for i,(m,t) in enumerate(((base,'BEFORE: best combo (soft p=2, M=30), dipole removed'),(totB,'AFTER: + packed clusters in CoE caps / voids in band (rec-era contrasts)'),(corrB,'body correction alone'))):
    hp.mollview(m/1e-5,fig=fig.number,sub=(2,2,i+1),cmap='RdBu_r',title=t,unit='dT/T ×1e5')
    for (L,B_),s,c in mk:hp.projscatter(L,B_,lonlat=True,marker='*',s=150,c=c,edgecolors='k');hp.projtext(L,B_+6,s,lonlat=True,fontsize=9)
ax=fig.add_subplot(2,2,4);l=np.arange(1,7);w=.2
for j,(y,lab) in enumerate(((b[1:],'before'),(pB[1:],'after (packed)'),(pA[1:],'after (full-coverage bound)'),([Pl[k] for k in l],'Planck'))):ax.bar(l+(j-1.5)*w,np.maximum(y,1e-3),w,label=lab)
ax.set_yscale('log');ax.set_xlabel('l');ax.set_ylabel('per-l rms (1e-5)');ax.legend();plt.savefig('void_cluster_cancel.png',dpi=100,bbox_inches='tight')
