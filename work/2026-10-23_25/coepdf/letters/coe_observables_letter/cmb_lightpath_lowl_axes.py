# Quadrupole/octopole axes (max angular-momentum dispersion, de Oliveira-Costa et al. 2004) of the monopole+dipole-removed map vs axis of evil (260,60).
import numpy as np,healpy as hp,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
m=np.load('cmb_lightpath_maps.npz')['tot_nodip'];alm=hp.map2alm(m,lmax=3,iter=3)
def full(l):
    a=np.zeros(2*l+1,complex)
    for mm in range(l+1):
        v=alm[hp.Alm.getidx(3,l,mm)];a[l+mm]=v
        if mm:a[l-mm]=(-1)**mm*np.conj(v)
    return a
def axis(l):
    ms=np.arange(-l,l+1);Lz=np.diag(ms).astype(complex);Lp=np.zeros((2*l+1,)*2,complex)
    for i,mm in enumerate(ms[:-1]):Lp[i+1,i]=np.sqrt(l*(l+1)-mm*(mm+1))
    Lx=(Lp+Lp.conj().T)/2;Ly=(Lp-Lp.conj().T)/2j;a=full(l);Ls=[Lx,Ly,Lz]
    A=np.array([[np.real(np.vdot(Li@a,Lj@a)) for Lj in Ls] for Li in Ls]);w,v=np.linalg.eigh(A);n=v[:,-1]
    lon,lat=hp.vec2ang(n,lonlat=True);return n,float(lon[0]),float(lat[0])
def ang(a,b):return float(np.degrees(np.arccos(abs(np.clip(a@b,-1,1)))))
aoe=np.array(hp.ang2vec(260,60,lonlat=True));coe=np.array(hp.ang2vec(180,-15,lonlat=True))
R={}
for l in(2,3):
    n,lo,la=axis(l);R[l]=dict(l=lo,b=la,angle_to_AoE=ang(n,aoe),angle_to_CoE_axis=ang(n,coe))
R['CoE_axis_to_AoE']=ang(coe,aoe);json.dump(R,open('cmb_lightpath_lowl_axes.json','w'),indent=1);print(R)
fig=plt.figure(figsize=(9,5));hp.mollview(m*2.7255e6,fig=fig.number,title='Monopole+dipole removed total [uK] with l=2,3 axes',cmap='coolwarm')
for l,c in((2,'g'),(3,'m')):
    hp.projscatter(R[l]['l'],R[l]['b'],lonlat=True,c=c,s=90,marker='o',label=f'l={l} axis');hp.projscatter((R[l]['l']+180)%360,-R[l]['b'],lonlat=True,c=c,s=90,marker='o')
hp.projscatter(260,60,lonlat=True,c='k',marker='x',s=100,label='AoE (260,60)');hp.projscatter(180,-15,lonlat=True,c='purple',marker='*',s=120,label='CoE')
plt.legend(loc='lower right',fontsize=7);plt.savefig('cmb_lightpath_lowl_axes.png',dpi=120,bbox_inches='tight')
