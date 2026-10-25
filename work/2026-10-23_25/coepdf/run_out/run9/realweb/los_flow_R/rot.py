import numpy as np,healpy as hp,itertools
sm=np.load('/workspace/redweb2/sm_red.npy');LM=8
th0,ph0=np.radians(105),np.radians(180);v0=hp.ang2vec(th0,ph0);V=np.array(hp.pix2vec(32,np.arange(12288))).T
tst=hp.map2alm(V@v0,lmax=LM)
for a,b,c in itertools.product([0,ph0,-ph0],[th0,-th0],[0,ph0,-ph0]):
  t=tst.copy();hp.rotate_alm(t,a,b,c)
  if abs(t[hp.Alm.getidx(LM,1,1)])<1e-6:print('OK',a,b,c);good=(a,b,c);break
# find smooth red axis: maximize m=0 power at l=3 over directions
a=hp.map2alm(hp.remove_dipole(sm),lmax=LM);cl=hp.alm2cl(a);print('smooth red C_l ratio l3 share:',cl[3]/cl[2:].sum())
best=(0,None)
for p in range(768):
  th,ph=hp.pix2ang(8,p);t=a.copy();hp.rotate_alm(t,-ph,-th,0) if True else 0
  f=abs(t[hp.Alm.getidx(LM,3,0)])**2/cl[3]/7
  if f>best[0]:best=(f,(90-np.degrees(th),np.degrees(ph)))
print('smooth red l3 m0 fraction best axis (b,l):',best)
t=a.copy();hp.rotate_alm(t,*good);P=[abs(t[hp.Alm.getidx(LM,3,m)])**2*(1 if m==0 else 2) for m in range(4)];print('CoE frame l3 m-fractions',np.round(np.array(P)/sum(P),3))
np.save('rotgood.npy',np.array(good))
