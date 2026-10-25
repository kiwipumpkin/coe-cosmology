import numpy as np,healpy as hp
from core import LM,l
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
Pl=np.loadtxt(B+'/planckTT.txt');ClP=np.zeros(LM+1)
for r in Pl:
    if int(r[0])<=LM:ClP[int(r[0])]=r[1]*2*np.pi/(r[0]*(r[0]+1))
NA=32;AX=np.array(hp.pix2vec(NA,np.arange(hp.nside2npix(NA)))).T
def obs_grid(fn):
    a=hp.read_alm(fn)*1e6;Ga=np.zeros((LM+1,len(AX)));Cl=hp.alm2cl(a)
    for L in range(2,LM+1):
        f=np.zeros(LM+1);f[L]=1;Ga[L]=np.sqrt(4*np.pi/(2*L+1))*hp.alm2map(hp.almxfl(a,f),NA)
    return Ga,Cl  # Ga[l,axis]=a_l0 (muK) with pole at that axis
