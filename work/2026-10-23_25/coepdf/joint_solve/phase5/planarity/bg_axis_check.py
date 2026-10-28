# Locate the true symmetry axis of the CoE-only background (min-<Lz^2> eigenvector) and dipole direction
import numpy as np,healpy as hp
exec(open('planarity.py').read().split('def stats')[0])
B=d['B'];alm=hp.map2alm(B,lmax=LM)
for l in (1,2,3):
    if l==1:
        a=hp.map2alm(B,lmax=1);v=np.array([np.real(a[hp.Alm.getidx(1,1,1)])*-np.sqrt(2),np.imag(a[hp.Alm.getidx(1,1,1)])*np.sqrt(2),np.real(a[hp.Alm.getidx(1,1,0)])])
        print('dipole dir (l,b)',hp.vec2ang(v/np.linalg.norm(v),lonlat=True));continue
    v=full(alm,l);L=LL[l];M=np.array([[np.real(v.conj()@L[i]@L[j]@v) for j in range(3)] for i in range(3)]);M=(M+M.T)/2
    w,e=np.linalg.eigh(M);print(l,'eig',w/np.real(v.conj()@v),'min-axis (l,b)',hp.vec2ang(e[:,0],lonlat=True))
