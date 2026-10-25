# web peculiar velocity field (unit mass scale s=1) on observer-centred grid: galactic-plane dirs (1 deg) and healpix nside16, distances DG
import os,sys,numpy as np,healpy as hp
A='/workspace/coe-cosmology-gh/draft/current/analyses';sys.path.insert(0,A+'/void_mass_network');os.chdir(A+'/shell_ensemble')
import vmn_gen as G
z=np.load('/workspace/redweb2/densweb/densweb_s1.npz');X,M=z['X'],z['M']   # Gpc (observer-centred, galactic xyz), kg
print('bodies',len(M),'sumM',M.sum(),flush=True)
GN=6.674e-11;GPC=3.0857e25;H0=70e3/3.0857e22;Om=0.31;f=Om**0.55;Fv=2*f/(3*H0*Om)   # linear theory v=2f g/(3 H Om)
DG=np.geomspace(0.005,2.0,40)
lp=np.radians(np.arange(360.));NP=np.c_[np.cos(lp),np.sin(lp),0*lp];NH=np.array(hp.pix2vec(16,np.arange(hp.nside2npix(16)))).T
def g(P):
    out=np.zeros((len(P),3))
    for i in range(0,len(M),20000):
        d=X[None,i:i+20000]-P[:,None];r=np.linalg.norm(d,axis=2);out+=np.einsum('pb,pbk->pk',M[i:i+20000]/r**3,d)
    return out*GN/GPC**2   # m/s^2
g0=g(np.zeros((1,3)))[0]
def vlos(N):
    V=np.empty((len(N),len(DG)))
    for k,dd in enumerate(DG):
        P=N*dd;gg=np.vstack([g(P[j:j+500]) for j in range(0,len(P),500)])
        V[:,k]=Fv*((gg-g0)*N).sum(1)/1e3   # km/s, relative to observer
        print(k,dd,np.abs(V[:,k]).max(),flush=True)
    return V
np.savez('/workspace/redweb2/densweb/web_v_s1.npz',DG=DG,Vp=vlos(NP),Vh=vlos(NH),v0=Fv*g0/1e3,n=len(M))
