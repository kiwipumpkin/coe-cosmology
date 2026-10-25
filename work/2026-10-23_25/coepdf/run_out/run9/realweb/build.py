# realistic web: groups of independent Poisson-thinned bodies (halos PS 1e14-2e15 Msun, SvdW compensated voids 10-50 Mpc)
# painted (CIC + Gaussian smoothing at cell scale) on two grids; FFT gravity; tidal sky signal on spheres R
import os,sys,json,time,numpy as np,healpy as hp,scipy.fft as F
from scipy.ndimage import map_coordinates
from cosmo import *
OUT='/workspace/redweb2/realweb/';NG=100;RS=[0.5,0.9,2.0,4.0,8.0,14.0]
MSUN=1.989e30;KG=7.7e51/0.0039;RV=14.3
ns_=32;V=np.array(hp.pix2vec(ns_,np.arange(hp.nside2npix(ns_)))).T
CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937);aC=np.array([1,0,np.tan(np.radians(15))]);aC/=np.linalg.norm(aC);PC=42.0*aC
shp=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
GR=[dict(name='fine',L=20.0,N=384,rin=0.0,rout=5.0,seed=11),dict(name='coarse',L=60.0,N=384,rin=5.0,rout=15.0,seed=12)]
Mb=np.geomspace(1e14,2e15,41);nM,sM=ps_dndlnM(Mb);bM=1+(dc**2/sM**2-1)/dc   # PS peak-background bias
Rv=np.geomspace(10,50,41);nV=svdw_dndlnR(Rv);BV=-1.5;DV=-0.85
def setup(g):
    N,L=g['N'],g['L'];dx=L/N;x=(np.arange(N)+0.5)*dx-L/2
    k1=2*np.pi*F.fftfreq(N,dx);k3=2*np.pi*F.rfftfreq(N,dx)
    KX,KY,KZ=np.meshgrid(k1,k1,k3,indexing='ij',sparse=True);K2=KX**2+KY**2+KZ**2
    rng=np.random.default_rng(g['seed']);w=rng.standard_normal((N,N,N),dtype=np.float32)
    kM=np.sqrt(K2)/1000.0   # 1/Mpc
    Pk=np.where(K2>0,P(np.maximum(kM,1e-6)),0)/(dx*1000)**3
    d=F.irfftn(F.rfftn(w,workers=8)*np.sqrt(Pk).astype(np.float32)*np.exp(-K2*dx*dx/2),s=w.shape,workers=8).astype(np.float32);del w
    X,Y,Z=np.meshgrid(x,x,x,indexing='ij',sparse=True);r=np.sqrt(X**2+Y**2+Z**2)
    mask=(r>=g['rin'])&(r<g['rout'])
    mod=shp(np.sqrt((X-PC[0])**2+(Y-PC[1])**2+(Z-PC[2])**2)/42).astype(np.float32);mod/=mod[mask].mean()
    return dict(dx=dx,x0=-L/2,KX=KX,KY=KY,KZ=KZ,K2=K2,d=d,mask=mask,mod=mod,N=N,L=L)
def cic(G,P_,m):
    N=G['N'];u=(P_-G['x0'])/G['dx']-0.5;i=np.floor(u).astype(np.int64);f=u-i;out=np.zeros(N**3)
    for a in (0,1):
      for b in (0,1):
        for c in (0,1):
          wgt=m*np.abs(1-a-f[:,0])*np.abs(1-b-f[:,1])*np.abs(1-c-f[:,2])
          ii=((i[:,0]+a)%N)*N*N+((i[:,1]+b)%N)*N+((i[:,2]+c)%N);out+=np.bincount(ii,wgt,N**3)
    return out.reshape(N,N,N)
def bodies(G,rng,frac):
    dV=(G['dx']*1000)**3;mk=G['mask'];d=G['d'][mk];mod=G['mod'][mk];idx=np.flatnonzero(mk)
    out=[];N=G['N']
    def place(cnt):
        c=np.repeat(idx,cnt);ii=np.stack(np.unravel_index(c,(N,N,N)),1);return G['x0']+(ii+rng.random(ii.shape))*G['dx']
    # halos, by mass bin
    dl=np.diff(np.log(Mb));Hm=[];Hp=[]
    for j in range(len(Mb)-1):
        lam=0.5*(nM[j]+nM[j+1])*dl[j]*dV*np.maximum(0,1+0.5*(bM[j]+bM[j+1])*d)*mod*frac;cnt=rng.poisson(lam)
        Hp.append(place(cnt));Hm.append(np.exp(rng.uniform(np.log(Mb[j]),np.log(Mb[j+1]),cnt.sum())))
    Hp=np.vstack(Hp);Hm=np.concatenate(Hm)
    dl=np.diff(np.log(Rv));Vp=[];Vr=[]
    for j in range(len(Rv)-1):
        lam=0.5*(nV[j]+nV[j+1])*dl[j]*dV*np.maximum(0,1+BV*d)*mod*frac;cnt=rng.poisson(lam)
        Vp.append(place(cnt));Vr.append(np.exp(rng.uniform(np.log(Rv[j]),np.log(Rv[j+1]),cnt.sum())))
    Vp=np.vstack(Vp);Vr=np.concatenate(Vr);Md=-DV*rho_m*4/3*np.pi*Vr**3   # deficit mass (Msun), local mean density
    return Hp,Hm,Vp,Vr,Md
def paint(G,rng,Hp,Hm,Vp,Vr,Md):
    rho=cic(G,Hp,Hm)+cic(G,Vp,-Md)
    for q in range(8):   # compensating wall: 8 random points on shell at 1.25 Rv
        n=rng.normal(size=Vp.shape);n/=np.linalg.norm(n,axis=1)[:,None];rho+=cic(G,Vp+n*(1.25*Vr/1000)[:,None],Md/8)
    rho[~G['mask']]=0;rho[G['mask']]-=rho[G['mask']].mean()   # neutral against uniform background in the region
    return rho/G['dx']**3   # Msun/Gpc^3
def accel(G,rho):
    rk=F.rfftn(rho.astype(np.float32),workers=8)*np.exp(-G['K2']*G['dx']**2/2)   # finite extent: Gaussian at cell scale
    with np.errstate(divide='ignore',invalid='ignore'):ph=np.where(G['K2']>0,-4*np.pi*rk/G['K2'],0)   # G=1, Phi
    return [F.irfftn(-1j*Kc*ph,s=rho.shape,workers=8) for Kc in (G['KX'],G['KY'],G['KZ'])]   # a=-grad Phi, Msun/Gpc^2
def skyfrom(G,a):
    pts=np.vstack([np.zeros((1,3))]+[R*V for R in RS]);u=((pts-G['x0'])/G['dx']-0.5).T
    A=np.stack([map_coordinates(c,u,order=1,mode='wrap') for c in a],1);a0=A[0];out=[]
    for i,R in enumerate(RS):Ai=A[1+i*len(V):1+(i+1)*len(V)];s=((Ai-a0)*V).sum(1);out.append(s-s.mean())
    return np.array(out)*MSUN/KG*RV**2   # same units as sky(): sum (M/KG) n.[d/|d|^3-y/|y|^3], lengths in RV
if __name__=='__main__':
    ck=OUT+'groups.npz';done=dict(np.load(ck)) if os.path.exists(ck) else {}
    Gs=[setup(g) for g in GR];print('setup done',[float(G['d'][G['mask']].std()) for G in Gs],flush=True)
    for G,g in zip(Gs,GR):np.save(OUT+f"delta_{g['name']}_slice.npy",G['d'][:,:,G['N']//2])
    S=done.get('S',np.zeros((0,2,len(RS),len(V))));stats=list(done.get('stats',np.zeros((0,2,4))))
    for k in range(len(S),NG):
        t=time.time();rng=np.random.default_rng(1000+k);row=[];st=[]
        for G,g in zip(Gs,GR):
            B=bodies(G,rng,1/NG);rho=paint(G,rng,*B);row.append(skyfrom(G,accel(G,rho)));st.append([len(B[1]),B[1].sum(),len(B[3]),B[4].sum()])
            if k==0:np.savez(OUT+f"bodies_g0_{g['name']}.npz",Hp=B[0],Hm=B[1],Vp=B[2],Vr=B[3]);np.save(OUT+f"rho_g0_{g['name']}_slice.npy",rho[:,:,G['N']//2])
        S=np.concatenate([S,np.array(row)[None]]);stats.append(st);np.savez(ck,S=S,stats=np.array(stats))
        print(k,'%.0fs'%(time.time()-t),st,flush=True)
