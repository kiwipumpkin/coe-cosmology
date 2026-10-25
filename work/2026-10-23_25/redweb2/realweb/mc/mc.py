# Monte Carlo of full realistic universes (natural amplitude, strength 1): new GRF + Poisson halo field per universe.
# Fast variant of build.py: halos deposited per cell (Poisson counts in 8 mass bins, NGP at cell, Gaussian smoothing at cell scale);
# compensated voids omitted (net ~0 at >=cell scale; checked against the 100-group sweep via web rms). N=256 grids.
import os,sys,time,json,numpy as np,scipy.fft as F
sys.path.insert(0,'/workspace/redweb2/realweb');import build as B
from cosmo import *
OUT='/workspace/redweb2/realweb/mc/';NU=int(sys.argv[1]) if len(sys.argv)>1 else 1000
GR=[dict(name='fine',L=20.0,N=256,rin=0.0,rout=5.0),dict(name='coarse',L=60.0,N=256,rin=5.0,rout=15.0)]
Me=np.geomspace(1e14,2e15,9);Mf=np.geomspace(1e14,2e15,400);nf,sf=ps_dndlnM(Mf);bf=1+(dc**2/sf**2-1)/dc
bins=[];lnM=np.log(Mf)
for j in range(8):
    s=(Mf>=Me[j])&(Mf<=Me[j+1]);n=np.trapezoid(nf[s],lnM[s]);bins.append((n,np.trapezoid(nf[s]*Mf[s],lnM[s])/n,np.trapezoid(nf[s]*bf[s],lnM[s])/n))
def setup(g):
    N,L=g['N'],g['L'];dx=L/N;x=(np.arange(N)+0.5)*dx-L/2;k1=2*np.pi*F.fftfreq(N,dx);k3=2*np.pi*F.rfftfreq(N,dx)
    KX,KY,KZ=np.meshgrid(k1,k1,k3,indexing='ij',sparse=True);K2=KX**2+KY**2+KZ**2
    X,Y,Z=np.meshgrid(x,x,x,indexing='ij',sparse=True);r=np.sqrt(X**2+Y**2+Z**2);mask=(r>=g['rin'])&(r<g['rout'])
    mod=B.shp(np.sqrt((X-B.PC[0])**2+(Y-B.PC[1])**2+(Z-B.PC[2])**2)/42).astype(np.float32);mod/=mod[mask].mean()
    Pk=(np.where(K2>0,P(np.maximum(np.sqrt(K2)/1000,1e-6)),0)/(dx*1000)**3);amp=(np.sqrt(Pk)*np.exp(-K2*dx*dx/2)).astype(np.float32)
    return dict(dx=dx,x0=-L/2,KX=KX,KY=KY,KZ=KZ,K2=K2,mask=mask,mod=mod[mask],amp=amp,N=N)
Gs=[setup(g) for g in GR]
def universe(seed):
    rng=np.random.default_rng(seed);out=[];nh=0
    for G in Gs:
        N=G['N'];w=rng.standard_normal((N,N,N),dtype=np.float32);d=F.irfftn(F.rfftn(w,workers=8)*G['amp'],s=w.shape,workers=8)[G['mask']];del w
        dV=(G['dx']*1000)**3;rho=np.zeros(int(G['mask'].sum()))
        for n,Mm,b in bins:c=rng.poisson(n*dV*np.maximum(0,1+b*d)*G['mod']);rho+=c*Mm;nh+=c.sum()
        rho-=rho.mean();R=np.zeros((N,N,N),np.float32);R[G['mask']]=rho/G['dx']**3
        out.append(B.skyfrom(G,B.accel(G,R)))
    return out[0]+out[1],nh
if __name__=='__main__':
    ck=OUT+'mc_sky.npy';S=list(np.load(ck)) if os.path.exists(ck) else []
    for i in range(len(S),NU):
        t=time.time();s,nh=universe(500000+i);S.append(s.astype(np.float32))
        if i%5==4 or i==NU-1:np.save(ck,np.array(S))
        print(i,'%.1fs'%(time.time()-t),nh,flush=True)
    np.save(ck,np.array(S))
