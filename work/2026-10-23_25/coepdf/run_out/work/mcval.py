# Monte Carlo validation: Gaussian potential field x W(|r|) on the fine grid, pushed through los_flow_R comp() at R=4 Gpc (grow, full-velocity moving lens).
import sys,json,time,numpy as np,healpy as hp,scipy.fft as Fft
sys.argv=['x'];exec(open('/workspace/redweb2/realweb/los_flow_R/flowR.py').read().split('\nif __name__')[0])
sys.path.insert(0,'/workspace/redweb2/stochfield');from fwd import Gc as G_,rho as RHO
NSV=int(sys.argv_[1]) if False else 64
V=np.array(hp.pix2vec(NSV,np.arange(hp.nside2npix(NSV)))).T
p=json.load(open('/workspace/redweb2/stochfield/mc_params.json'));nreal=int(p.pop('nreal'))
G=Gf;N=G['N'];K=np.sqrt(G['K2']);xs=(np.arange(N)+0.5)*G['dx']+G['x0'];X,Y,Z=np.meshgrid(xs,xs,xs,indexing='ij',sparse=True);r=np.sqrt(X**2+Y**2+Z**2)
Wr=np.exp(-0.5*((r-p['chiw'])/p['sw'])**2).astype(np.float32);del X,Y,Z,r
with np.errstate(divide='ignore',invalid='ignore'):
    Pphi=np.where(K>0,(4*np.pi*G_*RHO)**2*10**p['lgA']*K**p['n']*np.exp(-(K/p['kc'])**2)/K**4,0)
Lbox=N*G['dx'];amp=np.sqrt(Pphi*N**3/Lbox**3).astype(np.float32);del Pphi
out=[]
for i in range(nreal):
    t=time.time();rng=np.random.default_rng(9000+i);w=rng.standard_normal((N,N,N),dtype=np.float32)
    ph=Fft.irfftn(Fft.rfftn(w,workers=8)*amp,s=w.shape,workers=8).astype(np.float32)*Wr;del w
    pk=Fft.rfftn(ph,workers=8);ff=[ph]+[Fft.irfftn(1j*G[q]*pk,s=ph.shape,workers=8).astype(np.float32) for q in ('KX','KY','KZ')];del pk
    fq=[np.zeros((Gq['N'],)*3,np.float32) for _ in range(4)]
    o=comp(fq,ff,4.0);del ff
    sc=o['sw_grow']+o['isw_grow'];ml=o['mlf_grow']
    cl=lambda m:hp.anafast(m,lmax=100)
    ll=np.arange(101);f=ll*(ll+1)/(2*np.pi)*2.7255e6**2
    out.append(np.stack([cl(sc)*f,cl(ml)*f,cl(sc+ml)*f]));print(i,'%.0fs'%(time.time()-t),out[-1][2][[2,10,30]],flush=True)
    np.save('/workspace/redweb2/stochfield/mc_val.npy',np.array(out))
    if i==0:np.save('/workspace/redweb2/stochfield/mc_sky0.npy',sc+ml)
