# Line-of-sight CMB for CoE red universe + MC web (same seeds as mc/mc.py). Components saved separately so strength E can be applied afterwards.
import os,sys,time,numpy as np,healpy as hp,scipy.fft as F
from scipy.ndimage import map_coordinates
sys.path.insert(0,'/workspace/redweb2/realweb');sys.path.insert(0,'/workspace/redweb2/realweb/mc')
import mc as M;import build as B;from cosmo import rho_m
OUT='/workspace/redweb2/realweb/los/';NS=32;NPIX=hp.nside2npix(NS);V=np.array(hp.pix2vec(NS,np.arange(NPIX))).T
Gc=4.30e-12;c=2.998e5;H0=70e3          # G in Gpc (km/s)^2/Msun, c km/s, H0 km/s/Gpc
XS=14.0;zrec=1090.;k_=(1-1/np.sqrt(1+zrec))
chi=np.arange(0.5,XS,0.025)+0.0125;dchi=0.025
a_of=lambda x:(1-x/XS*k_)**2              # EdS-shaped a(chi) rescaled so chi=14 Gpc <-> z=1090 (assumption)
aa=a_of(chi);HH=H0*aa**-1.5
def grids(seed):
    # exact replay of mc.universe density fields (same rng order), returns per-grid rho (Msun/Gpc^3, mean removed in mask)
    rng=np.random.default_rng(seed);out=[]
    for G in M.Gs:
        N=G['N'];w=rng.standard_normal((N,N,N),dtype=np.float32);d=F.irfftn(F.rfftn(w,workers=8)*G['amp'],s=w.shape,workers=8)[G['mask']];del w
        dV=(G['dx']*1000)**3;rho=np.zeros(int(G['mask'].sum()))
        for n,Mm,b in M.bins:cc=rng.poisson(n*dV*np.maximum(0,1+b*d)*G['mod']);rho+=cc*Mm
        rho-=rho.mean();R=np.zeros((N,N,N),np.float32);R[G['mask']]=rho/G['dx']**3;out.append(R)
    return out
Gf,Gq=M.Gs
xf=(np.arange(Gf['N'])+0.5)*Gf['dx']-10;
def tocoarse(Rf):  # deposit fine-grid mass onto coarse grid (NGP)
    N=Gq['N'];ii=np.nonzero(Rf);P=np.stack([xf[i] for i in ii],1);m=Rf[ii]*Gf['dx']**3
    j=np.floor((P+30)/Gq['dx']).astype(int);out=np.zeros(N**3);np.add.at(out,(j[:,0]*N+j[:,1])*N+j[:,2],m);return out.reshape(N,N,N)/Gq['dx']**3
def phi(G,rho,smooth):
    rk=F.rfftn(rho.astype(np.float32),workers=8)*np.exp(-G['K2']*smooth**2/2)
    with np.errstate(divide='ignore',invalid='ignore'):pk=np.where(G['K2']>0,-4*np.pi*Gc*rk/G['K2'],0)
    return [F.irfftn(q,s=rho.shape,workers=8).astype(np.float32) for q in (pk,1j*G['KX']*pk,1j*G['KY']*pk,1j*G['KZ']*pk)]  # Phi,(km/s)^2 and grad /Gpc
def samp(G,fs,P):
    u=((P-G['x0'])/G['dx']-0.5).T;return np.stack([map_coordinates(f,u,order=1,mode='wrap') for f in fs],-1)
# CoE background potential: rho_bar*(mod-1) on coarse grid, whole box (red density shape about PC at 42 Gpc)
xq=(np.arange(Gq['N'])+0.5)*Gq['dx']-30;X,Y,Z=np.meshgrid(xq,xq,xq,indexing='ij',sparse=True)
modq=B.shp(np.sqrt((X-B.PC[0])**2+(Y-B.PC[1])**2+(Z-B.PC[2])**2)/42);rr=np.sqrt(X**2+Y**2+Z**2)
rhoC=(rho_m*1e9*(modq/modq[rr<15].mean()-1)*(rr<15)).astype(np.float32);del X,Y,Z,modq,rr
PC_=phi(Gq,rhoC,Gq['dx'])
def path(fsq,fsf,G0):
    """returns per pixel: Phi at source, RS integrals (web.web, web.coe, coe.coe), lensing deflection (web, coe) as 3-vectors"""
    nC=len(chi);out={k:np.zeros(NPIX) for k in ('rs_ww','rs_wc','rs_cc')};aw=np.zeros((NPIX,3));ac=np.zeros((NPIX,3))
    for i0 in range(0,nC,40):
        cs=chi[i0:i0+40];P=(V[:,None,:]*cs[None,:,None]).reshape(-1,3)
        w=samp(Gq,fsq,P)
        inf=np.linalg.norm(P,axis=1)<9.0
        if inf.any():w[inf]+=samp(Gf,fsf,P[inf])
        cC=samp(Gq,PC_,P);w=w.reshape(NPIX,len(cs),4);cC=cC.reshape(NPIX,len(cs),4)
        gw,gc=w[...,1:],cC[...,1:];h=HH[i0:i0+40];a_=aa[i0:i0+40]
        K=(2/c**3)*(2/(3*a_*h))*dchi        # dT/T = (2/c^3) int (2/(3aH)) |grad Phi|^2 dchi  (moving-lens / 2nd-order Phi-dot, EdS f=1)
        out['rs_ww']+=( (gw**2).sum(-1)*K).sum(1);out['rs_wc']+=(2*(gw*gc).sum(-1)*K).sum(1);out['rs_cc']+=((gc**2).sum(-1)*K).sum(1)
        L=(-2/c**2*(XS-cs)/XS*dchi)[None,:,None]   # alpha = -2 int (chi_s-chi)/chi_s grad_perp Phi/c^2 dchi
        for g,acc in ((gw,aw),(gc,ac)):
            gp=g-(g*V[:,None,:]).sum(-1,keepdims=True)*V[:,None,:];acc+=(gp*L).sum(1)
    Ps=XS*V;sw=samp(Gq,fsq,Ps)[:,0]/c**2/3   # SW: dT/T = Phi_src/3c^2 (web only; red CoE SW already in sm_red)
    out.update(sw=sw,aw=aw,ac=ac);return out
if __name__=='__main__':
    ids=np.random.default_rng(7).choice(1000,int(sys.argv[1]) if len(sys.argv)>1 else 60,replace=False)
    os.makedirs(OUT+'comp',exist_ok=True)
    for k in ids:
        f=OUT+f'comp/u{k}.npz'
        if os.path.exists(f):continue
        t=time.time();Rf,Rq=grids(500000+int(k))
        Rq=Rq+tocoarse(Rf)
        fq=phi(Gq,Rq,Gq['dx']);ff_=phi(Gf,Rf,Gf['dx']);fs=phi(Gf,Rf,Gq['dx']);ff=[a-b for a,b in zip(ff_,fs)];del ff_,fs
        o=path(fq,ff,None);np.savez(f,**o);print(k,'%.1fs'%(time.time()-t),{q:float(np.std(v)) for q,v in o.items()},flush=True)
