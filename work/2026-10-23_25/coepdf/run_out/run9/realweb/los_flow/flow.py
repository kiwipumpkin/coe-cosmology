# Light-path CMB with the web riding the red LTB (CoE) flow. Lagrangian (today-label) coordinates; photon paths from the red model's own geodesics.
import os,sys,time,numpy as np,healpy as hp,scipy.fft as F
from scipy.ndimage import map_coordinates
sys.argv=[sys.argv[0]];sys.path.insert(0,'/workspace/redweb2/realweb/los')
exec(open('/workspace/redweb2/realweb/los/los.py').read().split("if __name__=='__main__':")[0])   # grids(), tocoarse(), phi(), samp(), Gf,Gq, V, c, Gc
O='/workspace/redweb2/realweb/los_flow/';L=np.load(O+'ltb.npz');Dz=np.load(O+'deform.npz')
lT=np.log(L['T']);Ag=L['Agrow'];Af=L['Afix']
Afun={'grow':lambda t:np.interp(np.log(t),lT,Ag),'fixed':lambda t:np.interp(np.log(t),lT,Af)}
# deformation weight w(t): differential stretch eps(t,x)=G(t,x)/x-1, normalized to 1 at t_rec (shape-preserving approx), using x=0.65
xs,Ts,GG=Dz['xs'],Dz['Ts'],Dz['G'];ept=np.interp(0.65,xs,GG[0]/xs)-1;wv=np.array([np.interp(0.65,xs,g/xs)-1 for g in GG])/ept
wfun=lambda t:np.interp(np.log(t),np.log(Ts),wv)
u0=-B.aC/np.linalg.norm(B.aC)   # outward radial from CoE at our position (CoE is at +42 aC from us)
cps=np.clip(V@u0,-1,1);psi=np.degrees(np.arccos(cps));e=V-cps[:,None]*u0;en=np.linalg.norm(e,axis=1);e=np.where(en[:,None]>1e-9,e/np.maximum(en,1e-12)[:,None],np.cross(u0,[0,0,1])/np.linalg.norm(np.cross(u0,[0,0,1])))
PS=L['PS'];PA=L['paths']   # (37,300,4): lambda,t,x,phi
j=np.minimum((psi/5).astype(int),35);fr=(psi-PS[j])/5
pth=PA[j]*(1-fr)[:,None,None]+PA[j+1]*fr[:,None,None]   # (npix,300,4)
lam,tt,xx,ph=[pth[...,k] for k in range(4)]
XP=42*(xx[...,None]*(np.cos(ph)[...,None]*u0+np.sin(ph)[...,None]*e[:,None,:])-u0)   # Lagrangian position rel. us, Gpc
XP=XP[:,::-1];tt=tt[:,::-1];lam=lam[:,::-1]   # order source -> us (forward in time)
Xm=0.5*(XP[:,1:]+XP[:,:-1]);tm=0.5*(tt[:,1:]+tt[:,:-1]);dX=np.linalg.norm(XP[:,1:]-XP[:,:-1],axis=-1);lev=(lam[:,1:]+lam[:,:-1])/2/lam[:,:1]  # lever (1-l/ls): lam=0 at us
print('max |X| on light cone (Gpc):',float(np.linalg.norm(XP,axis=-1).max()),flush=True)
def dphi_deform(rho):
    # Phi change from displacing coarse-grid mass by psi(X)=(P-PC)(g(x)/x-1) at t_rec: d rho=-div(rho psi)
    N=Gq['N'];x=(np.arange(N)+0.5)*Gq['dx']-30;out=None
    Xg,Yg,Zg=np.meshgrid(x,x,x,indexing='ij',sparse=True);rel=[Xg-B.PC[0],Yg-B.PC[1],Zg-B.PC[2]];xx_=np.sqrt(rel[0]**2+rel[1]**2+rel[2]**2)/42
    s=(np.interp(xx_,xs,GG[0]/xs)-np.interp(1.0,xs,GG[0]/xs)).astype(np.float32)   # relative to observer shell
    dk=0
    for Kc,r_ in zip((Gq['KX'],Gq['KY'],Gq['KZ']),rel):dk=dk+1j*Kc*F.rfftn((rho*r_*s).astype(np.float32),workers=8)
    rk=-dk*np.exp(-Gq['K2']*Gq['dx']**2/2)
    with np.errstate(divide='ignore',invalid='ignore'):pk=np.where(Gq['K2']>0,-4*np.pi*Gc*rk/Gq['K2'],0)
    return [F.irfftn(q,s=rho.shape,workers=8).astype(np.float32) for q in (pk,1j*Gq['KX']*pk,1j*Gq['KY']*pk,1j*Gq['KZ']*pk)]
def field(fq,ff,P):
    w=samp(Gq,fq,P);inf=np.linalg.norm(P,axis=1)<9.0
    if inf.any():w[inf]+=samp(Gf,ff,P[inf])
    return w
if __name__=='__main__':
    ids=sorted(int(f[1:-4]) for f in os.listdir('/workspace/redweb2/realweb/los/comp'));os.makedirs(O+'comp',exist_ok=True)
    for k in ids:
        f=O+f'comp/u{k}.npz'
        if os.path.exists(f):continue
        t=time.time();Rf,Rq=grids(500000+k);Rq=Rq+tocoarse(Rf)
        fq=phi(Gq,Rq,Gq['dx']);ff_=phi(Gf,Rf,Gf['dx']);fs=phi(Gf,Rf,Gq['dx']);ff=[a-b for a,b in zip(ff_,fs)];del ff_,fs
        dd=dphi_deform(Rq)
        n=XP.shape[0];P0=field(fq,ff,Xm.reshape(-1,3)).reshape(n,-1,4);Pd=samp(Gq,dd,Xm.reshape(-1,3)).reshape(n,-1,4)
        S0=field(fq,ff,XP[:,0]);Sd=samp(Gq,dd,XP[:,0]);out={}
        g=P0[...,1:]+0*Pd[...,1:];gp=g-(g*V[:,None,:]).sum(-1,keepdims=True)*V[:,None,:]
        for nm,A in Afun.items():
            Ai=A(tt);wi=wfun(tt);Am=A(tm);wm=wfun(tm)
            out['sw_'+nm]=Ai[:,0]*(S0[:,0]+wi[:,0]*Sd[:,0])/(3*c**2)
            dAt=np.diff(Ai,axis=1);dw=np.diff(wi,axis=1)
            out['isw_'+nm]=(2/c**2)*(dAt*(P0[...,0]+wm*Pd[...,0])).sum(1)          # from dilution/growth of web potential (A-dot)
            out['iswdef_'+nm]=(2/c**2)*(Am*dw*Pd[...,0]).sum(1)                   # from anisotropic CoE stretching of the web
            out['al_'+nm]=(-2/c**2*(Am*lev*dX)[...,None]*gp).sum(1)
        np.savez(f,**out);print(k,'%.1fs'%(time.time()-t),{q:float(np.std(v)) for q,v in out.items()},flush=True)
