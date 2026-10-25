import os,sys,numpy as np,healpy as hp
from scipy.interpolate import CubicSpline
A='/workspace/coe-cosmology-gh/draft/current/analyses'
_cwd=os.getcwd();g={'__file__':A+'/shell_ensemble/ensemble.py','__name__':'imp'};os.chdir(A+'/shell_ensemble')
src=open(g['__file__']).read();exec(src[:src.index('# ---------- beta=0')],g);os.chdir(_cwd)
KG,RV=g['KG'],g['RV']
T0=2.7255e6;LM=30;ns=32;l=np.arange(LM+1)
DG=42.0  # CoE distance in Gpc as in redweb2/run.py
def mtot_web():
    f=os.path.join(os.path.dirname(__file__),'mtot_web.npy')
    if os.path.exists(f):return float(np.load(f))
    sys.path.insert(0,A+'/void_mass_network');import vmn_gen as G
    web=G.generate(Mmin=3e15,seed=1000,comp_mode='lcdm');X,M=G.kernel_bodies(web);np.save(f,M.sum());return float(M.sum())
# pole frame: z = anti-CoE axis aG; CoE at -DG*z
V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;PSI=np.degrees(np.arccos(np.clip(V[:,2],-1,1)))
def terms(o,C=None,shape=None,DRV=2.937,nr=160,nmu=96,nph=48):
    o=o[np.argsort(o[:,0])]
    sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(PSI));sm=sm/sm.mean()-1
    # emission distance per psi (meridian, 1 deg)
    ps=np.radians(np.arange(0,180.01,1.0));wf=CubicSpline(o[:,0],o[:,4]);pf=CubicSpline(o[:,0],o[:,5])
    z=np.array([0,0,1.]);n=np.stack([np.sin(ps),0*ps,np.cos(ps)],1);t=np.stack([np.cos(ps),0*ps,-np.sin(ps)],1)
    w=wf(np.degrees(ps));ph=pf(np.degrees(ps))
    E=DG*(-z[None]+w[:,None]*(np.cos(ph)[:,None]*z[None]+np.sin(ph)[:,None]*t));dist=np.linalg.norm(E,axis=1)
    # smooth spherical mass about CoE, density shape(x) truncated at Rb (famP convention)
    Rb=1.5*(1+1/DRV);x=(np.arange(nr)+.5)/nr*Rb;mu=(np.arange(nmu)+.5)/nmu*2-1;pa=(np.arange(nph)+.5)/nph*2*np.pi
    sx=shape(x)*x**2;Xg,Mg,Pg=np.meshgrid(x,mu,pa,indexing='ij');Wt=np.broadcast_to(sx[:,None,None],Xg.shape)
    sn=np.sqrt(1-Mg**2);P=DG*np.stack([Xg*sn*np.cos(Pg),Xg*sn*np.sin(Pg),Xg*Mg-1],-1).reshape(-1,3);Wt=Wt.reshape(-1).copy()
    keep=np.linalg.norm(P,axis=1)>1.0;P=P[keep];Wt=Wt[keep];Mb=Wt/Wt.sum()*mtot_web()
    Y=P/RV;D3=np.linalg.norm(Y,axis=1)**3;rho=dist/RV;out=np.zeros(len(ps))
    for i in range(0,len(Mb),20000):
        y=Y[i:i+20000];d=y[:,None,:]-rho[None,:,None]*n[None];r3=np.linalg.norm(d,axis=2)**3
        ar=np.einsum('bpk,pk->bp',d,n)/r3-(y@n.T)/D3[i:i+20000,None];out+=(Mb[i:i+20000]/KG)@ar
    gm=CubicSpline(np.degrees(ps),out)(PSI);gm=gm-gm.mean()
    return sm,gm,dist
def al0(mp):
    a=hp.map2alm(mp,lmax=LM,iter=3);return np.array([a[hp.Alm.getidx(LM,L,0)].real for L in l])*T0
def Dl_al0(a):return l*(l+1)*a**2/(2*l+1)/(2*np.pi)
def Dl_map(mp):return l*(l+1)*hp.anafast(np.asarray(hp.remove_dipole(mp)),lmax=LM)*T0**2/(2*np.pi)
def shape_poly(C,DRV):
    Rb=1.5*(1+1/DRV);return lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(C)))
def shape_slc(SL=0.041037,p=0.035):return lambda x:np.abs(x)**-SL*np.exp(-p*(np.abs(x)-1)**2)
