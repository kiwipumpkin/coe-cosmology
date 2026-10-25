import os,sys,json,numpy as np,healpy as hp
from scipy.interpolate import CubicSpline
A='/workspace/coe-cosmology-gh/draft/current/analyses'
g={'__file__':A+'/shell_ensemble/ensemble.py','__name__':'imp'};os.chdir(A+'/shell_ensemble')
src=open(g['__file__']).read();exec(src[:src.index('# ---------- beta=0')],g)
field3,th,KG,RV,FLOOR,phi=[g[k] for k in('field3','th','KG','RV','FLOOR','phi')]
sys.path.insert(0,A+'/void_mass_network');sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import vmn_gen as G;MSUN=1.989e30
ELL=np.radians(phi+180);Ch=np.array([np.cos(ELL),np.sin(ELL)]);L_G=922.65
FLOC=-0.9/KG*(RV/L_G)*(np.cos(th)*Ch[0]+np.sin(th)*Ch[1]);FLOC-=FLOC.mean()
web=G.generate(Mmin=3e15,seed=2000,comp_mode='lcdm');X,M=G.kernel_bodies(web)
pp=web['params'];vv=np.linalg.norm(web['S'],axis=1)<pp['R_out'];nv=int(vv.sum());nc=len(web['MC'])
gid=np.concatenate([np.arange(nv),np.arange(nv),nv+np.arange(nc)])  # void id shared by its -D and +W body; clusters own ids
far=np.linalg.norm(np.vstack([web['cen'][vv],web['cen'][vv],web['XC']]),axis=1)>pp['gap_us'];gid=gid[far];assert len(gid)==len(M)
import fit as _f0
CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937);aC=np.array([1,0,np.tan(np.radians(15))]);aC/=np.linalg.norm(aC);PC=42.0*aC
shp=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
PU=np.vstack([web['cen'][vv],web['XC']])/1000.0;xu=np.linalg.norm(PU-PC,axis=1)/42.0;pu=shp(xu)/shp(xu).max()
ku=np.random.default_rng(777).random(len(pu))<pu;np.savez('/workspace/redweb2/thin4_keep.npz',ku=ku,pu=pu,xu=xu,nv=nv)
dk=ku[gid];X,M=X[dk],M[dk];gid=gid[dk];print('density keep',ku[:nv].mean(),ku[nv:].mean(),pu.min(),pu.max(),flush=True)
XF,MF=X.copy(),M.copy()
P=M@field3(*[np.arctan2((X/RV)[:,1],(X/RV)[:,0]),np.arcsin((X[:,2]/np.linalg.norm(X,axis=1))),np.linalg.norm(X,axis=1)/RV])/KG if False else None
P=np.zeros(360)
for i in range(0,len(M),20000):
  Y=X[i:i+20000]/RV;D=np.linalg.norm(Y,axis=1);P+=M[i:i+20000]@field3(np.arctan2(Y[:,1],Y[:,0]),np.arcsin(Y[:,2]/D),D)/KG
P-=M.sum()*FLOC;ver=dict(n_bodies=len(M),pattern_rms_over_floor=float(np.sqrt(np.mean(P**2))/FLOOR),repo_value=149.87782489337775,repo_n_bodies=272175)
print(ver,flush=True)
# full-sky generalisation of field3: ring point at rho*n (RV units) for every sky direction n
CON='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained'
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
ns=32;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T
def sky(rho):
  Y=X/RV;out=np.zeros(len(V));D3=np.linalg.norm(Y,axis=1)**3
  for i in range(0,len(M),4000):
    y=Y[i:i+4000];d=y[:,None,:]-rho[None,:,None]*V[None];r3=np.linalg.norm(d,axis=2)**3
    ar=np.einsum('bpk,pk->bp',d,V)/r3-(y@V.T)/D3[i:i+4000,None];out+=(M[i:i+4000]/KG)@ar
  return out-out.mean()
# ring check of full-sky kernel in the plane b=0
o=np.load(CON+'/runs_poly/f20000_p_-0.041_-0.035_0.1353_0.npz')['o'];o=o[np.argsort(o[:,0])]
sys.path.insert(0,B+'/../../rotfit');os.chdir(B);import fit
aG=fit.vec(0,15);psi=np.degrees(np.arccos(np.clip(V@aG,-1,1)))
sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));sm=sm/sm.mean()-1
wf=CubicSpline(o[:,0],o[:,4]);pf=CubicSpline(o[:,0],o[:,5])
t=V-(V@aG)[:,None]*aG;tn=np.linalg.norm(t,axis=1);t=np.where(tn[:,None]>1e-9,t/np.maximum(tn,1e-12)[:,None],0)
w=wf(psi);ph=pf(psi);E=42.0*(-aG[None,:]+w[:,None]*(np.cos(ph)[:,None]*aG[None,:]+np.sin(ph)[:,None]*t));dist=np.linalg.norm(E,axis=1)
T0=2.7255e6;l=np.arange(31)
def Dl(x):return (l*(l+1)*hp.anafast(np.asarray(hp.remove_dipole(x)),lmax=30)*T0**2/(2*np.pi)).tolist()
Pl=np.loadtxt(B+'/planckTT.txt')
def chi(D):
  c=0
  for r in Pl:
    L=int(r[0])
    if 2<=L<=30:e=r[3] if D[L]>r[1] else r[2];c+=((D[L]-r[1])/e)**2
  return float(c)
res=dict(geom=dict(RV=float(RV),Xmin=[float(v) for v in XF.min(0)],Xmax=[float(v) for v in XF.max(0)],rnorm_minmax=[float(np.linalg.norm(XF,axis=1).min()),float(np.linalg.norm(XF,axis=1).max())],R_out=float(pp['R_out']),gap_us=float(pp['gap_us'])),verify=ver,n_voids=nv,n_clusters=nc,red_surface_Gpc=[float(dist.min()),float(np.median(dist)),float(dist.max())])
res['smooth']=dict(Dl=Dl(sm),chi2=chi(Dl(sm)))
rho=np.full(len(V),4.0/RV);X,M=XF,MF;wfull=sky(rho);res['full_web_only_rms']=float(wfull.std())
np.save("/workspace/redweb2/wfull4.npy",wfull);np.save("/workspace/redweb2/sm_red.npy",sm)
for s_ in(.01,.0175,.02,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0):d=Dl(sm+s_*wfull);res[f'scaled_{s_}']=dict(Dl=d,chi2=chi(d))
print({k:v['chi2'] for k,v in res.items() if k.startswith('scaled')},flush=True)
rng=np.random.default_rng(20261021);FS=[.01,.02,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0];res['thin']={str(f):[] for f in FS}
for k in range(5):
  uv=rng.random(nv);uc=rng.random(nc);u=np.concatenate([uv,uc])[gid];parts=[]
  for j in range(len(FS)):
    lo=0 if j==0 else FS[j-1];sel=(u>=lo)&(u<FS[j]) if j<len(FS)-1 else (u>=lo);X,M=XF[sel],MF[sel];parts.append(sky(rho))
  for j,f in enumerate(FS):
    wm=np.sum(parts[:j+1],axis=0);d=Dl(sm+wm)
    res['thin'][str(f)].append(dict(Dl=d,chi2=chi(d),n_voids=int((uv<f).sum()),n_clusters=int((uc<f).sum()),net_mass_Msun=float(MF[u<f].sum()/MSUN),rms=float(wm.std()),seed_draw=k))
    print(f,k,res['thin'][str(f)][-1]['chi2'],flush=True)
  json.dump(res,open('/workspace/redweb2/thin4_results.json','w'))
