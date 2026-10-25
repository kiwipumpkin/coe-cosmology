import os,sys,json,numpy as np,healpy as hp
from scipy.interpolate import CubicSpline
A='/workspace/coe-cosmology-gh/draft/current/analyses'
g={'__file__':A+'/shell_ensemble/ensemble.py','__name__':'imp'};os.chdir(A+'/shell_ensemble')
src=open(g['__file__']).read();exec(src[:src.index('# ---------- beta=0')],g)
field3,th,KG,RV,FLOOR,phi=[g[k] for k in('field3','th','KG','RV','FLOOR','phi')]
sys.path.insert(0,A+'/void_mass_network');import vmn_gen as G;MSUN=1.989e30
ELL=np.radians(phi+180);Ch=np.array([np.cos(ELL),np.sin(ELL)]);L_G=922.65
FLOC=-0.9/KG*(RV/L_G)*(np.cos(th)*Ch[0]+np.sin(th)*Ch[1]);FLOC-=FLOC.mean()
web=G.generate(Mmin=3e15,seed=2000,comp_mode='lcdm');X,M=G.kernel_bodies(web)
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
res=dict(verify=ver,red_surface_Gpc=[float(dist.min()),float(np.median(dist)),float(dist.max())])
res['smooth']=dict(Dl=Dl(sm),chi2=chi(Dl(sm)))
for tag,rho in [('R0.9',np.full(len(V),0.9/RV))]:
  wm=sky(rho);res[tag+'_web_only']=dict(Dl=Dl(wm),rms=float(wm.std()))
  for s in(.0175,):d=Dl(sm+s*wm);res[f'{tag}_{s}']=dict(Dl=d,chi2=chi(d))
  print(tag,{k:v.get('chi2') for k,v in res.items() if k.startswith(tag) and isinstance(v,dict)},flush=True)
json.dump(res,open(CON+'/red_web_s2000_09_175_results.json','w'),indent=1)
