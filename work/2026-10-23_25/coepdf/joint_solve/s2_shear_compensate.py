# Step 1 (tidal shear): place random voids/clusters inside the R=4 Gpc CMB sphere and solve their masses so that the
# radial tidal-shear sky (same kernel/units as realweb/build.py skyfrom: MSUN/KG*RV^2 * n.[d/|d|^3 - y/|y|^3], mean removed)
# cancels the deterministic CoE background (hl_R4_ns64: sm + rs_cc) for l=1..LC. Bounded least squares; seeds fixed.
import sys,json,numpy as np,healpy as hp
from scipy.optimize import lsq_linear
MSUN=1.989e30;KG=7.7e51/0.0039;RV=14.3;R=4.0;NS=64;T0=2.7255e6
seed=int(sys.argv[1]) if len(sys.argv)>1 else 2000;NB=int(sys.argv[2]) if len(sys.argv)>2 else 3000;LC=int(sys.argv[3]) if len(sys.argv)>3 else 12
MMAX=float(sys.argv[4]) if len(sys.argv)>4 else 1e17;RMAX=float(sys.argv[5]) if len(sys.argv)>5 else 3.0;EPS=0.1
z=np.load('/workspace/redweb2/realweb/los_flow_R/hl_R4_ns64.npz');B=z['sm']+z['rs_cc']
n=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T
rng=np.random.default_rng(seed);u=rng.normal(size=(NB,3));u/=np.linalg.norm(u,axis=1)[:,None]
r=(rng.uniform(0.3**3,RMAX**3,NB))**(1/3);y=u*r[:,None]
def kern(y):
    d=y[None,:]-R*n;dd=(np.sum(d*d,1)+EPS**2)**1.5;k=np.sum(n*d,1)/dd-n@y/np.linalg.norm(y)**3;return (k-k.mean())*MSUN/KG*RV**2
K=np.array([kern(v) for v in y]).T          # npix x NB, per Msun
LMX=3*NS-1
def alm(m):return hp.map2alm(m,lmax=LMX,iter=1)
ll,mm=hp.Alm.getlm(LMX);sel=(ll>=1)&(ll<=LC);w=np.where(mm[sel]==0,1.,np.sqrt(2.))
def vec(a):a=a[sel]*w;return np.r_[a.real,a.imag]
Ab=vec(alm(B));Ak=np.array([vec(alm(K[:,i])) for i in range(NB)]).T
sc=1e16;lam=float(sys.argv[6]) if len(sys.argv)>6 else 1e-3
# also penalise power the bodies add at l=LC+1..40 (keeps TT clean)
sel2=(ll>LC)&(ll<=40);w2=np.where(mm[sel2]==0,1.,np.sqrt(2.))
def vec2(a):a=a[sel2]*w2;return np.r_[a.real,a.imag]
Hb=vec2(alm(B));Hk=np.array([vec2(alm(K[:,i])) for i in range(NB)]).T
Amat=np.vstack([Ak*sc,0.3*Hk*sc,np.sqrt(lam)*np.abs(Ab).max()*np.eye(NB)]);bvec=np.r_[-Ab,-0.3*Hb,np.zeros(NB)]
res=lsq_linear(Amat,bvec,bounds=(-MMAX/sc,MMAX/sc),lsmr_tol='auto',max_iter=5000)
M=res.x*sc;C=K@M;Tot=B+C
def amp(m,L=40):cl=hp.anafast(m,lmax=L);return np.sqrt(cl*(2*np.arange(L+1)+1)/(4*np.pi))   # rms dT/T per l
out=dict(seed=seed,NB=NB,LC=LC,MMAX=MMAX,RMAX=RMAX,R=R,eps=EPS,lam=lam,
 n_clusters=int((M>0.01*MMAX).sum()),n_voids=int((M<-0.01*MMAX).sum()),M_pos=float(M[M>0].sum()),M_neg=float(M[M<0].sum()),
 rms_bg=float(B.std()),rms_tot=float(Tot.std()),max_abs_tot=float(np.abs(Tot-Tot.mean()).max()),
 dipole_bg=float(np.linalg.norm(hp.fit_dipole(B)[1])),dipole_tot=float(np.linalg.norm(hp.fit_dipole(Tot)[1])),
 amp_bg=amp(B).tolist(),amp_tot=amp(Tot).tolist())
tag=f's{seed}_N{NB}_L{LC}_M{MMAX:.0e}_r{RMAX}_lam{lam:.0e}'
np.savez(f's2_{tag}.npz',y=y,M=M,C=C,Tot=Tot,B=B);json.dump(out,open(f's2_{tag}.json','w'),indent=1)
print(tag,{k:v for k,v in out.items() if not k.startswith('amp')});print('amp_tot l1-6',np.array(out['amp_tot'][1:7]))
