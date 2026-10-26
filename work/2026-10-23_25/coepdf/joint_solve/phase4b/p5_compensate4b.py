# Phase4b (NEW) copy of phase4 script, outputs to phase4b.
# WQ: extra weight on l=2 rows (quadrupole-targeted); WH: weight on l=LC+1..40 penalty; cancel_err_* = |bg+bodies+field-l2-target|.
# Phase3 (NEW) copy of phase2 p5, output to phase3.
# Phase2 A4/B2: tidal-shear mass compensation at a chosen CMB-shell radius R, on the rebuilt background map, optionally ALSO cancelling
# a fraction (1-f) of the stochastic field's own l=2 a_lm on an ACTUAL seed-42 sky (a_lm drawn with rng(42), per-l power set to the
# seed-42 one-sky D_l used throughout: Df * chi2_{2l+1}/(2l+1) with rng(42)).  Same kernel/units as s2 (build.py skyfrom), geometry scaled with R:
# bodies at r in [0.075R, RMAXF*R], softening 0.025R.  Bounded least squares; l=1..LC targeted, l=LC+1..40 penalised (weight 0.3).
import sys,os,json,argparse,numpy as np,healpy as hp
from scipy.optimize import lsq_linear
ap=argparse.ArgumentParser();ap.add_argument('--bg');ap.add_argument('--R',type=float);ap.add_argument('--seed',type=int,default=2000);ap.add_argument('--NB',type=int,default=3000)
ap.add_argument('--LC',type=int,default=12);ap.add_argument('--MMAX',type=float);ap.add_argument('--RMAXF',type=float,default=0.875);ap.add_argument('--lam',type=float,default=1e-5)
ap.add_argument('--field',default='');ap.add_argument('--f',type=float,default=1.0);ap.add_argument('--tag',default='');ap.add_argument('--WQ',type=float,default=1.0);ap.add_argument('--WH',type=float,default=0.3)
a=ap.parse_args()
MSUN=1.989e30;KG=7.7e51/0.0039;RV=14.3;NS=64;T0=2.7255e6;R=a.R;EPS=0.025*R;O='/workspace/coepdf/joint_solve/phase4b/'
B=np.load(a.bg);n=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T;LMX=3*NS-1;ll,mm=hp.Alm.getlm(LMX)
Fmap=np.zeros_like(B);Ftarget=np.zeros_like(B)
if a.field:   # actual seed-42 field sky at l<=LMX (Delta T/T)
    Df=np.load(a.field)['Df'];L=np.arange(len(Df));g=np.random.default_rng(42).chisquare(2*np.maximum(L,1)+1)/(2*np.maximum(L,1)+1)
    Cl=np.zeros(LMX+1);l_=np.arange(2,LMX+1);Cl[2:]=Df[2:LMX+1]*g[2:LMX+1]*2*np.pi/(l_*(l_+1))/T0**2
    r=np.random.default_rng(42);al=(r.normal(size=len(ll))+1j*r.normal(size=len(ll)))/np.sqrt(2);al[mm==0]=al[mm==0].real*np.sqrt(2)
    for L_ in range(LMX+1):   # set each l's realised power exactly to the seed-42 one-sky value
        s=ll==L_;p=(np.abs(al[s][0])**2+2*np.sum(np.abs(al[s][1:])**2))/(2*L_+1)
        al[s]*=np.sqrt(Cl[L_]/p) if p>0 else 0
    Fmap=hp.alm2map(al,NS,lmax=LMX);a2=al*(ll==2);Ftarget=(1-a.f)*hp.alm2map(a2,NS,lmax=LMX)   # part of field l=2 to cancel
rng=np.random.default_rng(a.seed);u=rng.normal(size=(a.NB,3));u/=np.linalg.norm(u,axis=1)[:,None]
r=(rng.uniform((0.075*R)**3,(a.RMAXF*R)**3,a.NB))**(1/3);y=u*r[:,None]
def kern(y):
    d=y[None,:]-R*n;dd=(np.sum(d*d,1)+EPS**2)**1.5;k=np.sum(n*d,1)/dd-n@y/np.linalg.norm(y)**3;return (k-k.mean())*MSUN/KG*RV**2
K=np.array([kern(v) for v in y]).T
def alm(m):return hp.map2alm(m,lmax=LMX,iter=1)
sel=(ll>=1)&(ll<=a.LC);w=np.where(mm[sel]==0,1.,np.sqrt(2.))*np.where(ll[sel]==2,a.WQ,1.);sel2=(ll>a.LC)&(ll<=40);w2=np.where(mm[sel2]==0,1.,np.sqrt(2.))
def v1(x):x=x[sel]*w;return np.r_[x.real,x.imag]
def v2(x):x=x[sel2]*w2;return np.r_[x.real,x.imag]
Tg=B+Ftarget;At=alm(Tg);Ab=v1(At);Hb=v2(At);Ak=np.empty((len(Ab),a.NB));Hk=np.empty((len(Hb),a.NB))
for i in range(a.NB):x=alm(K[:,i]);Ak[:,i]=v1(x);Hk[:,i]=v2(x)   # memory-light: no stored alm list
K=K.astype(np.float32)
sc=a.MMAX;Am=np.vstack([Ak*sc,a.WH*Hk*sc,np.sqrt(a.lam)*np.abs(Ab).max()*np.eye(a.NB)]);bv=np.r_[-Ab,-a.WH*Hb,np.zeros(a.NB)]
res=lsq_linear(Am,bv,bounds=(-1,1),lsmr_tol='auto',max_iter=5000);M=res.x*sc;Cm=K.astype(float)@M
Resid=B+Cm;Sky=Resid+Fmap   # deterministic residual (bg+bodies) and the full seed-42 low-l sky
def amp(m,L=40):cl=hp.anafast(m,lmax=L);return np.sqrt(cl*(2*np.arange(L+1)+1)/(4*np.pi))
def Dl(m,L=LMX):l=np.arange(L+1);return l*(l+1)/(2*np.pi)*hp.anafast(hp.remove_dipole(m),lmax=L)*T0**2
out=dict(vars(a),EPS=EPS,n_clusters=int((M>0.01*sc).sum()),n_voids=int((M<-0.01*sc).sum()),M_pos=float(M[M>0].sum()),M_neg=float(M[M<0].sum()),
  rms_bg=float(B.std()),dipole_bg=float(np.linalg.norm(hp.fit_dipole(B)[1])),rms_resid=float(Resid.std()),max_resid=float(np.abs(Resid-Resid.mean()).max()),
  dipole_resid=float(np.linalg.norm(hp.fit_dipole(Resid)[1])),amp_bg=amp(B).tolist(),amp_resid=amp(Resid).tolist(),amp_field=amp(Fmap).tolist(),amp_sky=amp(Sky).tolist(),
  cancel_err_max=float(np.abs((Resid+Ftarget)-(Resid+Ftarget).mean()).max()),cancel_err_rms=float((Resid+Ftarget).std()),field_l2_target_max=float(np.abs(Ftarget).max()),frac_at_bound=float(np.mean(np.abs(res.x)>0.999)))
if a.field:out.update(D2_field=float(Dl(Fmap)[2]),D2_sky=float(Dl(Sky)[2]),D2_target=float(a.f**2*Dl(Fmap)[2]))
tag=a.tag or f'R{R:g}_s{a.seed}_L{a.LC}_M{a.MMAX:.0e}_f{a.f:g}'
np.savez(O+f'p5_{tag}.npz',y=y,M=M,C=Cm,Resid=Resid,B=B,Fmap=Fmap,Sky=Sky,Dl_resid=Dl(Resid),Dl_sky=Dl(Sky),Dl_field=Dl(Fmap))
json.dump(out,open(O+f'p5_{tag}.json','w'),indent=1)
print(tag,{k:v for k,v in out.items() if not k.startswith('amp')},flush=True);print('amp_resid l1-6',np.array(out['amp_resid'][1:7]))
