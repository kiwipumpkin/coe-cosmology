import sys,json,numpy as np,healpy as hp;sys.path.insert(0,'../../../rotfit')
import fit
from scipy.interpolate import CubicSpline
from scipy.special import eval_legendre as P
from scipy.optimize import least_squares
LM=12;nside=32;npix=hp.nside2npix(nside);V=np.array(hp.pix2vec(nside,np.arange(npix))).T
PD=3362/2.7255e6;PRMS=PD/np.sqrt(3)
# real orthonormal SH basis l=1..12 built numerically from healpy
idx=[];B=[]
for l in range(1,LM+1):
    for m in range(0,l+1):
        for part in ((1,) if m==0 else (1,1j)):
            a=np.zeros(hp.Alm.getsize(LM),complex);a[hp.Alm.getidx(LM,l,m)]=part if m==0 else part/np.sqrt(2)*1
            mp=hp.alm2map(a,nside,lmax=LM);mp/=np.sqrt(np.mean(mp**2)*4*np.pi);B.append(mp);idx.append(l)
B=np.array(B);L=np.array(idx);K=len(L)
# smooth target: shrink_rec b_1100 (mass b, CatWISE slope, 1+z~1100), axis anti-CoE (0,+15)
o=np.load('../../../shrink_rec/b_1100_k0.npz')['o'];aC=fit.vec(0,15)
lnz=CubicSpline(o[:,0],np.log(o[:,1]))(np.degrees(np.arccos(np.clip(V@aC,-1,1))));m=np.exp(-lnz);m=m/m.mean()-1
cs=B@m*4*np.pi/npix
def met(c):return np.sqrt(3/(4*np.pi)*np.sum(c[L==1]**2))/PD,np.sqrt(np.sum(c**2)/(4*np.pi))/PRMS
print('smooth: dipole x',met(cs)[0],'rms(l1-12) x',met(cs)[1])
DLSS=14.3  # Gpc, Chris's ring/visible radius (assumption: last-scattering surface there)
# body families; ranges from Chris's seed-3 laws
TYP={'mass':dict(dl=30,dh=300,rmax=0.30),'void':dict(dl=-0.45,dh=-0.25,rmax=2.14)}
def bl(alpha):c=np.cos(alpha);return np.array([2*np.pi*(P(l-1,c)-P(l+1,c))/(2*l+1) for l in L])
def Yat(th,ph):
    pp,ww=hp.get_interp_weights(nside,th,ph);return B[:,pp]@ww
def body(p,t):   # p=(theta,phi,t_size 0..1,t_delta 0..1); cap where the body cuts the surface
    T=TYP[t];al=p[2]*T['rmax']/DLSS;d=T['dl']+p[3]*(T['dh']-T['dl']);amp=np.log1p(d)/3
    return amp*bl(al)*Yat(p[0],p[1])
def total(P_,types):
    c=cs.copy()
    for p,t in zip(P_,types):c=c+body(p,t)
    return c
def obj(c,which):return lrms(c,[3,4,5])/1e-5 if which=='m345' else (met(c)[0] if which=='dip' else met(c)[1])
SEL=np.isin(L,[3,4,5])
def lrms(c,ls):return np.sqrt(np.sum(c[np.isin(L,ls)]**2)/(4*np.pi))
def vecr(c,which):
    if which=='m345':return c[SEL]/np.sqrt(4*np.pi)/1e-5
    return (np.sqrt(3/(4*np.pi))*c[L==1]/PD) if which=='dip' else c/np.sqrt(4*np.pi)/PRMS
starts=np.array(hp.pix2ang(2,np.arange(hp.nside2npix(2)))).T
info={ 'dipole_xPlanck':met(cs)[0],'dipole_dTT':met(cs)[0]*PD,'v_kms':met(cs)[0]*PD*299792.458,**{f'l{l}_rms_xPlanckRMS':lrms(cs,[l])/PRMS for l in (2,3,4,5)},**{f'l{l}_rms_uK':lrms(cs,[l])*2.7255e6 for l in (2,3,4,5)},'l345_rms_xPlanckRMS':lrms(cs,[3,4,5])/PRMS,'l345_rms_uK':lrms(cs,[3,4,5])*2.7255e6,'PlanckRMS_uK':PRMS*2.7255e6};print(info);json.dump(info,open('start.json','w'),indent=1)
def run(fam,which,NMAX):
    Ps=[];ty=[];curve=[(0,obj(cs,which))]
    for n in range(1,NMAX+1):
        best=None
        for t in (['mass','void'] if fam=='both' else [fam]):
            for th,ph in starts:
                for ts in (0.7,):
                    lo=np.array([0,-4*np.pi,0,0.]);hi=np.array([np.pi,4*np.pi,1,1.]);x0=np.clip([th,ph,ts,0.5],lo+1e-6,hi-1e-6)
                    f=lambda x:obj(total(Ps+[x],ty+[t]),which)
                    r=least_squares(lambda x:vecr(total(Ps+[x],ty+[t]),which),x0,bounds=(lo,hi),max_nfev=60)
                    if best is None or r.cost<best[0]:best=(r.cost,list(r.x),t)
        Ps.append(best[1]);ty.append(best[2])
        # joint refine
        x0=np.concatenate(Ps);nb=len(Ps)
        lo=np.tile([0,-4*np.pi,0,0],nb);hi=np.tile([np.pi,4*np.pi,1,1],nb)
        x0=np.clip(x0,lo+1e-9,hi-1e-9)
        r=least_squares(lambda x:vecr(total(x.reshape(-1,4),ty),which),x0,bounds=(lo,hi),max_nfev=300)
        Ps=[list(v) for v in r.x.reshape(-1,4)];v=obj(total(Ps,ty),which);curve.append((n,v))
        print(fam,which,n,round(v,4),ty[-1],flush=True)
        if v<=0.8 or (len(curve)>6 and curve[-6][1]-v<0.01):break
    c=total(Ps,ty);return dict(curve=curve,types=ty,params=Ps,final_dipole=met(c)[0],final_rms=met(c)[1],final_l2_1e5=lrms(c,[2])/1e-5,final_l345_1e5=lrms(c,[3,4,5])/1e-5,first_N_1x=next((n for n,v in curve if v<=1),None))
fam,which,NM=sys.argv[1],sys.argv[2],int(sys.argv[3])
R=run(fam,which,NM);json.dump(R,open(f'res_{fam}_{which}.json','w'),indent=1,default=float)
