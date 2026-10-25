# zero-centre law rho ∝ (x/xp)^a exp(-(x/xp)^b), x=r/D. Same dist2 LTB build + densweb/mfit smooth plane fit (real z, no vp, no Lambda).
import os,sys,csv,json,numpy as np
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8')
_cwd=os.getcwd();os.chdir(B+'dipole_amp');import dist2 as Dd;os.chdir(_cwd)
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar
RR=14.30048287260589
def shape_f(a,b,xp):return lambda x:np.where(x>0,(np.maximum(np.abs(x),1e-300)/xp)**a*np.exp(-(np.abs(x)/xp)**b),0.)
def ltb(a,b,xp,frac,D):
    sh=shape_f(a,b,xp);Rb=1.5*(1+RR/D)
    X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=sh(X)*(X<=Rb)
    m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
    F=Dd.build(frac*0.5*Dd.u**2,mt);S,t0=F['S'],F['t0']
    lo,hi=(0.9,1.1) if D>=10 else (max(F['xmin'],2e-3),min(1+3.2/D,39.))
    x=np.linspace(lo,hi,4001) if D>=10 else np.geomspace(lo,hi,8001);return F,x,S(t0,x)[0],S(t0,x,dx=1)[0]
W_=np.load('/workspace/redweb2/densweb/web_v_s1.npz');DG=np.r_[W_['DG'],2.5,3.0]
C=299792.458
vec=lambda l,b:np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))*np.ones_like(np.asarray(l,float))])
ax=vec(0.,15.)
M=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
L=np.array([float(r['l_deg']) for r in M]);Y=np.array([float(r['A_over_Aall']) for r in M]);Sg=0.5*(np.array([float(r['A_err_low']) for r in M])+np.array([float(r['A_err_high']) for r in M]));Wt=1/Sg**2
Tc=list(csv.DictReader(open('/workspace/cosmology/migkas/migkas2020_tableC1.csv')));CZ=C*np.array([float(r['z']) for r in Tc]);CL=np.array([float(r['l']) for r in Tc])
BINS=[np.where(np.abs((CL-x+180)%360-180)<=45)[0] for x in L]
chi=lambda m:float(np.sum(Wt*(m-Y)**2));r2=lambda m:float(np.corrcoef(m,Y)[0,1]**2) if np.std(m)>0 else 0.
NP=vec(np.arange(360.),0*np.arange(360.)).T
def migkas(xg,Rg,Rtg,D):
    Rt1=np.interp(1,xg,Rtg);H1=Rt1/np.interp(1,xg,Rg);SC=70/H1*D*1000;O=-ax;SP=np.empty((360,len(DG)))
    for k,d in enumerate(DG/D):
        P=O[None]+d*NP;r=np.linalg.norm(P,axis=1);x=np.interp(r,Rg,xg);v=np.interp(x,xg,Rtg)[:,None]*P/r[:,None]
        SP[:,k]=((v-Rt1*O[None])*NP).sum(1)*SC
    G=np.empty((360,len(CZ)))
    for j,c in enumerate(SP):
        if np.any(np.diff(c)<=0) or c[0]>CZ.min() or c[-1]<CZ.max():return None
        G[j]=CZ/(np.exp(np.interp(np.log(CZ),np.log(c),np.log(DG)))*1e3)
    G=G/G.mean(0)
    def mp(Dl):
        out=np.empty(len(L))
        for b in range(len(L)):
            t=(L[b]-Dl)%360;i0=int(np.floor(t))%360;f=t-np.floor(t);out[b]=np.mean((1-f)*G[i0,BINS[b]]+f*G[(i0+1)%360,BINS[b]])
        return out
    sc=[chi(mp(l)) for l in range(360)];l0=float(np.argmin(sc))
    r=minimize_scalar(lambda l:chi(mp(l)),bounds=(l0-1,l0+1),method='bounded');m=mp(r.x)
    return dict(chi2=chi(m),pearson_r2=r2(m),CoE_l=float(r.x%360),A=m.tolist())
def run(a,b,xp,frac,D):
    try:
        F,x,R,Rt=ltb(a,b,xp,frac,D);o=migkas(x,R,Rt,D)
        if o is None:return dict(a=a,b=b,xp=xp,frac=frac,D=D,error='nonmono')
        o.update(a=a,b=b,xp=xp,frac=frac,D=D,cross=bool(F['shell_crossing']),ndrop=int(F['n_dropped_inner']));return o
    except Exception as e:return dict(a=a,b=b,xp=xp,frac=frac,D=D,error=repr(e)[:200])
if __name__=='__main__':
    import time;t=time.time()
    # sanity: red reproduction through same migkas() using lowz ltb table
    T=np.load('/workspace/redweb2/densweb/ltb_t0.npz');print('red check',migkas(T['x'],T['R'],T['Rt'],42.)['chi2'])
    print(run(2,2,0.8,20000,42),time.time()-t)
