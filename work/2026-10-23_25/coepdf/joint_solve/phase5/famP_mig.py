# Phase5 (NEW) copy of famP.py with fam 'mig': the Migkas-fit density rho ~ exp(-(x/LDX)^ETA), x=r/D, truncated at Rb=1.5(1+1/DR) (same as prof_core rb_of(DR,1.5)). Env LDX=L/D, ETA. Original famP.py untouched.
# Shrinking-radius (critical density) CMB trace, 1+z~1100, for a density family re-locked to the CatWISE log-slope -0.041037 at our distance.
# Same shape is used for the mass distribution (enclosed-mass table, truncated at Rb as before) and the plasma density.
import sys,os,json,numpy as np
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');os.chdir(B+'dipole_amp')
import dist2 as D
from scipy.optimize import brentq
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp,cumulative_trapezoid
fam,p,frac,out=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),os.path.abspath(sys.argv[4]) if sys.argv[4].startswith('/') else sys.argv[4]
SLOPE=0.04103732750534994;DR=float(os.environ.get("DRV",42/14.30048287260589));Rb=1.5*(1+1/DR)
def base(kind,x,p):
    x=np.abs(x)
    return {'soft':lambda:1/(1+x**p),'exp':lambda:np.exp(-x),'eta':lambda:np.exp(-x**p)}[kind]()
if fam=='curv':  # slope locked at x=1, far/core curvature free: ln rho = -SLOPE ln x - p (x-1)^2
    shape=lambda x:np.abs(x)**-SLOPE*np.exp(-p*(np.abs(x)-1)**2);L=None
elif fam=='cub':  # slope locked: ln rho = -SLOPE ln x - p (x-1)^2 - q (x-1)^3, q from env Q
    qq=float(os.environ.get('Q','0'));shape=lambda x:np.abs(x)**-SLOPE*np.exp(-p*(np.abs(x)-1)**2-qq*(np.minimum(np.abs(x),Rb)-1)**3);L=None
elif fam=='shell':  # uniform background + Gaussian shell at x=1 of log-amplitude p, width W (env): ln rho = p*exp(-((x-1)/W)^2/2)
    Wv=float(os.environ.get('W','0.05'));shape=lambda x:np.exp(p*np.exp(-0.5*((np.abs(x)-1)/Wv)**2));L=None
elif fam=='combo':  # thin Gaussian shell (p,W) + broad power-law slope SL: ln rho = p*G - SL ln x
    Wv=float(os.environ.get('W','0.03'));SLv=float(os.environ.get('SL','0'));shape=lambda x:np.abs(x)**-SLv*np.exp(p*np.exp(-0.5*((np.abs(x)-1)/Wv)**2));L=None
elif fam=='hat':  # uniform background + odd (derivative-of-Gaussian) feature at x=1: local gradient p inside width W, flat outside
    Wv=float(os.environ.get('W','0.05'));shape=lambda x:np.exp(-p*(np.abs(x)-1)*np.exp(-0.5*((np.abs(x)-1)/Wv)**2));L=None
elif fam=='slc':  # FREE near-slope SL (env SL) + curvature p: ln rho = -SL ln x - p (x-1)^2
    SLv=float(os.environ.get('SL','0'));shape=lambda x:np.abs(x)**-SLv*np.exp(-p*(np.abs(x)-1)**2);L=None
elif fam=='outer':  # slope locked at x=1, mass rises beyond observer: ln rho = -SLOPE ln x + p*max(x-1,0)^2
    shape=lambda x:np.abs(x)**-SLOPE*np.exp(p*np.maximum(np.abs(x)-1,0)**2);L=None
elif fam=='slope':  # unlocked power law rho ~ x^-p (p=0.041 is the lock; p<=0 is outer-weighted)
    shape=lambda x:np.abs(x)**-p;L=None
elif fam=='core':  # slope locked at x=1 (core term ~e^-10 there), add central overdensity p
    shape=lambda x:np.abs(x)**-SLOPE*(1+p*np.exp(-np.abs(x)/0.05));L=None
elif fam=='poly':  # smooth polynomial: ln rho = sum_k C_k (x-1)^k, k=1..4 (env C1..C4); p unused
    CC=[float(os.environ.get('C%d'%k,'0')) for k in (1,2,3,4)]
    shape=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)));L=None
elif fam=='hyb':  # RED poly for x<=XS, tanh blend (width ~0.02) to REF slc 0.041037/0.035 renormalized at XS
    CC=[-0.041,-0.035,0.1353,0.0];XS=float(os.environ['XS'])
    lr=lambda x:sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC))
    lf=lambda x:-SLOPE*np.log(np.abs(x))-0.035*(np.abs(x)-1)**2
    off=lr(XS)-lf(XS)
    shape=lambda x:np.exp((1-0.5*(1+np.tanh((np.abs(x)-XS)/0.01)))*lr(x)+0.5*(1+np.tanh((np.abs(x)-XS)/0.01))*(lf(x)+off));L=None
elif fam=='mig':  # Migkas profile (phase5)
    LDX=float(os.environ['LDX']);ETA=float(os.environ['ETA']);shape=lambda x:np.exp(-(np.abs(x)/LDX)**ETA);L=None
elif fam=='pow': shape=lambda x:np.abs(x)**-SLOPE if True else 0;L=None
else:
    ls=lambda L:-(np.log(base(fam,1.001/L,p))-np.log(base(fam,0.999/L,p)))/0.002-SLOPE
    g=np.geomspace(0.05,1e5,600);v=np.array([ls(L) for L in g]);i=np.where(np.sign(v[:-1])!=np.sign(v[1:]))[0][-1]
    L=brentq(ls,g[i],g[i+1]);shape=lambda x:base(fam,x/L,p)
X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=np.where(X>0,shape(np.maximum(X,1e-12)),shape(1e-12))*(X<=Rb)
m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
F=D.build(frac*0.5*D.u**2,mt);S,t0=F['S'],F['t0']
R1=S(t0,1.)[0,0];Rp=S(t0,1.,dy=1)[0,0];Rt=S(t0,1.,dx=1)[0,0];Rtp=S(t0,1.,dx=1,dy=1)[0,0];ratio=(Rtp/Rp)/(Rt/R1)
def rho(t,x):
    R=S(t,x)[0,0];Rp=S(t,x,dy=1)[0,0];return shape(x)*x*x/(R*R*Rp)
tus=brentq(lambda t:S(t,1.0)[0,0]-1/1100.,F['T'][0]*1.001,t0);rc=rho(tus,1.0)
# FIX (phase5, Oct 26): te(x) was tabulated only on x in [0.6,1.6] and LINEARLY EXTRAPOLATED, but rays hit at x~0.2-2.6.
# Harmless for the flat poly te(x); garbage for the Migkas profile (te varies ~4x). Now: tabulate, trace, and if any hit lies
# outside the table, widen the table to cover it and retrace; abort loudly if it cannot be covered.
def tetab(XX):
    if os.environ.get('SURF','rho')=='t':  # (NEW Oct 26 PM) constant proper time since the (common) bang: t_rec = time our shell had a=1/1100 (or env TREC)
        return np.full(len(XX),float(os.environ.get('TREC',tus)))
    if os.environ.get('SURF','rho')=='a':  # recombination where each shell's own expansion since then is 1100 (T~1/a)
        return np.array([brentq(lambda t:np.log(S(t,x)[0,0]/S(t0,x)[0,0]*1100.),F['T'][0]*1.0001,t0,xtol=1e-16,rtol=1e-14) for x in XX])
    return np.array([brentq(lambda t:np.log(rho(t,x)/rc),F['T'][0]*1.0001,t0,xtol=1e-16,rtol=1e-14) for x in XX])
def AB(t,w):
    Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
    return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B_,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B_,0.5*At/A**2*pw**2+0.5*Bt/B_**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B_**2*pp**2,0]
A0=AB(t0,1.0)[1];B0=AB(t0,1.0)[4]
lo,hi=0.6,1.6
for it in range(6):
    XX=np.linspace(lo,hi,max(201,int(200*(hi-lo))+1));te=tetab(XX)
    if not np.all(np.isfinite(te)):sys.exit('te not finite')
    tf=interp1d(XX,te,kind='cubic',bounds_error=False,fill_value=np.nan)
    hit=lambda l,y:y[0]-float(tf(np.clip(y[1],lo,hi)));hit.terminal=True
    o=[]
    for ps in (np.array([float(v) for v in os.environ['PSI'].split(',')]) if os.environ.get('PSI') else np.arange(0,180.001,5)):  # PSI env (NEW): custom ray list for diagnostics
        pr=np.radians(ps);y0=[t0,1.0,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
        s=solve_ivp(rhs,(0,-200),y0,events=[hit],rtol=1e-10,atol=1e-13)
        if len(s.y_events[0])==0:sys.exit(f'ray psi={ps} never reached recombination')
        y=s.y_events[0][0];o.append((ps,-y[3],y[0],t0-y[0],y[1],y[2],S(y[0],y[1])[0,0]))
    w=np.array([q[4] for q in o]);print('iter',it,'table',lo,hi,'hits',w.min(),w.max(),flush=True)
    if w.min()>=lo+1e-9 and w.max()<=hi-1e-9:break
    lo=max(1e-3,min(lo,0.8*w.min()));hi=max(hi,1.2*w.max())
else:sys.exit('could not cover ray hits with te table')
O_=np.array(o)
if not np.all(np.isfinite(O_)):sys.exit('non-finite trace')
np.savez(out,o=np.array(o),fam=fam,p=p,frac=frac,L=L if L else np.nan,ratio=ratio,q=float(os.environ.get('Q','0')),sl=float(os.environ.get('SL','nan')),W=float(os.environ.get('W','nan')),XS=float(os.environ.get('XS','nan')),C=np.array([float(os.environ.get('C%d'%k,'0')) for k in (1,2,3,4)]),te=te,XX=XX)
print(out,'L/D',L,'ratio',ratio,'1+z',o[0][1],o[-1][1],flush=True)
