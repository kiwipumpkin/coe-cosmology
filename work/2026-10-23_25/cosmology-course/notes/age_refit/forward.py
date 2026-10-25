"""forward.py  CoE forward-integration refit (Chris's rules, 2026-10-15).
ASSUMPTIONS (shared by Migkas and CMB parts; one universe):
 A1 Spherical CoE: point mass M (kg, free, independent of D) + matter rho_i(r0) (law below) at start.
 A2 Radial Lagrangian shells, Newtonian: d2r/dtau2 = -G Menc(r0)/r^2, U=dr/dtau (proper speed, no cap, v>c allowed);
    beta=U/sqrt(c^2+U^2), gamma=sqrt(1+U^2/c^2). Common integration time tau = observer proper time (as unified2).
 A3 Launch: U_i(r0) = V0 (r0/r_i)^p at tau=0 (start); observer shell starts at r_i (start distance).
 A4 Temperature per shell: T = T_i(r0) (rho/rho_i)^(1/3), rho ~ 1/(r^2 dr/dr0); T_i(r0)=T_start (rho_i(r0)/rho_i(r_i))^(1/3), T_start=3e4 K at observer.
 A5 'Today' = when observer's local T = 2.725 K. D_today, H_today(dir), age = outputs.
 A6 Migkas: clusters on rings R=c z/Hbar (z=0.01..0.32, median .075), today's field, A=|v(x)-v_o|/R, mean-normalised; theta=0 away from CoE.
 A7 CMB: photons straight at c in CoE frame; trace back per direction to where emitting matter T=3000K; S=Dmotion*exp((Phi_e-Phi_o)/c^2) via gravshift.grav_factor.
 Exact symmetry r->lam r, tau->lam tau, M->lam M, rho->rho/lam^2 (U unchanged) leaves every dimensionless output unchanged,
 so r_i fixed=1 Gpc in the fit, then the solution is rescaled so Hbar_today=70 (presentation, not fitted).
"""
import sys,os,json,numpy as np,warnings;sys.path.insert(0,'/workspace/cosmology-course/notes/age_refit');from gravshift import grav_factor, gr_S;warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/cosmology/papers/spreadsheet_model/bubble'); _cwd=os.getcwd(); os.chdir('/workspace/cosmology/papers/spreadsheet_model/bubble')
from bubble_core import IMIG,MA,MHI,MLO,TH1,shifted
os.chdir(_cwd); sys.path.insert(0,'/workspace/cosmology-course/notes'); import importlib.util as iu
sp=iu.spec_from_file_location('cc','/workspace/cosmology-course/notes/cmb-cancel.py');
Y=MA; W=1/(0.5*(MLO+MHI))**2; CHI0=45.13662301218721
C=0.306601; GU=6.674e-11*(3.156e16)**2/(3.0857e25)**3; KMSMPC=977.8  # c Gpc/Gyr; G per kg in Gpc^3/Gyr^2; 1/Gyr->km/s/Mpc
ZMIG=np.array([0.01,0.025,0.04,0.055,0.075,0.10,0.14,0.20,0.32]); TST=3e4; TREC=3000.; TNOW=2.725
FAIL=[None]; GRMODE=os.environ.get('GRZ','newton')  # 'gr' = exact Schwarzschild (Cosmology), r_s from M_enc at observer
TMODE=os.environ.get('TMODE','uniform')  # A4 variant: 'uniform' = all matter at T_start at tau=0 (default); 'adiabatic' = T_i ~ rho_i^(1/3)
R0=np.geomspace(1e-3,3e4,900) if not os.environ.get('R0DENSE') else np.unique(np.r_[np.geomspace(1e-3,3e4,700),np.geomspace(0.2,5,1600)])  # dense near observer label 1 for GR rays
def shape(kind,x,s):
    if kind=='exp': return np.exp(-x)
    if kind=='eta': return np.exp(-x**s)
    if kind=='pl': return (1+x)**(-s)
    if kind=='soft': return 1/(1+x**s)
def setup(P):
    k=P['kind']; rho=P['rhoa']*shape(k,R0/P['L'],P.get('s')); 
    dm=4*np.pi*rho*R0**2*np.gradient(R0); Menc=P['M']+np.cumsum(dm)-dm/2
    return rho,dm,Menc
def Tof(Ti,r): return Ti*((R0/r)**2*np.gradient(R0)/np.gradient(r))**(1/3)
def integrate(P,keep=False):
    rho,dm,Me=setup(P); GM=GU*Me; io=np.searchsorted(R0,1.0)
    r=R0.copy(); U=P['V0']*R0**P['p']; t=0.0; Ti=TST*(rho/rho[io])**(1/3) if TMODE=='adiabatic' else np.full_like(R0,TST); T=Ti.copy()
    H=dict(Ts=[],rs=[],Us=[],ts=[])
    def f(r,U): return U,-GM/r**2
    def rec(T,r,U,t):
        if keep: H['Ts'].append(T);H['rs'].append(r);H['Us'].append(U);H['ts'].append(t)
    rec(T,r,U,t)
    for n in range(40000):
        dt=0.01*r[io]/max(abs(U[io]),1e-30)
        k1=f(r,U);k2=f(r+dt/2*k1[0],U+dt/2*k1[1]);k3=f(r+dt/2*k2[0],U+dt/2*k2[1]);k4=f(r+dt*k3[0],U+dt*k3[1])
        rn=r+dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); Un=U+dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        dead=~(rn>1e-6*R0)|~np.isfinite(rn); rn=np.where(dead,1e-6*R0,rn); Un=np.where(dead,0.0,Un)   # shells that fall into the CoE are absorbed (mass already in Menc)
        if dead[io:].any() or Un[io]<=0:
            FAIL[0]=('dead' if dead[io:].any() else 'recollapse',n); return None
        Tn=Tof(Ti,rn)
        if Tn[io]<=TNOW:
            w=np.log(T[io]/TNOW)/np.log(T[io]/Tn[io]); r=r+w*(rn-r);U=U+w*(Un-U);t+=w*dt;T=Tof(Ti,r); rec(T,r,U,t)
            o=dict(r=r,U=U,t=t,T=T,GM=GM,dm=dm,io=io,Ti=Ti,cross=bool(np.any(np.diff(r[io//2:])<=0)),ncollapsed=int(np.sum(r<=1.0001e-6*R0)),nstep=n)
            if keep: o.update({k:np.array(v) for k,v in H.items()})
            return o
        r,U,T,t=rn,Un,Tn,t+dt; rec(T,r,U,t)
    return None
def local(o):
    r,U,io,GM=o['r'],o['U'],o['io'],o['GM']; D=r[io]; Uo=U[io]
    lr=np.log(r); sl=slice(io-6,io+7); x=r[sl]-D; c=np.polyfit(x,U[sl],3)[::-1]   # U(r)=c0+c1 x+c2 x^2+c3 x^3
    Hpar=c[1]; Hper=Uo/D; Hbar=Hpar/3+2*Hper/3
    q=GM[io]*D/Uo**2/D*1.0; q=GM[io]/D**2*D/Uo**2      # q=-rdd r/U^2
    j=2*GM[io]/D**3*D**2/Uo**3*Uo                     # j = rddd r^2/U^3, rddd=2GM U/r^3
    return dict(D=D,Uo=Uo,Hpar=Hpar,Hper=Hper,Hbar=Hbar,U1=c[1],U2=2*c[2],U3=6*c[3],q=float(q),j=float(j),
                grav_over_kin=float(2*GM[io]/D/Uo**2),grav_over_kin_pointmass=float(2*GU*0+0))
def aflow(Ufun,D,Uo,Hbar,ring=ZMIG):
    th=np.radians(TH1); afs=[]
    for z in ring:
        R=C*z/Hbar; x=D+R*np.cos(th); y=R*np.sin(th); rr=np.hypot(x,y); u=Ufun(rr)
        dvx=u*x/rr-Uo; dvy=u*y/rr; a=np.hypot(dvx,dvy)/R; afs.append(a/a.mean())
    af=np.mean(afs,0); return af/af.mean()
def chi2_l0(af):
    L=np.arange(0,360,0.5); ch=[np.sum(W*(shifted(af,l)[IMIG]-Y)**2) for l in L]; i=int(np.argmin(ch))
    from scipy.optimize import minimize_scalar
    m=minimize_scalar(lambda l:np.sum(W*(shifted(af,l)[IMIG]-Y)**2),bounds=(L[i]-1,L[i]+1),method='bounded')
    return float(m.fun),float(m.x%360)
def migkas_forward(o):
    L=local(o); r,U=o['r'],o['U']; f=lambda x:np.interp(x,r,U)
    return aflow(f,L['D'],L['Uo'],L['Hbar']),L
def jerk_af(D,Hpar,Hper,k2,k3,Hbar=70/KMSMPC):
    Uo=Hper*D; f=lambda rr:Uo+Hpar*(rr-D)*(1+k2*(rr-D)/D/2+k3*((rr-D)/D)**2/6)
    return aflow(f,D,Uo,Hbar)
def params_from(kind,lM,G_i,p,lf,lL,s=None):
    """lM=log10 M[kg]; G_i=2 G Menc(r_i)/(r_i V0^2) (initial gravity/kinetic of observer shell); lf=log10 (matter mass inside r_i / M)."""
    M=10.0**lM; L=10**lL; P=dict(kind=kind,M=M,L=L,p=p,s=s,rhoa=1.0)
    rho,dm,Me=setup(P); io=np.searchsorted(R0,1.0); ms=np.sum(dm[:io])+dm[io]/2
    P['rhoa']=10**lf*M/ms if ms>0 else 0; rho,dm,Me=setup(P)
    P['V0']=np.sqrt(2*GU*Me[io]/(G_i*R0[io])); P['G_i']=G_i; P['f_matter']=10**lf; return P
def cmb(o,ell=np.arange(360.)):
    rs,Ts,Us,ts=o['rs'],o['Ts'],o['Us'],o['ts']; io=o['io']; D=o['r'][io]; Uo=o['U'][io]; tn=o['t']
    th=np.radians(ell); n=len(ell); found=np.zeros(n,bool); out={k:np.full(n,np.nan) for k in ('te','re','Ue','r0e','x','y','phie')}
    GM=o['GM']; dm=o['dm']
    def phi(r):  # potential of every shell at that snapshot
        inner=-GM/r; outer=-GU*(np.cumsum((dm/r)[::-1])[::-1]-dm/r); return inner+outer
    Tprev=None
    for k in range(len(ts)-1,-1,-1):
        Lb=C*(tn-ts[k]); x=D+Lb*np.cos(th); y=Lb*np.sin(th); rad=np.hypot(x,y)
        live=rs[k]>1.1e-6*R0; rr=rs[k][live]
        T=np.interp(rad,rr,Ts[k][live],left=np.nan,right=np.nan)
        hit=(~found)&(T>=TREC)
        if hit.any():
            for kk in np.where(hit)[0]:
                if Tprev is not None and np.isfinite(Tprev[kk]):
                    w=np.log(TREC/Tprev[kk])/np.log(T[kk]/Tprev[kk])
                else: w=0.0
                kp=min(k+1,len(ts)-1); tt=ts[kp]+w*(ts[k]-ts[kp]); Lb2=C*(tn-tt)
                xe,ye=D+Lb2*np.cos(th[kk]),Lb2*np.sin(th[kk]); re=np.hypot(xe,ye)
                # shell state interpolated between snapshots kp,k
                Ue=(1-w)*np.interp(re,rs[kp][live],Us[kp][live])+w*np.interp(re,rr,Us[k][live])
                ph=(1-w)*np.interp(re,rs[kp][live],phi(rs[kp])[live])+w*np.interp(re,rr,phi(rs[k])[live])
                for nm,v in zip(('te','re','Ue','x','y','phie','r0e'),(tt,re,Ue,xe,ye,ph,np.interp(re,rr,R0[live]))): out[nm][kk]=v
            found|=hit
        Tprev=T
        if found.all(): break
    out['found']=found; uo=Uo/C; ue=out['Ue']/C
    bo=np.array([uo/np.sqrt(1+uo**2),0.]); go=np.sqrt(1+uo**2)
    nx,ny=D-out['x'],-out['y']; nn=np.hypot(nx,ny); nx/=nn; ny/=nn   # propagation direction source->observer
    be=ue/np.sqrt(1+ue**2); ge=np.sqrt(1+ue**2); ex,ey=out['x']/out['re'],out['y']/out['re']
    Dm=go*(1-bo[0]*nx)/(ge*(1-be*(ex*nx+ey*ny)))
    Dg=grav_factor(out['phie'],phi(o['r'])[io],C**2)   # FIXED 2026-10-15 (was Phi_o-Phi_e)
    rs=2*GM[io]/C**2; S=Dm*Dg
    if GRMODE=='gr': S=gr_S(rs,out['x'],out['y'],out['Ue']/C,D,0.0,Uo/C,nx,ny)
    out.update(S=S,Dm=Dm,Dg=Dg,Tobs=TREC*S,D=D,Uo=Uo,tn=tn,rs=rs,rs_pm=2*GU*o['GM'][0]/GU*GU/C**2 if False else 2*GM[0]/C**2); return out
def dt_ls_formula(D,U,Rrec):
    """Closed form, LS sphere of fixed radius Rrec about the CoE (D<Rrec), observer at D moving outward with proper speed U.
    CoE-frame: dt = t_e(0deg) - t_e(180deg) = 2D/c.  Observer rest frame (Lorentz, gamma=sqrt(1+U^2/c^2), gamma*beta=U/c,
    events separated by dx=+2Rrec along the axis):  dt' = (2/c)[ D sqrt(1+U^2/c^2) - U Rrec/c ]."""
    return 2*D/C, (2/C)*(D*np.sqrt(1+(U/C)**2)-U*Rrec/C)
def cmb_rrec(o,Rrec,ell=np.arange(360.)):
    """Chris's picture: an element recombines when its own r crosses fixed R_rec; plasma = |x|<R_rec at all times
    (radius monotone in label).  Past light cone of observer (at D today) vs the cylinder r=R_rec: per direction the
    first (near) intersection  L(psi)=D cos psi - sqrt(R_rec^2 - D^2 sin^2 psi), psi = angle from the CoE direction
    (psi=180-theta).  Directions with D sin psi > R_rec never meet plasma (transparent back to tau=0)."""
    io=o['io']; D=o['r'][io]; Uo=o['U'][io]; tn=o['t']; th=np.radians(ell); psi=np.pi-th
    disc=Rrec**2-(D*np.sin(psi))**2; hit=(disc>=0)&(np.cos(psi)>0)
    L=np.where(hit,D*np.cos(psi)-np.sqrt(np.clip(disc,0,None)),np.nan); tLS=tn-L/C
    rs,Us,ts=o['rs'],o['Us'],o['ts']; Ue=np.full_like(L,np.nan)
    for i in np.where(hit)[0]:
        k=np.searchsorted(ts,tLS[i]); k=min(max(k,1),len(ts)-1); w=(tLS[i]-ts[k-1])/(ts[k]-ts[k-1])
        Ue[i]=(1-w)*np.interp(Rrec,rs[k-1],Us[k-1])+w*np.interp(Rrec,rs[k],Us[k])
    xr=np.where(hit,D+L*np.cos(th),np.nan); yr=np.where(hit,L*np.sin(th),np.nan)
    nx,ny=-np.cos(th),-np.sin(th); uo=Uo/C; ue=Ue/C
    Dm=np.sqrt(1+uo**2)*(1-uo/np.sqrt(1+uo**2)*nx)/(np.sqrt(1+ue**2)*(1-ue/np.sqrt(1+ue**2)*(xr*nx+yr*ny)/Rrec))
    tx=np.interp(Rrec,o['rs'][:,io],ts)   # observer's own crossing of R_rec
    psic=np.degrees(np.arcsin(min(1,Rrec/D)))
    return dict(hit=hit,tLS=tLS,rLS=np.where(hit,Rrec,np.nan),L=L,Ue=Ue,Dm=Dm,t_obs_cross=tx,psi_crit_deg=psic,
                R_critical_dist_from_us=D-Rrec,
                gap_180_integrated=float(tLS[180]-tx),
                gap_formula_constU=float((D-Rrec)*(1/Uo-1/C)),   # (D-Rrec)/U_o - (D-Rrec)/c, U_o today
                gap_formula_int=float(np.trapezoid(1/np.interp(np.linspace(Rrec,D,400),o['rs'][:,io],o['Us'][:,io]),np.linspace(Rrec,D,400))-(D-Rrec)/C))
def gap_dr(o,Rrec,dr):
    """time between observer crossing R_rec and the element initially dr closer to CoE crossing: integrated vs dr/U_o."""
    io=o['io']; rs,ts=o['rs'],o['ts']; tx=np.interp(Rrec,rs[:,io],ts)
    k=np.searchsorted(ts,tx); lab=np.interp(Rrec-dr,rs[k],R0); j=np.searchsorted(R0,lab)
    tj=np.interp(Rrec,rs[:,j],ts) if rs[-1,j]>Rrec else np.nan
    Ux=np.interp(tx,ts,o['Us'][:,io]); return float(tj-tx),float(dr/Ux)
def cmb_sr_local(o,q,ell=np.arange(360.)):
    """LOCAL-SR CMB, cmb-cancel.py structure. Emitter for emission-angle theta = matter at offset R=q*D_today from OUR position
    at emission time t_e (theta=0 away from CoE), on our past light cone: |x_o(t_n)-x_o(t_e)-R n_theta| = c (t_n-t_e).
    Emitter velocity = local flow U(r_e,t_e) (radial), observer U_o, Phi from same mass distribution.
    S = gamma_o(1-b_o.n)/[gamma_e(1-b_e.n)] * exp[(Phi_e-Phi_o)/c^2]  (T = T0*S)."""
    rs,Us,ts=o['rs'],o['Us'],o['ts']; io=o['io']; D=o['r'][io]; tn=o['t']; R=q*D; th=np.radians(ell)
    GM,dm=o['GM'],o['dm']
    def phi(r): return -GM/r-GU*(np.cumsum((dm/r)[::-1])[::-1]-dm/r)
    xo=rs[:,io]; ex=np.cos(th); ey=np.sin(th)
    # g_k(theta) = |x_o(tn)-x_o(tk)-R n| - c(tn-tk); find sign change scanning back
    X=xo[-1]-xo[:,None]-R*ex[None,:]; Yy=-R*ey[None,:]; g=np.hypot(X,Yy)-C*(tn-ts)[:,None]
    out=dict(te=np.full(len(th),np.nan),S=np.full(len(th),np.nan),Dm=None,Dg=None,re=np.full(len(th),np.nan))
    Dm=np.full(len(th),np.nan); Dg=np.full(len(th),np.nan)
    for i in range(len(th)):
        s=np.where(np.sign(g[1:,i])!=np.sign(g[:-1,i]))[0]
        if len(s)==0: continue
        k=s[-1]; w=g[k,i]/(g[k,i]-g[k+1,i]); te=ts[k]+w*(ts[k+1]-ts[k])
        xoe=xo[k]+w*(xo[k+1]-xo[k]); xe,ye=xoe+R*ex[i],R*ey[i]; re=np.hypot(xe,ye)
        live=rs[k]>1.1e-6*R0
        Ue=(1-w)*np.interp(re,rs[k][live],Us[k][live])+w*np.interp(re,rs[k+1][live],Us[k+1][live])
        Pe=(1-w)*np.interp(re,rs[k][live],phi(rs[k])[live])+w*np.interp(re,rs[k+1][live],phi(rs[k+1])[live])
        nx,ny=xo[-1]-xe,-ye; nn=np.hypot(nx,ny); nx/=nn; ny/=nn
        uo=o['U'][io]/C; ue=Ue/C; bo=uo/np.sqrt(1+uo**2); be=ue/np.sqrt(1+ue**2)
        Dm[i]=np.sqrt(1+uo**2)*(1-bo*nx)/(np.sqrt(1+ue**2)*(1-be*(xe*nx+ye*ny)/re))
        Dg[i]=grav_factor(Pe,phi(o['r'])[io],C**2); out['te'][i]=te; out['re'][i]=re
        if GRMODE=='gr': Dm[i]=gr_S(2*GM[io]/C**2,np.array([xe]),np.array([ye]),np.array([Ue/C]),xo[-1],0.0,o['U'][io]/C,np.array([nx]),np.array([ny]))[0]; Dg[i]=1.0
    out.update(S=Dm*Dg,Dm=Dm,Dg=Dg,D=D,R=R,rs=2*GM[io]/C**2,rs_pm=2*GM[0]/C**2); return out

# ---------------- light-cone (corrected) Migkas, backward integration -----------------
def _state_at(o,t):
    ts=o['ts']; k=int(np.clip(np.searchsorted(ts,t)-1,0,len(ts)-2)); w=(t-ts[k])/(ts[k+1]-ts[k])
    return (1-w)*o['rs'][k]+w*o['rs'][k+1],(1-w)*o['Us'][k]+w*o['Us'][k+1]
def _phi(o,r):
    GM,dm=o['GM'],o['dm']; return -GM/r-GU*(np.cumsum((dm/r)[::-1])[::-1]-dm/r)
def z_lightcone(o,ell=np.arange(360.),nD=240):
    """1+z(theta,Delta) along the past light cone: photon arrives at observer (today, x_o=(D,0)) from direction n_theta;
    emitter at x=x_o+n c Delta at time t_n-Delta, moving with local flow U(|x|,t). 1+z = 1/[Dmotion*grav_factor]."""
    io=o['io']; D=o['r'][io]; Uo=o['U'][io]; tn=o['t']; th=np.radians(ell); nx,ny=np.cos(th),np.sin(th)
    Dl=np.geomspace(1e-5,0.6,nD)*tn; Z=np.full((nD,len(th)),np.nan)
    uo=Uo/C; bo=uo/np.sqrt(1+uo**2); go=np.sqrt(1+uo**2); pho=_phi(o,o['r'])[io]
    for j,dl in enumerate(Dl):
        r,U=_state_at(o,tn-dl); live=r>1.1e-6*R0; x=D+C*dl*nx; y=C*dl*ny; re=np.hypot(x,y)
        if np.any(np.diff(r[live])<=0): order=np.argsort(r[live])
        else: order=slice(None)
        rr=r[live][order]; Ue=np.interp(re,rr,U[live][order])/C; ph=np.interp(re,rr,_phi(o,r)[live][order])
        be=Ue/np.sqrt(1+Ue**2); ge=np.sqrt(1+Ue**2); kx,ky=-nx,-ny   # propagation direction
        Dm=go*(1-bo*kx)/(ge*(1-be*(x*kx+y*ky)/re)); S=Dm*grav_factor(ph,pho,C**2); Z[j]=1/S-1
    return Dl,Z
def migkas_lightcone(o,zs=ZMIG):
    """For each direction & ring redshift z_k: find light-travel lookback Delta_k(theta) with z(Delta)=z_k on the past light cone,
    apparent Hubble rate H_app = c z_k/(c Delta_k) ; A/A_all ~ H_app (per ring mean-normalised, averaged over rings)."""
    Dl,Z=z_lightcone(o); afs=[]
    for zk in zs:
        d=np.full(Z.shape[1],np.nan)
        for i in range(Z.shape[1]):
            zi=Z[:,i]; s=np.where((zi[:-1]-zk)*(zi[1:]-zk)<=0)[0]
            if len(s): k=s[0]; w=(zk-zi[k])/(zi[k+1]-zi[k]); d[i]=Dl[k]+w*(Dl[k+1]-Dl[k])
        if not np.all(np.isfinite(d)): return None
        h=zk/d; afs.append(h/h.mean())
    af=np.mean(afs,0); return af/af.mean()
def integrate_back(P,keep=True):
    """Backward from TODAY: labels = today's radii R0 (observer at R0=1), U_today(r)=U_o (r/D)^p, density today rho(r).
    Integrate d2r/dtau2=-G Menc/r^2 backward until the observer shell reaches r<=1e-3 D (Big-Bang-like origin) or turns."""
    rho,dm,Me=setup(P); GM=GU*Me; io=np.searchsorted(R0,1.0); D=R0[io]
    r=R0.copy(); U=P['Uo']*(R0/D)**P['p']; t=0.0; H=dict(rs=[r],Us=[U],ts=[0.0])
    f=lambda r,U:(U,-GM/r**2)
    for n in range(40000):
        dt=-0.01*r[io]/max(abs(U[io]),1e-30)
        k1=f(r,U);k2=f(r+dt/2*k1[0],U+dt/2*k1[1]);k3=f(r+dt/2*k2[0],U+dt/2*k2[1]);k4=f(r+dt*k3[0],U+dt*k3[1])
        r=r+dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); U=U+dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]); t+=dt
        dead=~(r>1e-6*R0); r=np.where(dead,1e-6*R0,r); U=np.where(dead,0,U)
        H['rs'].append(r);H['Us'].append(U);H['ts'].append(t)
        if r[io]<=1e-3*D or U[io]<=0 or not np.isfinite(r[io]): break
    ts=np.array(H['ts'])[::-1]; age=-ts[0]; ts=ts+age
    rs=np.array(H['rs'])[::-1]; Us=np.array(H['Us'])[::-1]
    # Kepler bang time per shell (time since r=0), constant Menc: consistency of a common origin
    E=0.5*R0*0+0.5*(P['Uo']*(R0/D)**P['p'])**2-GM/R0; mu=GM; tb=np.full_like(R0,np.nan)
    a=mu/(2*np.abs(E)); 
    hp=E>0; ch=1+R0[hp]/a[hp]; eta=np.arccosh(ch); tb[hp]=np.sqrt(a[hp]**3/mu[hp])*(np.sinh(eta)-eta)
    hm=(E<0)&(R0<2*a); c_=1-R0[hm]/a[hm]; eta=np.arccos(c_); tb[hm]=np.sqrt(a[hm]**3/mu[hm])*(eta-np.sin(eta))
    o=dict(r=rs[-1],U=Us[-1],t=age,GM=GM,dm=dm,io=io,rs=rs,Us=Us,ts=ts,cross=bool(np.any(np.diff(rs[-1][io//2:])<=0)),
           ncollapsed=0,reached_origin=bool(rs[0][io]<=1.01e-3*D),age_kepler=float(tb[io]),tb=tb,E=E)
    return o
def cmb_uniform(o,a_rec=1/1100.,ell=np.arange(360.)):
    """Reference 'uniform simultaneous surface': all directions recombine at ONE time t_rec, defined by the observer
    shell's scale r_o(t_rec)/D_today = a_rec. Emitter in direction theta = local-flow matter at x=x_o+n c (t_n-t_rec).
    nu_obs/nu_emit = SR Doppler (proper speeds) * shared gravity factor (GRZ=gr -> exact Schwarzschild)."""
    io=o['io']; D=o['r'][io]; Uo=o['U'][io]; tn=o['t']; ts=o['ts']; ro_t=o['rs'][:,io]
    trec=float(np.interp(a_rec*D,ro_t,ts)); dl=tn-trec; th=np.radians(ell); nx,ny=np.cos(th),np.sin(th)
    r,U=_state_at(o,trec); live=r>1.1e-6*R0; x=D+C*dl*nx; y=C*dl*ny; re=np.hypot(x,y)
    inside=(re>=r[live].min())&(re<=r[live].max())
    o2=np.argsort(r[live]); rr=r[live][o2]; Ue=np.interp(re,rr,U[live][o2]); ph=np.interp(re,rr,_phi(o,r)[live][o2])
    uo=Uo/C; ue=Ue/C; bo=uo/np.sqrt(1+uo**2); be=ue/np.sqrt(1+ue**2); kx,ky=-nx,-ny
    Dm=np.sqrt(1+uo**2)*(1-bo*kx)/(np.sqrt(1+ue**2)*(1-be*(x*kx+y*ky)/re)); S=Dm*grav_factor(ph,_phi(o,o['r'])[io],C**2)
    if GRMODE=='gr': S=gr_S(2*o['GM'][io]/C**2,x,y,ue,D,0.0,uo,kx,ky)
    S=np.where(inside,S,np.nan)
    return dict(S=S,trec=trec,lookback=dl,r_emit=re,inside=inside,mean_z=float(np.nanmean(1/S-1)))
