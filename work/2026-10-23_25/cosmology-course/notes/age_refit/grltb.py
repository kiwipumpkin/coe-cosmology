"""Exact-GR (LTB) light propagation in the CoE universes.
Dynamics: LTB dust + central point mass. Areal radius R(t,r0), comoving label r0, synchronous proper time t.
  (dR/dt)^2 = 2 m(r0)/R + 2E(r0)   [c=1, m=G M_enc/c^2]  -- EXACT GR, identical in form to the Newtonian shell eq.
  so forward.integrate / integrate_back already give exact LTB R(t,r0) (U = dR/dtau, may exceed c; no horizon issue
  in comoving coordinates, white-hole/expanding branch is just Rdot>0 inside R<2m).
Metric (equatorial plane): ds^2 = -dt^2 + A^2 dr0^2 + R^2 dphi^2,  A = R'/sqrt(1+2E).
Null geodesics (radial AND non-radial), parameter t, k=(k^t,k^r,k^phi), L=R^2 k^phi conserved:
  dk^t/dl = -A Adot (k^r)^2 - R Rdot (k^phi)^2
  dk^r/dl = -2 (Adot/A) k^t k^r - (A'/A)(k^r)^2 + (R R'/A^2)(k^phi)^2
Observer and emitters are comoving dust (u=d/dt) => 1+z = k^t_emit / k^t_obs  (all Doppler+gravity+expansion, exact).
Emitters reachable = where the backward-integrated null ray actually goes (no straight-line assumption)."""
import numpy as np, forward as F
from scipy.interpolate import RegularGridInterpolator as RGI
C=F.C
def fields(o):
    ts=o['ts']*C; R=o['rs']; U=o['Us']/C; lab=F.R0; m=o['GM']/C**2
    E=0.5*U[-1]**2-m/R[-1]                     # conserved per shell (c=1)
    W2=1+2*E; ok_shell=W2>0; W=np.sqrt(np.clip(W2,1e-30,None))
    x=np.log(lab); Rp=np.gradient(R,x,axis=1)/lab; Up=np.gradient(U,x,axis=1)/lab
    A=Rp/W; Ad=Up/W; Ap=np.gradient(A,x,axis=1)/lab
    dead=~(R>1.1e-6*lab)
    rho=np.gradient(m,x)/lab/np.clip(R**2*Rp,1e-300,None)
    kw=dict(bounds_error=False,fill_value=np.nan)
    f={k:RGI((ts,x),v,**kw) for k,v in dict(R=R,U=U,A=A,Ad=Ad,Ap=Ap,Rp=Rp,dead=dead.astype(float),rho=rho).items()}
    f.update(ts=ts,x=x,io=o['io'],tn=ts[-1],ok_shell=ok_shell,cross=bool(np.any(Rp[-1]<=0)))
    return f
def trace(f,ell=np.arange(0,360,2.),nstep=500,tmin_frac=1e-4,lbmax_frac=1-1e-7):
    """Backward null rays from observer (label io, t_n) arriving from direction theta (0=away from CoE).
    Returns lookback grid and per-ray k^t (=> 1+z), label, phi, alive."""
    io=f['io']; tn=f['tn']; x0=f['x'][io]; th=np.radians(ell); n=len(th)
    def ev(t,x):
        p=np.c_[np.full_like(x,t),x]; return {k:f[k](p) for k in ('R','U','A','Ad','Ap','Rp','dead')}
    g=ev(tn,np.full(n,x0)); Ro,Ao=g['R'],g['A']
    kt=np.ones(n); kr=-np.cos(th)/Ao; L=-Ro*np.sin(th)*1.0   # k^phi = L/R^2 ; physical components -cos,-sin
    lam=np.exp(x0)*0+0; X=np.full(n,x0); phi=np.zeros(n)
    span=lbmax_frac*(tn-f['ts'][0]); t0s=f['ts'][0]
    if lbmax_frac<0.9999: lb=np.concatenate([[0],np.geomspace(tmin_frac*tn*1e-2,span,nstep)])   # Migkas (low z): resolve near today
    else:   # CMB: resolve BOTH near today and near the start (log in t), FIXED 2026-10-15
        early=tn-(t0s+np.geomspace(max(1e-9*tn,(f['ts'][1]-t0s)*1e-3),0.5*(tn-t0s),nstep//2))
        lb=np.unique(np.r_[0,np.geomspace(tmin_frac*tn*1e-2,0.5*(tn-t0s),nstep//2),early[early<=span]])
    KR=[kr.copy()];KT=[kt.copy()];XX=[X.copy()];PH=[phi.copy()];alive=np.ones(n,bool);AL=[alive.copy()]
    def rhs(t,X,kt,kr):
        lab=np.exp(X); g=ev(t,X); R=g['R'];A=g['A'];Ad=g['Ad'];Ap=g['Ap'];Rp=g['Rp'];U=g['U']
        kph=L/R**2
        dkt=-A*Ad*kr**2-R*U*kph**2
        dkr=-2*(Ad/A)*kt*kr-(Ap/A)*kr**2+(R*Rp/A**2)*kph**2
        # parameter t: d/dt = (1/kt) d/dl ; X=ln r0 => dX/dt = k^r/(kt*lab)
        return kr/(kt*lab), kph/kt, dkt/kt, dkr/kt, g['dead']
    for j in range(1,len(lb)):
        t0=tn-lb[j-1]; h=-(lb[j]-lb[j-1])
        a1=rhs(t0,X,kt,kr); a2=rhs(t0+h/2,X+h/2*a1[0],kt+h/2*a1[2],kr+h/2*a1[3])
        a3=rhs(t0+h/2,X+h/2*a2[0],kt+h/2*a2[2],kr+h/2*a2[3]); a4=rhs(t0+h,X+h*a3[0],kt+h*a3[2],kr+h*a3[3])
        X=X+h/6*(a1[0]+2*a2[0]+2*a3[0]+a4[0]); phi=phi+h/6*(a1[1]+2*a2[1]+2*a3[1]+a4[1])
        kt=kt+h/6*(a1[2]+2*a2[2]+2*a3[2]+a4[2]); kr=kr+h/6*(a1[3]+2*a2[3]+2*a3[3]+a4[3])
        bad=~np.isfinite(X)|~np.isfinite(kt)|(a4[4]>0.5)|(kt<=0)
        alive&=~bad; X=np.where(alive,X,np.nan); kt=np.where(alive,kt,np.nan)
        KR.append(np.where(alive,kr,np.nan));KT.append(kt.copy());XX.append(X.copy());PH.append(phi.copy());AL.append(alive.copy())
    return dict(L=L,kr=np.array(KR),lb=lb,kt=np.array(KT),X=np.array(XX),phi=np.array(PH),alive=np.array(AL),ell=ell)
def migkas(tr,zs=F.ZMIG):
    Z=tr['kt']-1; lb=tr['lb']; afs=[]
    for zk in zs:
        d=np.full(Z.shape[1],np.nan)
        for i in range(Z.shape[1]):
            zi=Z[:,i]; s=np.where(np.isfinite(zi[:-1])&np.isfinite(zi[1:])&((zi[:-1]-zk)*(zi[1:]-zk)<=0))[0]
            if len(s): k=s[0]; w=(zk-zi[k])/(zi[k+1]-zi[k]); d[i]=lb[k]+w*(lb[k+1]-lb[k])
        if not np.all(np.isfinite(d)): return None
        h=zk/d; afs.append(h/h.mean())
    af=np.mean(afs,0); af=np.interp(np.arange(360.),np.append(tr['ell'],360),np.append(af,af[0])); return af/af.mean()
def cmb_uniform(f,tr,o,a_rec=1/1100.):
    io=o['io']; D=o['r'][io]; trec=np.interp(a_rec*D,o['rs'][:,io],o['ts'])*C; dl=f['tn']-trec
    if dl>tr['lb'][-1]: return np.full(360,np.nan),trec/C
    kt=np.array([np.interp(dl,tr['lb'],tr['kt'][:,i]) if tr['alive'][:,i][np.searchsorted(tr['lb'],dl)-1 if np.searchsorted(tr['lb'],dl)>0 else 0] else np.nan for i in range(tr['kt'].shape[1])])
    S=1/kt; return np.interp(np.arange(360.),np.append(tr['ell'],360),np.append(S,S[0])),trec/C
def cmb_Ttraced(f,tr,Trec_ratio=3000/2.725):
    """direction-dependent LS: dust T ~ rho^(1/3) normalised to 2.725 K at observer today; LS where T reaches 3000 K."""
    io=f['io']; rho_o=f['rho']([[f['tn'],f['x'][io]]])[0]; n=tr['kt'].shape[1]; S=np.full(n,np.nan); tls=np.full(n,np.nan)
    for i in range(n):
        al=tr['alive'][:,i]; t=f['tn']-tr['lb'][al]; X=tr['X'][al,i]
        if len(t)<2: continue
        Tr=(f['rho'](np.c_[t,X])/rho_o)**(1/3); s=np.where(np.isfinite(Tr)&(Tr>=Trec_ratio))[0]
        if len(s): k=s[0]; S[i]=1/tr['kt'][al][k,i]; tls[i]=t[k]/C
    sf=lambda v: np.interp(np.arange(360.),np.append(tr['ell'],360),np.append(v,v[0]))
    return sf(S),sf(tls)
def _fgrid(trec,tn,n):
    a=max(trec,1e-9*tn); mid=0.5*(a+tn)
    t1=np.geomspace(a,mid,n//2); t2=tn-np.geomspace(tn-mid,1e-7*tn,n//2)
    return np.unique(np.r_[trec,t1,t2,tn])
def fan_forward(f,trec,nlab=240,npsi=181,nstep=500):
    """Forward map from the recombination slice t=trec: by spherical symmetry a photon's history depends only on
    (emitter label r0e, emission angle psi w.r.t. outward radial in the emitter's comoving frame); phi just rotates.
    Integrate forward to t_n; rays with final label == observer label are on our past light cone today
    (rotate phi so they land on us). Returns emission label, psi, arrival sky angle theta (0=from away-CoE), z, phi offset."""
    io=f['io']; x0=f['x'][io]; tn=f['tn']
    xs=f['x']; alive_lab=np.isfinite(f['R'](np.c_[np.full_like(xs,trec),xs]))&(f['dead'](np.c_[np.full_like(xs,trec),xs])<0.5)
    lab=np.linspace(xs[alive_lab].min(),xs[alive_lab].max(),nlab) if alive_lab.any() else np.array([])
    if not len(lab): return None
    psi=np.radians(np.linspace(0,180,npsi)); Xe,Pe=np.meshgrid(lab,psi,indexing='ij'); Xe=Xe.ravel(); Pe=Pe.ravel(); n=len(Xe)
    def ev(t,x):
        p=np.c_[np.full_like(x,t),x]; return {k:f[k](p) for k in ('R','U','A','Ad','Ap','Rp','dead')}
    g=ev(trec,Xe); kt=np.ones(n); kr=np.cos(Pe)/g['A']; L=g['R']*np.sin(Pe); X=Xe.copy(); phi=np.zeros(n); alive=np.isfinite(kr)
    ts=_fgrid(trec,tn,nstep); prevX=X.copy()
    def rhs(t,X,kt,kr):
        lab_=np.exp(X); g=ev(t,X); R=g['R'];A=g['A'];Ad=g['Ad'];Ap=g['Ap'];Rp=g['Rp'];U=g['U']; kph=L/R**2
        return kr/(kt*lab_), kph/kt, (-A*Ad*kr**2-R*U*kph**2)/kt, (-2*(Ad/A)*kt*kr-(Ap/A)*kr**2+(R*Rp/A**2)*kph**2)/kt, g['dead']
    for j in range(len(ts)-1):
        t0=ts[j]; h=ts[j+1]-ts[j]
        a1=rhs(t0,X,kt,kr); a2=rhs(t0+h/2,X+h/2*a1[0],kt+h/2*a1[2],kr+h/2*a1[3]); a3=rhs(t0+h/2,X+h/2*a2[0],kt+h/2*a2[2],kr+h/2*a2[3]); a4=rhs(t0+h,X+h*a3[0],kt+h*a3[2],kr+h*a3[3])
        X=X+h/6*(a1[0]+2*a2[0]+2*a3[0]+a4[0]); phi=phi+h/6*(a1[1]+2*a2[1]+2*a3[1]+a4[1]); kt=kt+h/6*(a1[2]+2*a2[2]+2*a3[2]+a4[2]); kr=kr+h/6*(a1[3]+2*a2[3]+2*a3[3]+a4[3])
        alive&=np.isfinite(X)&np.isfinite(kt)&(a4[4]<0.5)&(kt>0)
    g=ev(tn,np.where(alive,X,x0))
    # arrival local direction: photon propagation (k^rhat,k^phihat)=(A kr, R kphi)/kt ; source direction = -that
    kp=(g['A']*kr)/kt; kq=(g['R']*L/g['R']**2)/kt
    res=dict(Xe=Xe.reshape(nlab,npsi),psi=Pe.reshape(nlab,npsi),Xf=np.where(alive,X,np.nan).reshape(nlab,npsi),
             z=(kt-1).reshape(nlab,npsi),theta=np.degrees(np.arctan2(-kq,-kp)%(2*np.pi)).reshape(nlab,npsi),dphi=phi.reshape(nlab,npsi),x0=x0)
    # reachable: per emitter label, psi where Xf crosses x0
    hits=[]
    for i in range(nlab):
        xf=res['Xf'][i]-x0; s=np.where(np.isfinite(xf[:-1])&np.isfinite(xf[1:])&(xf[:-1]*xf[1:]<=0))[0]
        for k in s:
            w=xf[k]/(xf[k]-xf[k+1]); hits.append(dict(lab=float(np.exp(lab[i])),psi=float(np.degrees(Pe.reshape(nlab,npsi)[i,k]+w*(psi[1]-psi[0]))),
                 z=float(res['z'][i,k]+w*(res['z'][i,k+1]-res['z'][i,k])),theta=float(res['theta'][i,k]),emit_phi=float(-res['dphi'][i,k])))
    res['hits']=hits; return res

def forward_check(f,tr,trec,nstep=800):
    """Backward trace = guess; re-integrate each surviving ray FORWARD from E_rec (t=trec) with its (r0,k) and check it
    lands on the observer today (label and arrival angle). Returns per-ray miss in label and angle, and z both ways."""
    tn=f['tn']; dl=tn-trec; j=int(np.searchsorted(tr['lb'],dl)); j=min(j,len(tr['lb'])-1)
    w=(dl-tr['lb'][j-1])/(tr['lb'][j]-tr['lb'][j-1]); I=lambda a:(1-w)*a[j-1]+w*a[j]
    X=I(tr['X']); kt=I(tr['kt']); kr=I(tr['kr']); L=tr['L'].copy(); ok=np.isfinite(X)&np.isfinite(kt)&np.isfinite(kr)
    X=np.where(ok,X,f['x'][f['io']]); kt=np.where(ok,kt,1.0); kr=np.where(ok,kr,0.0)
    kt0=kt.copy(); kt=kt/kt0; kr=kr/kt0; L=L/kt0      # renormalise k^t=1 at emission (affine rescale)
    def ev(t,x):
        p=np.c_[np.full_like(x,t),x]; return {k:f[k](p) for k in ('R','U','A','Ad','Ap','Rp','dead')}
    def rhs(t,X,kt,kr):
        lab_=np.exp(X); g=ev(t,X); R=g['R'];A=g['A'];Ad=g['Ad'];Ap=g['Ap'];Rp=g['Rp'];U=g['U']; kph=L/R**2
        return kr/(kt*lab_), kph/kt, (-A*Ad*kr**2-R*U*kph**2)/kt, (-2*(Ad/A)*kt*kr-(Ap/A)*kr**2+(R*Rp/A**2)*kph**2)/kt
    # sign: backward trace k is future-pointing propagation; integrate forward in t
    ts=_fgrid(trec,tn,nstep)
    phi=np.zeros_like(X)
    for i in range(len(ts)-1):
        t0=ts[i]; h=ts[i+1]-t0
        a1=rhs(t0,X,kt,kr); a2=rhs(t0+h/2,X+h/2*a1[0],kt+h/2*a1[2],kr+h/2*a1[3]); a3=rhs(t0+h/2,X+h/2*a2[0],kt+h/2*a2[2],kr+h/2*a2[3]); a4=rhs(t0+h,X+h*a3[0],kt+h*a3[2],kr+h*a3[3])
        X=X+h/6*(a1[0]+2*a2[0]+2*a3[0]+a4[0]); phi=phi+h/6*(a1[1]+2*a2[1]+2*a3[1]+a4[1]); kt=kt+h/6*(a1[2]+2*a2[2]+2*a3[2]+a4[2]); kr=kr+h/6*(a1[3]+2*a2[3]+2*a3[3]+a4[3])
    g=ev(tn,X); kp=g['A']*kr/kt; kq=g['R']*(L/g['R']**2)/kt; th_arr=np.degrees(np.arctan2(-kq,-kp))%360
    z_fwd=1/kt-1   # nu_obs/nu_emit = kt_obs/kt_emit with kt_emit=1  -> 1+z = 1/kt_obs
    return dict(ok=ok,label_miss=np.where(ok,np.exp(X)/np.exp(f['x'][f['io']])-1,np.nan),theta_arr=np.where(ok,th_arr,np.nan),
                z_back=np.where(ok,I(tr['kt'])-1,np.nan),z_fwd=np.where(ok,z_fwd,np.nan),emit_label=np.where(ok,np.exp(I(tr['X'])),np.nan),
                emit_phi=np.where(ok,I(tr['phi']),np.nan))
