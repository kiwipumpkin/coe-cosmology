"""Light-cone CoE model (Chris's confirmed construction, 2026-10-15).  Exact GR: LTB dust + central point mass, c=1, Gpc.
ASSUMPTIONS ADDED (listed for the report):
 L1 Simultaneous bang: every shell starts at t=0, R=0 (common origin)  -> tB(x)=0.
 L2 Shell energy E(r0)=0.5*(V0*r0^p)^2 (>0, hyperbolic; V0*r0^p is each shell's asymptotic speed in c units). r0 = mass label.
 L3 Observer = comoving dust element at label r0=1 (scale symmetry; results rescaled so the observer's mean H today = 70).
 L4 Element temperature T = K*(rho/rho_obs,today)^(1/3) (adiabatic photon gas, uniform photon/baryon ratio). Recombination at 3000 K.
 L5 Emission point per direction: earliest point on the past null geodesic (t_e>=t*) whose matter is recombined at t_e.
    If the ray reaches t* before meeting plasma, emitter is the point at t_e=t* (Chris: nothing emitted before our crossing).
 L6 Migkas statistic: clusters on our past light cone; per direction lookback distance d(z) where the ray's z=zk; A=z/d, mean-normalised.
"""
import sys,os,numpy as np,warnings;warnings.filterwarnings('ignore')
sys.path[:0]=['/workspace/cosmology-course/notes/age_refit','/workspace/cosmology/papers/spreadsheet_model/bubble']; _c=os.getcwd(); os.chdir('/workspace/cosmology-course/notes/age_refit')
import forward as F, ltb2 as T
import importlib.util as iu
sp=iu.spec_from_file_location('cc','/workspace/cosmology-course/notes/cmb-cancel.py'); CC=iu.module_from_spec(sp); sp.loader.exec_module(CC)
os.chdir(_c)
from scipy.integrate import solve_ivp
from scipy.stats import norm,chi2 as C2
def sig(chi,k): d=F.CHI0-chi; return float(norm.isf(C2.sf(d,k))) if d>0 else 0.0
class LC(T.LTB):
    def __init__(s,P):
        P=dict(P); P.setdefault('V0',1); super().__init__(P,'fwd'); s.P=P
    def E(s,x): return 0.5*(s.P['V0']*np.exp(x)**s.P['p'])**2
    def tB(s,x): return 0.0*np.asarray(x)
def make(law,z):
    """z = [log10 m_c (Gpc), log10 rhoa(kg/Gpc^3 rel), log10 L, s?, log10 V0, p, log10 tn]"""
    d=dict(zip(NAMES[law],z))
    P=dict(kind=law,M=10**d['lM']*F.C**2/F.GU,rhoa=10**d['lrho']*F.C**2/F.GU,L=10**d['lL'],s=d.get('s'),V0=10**d['lV'],p=d['p'])
    return LC(P),10**d['ltn']
NAMES={k:['lM','lrho','lL']+(['s'] if k in('eta','pl','soft') else [])+['lV','p','ltn'] for k in ('exp','eta','pl','soft')}
def obs(L,tn):
    X=np.array([0.0]); R,U,Rp,A,Ad,Ap=[v[0] for v in L.geo(tn,X)]
    Hr=Ad/A; Ht=U/R; return dict(R=R,U=U,Hpar=Hr,Hper=Ht,Hbar=(Hr+2*Ht)/3,rho=L.rho_rel(tn,X)[0])
class _Budget(Exception): pass
def ray(L,tn,th,zmax=None,tmin=None,rtol=1e-9,budget=40000,evx=None):
    """past null geodesic from observer (tn, r0=1) arriving from sky angle th (0 = away from CoE). Integrated in chunks so that a
    stalled solver keeps its earlier samples and is FLAGGED (stop='integrator_stopped...'), never scored as an emitter."""
    o=obs(L,tn); Lang=-o['R']*np.sin(th); f0=T.rhs(L,Lang); y=np.array([0.,0.,0.,-np.cos(th)]); t=tn
    t1=tmin if tmin is not None else 1e-6*tn
    ec=lambda t,y: y[0]-(L.xmin+1.0); ec.terminal=True; ec.direction=-1
    ez=(lambda t,y: y[2]-np.log1p(zmax)) if zmax else (lambda t,y:1.); ez.terminal=True
    ex=lambda t,y: y[0]-(L.xmax-1.5); ex.terminal=True
    evs=(ec,ez,ex)+((evx,) if evx else ())
    cnt=[0]
    def fr(t,y):
        cnt[0]+=1
        if cnt[0]>budget: raise _Budget()
        return f0(t,y)
    out=[]; nc=0; stop='reached_tmin'; edges=np.geomspace(tn,t1,41)[1:]; j=0
    while j<len(edges):
        try: s=solve_ivp(fr,(t,edges[j]),y,method='LSODA',rtol=rtol,atol=1e-12,events=evs,dense_output=True)
        except _Budget: stop='integrator_stopped(step budget)'; break
        except Exception: stop='integrator_exception'; break
        fin=np.isfinite(s.y).all(0); tt=s.t[fin]
        if len(tt)>1:
            q=np.unique(np.concatenate([np.linspace(tt[k],tt[k+1],5) for k in range(len(tt)-1)]))[::-1]
            Y=np.array([s.sol(v) for v in q]).T; out.append(np.vstack([q,Y,np.full(len(q),nc)]))
        if not fin.all() or s.status<0: stop='integrator_stopped(nonfinite/fail)'; break
        if s.status==1:
            if len(s.t_events[0]):
                if nc>=3: stop='integrator_stopped(>3 centre passages)'; break
                t=s.t_events[0][0]; y=s.sol(t).copy(); y[1]+=np.pi; y[3]=-y[3]; nc+=1; t+=(edges[j]-t)*1e-9; continue
            stop='zmax' if len(s.t_events[1]) else ('circle_hit' if evx and len(s.t_events[3]) else 'left_label_grid'); break
        t=s.t[-1]; y=s.y[:,-1]; j+=1
    D=np.hstack(out) if out else np.zeros((6,0))
    return dict(t=D[0],x=D[1],phi=D[2],lnz1=D[3],nc=D[5],stop=stop,th=th)
def migkas(L,tn,ells=np.arange(0,360,6.)):
    afs=np.zeros((len(F.ZMIG),len(ells))); rays=[]
    for i,e in enumerate(ells):
        r=ray(L,tn,np.radians(e),zmax=F.ZMIG[-1]*1.05,rtol=1e-8,budget=15000); rays.append(r)
        z=np.expm1(r['lnz1']); lb=tn-r['t']
        if len(z)<3 or z.max()<F.ZMIG[-1] or np.any(np.diff(z)<=0): return None
        d=np.interp(F.ZMIG,z,lb); afs[:,i]=F.ZMIG/d
    afs/=afs.mean(1,keepdims=True); af=afs.mean(0)
    af=np.interp(np.arange(360.),np.append(ells,360),np.append(af,af[0])); return af/af.mean()
from scipy.optimize import brentq
LRMIN=np.log(1e-12)
def lnJ(L,t,x,dx=1e-4):
    x=np.atleast_1d(x).astype(float); t=np.broadcast_to(t,x.shape).astype(float)
    R=L.R(t,x); Rp=(L.R(t,x+dx)-L.R(t,x-dx))/(2*dx)/np.exp(x); return np.log(R**2*Rp/np.exp(x)**2)
def cmb(L,tn,ells=np.arange(0,360,6.),tmode='A',RD=None,budget=40000):
    """tmode A: every element starts equally hot at the common origin and cools with its OWN volume expansion, T=K*(J_obs,today/J)^(1/3),
    J=R^2 R'/r0^2 (forward.py A4 convention).  tmode B: uniform photon/baryon ratio, T=K*(rho/rho_obs,today)^(1/3) (cmb4 convention)."""
    o=obs(L,tn); X0=np.array([0.]); lJo=lnJ(L,tn,X0)[0]
    if RD is None:
      RD=[]
      for e in ells:
        r=ray(L,tn,np.radians(e),rtol=1e-9,budget=budget)
        r['lr']=np.log(L.rho_rel(r['t'],r['x'])/o['rho']) if len(r['t']) else np.array([]); r['lJ']=lJo-lnJ(L,r['t'],r['x']) if len(r['t']) else np.array([]); RD.append(r)
    def tstar(K):
        f=lambda t: (np.log(L.rho_rel(t,X0)[0]/o['rho']) if tmode=='B' else lJo-lnJ(L,t,X0)[0])/3+np.log(K/3000.)
        return brentq(f,1e-9*tn,tn) if f(1e-9*tn)>0>f(tn) else None
    def emit(K,ts):
        rows=[]
        for r in RD:
            LT=r['lr'] if tmode=='B' else r['lJ']; Tel=K*np.exp(np.where(np.isfinite(LT),LT,-np.inf)/3)
            row=dict(ell=float(np.degrees(r['th'])),stop=r['stop'])
            if len(r['t'])<2: row['status']='stopped:'+r['stop']; rows.append(row); continue
            ok=r['t']>=ts; hot=(Tel>=3000)&ok&np.isfinite(r['lr'])&(r['lr']>LRMIN)
            if hot.any():
                k=int(np.argmax(hot)); 
                if k==0: row['status']='observer_in_plasma'; rows.append(row); continue
                w=(np.log(3000)-np.log(Tel[k-1]))/(np.log(Tel[k])-np.log(Tel[k-1])); g=lambda a:a[k-1]+w*(a[k]-a[k-1])
                row.update(status='plasma_edge',t_e=g(r['t']),lnz1=g(r['lnz1']),x_e=g(r['x']),phi=g(r['phi']),nc=int(r['nc'][k]),T_emit=3000.)
            elif r['t'][-1]<=ts:
                k=int(np.argmax(r['t']<ts)); w=(ts-r['t'][k-1])/(r['t'][k]-r['t'][k-1]); g=lambda a:a[k-1]+w*(a[k]-a[k-1])
                if np.isfinite(g(r['lr'])) and g(r['lr'])>LRMIN:
                    row.update(status='tstar_limit',t_e=ts,lnz1=g(r['lnz1']),x_e=g(r['x']),phi=g(r['phi']),nc=int(r['nc'][k]),T_emit=float(g(Tel)))
                else:   # at t* the ray is outside the matter; last matter point (t>=t*) along the ray is the emitter
                    mt=np.where(ok&np.isfinite(r['lr'])&(r['lr']>LRMIN))[0]
                    if len(mt): k=mt[-1]; row.update(status='matter_edge',t_e=r['t'][k],lnz1=r['lnz1'][k],x_e=r['x'][k],phi=r['phi'][k],nc=int(r['nc'][k]),T_emit=float(Tel[k]))
                    else: row['status']='no_matter_on_ray'
            elif not (np.isfinite(r['lr'][-1]) and r['lr'][-1]>LRMIN):   # ray ran out of the matter before t*: last matter point emits
                mt=np.where(ok&np.isfinite(r['lr'])&(r['lr']>LRMIN))[0]
                if len(mt): k=mt[-1]; row.update(status='matter_edge',t_e=r['t'][k],lnz1=r['lnz1'][k],x_e=r['x'][k],phi=r['phi'][k],nc=int(r['nc'][k]),T_emit=float(Tel[k]))
                else: row['status']='no_matter_on_ray'
            else: row['status']='stopped:'+r['stop']
            if 'lnz1' in row: row['T_obs']=row['T_emit']*np.exp(-row['lnz1'])
            rows.append(row)
        return rows
    def mean(K):
        ts=tstar(K)
        if ts is None: return np.nan,None
        rows=emit(K,ts); Tt=[q['T_obs'] for q in rows if 'T_obs' in q]; return (np.mean(Tt) if Tt else np.nan),rows
    lk=np.linspace(np.log(3000.)-25,np.log(3000.)-1e-6,120); mv=np.array([mean(np.exp(v))[0] for v in lk])
    good=np.isfinite(mv); res=dict(K=None)
    sgn=np.where(good[:-1]&good[1:]&((mv[:-1]-2.725)*(mv[1:]-2.725)<=0))[0]
    if not len(sgn): res['K_scan']=list(zip(np.exp(lk[good]).tolist(),mv[good].tolist())); return res,RD
    i=sgn[0]; K=np.exp(brentq(lambda v: mean(np.exp(v))[0]-2.725,lk[i],lk[i+1]))
    _,rows=mean(K); res.update(K=K,tstar=tstar(K),rows=rows,tmode=tmode); return res,RD
def score(rows,l0,ells):
    ok=[q for q in rows if 'T_obs' in q]; n=len(ok); out=dict(n_emit=n,n_total=len(rows))
    out['status_counts']={s:sum(q['status']==s for q in rows) for s in set(q['status'] for q in rows)}
    out['beyond_CoE_emitters']=sum(q.get('nc',0)>0 for q in ok)
    Tt=np.array([q['T_obs'] for q in ok]); out['T_range']=[float(Tt.min()),float(Tt.max())]; out['T_mean']=float(Tt.mean())
    el=np.array([q['ell'] for q in ok]); gl=(el+l0)%360
    Gt=np.interp(gl,np.arange(361.),np.append(CC.G_FULL,CC.G_FULL[0])); th=np.radians(gl); Xp=np.c_[np.ones(n),np.cos(th),np.sin(th)]
    def rmd(v): c=np.linalg.lstsq(Xp,v,rcond=None)[0]; return v-Xp[:,1:]@c[1:]
    m=rmd(Tt/Tt.mean()-1); b=rmd(Gt/Gt.mean()-1); R=rmd(Gt*Tt/np.mean(Gt*Tt)-1); rms=lambda v: float(np.sqrt(np.mean(v**2)))
    out.update(dTT_rms_m1rem=rms(m),target_rms_m1rem=rms(b),resid_rms=rms(R),cancel_pct=100*(1-rms(R)/rms(b)),corr_with_Gm1=float(np.corrcoef(m,b)[0,1]),
       gal_l=gl.tolist(),dTT=m.tolist(),Gm1=b.tolist(),resid=R.tolist())
    te={q['ell']:q['t_e'] for q in ok}
    if 0.0 in te:
        far=max(te,key=lambda e: -abs(e-180)); out['order']=dict(t0=te[0.0],ell_coe_side=far,t_coe_side=te[far],zero_first=bool(te[0.0]<te[far]))
    return out
