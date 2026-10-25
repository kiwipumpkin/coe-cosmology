"""CMB v4: recombination from the RADIATION temperature carried by each element.
 Each dust element carries thermal radiation (tightly coupled before recombination) that cools with that element's own
 volume expansion: T_elem(t) = K * (rho_elem(t)/rho_obs,today)^{1/3}   (adiabatic photon gas, V ~ 1/rho, uniform photon/baryon ratio).
 Recombination: first point (backward along each null ray, i.e. latest emission) where T_elem = 3000 K.
 K is NOT assumed = 2.725 K; it is SOLVED so that the sky-mean observed temperature <3000/(1+z)> = 2.725 K (requirement T0).
usage: cmb4.py <law> <fwd|back> <gtc|le> [dstep]   env XOVR (param vector JSON) OUTTAG PAGE=1"""
import os,sys,json,time,numpy as np,warnings;warnings.filterwarnings('ignore')
law,mode,vm=sys.argv[1:4]; DST=float(sys.argv[4]) if len(sys.argv)>4 else 6.
os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
else: os.environ.pop('VMODE',None)
sys.argv=['x',law,mode]
exec(open('fit_lc.py').read().split("ck=f'fwd")[0])
import grltb as G, ltb2 as T
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
tag=f'{law}_{mode}_{vm}'
x=np.array(json.loads(os.environ['XOVR'])) if os.environ.get('XOVR') else np.array(json.load(open(f'gr/ck_{tag}.json'))['x'])
OT=os.environ.get('OUTTAG',tag); OUTD=os.environ.get('OUTD','gr2/v4')
os.makedirs(OUTD,exist_ok=True)
def bail(msg): json.dump(dict(x=list(x),tag=tag,failed=msg,names=nm),open(f'{OUTD}/{OT}.json','w')); print('FAIL',msg); sys.exit()
try: P,o=build(x)
except Exception as e: bail('build exception '+repr(e))
if o is None: bail('build rejected (unphysical/non-integrable)')
io=o['io']; xo=np.log(F.R0[io])
if mode=='fwd': L=T.LTB(P,'fwd'); tn=o['t']*F.C
else: L=T.LTB(P,'back',r_ref_obs=F.R0[io],t_ref=0.0); tn=0.0
X=np.array([xo]); Ro=L.R(tn,X)[0]; Uo=L.Rdot(tn,X)[0]; tbo=L.tB(X)[0]; rho_o=L.rho_rel(tn,X)[0]
TB_MIN=float(np.nanmin(L.tB(np.linspace(np.log(1e-3)+1,np.log(3e4)-1,3000))))+1e-9*abs(tn)
f=G.fields(o); tr=G.trace(f,np.arange(0,360,4.),nstep=320,lbmax_frac=0.999); af=G.migkas(tr); chi,l0=F.chi2_l0(af)
xi=np.linspace(np.log(1e-3)+1,np.log(3e4)-1,4000); inter=L.interior(xi); x_int=xi[inter].max() if inter.any() else -np.inf
lam=F.local(o)['Hbar']/H70
print(tag,OT,'chi2',chi,'l0',l0,'R_o',Ro,'U_o',Uo,flush=True)
def lnrho(t,y):
    v=L.rho_rel(t,np.array([y[0]]))[0]; return np.log(v/rho_o) if np.isfinite(v) and v>0 else (1e3 if np.isfinite(L.R(t,np.array([y[0]]))[0]) else -1e3)
def trace(th):
    """integrate the past null ray all the way back (no recombination event); return dense samples (t,x,phi,lnz1,lnrho)."""
    Lang=-Ro*np.sin(th); fr=T.rhs(L,Lang); y=np.array([xo,0.,0.,-np.cos(th)]); t=tn; t1=TB_MIN
    e2=lambda t,y: y[0]-(L.xmin+1.0); e2.terminal=True; e2.direction=-1
    e3=lambda t,y: y[0]-x_int-1e-3 if np.isfinite(x_int) else 1.0; e3.terminal=True; e3.direction=-1
    segs=[]; nc=0; stop='reached_TB_MIN'
    for _ in range(4):
        try: s=solve_ivp(fr,(t,t1),y,method='LSODA',rtol=1e-10,atol=1e-12,events=(e2,e3),dense_output=True)
        except ValueError:
            s=solve_ivp(fr,(t,t1),y,method='LSODA',rtol=1e-10,atol=1e-12,dense_output=True); s.status=0; s.t_events=[[],[]]
        fin=np.isfinite(s.y).all(0); tt=s.t[fin]
        # densify between solver steps
        if len(tt)>1:
            q=np.unique(np.concatenate([np.linspace(tt[k],tt[k+1],6) for k in range(len(tt)-1)])) [::-1]
            Y=np.array([s.sol(v) for v in q]).T; segs.append((q,Y,nc))
        if s.status==1 and len(s.t_events[0]):
            t=s.t_events[0][0]; y=s.sol(t).copy(); y[1]+=np.pi; y[3]=-y[3]; nc+=1; t+=(t1-t)*1e-12; continue
        if s.status==1 and len(s.t_events[1]): stop='trapped'
        elif not fin.all() or s.status<0: stop='integrator_stop'
        break
    tt=np.concatenate([g[0] for g in segs]); Y=np.concatenate([g[1] for g in segs],1); NC=np.concatenate([np.full(len(g[0]),g[2]) for g in segs])
    lr=np.array([lnrho(tt[k],Y[:,k]) for k in range(len(tt))])
    return dict(t=tt,x=Y[0],phi=Y[1],lnz1=Y[2],lr=lr,nc=NC,stop=stop)
CAUS=os.environ.get('CAUSTIC','opaque')
def recpoint(rd,lnth):
    lr=rd['lr'].copy()
    if CAUS=='transparent':
        m=lr<1e2; lr=lr[m]; rd={k:(v[m] if isinstance(v,np.ndarray) else v) for k,v in rd.items()}
    k=np.where((lr[:-1]<lnth)&(lr[1:]>=lnth))[0]
    if not len(k): return None
    k=k[0]; w=(lnth-lr[k])/(lr[k+1]-lr[k]) if np.isfinite(lr[k+1]-lr[k]) and lr[k+1]<1e2 else 0.0
    g=lambda a:a[k]+w*(a[k+1]-a[k])
    return dict(t=g(rd['t']),x=g(rd['x']),phi=g(rd['phi']),lnz1=g(rd['lnz1']),nc=int(rd['nc'][k]))
t0=time.time(); ells=np.arange(0,360,DST); RD=[]
for e in ells: RD.append(trace(np.radians(e)))
print('traced %d rays %.0fs'%(len(ells),time.time()-t0),flush=True)
def sky(K):
    lnth=3*np.log(3000/K); return [recpoint(r,lnth) for r in RD]
def meanT(K):
    p=sky(K); T_=[3000*np.exp(-q['lnz1']) for q in p if q is not None]; return (np.mean(T_) if T_ else np.nan), sum(q is not None for q in p)
# solve K: scan log K, then brentq on f=ln(meanT/2.725) restricted to K where every direction recombines (else max coverage)
Ks=np.logspace(-8,6,141); vals=[meanT(K) for K in Ks]; cov=np.array([v[1] for v in vals]); mT=np.array([v[0] for v in vals])
N=len(ells); full_=cov==cov.max(); fk=np.log(mT/2.725)
Ksol=None; note=''
idx=[i for i in range(len(Ks)-1) if np.isfinite(fk[i]) and np.isfinite(fk[i+1]) and fk[i]*fk[i+1]<=0]
idx=sorted(idx,key=lambda i:-min(cov[i],cov[i+1]))
if idx:
    i=idx[0]; Ksol=brentq(lambda K:np.log(meanT(K)[0]/2.725),Ks[i],Ks[i+1],xtol=1e-12*Ks[i]); note='solved'
else:
    ok=np.isfinite(fk)&full_
    if ok.any(): i=np.where(ok)[0][np.argmin(abs(fk[ok]))]; Ksol=Ks[i]; note=f'no root with max coverage; closest <T>={mT[i]:.4g} K'
out=dict(tag=tag,outtag=OT,names=nm,x=list(x),chi2=chi,sigma=sig(chi,len(nm)+1),l0=l0,Uo_c=float(Uo),D_Gpc=float(Ro*lam),age_obs_Gyr=float((tn-tbo)/F.C*lam),
  M=float(P.get('M',np.nan)),K_note=note,K_scan=dict(K=Ks.tolist(),meanT=mT.tolist(),coverage=cov.tolist()),n_dir=N,caustic_mode=CAUS,ray_zmax=[float(np.exp(np.nanmax(r['lnz1']))-1) for r in RD],ray_lnrho_max_finite=[float(r['lr'][(r['lr']<1e2)&(r['lr']>-1e2)].max()) if ((r['lr']<1e2)&(r['lr']>-1e2)).any() else None for r in RD],ray_caustic=[bool((r['lr']>=1e2).any()) for r in RD],ray_stop={s_:sum(r['stop']==s_ for r in RD) for s_ in set(r['stop'] for r in RD)})
out['migkas_model']=[float(v) for v in F.shifted(af,l0)]; out['migkas_data']=dict(l=[float(v) for v in F.IMIG],y=[float(v) for v in F.Y],err=[float(v) for v in 0.5*(F.MLO+F.MHI)])
rows=[]
if Ksol is not None:
    p=sky(Ksol); out['K_local_today_K']=float(Ksol); out['coverage']=int(sum(q is not None for q in p))
    for e,q,rd in zip(ells,p,RD):
        r=dict(ell=float(e),status='rec' if q else 'no_rec('+rd['stop']+')')
        if q:
            xx=np.array([q['x']]); r.update(z=float(np.exp(q['lnz1'])-1),T_obs_K=float(3000*np.exp(-q['lnz1'])),t_e_Gyr=float((q['t']-tbo)/F.C*lam),
              flight_Gyr=float((tn-q['t'])/F.C*lam),R_e_over_D=float(L.R(q['t'],xx)[0]/Ro),phi_e_deg=float(np.degrees(q['phi'])),ncross=q['nc'],U_e_c=float(L.Rdot(q['t'],xx)[0]))
        rows.append(r)
    out['rows']=rows; ok=[r for r in rows if r['status']=='rec']
    if ok:
        Tt=np.array([r['T_obs_K'] for r in ok]); out['T_obs_mean_K']=float(Tt.mean()); out['T_obs_range_K']=[float(Tt.min()),float(Tt.max())]
        te={r['ell']:r['t_e_Gyr'] for r in ok}
        if 0.0 in te and 180.0 in te: out['t_rec_0_minus_180_Gyr']=te[0.0]-te[180.0]; out['order_0_before_180']=bool(te[0.0]<te[180.0])

    # partial-sky score (only directions that recombine; non-recombining = shadow/no-solution directions excluded and reported)
    if ok and len(ok)>=5:
        el=np.array([r['ell'] for r in ok]); gl=(el+l0)%360; Tt=np.array([r['T_obs_K'] for r in ok]); Gt=np.interp(gl,np.arange(361.),np.append(CC.G_FULL,CC.G_FULL[0]))
        thg=np.radians(gl); Xp=np.c_[np.ones(len(gl)),np.cos(thg),np.sin(thg)]
        def rmd(v): c=np.linalg.lstsq(Xp,v,rcond=None)[0]; return v-Xp[:,1:]@c[1:]
        Rp=rmd(Gt*Tt/np.mean(Gt*Tt)-1); bp=rmd(Gt/Gt.mean()-1); mp=rmd(Tt/Tt.mean()-1)
        out['score_partial']=dict(n=len(ok),res_1e5=float(np.sqrt(np.mean(Rp**2))/1e-5),cancel=float(100*(1-np.sqrt(np.mean(Rp**2))/np.sqrt(np.mean(bp**2)))),
            model_dTT_rms_m1rem=float(np.sqrt(np.mean(mp**2))),target_rms_m1rem=float(np.sqrt(np.mean(bp**2))),gal_l=gl.tolist(),dTT=(Tt/Tt.mean()-1).tolist(),resid=Rp.tolist())
    if len(ok)==N:
        th=np.radians(np.arange(360.)); Xm=np.c_[np.ones(360),np.cos(th),np.sin(th)]
        def rmdip(v): c=np.linalg.lstsq(Xm,v,rcond=None)[0]; return v-Xm[:,1:]@c[1:]
        S=np.array([r['T_obs_K'] for r in ok])/2.725
        S360=np.exp(np.interp(np.arange(360.),np.append(ells,360),np.log(np.append(S,S[0]))))
        So=np.interp((np.arange(360.)-l0)%360,np.arange(361.),np.append(S360,S360[0]))   # model sky placed on galactic longitude (CoE at l0+180)
        Gt=CC.G_FULL; base=Gt/Gt.mean()-1; R=Gt*So/np.mean(Gt*So)-1
        # (A) previous convention: G*S residual, no dipole removal anywhere
        A_=dict(res_1e5=CC.rms(R)/1e-5,cancel=100*(1-CC.rms(R)/CC.rms(base)))
        # (B) observed-CMB convention: remove m=1 from the OBSERVED product (our kinetic dipole) and from the baseline alike
        Rb=rmdip(R); bb=rmdip(base); B_=dict(res_1e5=CC.rms(Rb)/1e-5,cancel=100*(1-CC.rms(Rb)/CC.rms(bb)))
        # (C) orientation free (phi free instead of Migkas l0), convention B
        best=(1e9,None)
        for ph in np.arange(0,360,1.):
            Sp=np.interp((np.arange(360.)-ph)%360,np.arange(361.),np.append(S360,S360[0])); r_=rmdip(Gt*Sp/np.mean(Gt*Sp)-1); v=CC.rms(r_)
            if v<best[0]: best=(v,ph)
        out['score']=dict(A_nodipremoval=A_,B_m1removed=B_,C_m1removed_phifree=dict(res_1e5=best[0]/1e-5,cancel=100*(1-best[0]/CC.rms(bb)),phi=best[1]),
             model_dTT_rms=float(np.std(So/So.mean())),model_dTT_rms_m1removed=float(np.std(rmdip(So/So.mean()-1))),target_rms_m1removed=float(CC.rms(bb)),
             modes_model={int(m['m']):float(m['A']) for m in CC.fourier(So/So.mean()-1,6)['modes']},modes_target={int(m['m']):float(m['A']) for m in CC.fourier(base,6)['modes']})
        out['_curves']=dict(dTT=(So/So.mean()-1).tolist(),resid_B=Rb.tolist(),target_m1rem=bb.tolist(),inv_target=(-bb).tolist())
json.dump(out,open(f'{OUTD}/{OT}.json','w'),indent=1,default=float)
print('DONE',OT,'chi2 %.2f'%chi,'K',out.get('K_local_today_K'),note,'cov',out.get('coverage'),'/',N,'T',out.get('T_obs_range_K'),out.get('score',{}).get('B_m1removed'),'%.0fs'%(time.time()-t0),flush=True)
