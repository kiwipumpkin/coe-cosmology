"""Observer-centred / common-origin rerun (exact LTB, ltb2 geodesics).
Change vs cmb4/ltb2 'fwd': there every shell starts at t=0 at its OWN radius R=r0 (extended start around the CoE),
which makes the bang time tB(x) = -T(x,r0) different for every shell.  Here tB(x) == 0 for ALL shells: every element
(us included) leaves R=0 at t=0, i.e. one common origin event.  Same m(x) (central M + profile) and same E(x) launch law.
Observer = same label x_o, placed at the time tn' where its areal distance to the CoE equals the old fit's R_o.
Migkas recomputed from ltb2 null rays in this spacetime (same estimator as grltb.migkas: h=z/lookback, ZMIG ring).
usage: oc.py <tag> '<x json>'   (6 deg rays)"""
import os,sys,json,numpy as np
os.chdir('/workspace/cosmology-course/notes/age_refit'); sys.path.insert(0,os.getcwd())
TAG=OTAG=sys.argv[1]; os.environ['XOVR']=sys.argv[2]; os.environ['OUTD']='obs_centred'; os.environ['OUTTAG']=TAG
sys.argv=['cmb4.py','eta','fwd','gtc','6']
src=open('cmb4.py').read(); head,tail=src.split("def lnrho")
exec(head)
chi_old,l0_old,Ro_old,Uo_old,tn_old,tbo_old=chi,l0,Ro,Uo,tn,tbo
import ltb2 as T
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
class LTB0(T.LTB):
    def tB(s,x): return np.zeros_like(np.atleast_1d(x),dtype=float)
L=LTB0(P,'fwd')
g=lambda t: L.R(t,X)[0]-Ro_old
tt=np.geomspace(1e-8*tn_old,50*tn_old,400); gv=np.array([g(t) for t in tt]); k=np.where(np.isfinite(gv[:-1])&(gv[:-1]<0)&(gv[1:]>=0))[0][0]
tn=brentq(g,tt[k],tt[k+1],xtol=1e-14*tn_old)
Ro=L.R(tn,X)[0]; Uo=L.Rdot(tn,X)[0]; tbo=0.0; rho_o=L.rho_rel(tn,X)[0]
TB_MIN=1e-9*tn
inter=L.interior(xi); x_int=xi[inter].max() if inter.any() else -np.inf
# Migkas in the new spacetime
def zray(th):
    Lang=-Ro*np.sin(th); fr=T.rhs(L,Lang); y=[xo,0.,0.,-np.cos(th)]
    ev=lambda t,y: y[2]-np.log(1.4); ev.terminal=True
    s=solve_ivp(fr,(tn,TB_MIN),y,method='LSODA',rtol=1e-10,atol=1e-12,events=ev,dense_output=True)
    q=np.linspace(tn,s.t[-1],600); Z=np.exp(np.array([s.sol(v)[2] for v in q]))-1; return tn-q,Z
ells_m=np.arange(0,360,6.); afs=[]
ZR=[zray(np.radians(e)) for e in ells_m]
ok=True
for zk in F.ZMIG:
    d=[]
    for lb,Z in ZR:
        s_=np.where((Z[:-1]-zk)*(Z[1:]-zk)<=0)[0]
        if not len(s_): ok=False; d.append(np.nan); continue
        j=s_[0]; w=(zk-Z[j])/(Z[j+1]-Z[j]); d.append(lb[j]+w*(lb[j+1]-lb[j]))
    h=zk/np.array(d); afs.append(h/np.nanmean(h))
af=np.mean(afs,0); af=np.interp(np.arange(360.),np.append(ells_m,360),np.append(af,af[0])); af/=af.mean()
chi,l0=F.chi2_l0(af)
print('new tn',tn,'old',tn_old,'Uo',Uo,'old',Uo_old,'chi2 new',chi,'old',chi_old,flush=True)
tail=tail.replace("out=dict(tag=tag","out=dict(chi2_oldspacetime=chi_old,l0_old=l0_old,tn_code=tn,tn_old_code=tn_old,tbo_old_code=float(tbo_old),Uo_old_c=float(Uo_old),tag=tag")
exec("def lnrho"+tail)
json.dump(dict(ZR=[[a.tolist(),b.tolist()] for a,b in ZR]),open(f'obs_centred/{OTAG}_zrays.json','w'))
# ---- Chris's construction (sync): t*=our own 3000 K crossing; today's CMB = past-light-cone points at t=t* ----
def sync(K):
    lr_us=lambda t: np.log(L.rho_rel(t,X)[0]/rho_o)
    lnth=3*np.log(3000/K)
    a_,b_=TB_MIN*10,tn
    if not (lr_us(a_)>lnth>lr_us(b_)): return np.nan,[]
    ts_=brentq(lambda t: lr_us(t)-lnth,a_,b_)
    Ts=[];rows_=[]
    for e,rd in zip(ells,RD):
        t_=rd['t']; 
        if not (t_.min()<=ts_): rows_.append(dict(ell=float(e),status='ray_never_reaches_t*')); continue
        o_=np.argsort(t_); ip=lambda a: float(np.interp(ts_,t_[o_],a[o_]))
        lz=ip(rd['lnz1']); lr_=ip(rd['lr']); Te=K*np.exp(lr_/3); To=Te*np.exp(-lz)
        rows_.append(dict(ell=float(e),z=float(np.exp(lz)-1),T_emit_K=float(Te),recombined=bool(Te<=3000),T_obs_K=float(To),
             R_e_over_D=float(L.R(ts_,np.array([ip(rd['x'])]))[0]/Ro),phi_e_deg=float(np.degrees(ip(rd['phi'])))))
    return ts_,rows_
def smean(K):
    ts_,r_=sync(K); v=[r['T_obs_K'] for r in r_ if 'T_obs_K' in r]; return np.mean(v) if v else np.nan
Kg=np.logspace(-6,4,81); mv=np.array([smean(k) for k in Kg]); f_=np.log(mv/2.725)
ii=[i for i in range(80) if np.isfinite(f_[i]*f_[i+1]) and f_[i]*f_[i+1]<=0]
syn=dict(note='no K root')
if ii:
    Ks_=brentq(lambda k: np.log(smean(k)/2.725),Kg[ii[0]],Kg[ii[0]+1]); ts_,r_=sync(Ks_)
    okr=[r for r in r_ if 'T_obs_K' in r]; Tt=np.array([r['T_obs_K'] for r in okr]); el=np.array([r['ell'] for r in okr])
    gl=(el+l0)%360; Gt=np.interp(gl,np.arange(361.),np.append(CC.G_FULL,CC.G_FULL[0])); th=np.radians(gl); Xp=np.c_[np.ones(len(gl)),np.cos(th),np.sin(th)]
    rmd=lambda v: v-Xp[:,1:]@np.linalg.lstsq(Xp,v,rcond=None)[0][1:]
    Rp=rmd(Gt*Tt/np.mean(Gt*Tt)-1); bp=rmd(Gt/Gt.mean()-1); mp=rmd(Tt/Tt.mean()-1)
    syn=dict(K=Ks_,t_star_Myr=ts_/F.C*lam*1e3,n=len(okr),n_not_recombined=sum(not r['recombined'] for r in okr),
      dTT_rms_m1rem=float(np.sqrt(np.mean(mp**2))),target=float(np.sqrt(np.mean(bp**2))),res_rms=float(np.sqrt(np.mean(Rp**2))),
      cancel=float(100*(1-np.sqrt(np.mean(Rp**2))/np.sqrt(np.mean(bp**2)))),corr=float(np.corrcoef(mp,bp)[0,1]),rows=r_)
    print('SYNC',{k:v for k,v in syn.items() if k!='rows'},flush=True)
J=json.load(open(f'obs_centred/{OTAG}.json')); J['sync_chris']=syn; json.dump(J,open(f'obs_centred/{OTAG}.json','w'),indent=1,default=float)
import pickle; pickle.dump(dict(RD=RD,ells=ells,tn=tn,Ro=Ro,lam=lam,K=J.get('K_local_today_K'),tstar=syn.get('t_star_Myr')),open(f'obs_centred/{OTAG}_rays.pkl','wb'))
