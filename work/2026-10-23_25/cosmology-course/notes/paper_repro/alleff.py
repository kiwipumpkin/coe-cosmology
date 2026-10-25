# CMB with ALL effects per direction, Chris's split-annulus picture. arg: 1 = fitted mass, 0 = no mass (same launch)
import numpy as np,json,sys
from scipy.optimize import brentq
MUSCALE=float(sys.argv[1])
s=open('cmbann.py').read().split("meanT=lambda")[0]
s=s.replace("exec(src)","src=src.replace(\"mu,MT,RG,S=rr['mu'],rr['MT'],pc.RG,rr['S']\",\"mu,MT,RG,S=rr['mu']*%r,rr['MT'],pc.RG,rr['S']\"); exec(src)"%MUSCALE)
exec(s)
b=np.arange(360.); br=np.radians(b); cb=np.cos(br)     # b=0: anti-CoE (l=15.2), b=180: toward CoE
D1=Dt(1.0); Vo=np.array([v_of_r(D1,1.0),0.0])          # observer velocity now (radial, +x away from CoE)
def GMr(r): return mu*np.interp(r,RG,MT)/r             # G M(<r)/r, code units (same acc as the code)
def te_of_q(q,Rr): # time layer q crosses R_rec
    j=np.searchsorted(QG,q); w=(q-QG[j-1])/(QG[j]-QG[j-1])
    return np.interp(Rr,XR[:,j-1]*(1-w)+XR[:,j]*w,T)
def calc(Rr):
    tr=np.interp(Rr,XR[:,i1],T); Rn=C*(1-tr)
    qi,qo=q_of_r(D1-Rn,1.0),q_of_r(D1+Rn,1.0)
    V=lambda a,c,t: r_of(c,t)**3-r_of(a,t)**3
    si=(V(qi,1.0,1.0)/V(qi,1.0,tr))**(1/3); so=(V(1.0,qo,1.0)/V(1.0,qo,tr))**(1/3)
    S=si*(1-cb)/2+so*(1+cb)/2
    n=np.c_[cb,np.sin(br)]                                  # sky direction (from us)
    P=np.array([D1,0])+Rn*n; rP=np.hypot(*P.T); uP=P/rP[:,None]   # emitting layer, where it is now
    q=np.array([q_of_r(x,1.0) for x in rP]); te=np.array([te_of_q(x,Rr) for x in q])
    ve=np.array([v_of_r(Rr,t) for t in te])                 # its speed when it crossed R_rec
    brel=(ve[:,None]*uP-Vo)/C                               # relative velocity / c
    blos=np.sum(brel*n,1)                                   # + = moving away from us along sight line
    b2=np.sum(brel**2,1); gam=np.where(b2<1,1/np.sqrt(np.clip(1-b2,1e-300,None)),np.nan)
    dop1=1+blos                                             # speed effect, first order (always defined)
    dopR=gam*(1+blos)                                       # with the relativistic time-slowing (gamma)
    xe=2*GMr(Rr)/C**2; xo=2*GMr(D1)/C**2                    # 2GM/(r c^2) at emitter / observer
    grav=np.sqrt(1-xo)/np.sqrt(1-xe) if (xe<1 and xo<1) else np.nan  # T_obs factor = sqrt(1-xe)/sqrt(1-xo) -> 1/grav
    delay_avail=1-te; delay_need=Rn/C*np.ones_like(te)      # code time units (t0)
    return dict(tr=tr,Rn=Rn,si=si,so=so,S=S,te=te,ve=ve,brel=brel,blos=blos,gam=gam,dop1=dop1,dopR=dopR,xe=xe,xo=xo,grav=grav,
                da=delay_avail,dn=delay_need)
def Tobs(Rr,which):
    c=calc(Rr); f=c['S']*(c['dopR'] if which=='full' else c['dop1'] if which=='dop1' else 1)
    if which=='full' and np.isfinite(c['grav']): f=f*c['grav']
    return 3000/f
def solve(which):
    g=lambda x: np.nanmean(Tobs(np.exp(x),which))-2.725
    return np.exp(brentq(g,np.log(0.0005/DG),np.log(2/DG),xtol=1e-10))
def stats(Tv,step):
    k=(b%step==0); Tv=Tv[k]; bb=br[k]; dT=(Tv-Tv.mean())/Tv.mean()
    A=np.c_[np.ones_like(bb),np.cos(bb),np.sin(bb)]; c,*_=np.linalg.lstsq(A,dT,rcond=None)
    return dict(full_amp_1e5=float((dT.max()-dT.min())/2/1e-5),dipole_1e5=float(np.hypot(c[1],c[2])/1e-5),
                hot_l=float((np.degrees(np.arctan2(c[2],c[1]))+15.2)%360),leftover_rms_1e5=float((dT-A@c).std()/1e-5))
out=dict(mu=float(mu),D_now_Gpc=float(D1*DG),obs_speed_c=float(Vo[0]/C))
R={}
for which in ('stretch','dop1','full'):
    try: Rr=solve(which)
    except Exception as e: out[which]=dict(error=repr(e)); continue
    c=calc(Rr); Tv=Tobs(Rr,which); R[which]=(Rr,Tv,c)
    d=dict(R_rec_Mpc=float(Rr*DG*1e3),t_rec_Myr=float(c['tr']*t0*1e3),mean_T=float(np.nanmean(Tv)),
        finite_dirs=int(np.isfinite(Tv).sum()),
        stretch_inner=float(c['si']),stretch_outer=float(c['so']),
        speed_factor_first_order=[float(c['dop1'].min()),float(c['dop1'].max())],
        gamma_defined_dirs=int(np.isfinite(c['gam']).sum()),
        rel_speed_c_toward_CoE=float(np.hypot(*c['brel'][180])),rel_speed_c_away_from_CoE=float(np.hypot(*c['brel'][0])),
        los_speed_c_toward_CoE=float(c['blos'][180]),los_speed_c_away=float(c['blos'][0]),
        emitter_speed_c=[float(c['ve'].min()/C),float(c['ve'].max()/C)],
        twoGM_rc2_emitter=float(c['xe']),twoGM_rc2_observer=float(c['xo']),gravity_factor=None if not np.isfinite(c['grav']) else float(c['grav']),
        delay_needed_Gyr=float(c['dn'][0]*t0),delay_available_Gyr=[float(c['da'].min()*t0),float(c['da'].max()*t0)],
        emit_time_Myr=[float(c['te'].min()*t0*1e3),float(c['te'].max()*t0*1e3)])
    if np.all(np.isfinite(Tv)): d['6deg']=stats(Tv,6); d['1deg']=stats(Tv,1)
    # size of each effect on the sky pattern (1e-5): its own half-swing
    for nm,f in (('stretch',c['S']),('speed_first_order',c['dop1']),('gamma',c['gam'])):
        f=np.asarray(f,float); d['swing_1e5_'+nm]=float((f.max()-f.min())/2/f.mean()/1e-5) if np.all(np.isfinite(f)) else None
    out[which]=d
assert abs(out['stretch']['R_rec_Mpc']-json.load(open('cmb_mass%g.json'%MUSCALE))['R_rec_Mpc'])<0.05  # matches earlier stretch-only run
np.save('alleff%g.npy'%MUSCALE,np.array([b]+[R[k][1] for k in ('stretch','dop1','full') if k in R]))
json.dump(out,open('alleff%g.json'%MUSCALE,'w'),indent=1); print(json.dumps(out,indent=1))
