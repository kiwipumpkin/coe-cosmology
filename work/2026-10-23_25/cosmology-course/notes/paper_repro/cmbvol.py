import numpy as np,json
src=open('steps.py').read(); exec(src[:src.index('import matplotlib')])
from scipy.integrate import solve_ivp
Dt=lambda t: r_of(1.0,t)
def rhs(t,y,n): r=np.hypot(*y); return (v_of_r(r,t)/r)*y + C*n
ev=lambda t,y,n: np.hypot(*y)-RREC; ev.terminal=True
dq=2e-3
def V(q,t): r=r_of(q,t); return r*r*(r_of(q+dq,t)-r_of(q-dq,t))      # V ∝ r^2 dr  (4π dropped)
PH=pc.phi_coe(np.array([1.0]),mu,MT)  # check callable
def phi(r): return float(pc.phi_coe(np.atleast_1d(r),mu,MT)[0])
def run(step):
    rows=[]
    for beta in np.arange(0,360,step):
        n=-np.array([np.cos(np.radians(beta)),np.sin(np.radians(beta))])
        s=solve_ivp(rhs,(1.0,T[5]),[Dt(1.0),0.0],args=(n,),events=ev,rtol=1e-10,atol=1e-13,max_step=2e-3)
        te=s.t_events[0][0]; xe=s.y_events[0][0]; qe=q_of_r(RREC,te)
        assert abs(r_of(qe,te)-RREC)<1e-8
        S=(V(qe,1.0)/V(qe,te))**(1/3)
        # emitter matter now (same layer, same direction) relative to us now
        u=xe/RREC; re=r_of(qe,1.0); Pe=re*u; Ve=v_of_r(re,1.0)*u; Vu=np.array([v_of_r(Dt(1.0),1.0),0])
        vr=(Ve-Vu)/C; los=(Pe-[Dt(1.0),0]); los/=np.hypot(*los); br=vr@los; b2=vr@vr
        dop=1+br; gam=1/np.sqrt(1-b2) if b2<1 else np.nan
        ge=1+2*phi(RREC)/C**2; go=1+2*phi(Dt(1.0))/C**2
        grav=np.sqrt(go/ge) if (ge>0 and go>0) else np.nan     # clock ratio, as letter's D_grav
        rows.append(dict(beta=beta,te=te,qe=qe,S=S,dop=dop,gam=gam,grav=grav,b=np.sqrt(b2),phiE_c2=phi(RREC)/C**2,phiO_c2=phi(Dt(1.0))/C**2))
    return rows
def stats(T,beta):
    b=np.radians(beta); A=np.c_[np.ones_like(b),np.cos(b),np.sin(b)]; c,*_=np.linalg.lstsq(A,T,rcond=None)
    res=T-A@c; return dict(mean=float(T.mean()),dip_amp=float(np.hypot(c[1],c[2])),dip_dir=float(np.degrees(np.arctan2(c[2],c[1]))),
        left_rms=float(res.std()),left_pp=float(res.max()-res.min()),left_rms_1e5=float(res.std()/T.mean()/1e-5))
out={}
for step in (1,6):
    R=run(step); b=np.array([r['beta'] for r in R]); g=lambda k: np.array([r[k] for r in R])
    S,dop,gam,grav=g('S'),g('dop'),g('gam'),g('grav')
    res={'terms':{k:[float(np.nanmin(v)),float(np.nanmean(v)),float(np.nanmax(v))] for k,v in dict(stretch=S,doppler=dop,gamma=gam,grav=grav,beta_rel=g('b')).items()},
         'phi_c2_emit':R[0]['phiE_c2'],'phi_c2_obs':R[0]['phiO_c2'],'te_Gyr':[float(g('te').min()*t0),float(g('te').max()*t0)]}
    for lab,z1 in (('full',S*dop*grav*gam),('no_grav_no_gamma',S*dop),('stretch_only',S)):
        Tt=3000/z1
        if np.all(np.isfinite(Tt)): res[lab]=stats(Tt,b)
        else: res[lab]='undefined (NaN terms)'
    lz=np.log(S*dop*grav*gam)
    if np.all(np.isfinite(lz)): res['share_lnz']=dict(stretch=float(np.mean(np.log(S))/lz.mean()),motion=float(np.mean(np.log(dop*gam))/lz.mean()),grav=float(np.mean(np.log(grav))/lz.mean()))
    out[f'{step}deg']=res
    if step==1: keep=(b,S,dop,gam,grav,g('qe'),g('te'))
out['needed_stretch']=3000/2.725
json.dump(out,open('cmbvol.json','w'),indent=1,default=float); print(json.dumps(out,indent=1,default=float))
np.save('cmbvol_keep.npy',np.array(keep))
