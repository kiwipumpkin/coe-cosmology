# CMB vs angle with relativistic Doppler on PECULIAR velocity only. arg: 1 mass, 0 no mass
import numpy as np,json,sys
from scipy.optimize import brentq
MUSCALE=float(sys.argv[1])
s=open('cmbann.py').read().split("meanT=lambda")[0]
s=s.replace("exec(src)","src=src.replace(\"mu,MT,RG,S=rr['mu'],rr['MT'],pc.RG,rr['S']\",\"mu,MT,RG,S=rr['mu']*%r,rr['MT'],pc.RG,rr['S']\"); exec(src)"%MUSCALE)
exec(s)
b=np.arange(360.); br=np.radians(b); cb=np.cos(br); n=np.c_[cb,np.sin(br)]   # b=0 anti-CoE (l=15.2)
D1=Dt(1.0)
def te_of_q(q,Rr):
    j=np.searchsorted(QG,q); w=(q-QG[j-1])/(QG[j]-QG[j-1]); return np.interp(Rr,XR[:,j-1]*(1-w)+XR[:,j]*w,T)
def Hchoice(kind,t,qi,qo):
    D=Dt(t); e=1e-4*D
    if kind=='local_dvdr_at_us': return (v_of_r(D+e,t)-v_of_r(D-e,t))/(2*e)
    if kind=='us_v_over_D':      return v_of_r(D,t)/D
    if kind=='annulus_volume':   # (1/3) dlnV/dt of our annulus D-R..D+R (same layers as the stretch)
        V=lambda tt: r_of(qo,tt)**3-r_of(qi,tt)**3; h=1e-4*t; return (np.log(V(t+h))-np.log(V(t-h)))/(2*h)/3
HK=('local_dvdr_at_us','us_v_over_D','annulus_volume')
def calc(Rr,kind):
    tr=np.interp(Rr,XR[:,i1],T); Rn=C*(1-tr); qi,qo=q_of_r(D1-Rn,1.0),q_of_r(D1+Rn,1.0)
    V=lambda a,c,t: r_of(c,t)**3-r_of(a,t)**3
    si=(V(qi,1.0,1.0)/V(qi,1.0,tr))**(1/3); so=(V(1.0,qo,1.0)/V(1.0,qo,tr))**(1/3); S=si*(1-cb)/2+so*(1+cb)/2
    P=np.array([D1,0])+Rn*n; rP=np.hypot(*P.T); uP=P/rP[:,None]
    te=np.array([te_of_q(q_of_r(x,1.0),Rr) for x in rP])
    xe=Rr*uP; ve=np.array([v_of_r(Rr,t) for t in te])[:,None]*uP          # emitter at emission
    xo=np.c_[[Dt(t) for t in te],0*te]; vo=np.c_[[v_of_r(Dt(t),t) for t in te],0*te]  # observer, same moment
    H=np.array([Hchoice(kind,t,qi,qo) for t in te])
    pe=(ve-vo-H[:,None]*(xe-xo))/C          # emitter peculiar velocity / c, relative to growth flow centred on us
    po=0.0                                   # observer has zero peculiar velocity in a flow centred on itself
    blos=np.sum(pe*n,1); b2=np.sum(pe**2,1)   # + = away from us along sight line
    ok=bool(np.all(b2<1)); gam=1/np.sqrt(np.where(b2<1,1-b2,np.nan))
    dop=gam*(1+blos)
    pe_=float(phi(Rr)/C**2); po_=float(phi(D1)/C**2) if mu>0 else 0.0
    if mu==0: pe_=0.0
    weak=max(abs(pe_),abs(po_))<0.1; grav=(1+pe_-po_) if weak else np.nan   # T_obs divided by this
    return dict(tr=tr,Rn=Rn,si=si,so=so,S=S,te=te,H=H,pe=pe,blos=blos,b=np.sqrt(b2),ok=ok,gam=gam,dop=dop,
                phie=pe_,phio=po_,weak=weak,grav=grav,dn=Rn/C,da=1-te)
def Tobs(Rr,kind,use_dop=True):
    c=calc(Rr,kind); f=c['S']*(c['dop'] if use_dop else 1)
    if c['weak']: f=f*c['grav']
    return 3000/f
def stats(Tv,step):
    k=(b%step==0); Tv=Tv[k]; bb=br[k]; dT=(Tv-Tv.mean())/Tv.mean()
    A=np.c_[np.ones_like(bb),np.cos(bb),np.sin(bb)]; c,*_=np.linalg.lstsq(A,dT,rcond=None)
    return dict(full_amp_1e5=float((dT.max()-dT.min())/2/1e-5),dipole_1e5=float(np.hypot(c[1],c[2])/1e-5),
                hot_l=float((np.degrees(np.arctan2(c[2],c[1]))+15.2)%360),leftover_rms_1e5=float((dT-A@c).std()/1e-5))
sw=lambda f:float((np.max(f)-np.min(f))/2/np.mean(f)/1e-5)
out=dict(mu=float(mu),D_now_Gpc=float(D1*DG)); curves=[b]
for kind in HK:
    c0=calc(json.load(open('cmb_mass%g.json'%MUSCALE))['R_rec_Mpc']/DG/1e3,kind)
    d=dict(peculiar_below_c_at_stretch_Rrec=c0['ok'],max_peculiar_c=float(np.nanmax(c0['b'])))
    if not c0['ok']: out[kind]=d; continue
    try: Rr=np.exp(brentq(lambda x: np.mean(Tobs(np.exp(x),kind))-2.725,np.log(0.0005/DG),np.log(2/DG),xtol=1e-10))
    except Exception as e: d['error']=repr(e); out[kind]=d; continue
    c=calc(Rr,kind); assert c['ok'],'peculiar speed >= c'; Tv=Tobs(Rr,kind); assert abs(Tv.mean()-2.725)<1e-6
    d.update(R_rec_Mpc=float(Rr*DG*1e3),t_rec_Myr=float(c['tr']*t0*1e3),growth_rate_per_Gyr=[float(c['H'].min()/t0),float(c['H'].max()/t0)],
      stretch_inner=float(c['si']),stretch_outer=float(c['so']),swing_1e5_stretch=sw(c['S']),
      peculiar_c_toward_CoE=float(c['b'][180]),peculiar_c_away_from_CoE=float(c['b'][0]),
      los_peculiar_c_toward_CoE=float(c['blos'][180]),los_peculiar_c_away=float(c['blos'][0]),max_peculiar_c=float(c['b'].max()),
      doppler_range=[float(c['dop'].min()),float(c['dop'].max())],swing_1e5_doppler=sw(c['dop']),
      phi_c2_emitter=c['phie'],phi_c2_observer=c['phio'],gravity_weak_field_ok=c['weak'],
      gravity_factor=None if not c['weak'] else float(c['grav']),
      delay_needed_Gyr=float(c['dn']*t0),delay_available_Gyr=[float(c['da'].min()*t0),float(c['da'].max()*t0)])
    d['6deg']=stats(Tv,6); d['1deg']=stats(Tv,1); out[kind]=d; curves.append(Tv)
    out[kind+'_curve_index']=len(curves)-1
np.save('pec%g.npy'%MUSCALE,np.array(curves)); json.dump(out,open('pec%g.json'%MUSCALE,'w'),indent=1); print(json.dumps(out,indent=1))
