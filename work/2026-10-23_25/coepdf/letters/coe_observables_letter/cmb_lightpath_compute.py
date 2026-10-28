# Predicted CMB dT/T for the CoE model by tracing light from the recombination surface (no PNG; writes cmb_lightpath.json + cmb_lightpath_maps.npz).
# Pieces reused unchanged (see JSON 'sources'):
#   Newtonian integrator + eta density law: /workspace/cosmology/papers/spreadsheet_model/bubble/prof_core.py (kind 'eta', Ob=0)
#   CatWISE-locked Migkas fit at D=34.96 Gpc: /workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json
#   CatWISE lock L(D,eta,d): /workspace/cosmology-course/notes/catwise-eta-density-migkas.py (L_from_delta_eta, d=RV=14.3005 Gpc, delta=0.031)
#   Recombination rule (f=T0/Trec=t_rec/t_now, R_see=c(t_now-t_rec)): cmb_reconcile/cycle7_recombination + cycle4/dtt_rec_track/compute_track.py
#   Doppler/static conventions: cycle4/lib4.py  (dT/T_kin=-n.(V_e-V_o)/c linear; dT/T_grav=(Phi_e-Phi_o)/c^2)
#   Orientation: anti-CoE at (l,b)=(0,+15) (blended figure, make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py); field is axisymmetric so roll/gradient do not enter.
import os,sys,json,numpy as np,healpy as hp,importlib.util as iu
HERE=os.path.dirname(os.path.abspath(__file__))
BUB='/workspace/cosmology/papers/spreadsheet_model/bubble';sys.path.insert(0,BUB);cwd=os.getcwd();os.chdir(BUB)
import prof_core as pc
from bubble_core import RKM,GPC,TH0
os.chdir(cwd)
sp=iu.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=iu.module_from_spec(sp);sp.loader.exec_module(cw);os.chdir(cwd)
FIT=json.load(open('/workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json'))['34.96']
D,eta,mc=FIT['D'],FIT['eta'],FIT['mc'];RV=float(RKM/GPC)
L=cw.L_from_delta_eta(D,RV,eta,cw.DELTA_USE);LD=L/D
rr=pc.run('eta',dict(DR=D/RV,LD=LD,eta=eta),mc,0.0)          # same call as the fit (sky_locked)
mu,MT,Rb,Ht=rr['mu'],rr['MT'],rr['Rb'],rr['Ht'];t0=Ht*TH0     # Gyr
CK=299792.458;cG=CK*3.15576e16/3.0857e22                    # c in Gpc/Gyr
vc=D/t0/cG                                                  # code velocity unit D/t0 in units of c
T0K,TREC=2.725,3000.0;f=T0K/TREC                            # cycle7 rule
Rsee=cG*t0*(1-f)
TAUI,NST=0.01,200
th=np.linspace(0,180,1441)                                  # angle about observer from anti-CoE direction (integrator convention)
def trace(R):
    rho=R/D;q=np.empty((len(th)+1,2));q[0]=(1,0);q[1:,0]=1+rho*np.cos(np.radians(th));q[1:,1]=rho*np.sin(np.radians(th))
    X,V,VR=pc.integrate(q,mu,pc.RG,MT,0.0,Rb,TAUI,NST)
    Xo,Vo=X[0],V[0];d=X[1:]-Xo;sep=np.hypot(d[:,0],d[:,1]);n=d/sep[:,None]
    th_obs=np.degrees(np.arctan2(n[:,1],n[:,0]))
    phi=lambda r:pc.phi_coe(np.maximum(r,1e-12),mu,MT)
    # model state at recombination: before tau_i the integrator is free coasting, x=q*tau, v=q (its initial condition)
    xe_rec=q[1:]*f;ve_rec=q[1:]
    kin_rec=-np.sum(n*(ve_rec-Vo),1)*vc
    kin_now=-np.sum(n*(V[1:]-Vo),1)*vc                       # cycle4 'today field' variant
    re=np.hypot(xe_rec[:,0],xe_rec[:,1]);ro=np.hypot(*Xo)
    grav_rec=(phi(re)-phi(ro))*vc**2
    grav_now=(phi(np.hypot(X[1:,0],X[1:,1]))-phi(ro))*vc**2
    # straight-ray line integral of the static field from the emission event to us: -int grad(Phi).dl = Phi_e-Phi_o (check)
    s=np.linspace(0,1,4001)[:,None,None];P=xe_rec[None]+(Xo[None,None,:]-xe_rec[None])*s
    pr=np.hypot(P[...,0],P[...,1]);pp=phi(pr.ravel()).reshape(pr.shape)
    line=(pp[0]-pp[-1])*vc**2
    # absolute-space light geometry at t_rec: direction and distance of the emission event as seen from us now
    de=xe_rec-Xo;dist_rec=np.hypot(de[:,0],de[:,1])*D;ang_rec=np.degrees(np.arccos(np.clip(-(de[:,0])/np.hypot(de[:,0],de[:,1]),-1,1)))
    return dict(rho=rho,th_obs=th_obs,sep_now=sep*D,kin_rec=kin_rec,kin_now=kin_now,grav_rec=grav_rec,grav_now=grav_now,line=line,
                dist_rec=dist_rec,ang_from_CoE_rec=ang_rec,Xo=Xo,Vo=Vo,vrmin=float(VR.min()))
NS=256
def lbvec(l,b):return np.array(hp.ang2vec(l,b,lonlat=True))
anti=lbvec(0.,15.);nv=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T
psi=np.degrees(np.arccos(np.clip(nv@anti,-1,1)))           # angle of each pixel from anti-CoE
def paint(tr,y):
    t=tr['th_obs'];o=np.argsort(t);return np.interp(psi,t[o],y[o])
out={};maps={}
for tag,R in (('Rsee',Rsee),('RV',RV)):
    tr=trace(R);mono=np.all(np.diff(tr['th_obs'])>0)
    geo=tr['kin_rec']+tr['grav_rec'];path=np.zeros_like(geo);tot=geo+path
    mg=paint(tr,geo);mp=paint(tr,path);mt=mg+mp
    m0=hp.remove_monopole(mt);m1=hp.remove_dipole(mt)
    out[tag]=dict(R_Gpc=R,rho=tr['rho'],th_obs_monotonic=bool(mono),vrmin=tr['vrmin'],
      path_len_today_Gpc=[float(tr['sep_now'].min()),float(tr['sep_now'].max())],
      emission_event_distance_Gpc=[float(tr['dist_rec'].min()),float(tr['dist_rec'].max())],
      emission_event_max_angle_from_CoE_deg=float(tr['ang_from_CoE_rec'].max()),
      light_travel_c_t0_Gpc=cG*t0,
      kin_rec_range=[float(tr['kin_rec'].min()),float(tr['kin_rec'].max())],grav_rec_range=[float(tr['grav_rec'].min()),float(tr['grav_rec'].max())],
      kin_now_range=[float(tr['kin_now'].min()),float(tr['kin_now'].max())],grav_now_range=[float(tr['grav_now'].min()),float(tr['grav_now'].max())],
      line_integral_minus_endpoint_maxabs=float(np.max(np.abs(tr['line']-tr['grav_rec']))),
      Vo_code=tr['Vo'].tolist(),Xo_code=tr['Xo'].tolist(),
      monopole_dTT=float(np.mean(mt)),rms_dTT_monopole_removed=float(np.std(m0)),rms_dTT_dipole_removed=float(np.std(m1)),
      rms_uK_monopole_removed=float(np.std(m0)*2.7255e6),rms_uK_dipole_removed=float(np.std(m1)*2.7255e6),
      rms_uK_variant_today_field_dipole_removed=float(np.std(hp.remove_dipole(paint(tr,tr['kin_now']+tr['grav_now'])))*2.7255e6),
      profile_theta_obs=tr['th_obs'][::8].tolist(),profile_kin=tr['kin_rec'][::8].tolist(),profile_grav=tr['grav_rec'][::8].tolist())
    if tag=='Rsee':maps=dict(geo=mg,path=mp,tot=mt,tot_nodip=m1)
LMAX=3*NS-1
cl=hp.anafast(hp.remove_monopole(maps['tot']),lmax=LMAX,iter=3)
cl1=hp.anafast(maps['tot_nodip'],lmax=LMAX,iter=3)
alm=hp.map2alm(maps['tot_nodip'],lmax=10,iter=3)
np.savez('cmb_lightpath_maps.npz',**maps,cl=cl,cl_nodip=cl1)
J=dict(sources=dict(integrator='prof_core.py integrate (RK4 in log t, taui=0.01, nsteps=200), kind eta, Ob=0',fit='p1_migkas_D34.96.json',
  lock='catwise-eta-density-migkas.py L_from_delta_eta(D,d=RV,eta,0.031)',recomb_rule='cycle7_recombination/README.md, cycle4/dtt_rec_track/compute_track.py',
  conventions='cycle4/lib4.py',orientation='anti-CoE (0,+15) from blended fig8 rows; axisymmetric field'),
 params=dict(D_Gpc=D,eta=eta,mc=mc,l0_fit_deg=FIT['l0'],L_Gpc=L,LD=LD,Ht=Ht,t0_Gyr=t0,mu=mu,Rb=Rb,v_unit_over_c=vc,f_rec=f,R_see_Gpc=Rsee,RV_Gpc=RV,nside=NS,lmax=LMAX,T0_uK=2.7255e6),
 not_in_model=['void (no density profile exists; 15 Mpc is only the Alnes-Amarzguioui literature bound in notes)','rings as 3D densities (only sky-pattern Gaussians in A_required_rings_*)','time evolution of Phi (prof_core CoE profile fixed in space) -> ISW/Rees-Sciama term identically 0'],
 runs=out)
json.dump(J,open('cmb_lightpath.json','w'),indent=1)
print(json.dumps({k:v for k,v in J['params'].items()},indent=0))
for k,v in out.items():print(k,{a:b for a,b in v.items() if not a.startswith('profile')})
