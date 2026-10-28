# v2: per-shell recombination in the time-dependent CoE density, and search for emission events on our past light cone.
# Sources: prof_core.py (eta law, acc/integrate scheme, phi_coe), p1_migkas_D34.96.json, catwise L_from_delta_eta,
#   cycle4/dtt_rec_track/README.md l.9-11 ("T ∝ 1/a ∝ 1/t for each shell"; T_rec=3000 K, T0=2.725 K).
# Per-shell T∝1/a with mass conservation (rho∝a^-3) is used as T = T0 (rho/rho_obs,now)^(1/3): one T(rho) law for all shells,
#   normalised to T0 at our shell today (note: no explicit T-rho line exists in his notes; this is the T∝1/a rule rewritten).
# Time-dependent profile: radial shell trajectories r(q,t) from the same acceleration and RK4/log-t scheme as prof_core.integrate
#   (free coasting x=q t before taui, as in the integrator); today's eta-law density S(r) is assigned to shells at r(q,t0) and
#   carried back by mass conservation: rho(q,t)=S(r0) r0^2 dr0/dq / (r^2 dr/dq).  Writes cmb_lightpath_v2.json / _v2.npz.
import os,sys,json,numpy as np,importlib.util as iu
HERE=os.path.dirname(os.path.abspath(__file__));BUB='/workspace/cosmology/papers/spreadsheet_model/bubble';sys.path.insert(0,BUB);cwd=os.getcwd();os.chdir(BUB)
import prof_core as pc
from bubble_core import RKM,GPC,TH0
sp=iu.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=iu.module_from_spec(sp);sp.loader.exec_module(cw);os.chdir(cwd)
FIT=json.load(open('/workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json'))['34.96'];D,eta,mc=FIT['D'],FIT['eta'],FIT['mc'];RV=float(RKM/GPC)
L=cw.L_from_delta_eta(D,RV,eta,cw.DELTA_USE);rr=pc.run('eta',dict(DR=D/RV,LD=L/D,eta=eta),mc,0.0)
mu,MT,Rb,Ht,S=rr['mu'],rr['MT'],rr['Rb'],rr['Ht'],rr['S'];t0=Ht*TH0;cG=299792.458*3.15576e16/3.0857e22;c=1/(D/t0/cG)  # c in D/t0
T0,TR=2.725,3000.0;TAUI,NST=0.01,200
q=np.geomspace(1e-3,3.0,3000)
def a_r(r):r=np.maximum(r,1e-12);return -mu*np.interp(r,pc.RG,MT)/r**2
s0=np.log(TAUI);h=-s0/NST;x=q*TAUI;v=q.copy();tt=[TAUI];XR=[x.copy()];VV=[v.copy()];s=s0
for n in range(NST):
    t1,t2,t4=np.exp(s),np.exp(s+h/2),np.exp(s+h)
    k1=t1*v;l1=t1*a_r(x);k2=t2*(v+h/2*l1);l2=t2*a_r(x+h/2*k1);k3=t2*(v+h/2*l2);l3=t2*a_r(x+h/2*k2);k4=t4*(v+h*l3);l4=t4*a_r(x+h*k3)
    x=x+h/6*(k1+2*k2+2*k3+k4);v=v+h/6*(l1+2*l2+2*l3+l4);s+=h;tt.append(np.exp(s));XR.append(x.copy());VV.append(v.copy())
tt=np.array(tt);XR=np.array(XR);VV=np.array(VV)
# observer shell: the integrator's observer starts at q=1 on the axis
qo=1.0;ro_now=np.interp(qo,q,XR[-1]);vo_now=np.interp(qo,q,VV[-1])
r0=XR[-1];dr0=np.gradient(r0,q);Mq=4*np.pi*np.interp(r0,pc.RG,S)*r0**2*dr0        # Lagrangian mass density per dq
def rho(i):r=XR[i];return Mq/(r**2*np.gradient(r,q))
rho_obs_now=np.interp(ro_now,pc.RG,S);rho_rec=rho_obs_now*(TR/T0)**3
shell_cross=bool(np.any(np.gradient(XR,q,axis=1)<=0))
# rho(q,t): coasting before taui => rho=rho(q,taui)(taui/t)^3 ; after: tabulated
R=np.array([rho(i) for i in range(len(tt))])
trec=np.full(len(q),np.nan);xrec=np.full(len(q),np.nan);vrec=np.full(len(q),np.nan)
for j in range(len(q)):
    if R[0,j]<rho_rec: trec[j]=TAUI*(R[0,j]/rho_rec)**(1/3);xrec[j]=q[j]*trec[j];vrec[j]=q[j]   # before taui (coasting)
    else:
        k=np.where(R[:,j]<rho_rec)[0]
        if len(k):k=k[0];w=(np.log(R[k-1,j])-np.log(rho_rec))/(np.log(R[k-1,j])-np.log(R[k,j]));trec[j]=np.exp((1-w)*np.log(tt[k-1])+w*np.log(tt[k]));xrec[j]=(1-w)*XR[k-1,j]+w*XR[k,j];vrec[j]=(1-w)*VV[k-1,j]+w*VV[k,j]
ok=np.isfinite(trec)
# light cone: for direction psi (angle from anti-CoE), point X_o + n s with s=c(1-t); need |X_o+n s| = r_front(t)
# front: emission radius as function of emission time from the shell family (r_front may be multivalued; use all shells)
psi=np.radians(np.linspace(0,180,361));best=[]
for p in psi:
    tq=trec[ok];sq=c*(1-tq);dist=np.sqrt(ro_now**2+sq**2+2*ro_now*sq*np.cos(p))   # distance from CoE of our cone point at time tq
    g=dist-xrec[ok];best.append(np.min(np.abs(g)))
best=np.array(best)
beta_e=np.abs(vrec[ok])/c;beta_o=abs(vo_now)/c
J=dict(sources=dict(fit='p1_migkas_D34.96.json',integrator='prof_core.py (acc/integrate scheme, eta law, phi_coe)',lock='catwise-eta-density-migkas.py L_from_delta_eta',
   T_rule='cycle4/dtt_rec_track/README.md l.9-11: T ∝ 1/a ∝ 1/t for each shell; T_rec=3000 K, T0=2.725 K (also cycle7_recombination/README.md l.3)',
   not_found='no explicit T-rho relation, no rho(r,t) formula and no per-shell t_rec code exist in his notes; the spreadsheet exponential is the same shape family as prof_core exp/eta'),
 params=dict(D=D,eta=eta,mc=mc,L_Gpc=L,Ht=Ht,t0_Gyr=t0,mu=mu,c_code=c,taui=TAUI,nsteps=NST,r_obs_now_code=float(ro_now),v_obs_now_code=float(vo_now),
   rho_rec_over_rho_obs_now=float((TR/T0)**3),shell_crossing=shell_cross,n_shells=len(q),n_shells_recombined=int(ok.sum())),
 trec=dict(t_rec_over_t0_min=float(np.nanmin(trec)),t_rec_over_t0_max=float(np.nanmax(trec)),t_rec_obs_shell=float(np.interp(qo,q[ok],trec[ok])),
   t_rec_Myr_range=[float(np.nanmin(trec)*t0*1e3),float(np.nanmax(trec)*t0*1e3)],x_rec_Gpc_range=[float(np.nanmin(xrec)*D),float(np.nanmax(xrec)*D)]),
 light_cone=dict(c_t0_Gpc=c*D,our_distance_from_CoE_Gpc=float(ro_now*D),min_over_directions_of_min_mismatch_Gpc=float(best.min()*D),
   max_over_directions_Gpc=float(best.max()*D),n_directions_with_solution=int(np.sum(best*D<1e-3))),
 doppler=dict(beta_emitter_range=[float(beta_e.min()),float(beta_e.max())],beta_observer=float(beta_o)))
json.dump(J,open('cmb_lightpath_v2.json','w'),indent=1)
np.savez('cmb_lightpath_v2.npz',q=q,trec=trec,xrec=xrec,vrec=vrec,psi=psi,mismatch=best,ro=ro_now,c=c,D=D,t0=t0,tt=tt,robs_t=np.interp(qo,q,XR.T[:,0]) if False else XR[:,np.searchsorted(q,qo)])
print(json.dumps(J,indent=1))
