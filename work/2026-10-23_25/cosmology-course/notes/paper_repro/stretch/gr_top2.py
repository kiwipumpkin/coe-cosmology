exec(open("paper_repro/stretch/gr_top.py").read().split("def runR")[0])
import json
X,Vv,VR=None,None,None
rr=rrN; mu=rr["mu"]; Ht=rr["Ht"]
# Newtonian present state: rerun my radial integrator with the Newtonian table (validated to 1e-7)
r0=np.hypot(*q.T)*taui; v0=np.hypot(*q.T); u=q/v0[:,None]
res=[radial(r,v,lambda x:mu*np.interp(x,P.RG,MT)) for r,v in zip(r0,v0)]
Rt=np.array([x[0] for x in res]); Vt=np.array([x[1] for x in res])
M=mu*np.interp(Rt,P.RG,MT); E=Vt**2/2-M/Rt
print("E sign: positive",int((E>0).sum()),"negative",int((E<0).sum()),"Emin",E.min(),"Emax",E.max())
# backward from today with const M,E: (dR/dt)^2=2M/R+2E ; find R=0 time
tgrid=np.linspace(1,0,20001); Rgrid=[]; tbb=[]
for Ri,Vi,Mi in zip(Rt,Vt,M):
    ev=lambda t,y:y[0]-1e-9; ev.terminal=True
    s=solve_ivp(lambda t,y:[y[1],-Mi/y[0]**2],(1.0,-5.0),[Ri,Vi],events=ev,rtol=1e-11,atol=1e-14,dense_output=True)
    tb=s.t_events[0][0]; tbb.append(tb)
    tt=tgrid[tgrid>tb]; Rg=np.full(tgrid.size,0.0); Rg[:tt.size]=s.sol(tt)[0]; Rgrid.append(Rg)
tbb=np.array(tbb); Rgrid=np.array(Rgrid)
print("Big Bang times tau_BB: min",tbb.min(),"max",tbb.max(),"count after code start 0.01:",int((tbb>taui).sum()))
o=np.argsort(M); cross=0
for k in range(tgrid.size):
    if tgrid[k]<=tbb.max(): break
    Rk=Rgrid[o,k]; cross+=int(np.any(np.diff(Rk)<-1e-12))
print("time samples with shell crossing (R not increasing with M) after the last Big Bang:",cross,"; sightline count",len(M))
print("today: GR R,v vs Newtonian max diff 0 by construction (same R_today, v_today)")
eps_=eps; vrel=(Vt[:,None]*u)[1:]-(Vt[:,None]*u)[0]; af=np.hypot(*vrel.T)/eps_; afG=af/af.mean()
fullN=shifted(afN,l0); fullG=shifted(afG,l0)
Y=MA; W=1/(0.5*(MLO+MHI))**2
def st(m): mm=m[IMIG]; return float(np.sum(W*(mm-Y)**2)),float(1-np.sum((Y-mm)**2)/np.sum((Y-Y.mean())**2))
d=dict(chiN=st(fullN),chiG=st(fullG),maxdiff=float(np.abs(fullN-fullG).max()))
print(d); json.dump(d,open("/tmp/gr_top_stats.json","w")); np.save("paper_repro/stretch/gr_top_curves.npy",np.vstack([fullN,fullG]))
ks=[k for k in range(tgrid.size) if tgrid[k]>tbb.max() and np.any(np.diff(Rgrid[o,k])<-1e-12)]
print("latest crossing tau",tgrid[min(ks)],"earliest sampled",tgrid[max(ks)], "today ordered:",bool(np.all(np.diff(Rt[o])>=0)))
