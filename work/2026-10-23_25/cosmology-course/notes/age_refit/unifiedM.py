import sys,json,os,numpy as np,warnings;warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import joint_vgtc as J
def Smodel_N(b1,q,mu,s,phi,p):
    th=np.radians(J.cc.ELL-phi); xe0=1+q*np.cos(th); ye0=q*np.sin(th); re0=np.hypot(xe0,ye0)
    ro,re=s,s*re0; vo2=J.v2(1,ro,b1,mu,p); ve2=J.v2(re0,re,b1,mu,p)
    if vo2<=0 or np.any(ve2<=0) or np.any(re0<1e-6): return None
    vo,ve=np.sqrt(vo2/(1+vo2)),np.sqrt(ve2/(1+ve2)); xe,ye=s*xe0,s*ye0; dx,dy=ro-xe,-ye
    a=vo**2-1;b=2*dx*vo;c=dx**2+dy**2; tf=(-b-np.sqrt(b*b-4*a*c))/(2*a); nx,ny=(dx+vo*tf)/tf,dy/tf
    Dm=np.sqrt(1+vo2)/np.sqrt(1+ve2)*(1-vo*nx)/(1-ve*(xe*nx+ye*ny)/re)
    return Dm if SRONLY else Dm*np.exp(mu/ro-mu/re)   # Newtonian-potential redshift: no horizon singularity
NS={};os.chdir('/workspace/cosmology-course/notes');exec(open('catwise-density-levelset-migkas.py').read().replace('if __name__ == "__main__":\n    main()','pass'),NS)
os.chdir('/workspace/cosmology-course/notes/age_refit')
pc=NS['pc'];RV=NS['RV'];IMIG=NS['IMIG'];Y=NS['Y'];W=NS['W'];TH1=NS['TH1'];shifted=NS['shifted']
C=299792.458;GPCm=3.0857e25;G=6.674e-11;GYR=3.156e16
LAW=sys.argv[1]; DFIX=float(sys.argv[2]) if len(sys.argv)>2 and sys.argv[2]!='-' else None; SRONLY=len(sys.argv)>3 and sys.argv[3]=='sr'; LOCK=LAW.endswith('_lock'); KIND=LAW.replace('_lock','')
SH={'exp':[],'soft':['alpha'],'eta':['eta'],'pl':['n'],'gauss':[],'cored':[]}[KIND]
names=['D','M','l0']+([] if LOCK else ['L'])+SH+['H','q','s','p']
def integ(P):
    shape={k:P[k] for k in SH}
    if KIND=='soft': shape['beta']=1.0
    if LOCK:
        o=NS['sky_locked'](KIND,P['D'],P['mc'],NS['D_LOCK'],shape)
        if o is None: return None
        _,rr,LD,L=o
    else:
        a=NS['_pc_args'](KIND,P['D']/RV,P['L']/P['D'],shape); rr=pc.run(KIND,a,P['mc'],0.0); L=P['L']
        if not np.all(np.isfinite(rr['af'])) or rr['vrmin']<=0: return None
    return rr,L
def ev(P):
    P=dict(P); P['mc']=2*G*P['M']/((P['H']*1e3/3.0857e22)**2*(P['D']*GPCm)**3)
    if P['D']<=RV*1.05 or P['M']<=0 or P['H']<=0 or P['q']<=0 or P['s']<=0: return None
    if not LOCK and P['L']<=0: return None
    try: o=integ(P)
    except Exception: return None
    if o is None: return None
    rr,L=o
    af=rr['af']; A=shifted(af/af.mean(),P['l0'])[IMIG]; chi=float(np.sum(W*(A-Y)**2))
    t0=rr['Ht']/P['H']*977.8*GYR            # s ; H = Ht/t0
    Du=P['D']*GPCm/1e3*1e3                    # m  (length unit)
    rE=rr['rE']*Du; vo=rr['vE']*Du/t0/(C*1e3) # v_o/c
    M=P['M']*np.interp(rr['rE'],pc.RG,rr['MT'])  # enclosed mass at Earth (free M = M(<D))
    s=P['s']; d1=rE/s; mu=G*M/((C*1e3)**2*d1)
    b2=vo**2-2*mu*(1/s-1)
    if b2<=0: return None
    b1=np.sqrt(b2)
    S=Smodel_N(b1,P['q'],mu,s,P['l0'],P['p'])
    if S is None: return None
    rms=J.cc.rms(J.cc.resid(J.cc.G_FULL,S))
    age=J.ageF(b1,mu,s)*977.8/P['H']
    th=np.radians(TH1);w=np.abs(np.sin(th));Hm=float(np.sum(af*w)/np.sum(w)/rr['Ht']*P['H'])
    return dict(chi2=chi,rms=rms,cancel=100*(1-rms/J.cc.rms(J.cc.resid(J.cc.G_FULL,np.ones(360)))),age=float(age),
      age_integrator=rr['Ht']/P['H']*977.8,Hmean=Hm,vo=float(vo),beta1=float(b1),mu=float(mu),M=float(M),D_Earth=rE/GPCm,L=float(L),
      mc=float(P['mc']),Rs_over_r=float(2*mu/s),horizon_flag=bool(2*mu/s>1),A=A.tolist(),R2=NS['r2_u'](Y,A))
def cost(z):
    P=dict(zip(names,z)); P['l0']%=360
    if DFIX: P['D']=np.log(DFIX)
    for k in ('D','M','L','H','q','s','alpha','eta','n'):
        if k in P: P[k]=np.exp(P[k])
    e=ev(P)
    if e is None: return 1e6
    return (10*max(0,e['chi2']-11.6)/45.14+e['chi2']/45.14+10*max(0,e['rms']/2.37e-3-1)**2
            +100*max(0,(13-e['age'])/0.5)**2+100*max(0,67-e['Hmean'],e['Hmean']-73)**2)
from scipy.optimize import differential_evolution
BX={'D':(np.log(15.02),np.log(2000)),'M':(80,140),'l0':(0,360),'L':(-10,8),'alpha':(-3,3),'eta':(-3,3),'n':(-3,6),'H':(0,np.log(1e4)),'q':(-12,3),'s':(-10,10),'p':(-20,20)}
B=[BX[n] for n in names]
r=differential_evolution(cost,B,seed=1,popsize=12,maxiter=60,tol=1e-8,polish=False,init='sobol')
print('DE',r.fun,flush=True)
best=minimize(cost,r.x,method='Nelder-Mead',options=dict(maxiter=2500,xatol=1e-7,fatol=1e-9,adaptive=True))
if r.fun<best.fun: best=r
P=dict(zip(names,best.x));P['l0']%=360
if DFIX: P['D']=np.log(DFIX)
for k in P:
    if k not in('l0','p'): P[k]=float(np.exp(P[k]))
edges=[n for n,v,(lo,hi) in zip(names,best.x,B) if n!='l0' and not (DFIX and n=='D') and (v<lo+1e-3*(hi-lo) or v>hi-1e-3*(hi-lo))]
e=ev(P);e.update(params=P,cost=float(best.fun),law=LAW,edges=edges,DFIX=DFIX,sronly=SRONLY)
SRONLY=not SRONLY; e['cost_other_gravmode']=float(cost(best.x)); SRONLY=not SRONLY
tag=LAW+('_D%g'%DFIX if DFIX else '')+('_sr' if SRONLY else '')
json.dump(e,open(f'unifiedM_{tag}.json','w'),indent=1,default=float)
print({k:v for k,v in e.items() if k!='A'})
