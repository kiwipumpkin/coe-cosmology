"""Unified CoE model: ONE Newtonian integrator (extended density law, shared launch exponent p) feeds BOTH
the Migkas A/A_all (|dv|/R at R=14.3 Gpc) and the CMB Doppler+gravity cancellation (emitters at comoving offset qc).
Velocities = integrator dr/dt treated as proper speed U (v>c allowed): beta=U/sqrt(c^2+U^2), gamma=sqrt(1+U^2/c^2).
Grav redshift: exp[(Phi_o-Phi_e)/c^2] from the SAME enclosed-mass potential. Age = t0 = Ht/H (=proper time since U=dr/dtau)."""
import sys,json,os,numpy as np,warnings;warnings.filterwarnings('ignore')
from scipy.optimize import minimize, differential_evolution
import joint_vgtc as J
NS={};os.chdir('/workspace/cosmology-course/notes');exec(open('catwise-density-levelset-migkas.py').read().replace('if __name__ == "__main__":\n    main()','pass'),NS)
os.chdir('/workspace/cosmology-course/notes/age_refit')
pc=NS['pc'];RV=NS['RV'];IMIG=NS['IMIG'];Y=NS['Y'];W=NS['W'];TH1=NS['TH1'];shifted=NS['shifted']
from bubble_core import geometry, rb_of
C=299792.458;GPCm=3.0857e25;G=6.674e-11;GYR=3.156e16;HU=1e3/3.0857e22
LAW=sys.argv[1]; LOCK=LAW.endswith('_lock'); KIND=LAW.replace('_lock','')
DFIX=float(sys.argv[2]) if len(sys.argv)>2 and sys.argv[2]!='-' else None
SRONLY=len(sys.argv)>3 and sys.argv[3]=='sr'
SH={'exp':[],'soft':['alpha'],'eta':['eta'],'pl':['n']}[KIND]
names=['D','M','l0']+([] if LOCK else ['L'])+SH+['H','qc','p']
ZMIG=np.array([0.01,0.025,0.04,0.055,0.075,0.10,0.14,0.20,0.32])  # Migkas+2020: z 0.003-0.45, median 0.075 (no per-cluster table on box)
RLOW=C*ZMIG/70/1e3   # Gpc, representative Migkas spread (no per-cluster z table used)
RING=RLOW if os.environ.get('RING','low')=='low' else np.array([RV])
G0=J.cc.G_FULL; S0=J.cc.rms(J.cc.resid(G0,np.ones(360)))
def model(P,sronly=SRONLY):
    D,M,H=P['D'],P['M'],P['H']
    if D<=0 or M<=0 or H<=0: return None
    if P.get('qc') is None: P=dict(P); P['qc']=C*1e3*(0.95/(H*HU))/(D*GPCm)  # last scattering: comoving offset = c t0 / D (light horizon)
    shape={k:P[k] for k in SH}
    if KIND=='soft': shape['beta']=1.0
    if LOCK:
        try: L=NS['L_from_delta'](KIND,D,NS['D_LOCK'],shape,NS['DELTA_USE'])
        except Exception: return None
    else: L=P['L']
    if not (np.isfinite(L) and L>0): return None
    a=NS['_pc_args'](KIND,D/RV,L/D,shape)
    DR=D/RV; qC,_=geometry(J.cc.ELL,1/P['qc'])
    rings=[geometry(TH1,D/Rk) for Rk in RING]
    q=np.vstack([qC[:1]]+[g[0][1:] for g in rings]+[qC[1:]]); nr=np.hypot(q[:,0],q[:,1]); q=q*(nr**(P['p']-1))[:,None]   # shared launch law v0 ~ r0^p
    Rb=rb_of(DR,1.5); S,MT=pc.tables(KIND,a,Rb); Ht=0.95
    mc=2*G*M/((H*HU)**2*(D*GPCm)**3)
    for _ in range(3):
        mu=mc*Ht**2/2
        X,V,VR=pc.integrate(q,mu,pc.RG,MT,0.0,Rb,0.01,200)
        if not np.all(np.isfinite(X)): return None
        Ht=np.hypot(*V[0])/np.hypot(*X[0])
    if VR.min()<=0: return None
    afs=[]
    for k,(gq,rr) in enumerate(rings):
        vr=V[1+360*k:1+360*(k+1)]-V[0]; a_=np.hypot(vr[:,0],vr[:,1])/rr; afs.append(a_/a_.mean())
    af=np.mean(afs,0); af_raw=np.mean([ (lambda vr: np.hypot(vr[:,0],vr[:,1])/rr)(V[1+360*k:1+360*(k+1)]-V[0]) for k,(gq,rr) in enumerate(rings)],0)
    n1=1+360*len(rings)
    A=shifted(af/af.mean(),P['l0'])[IMIG]; chi=float(np.sum(W*(A-Y)**2))
    t0=Ht/(H*HU)                         # s
    cu=C*1e3*t0/(D*GPCm)                 # c in integrator units
    xo,vo=X[0],V[0]; xe,ve=X[n1:],V[n1:]
    Uo=vo/cu; Ue=ve/cu                   # proper-speed vectors /c
    bo=Uo/np.sqrt(1+np.sum(Uo**2)); be=Ue/np.sqrt(1+np.sum(Ue**2,1))[:,None]
    d=xo-xe; dn=np.hypot(d[:,0],d[:,1])/cu   # light-units
    # observer moves during flight: solve |d + bo t| = t (t in light units)
    dd=d/cu; a2=np.sum(bo**2)-1; b2=2*(dd@bo); c2=np.sum(dd**2,1)
    tf=(-b2-np.sqrt(b2*b2-4*a2*c2))/(2*a2); n=(dd+np.outer(tf,bo))/tf[:,None]
    go=np.sqrt(1+np.sum(Uo**2)); ge=np.sqrt(1+np.sum(Ue**2,1))
    Dm=go*(1-n@bo)/(ge*(1-np.sum(be*n,1)))
    def phi(r): return float(pc.phi_coe(r,mu,MT)) if np.isscalar(r) else np.array([pc.phi_coe(x,mu,MT) for x in r])
    ro=np.hypot(*xo); re=np.hypot(xe[:,0],xe[:,1])
    Dg=1.0 if sronly else np.exp((phi(ro)-phi(re))/cu**2)
    Sx=Dm*Dg; Rres=J.cc.resid(G0,Sx); rms=J.cc.rms(Rres); fo=J.cc.fourier(Rres,6)
    menc=np.interp(ro,pc.RG,MT); Rs=2*mu*menc/ro/cu**2
    w=np.abs(np.sin(np.radians(TH1))); Hm=float(np.sum(af_raw*w)/np.sum(w)/Ht*H)
    return dict(chi2=chi,R2=NS['r2_u'](Y,A),rms=rms,cancel=100*(1-rms/S0),age=t0/GYR,Hmean=Hm,vo=float(np.hypot(*Uo)),
                D_Earth=float(ro*D),M=M,mc=float(mc),L=float(L),Rs_over_r=float(Rs),horizon_flag=bool(Rs>1),A=A.tolist(),modes={m['m']:m['A'] for m in fo['modes']},qc_used=float(P['qc']))
def unpack(z):
    P=dict(zip(names,z)); P['l0']%=360
    if DFIX: P['D']=np.log(DFIX)
    for k in P:
        if k not in('l0','p'): P[k]=float(np.exp(P[k]))
    return P
def cost(z,sronly=SRONLY):
    try: e=model(unpack(z),sronly)
    except Exception: e=None
    if e is None: return 1e6
    return (10*max(0,e['chi2']-11.6)/45.14+e['chi2']/45.14+10*max(0,e['rms']/2.37e-3-1)**2
            +100*max(0,(13-e['age'])/0.5)**2+100*max(0,67-e['Hmean'],e['Hmean']-73)**2)
if __name__=='__main__':
    BX={'D':(np.log(15.02),np.log(2000)),'M':(80,140),'l0':(0,360),'L':(-6,8),'alpha':(-3,3),'eta':(-3,3),'n':(-3,6),'H':(0,np.log(1e4)),'qc':(-8,3),'p':(-5,5)}
    B=[BX[n] for n in names]
    r=differential_evolution(cost,B,seed=1,popsize=8,maxiter=30,tol=1e-8,polish=False,init='sobol')
    print('DE',r.fun,flush=True)
    best=minimize(cost,r.x,method='Nelder-Mead',options=dict(maxiter=1000,xatol=1e-7,fatol=1e-9,adaptive=True))
    if r.fun<best.fun: best=r
    P=unpack(best.x); e=model(P)
    edges=[nm for nm,v,(lo,hi) in zip(names,best.x,B) if nm!='l0' and not(DFIX and nm=='D') and (v<lo+1e-3*(hi-lo) or v>hi-1e-3*(hi-lo))]
    e.update(params=P,cost=float(best.fun),law=LAW,edges=edges,DFIX=DFIX,sronly=SRONLY,cost_other_gravmode=float(cost(best.x,not SRONLY)))
    tag=LAW+('_D%g'%DFIX if DFIX else '')+('_sr' if SRONLY else '')
    json.dump(e,open(f'u2_{tag}.json','w'),indent=1,default=float)
    print({k:v for k,v in e.items() if k!='A'})
