"""CMB vs angle for CLASSICAL and SR frameworks with growth fraction f. usage: python3 cmb_flat.py classical|sr"""
import numpy as np,json,sys
from scipy.optimize import brentq
import core,dyn
fw=sys.argv[1]; REL=fw=='srrel'
if REL: BF=sys.argv[2]; F=json.load(open(f'migkas_srrel_{BF}.json'))['1deg']; fw='sr'; TAG=f'srrel_{BF}'
else: F=json.load(open(f'migkas_{fw}.json'))['1deg']; TAG=fw
p=F['params']
if REL: p['U']=F['U']
MT,LD,Rb=core.profile(p['DR'],p['alpha']); Dphys=p['DR']*core.RV; phiT=core.phi_table(MT)
QG=np.geomspace(1e-5,(F['qf'] if REL else 3.0),1500); NS=6000 if fw=='classical' else 3000
B=np.arange(0,181,1.0); nb=np.c_[np.cos(np.radians(B)),np.sin(np.radians(B))]; wgt=np.sin(np.radians(B)); wgt[[0,-1]]=np.sin(np.radians(0.5))/4  # area weights (sky symmetric about CoE axis)
FS=np.round(np.arange(0,1.0001,0.1),2)
def layers(mass):
    if fw=='classical':
        TI=1e-6; V0=dyn.cl_V0(p['mu'],MT,ti=TI,NS=NS); mu=p['mu']*mass   # launch moved earlier than paper's 0.01 so recombination-era light exists
        T,X,V=dyn.cl_integrate(QG,V0,mu,MT,NS=NS,ti=TI,keep=True)
        vo1=float(np.interp(1.0,*dyn.cl_integrate(QG,V0,p['mu'],MT,NS=NS,ti=TI)))  # mass case observer speed (D/t_obs)
        t0=vo1*core.TH0   # Gyr: observer v/D = H0 (paper convention), same time unit for both cases
        c=core.C_GPC_GYR*t0/Dphys; Phi=lambda r: mu*phiT(r)
        meta=dict(t0_Gyr=t0,c_code=c,V0=V0)
    else:
        U,K=p['U'],p['K']*mass; to=F['t_obs']
        T,X,P=dyn.sr_integrate(QG,U,K,MT,to,NS=NS,keep=True); V=P/np.sqrt(1+P*P); c=1.0
        Phi=lambda r: (K/2)*phiT(r); meta=dict(t_obs_c_units=to,t_obs_Gyr=to*Dphys/core.C_GPC_GYR)
    keep=np.argmax(V.min(0)>0); keep=max(keep,1)+1    # keep only layers that move outward at all times (Chris: everything moves outward); slower ones fall back
    meta['inner_outward_layer_q']=float(QG[keep]); X,V=X[:,keep:],V[:,keep:]; qg=QG[keep:]
    assert all(np.all(np.diff(X[n])>0) for n in range(0,len(T),50)), 'layers crossed'
    G=np.array([np.gradient(V[n],X[n])+2*V[n]/X[n] for n in range(len(T))])
    # observer = layer that is at the mass-case observer label
    if fw=='classical': qo=1.0
    else: qo=1.0
    Xo=float(np.interp(qo,qg,X[-1])); Vo=float(np.interp(qo,qg,V[-1]))
    return T,X,V,G,c,Phi,Xo,Vo,meta
def trace(f,T,X,V,G,c,Xo):
    n=len(B); x=np.tile([Xo,0.0],(n,1)); k=-nb; lnS=np.zeros(n); R=[np.hypot(*x.T)]; L=[lnS.copy()]; P=[x.copy()]; out=np.zeros(n,bool); O=[out.copy()]
    def vel(x,j):
        r=np.hypot(*x.T); v=np.interp(r,X[j],V[j]); g=np.interp(r,X[j],G[j]); return f*v[:,None]*x/r[:,None]+c*k, f*g/3, r>X[j][-1]
    for j in range(len(T)-1,0,-1):
        dt=T[j-1]-T[j]
        a1,s1,o1=vel(x,j); xm=x+dt*a1; a2,s2,o2=vel(xm,j-1)
        x=x+0.5*dt*(a1+a2); lnS=lnS-0.5*dt*(s1+s2); out|=o1|o2
        R.append(np.hypot(*x.T)); L.append(lnS.copy()); P.append(x.copy()); O.append(out.copy())
    global OUT; OUT=np.array(O); return np.array(R),np.array(L),np.array(P),out
def emission(Rr,Rp,Lp,Pp,T):
    Tb=T[::-1]; n=Rp.shape[1]; te=np.full(n,np.nan); lnS=np.full(n,np.nan); u=np.full((n,2),np.nan)
    for i in range(n):
        j=np.argmax(Rp[:,i]<Rr)
        if Rp[j,i]>=Rr: continue
        if OUT[j,i]: continue      # ray was outside the outermost layer (outside the universe's matter) between emission and us
        w=(Rp[j-1,i]-Rr)/(Rp[j-1,i]-Rp[j,i]); te[i]=Tb[j-1]+w*(Tb[j]-Tb[j-1]); lnS[i]=Lp[j-1,i]+w*(Lp[j,i]-Lp[j-1,i])
        xe=Pp[j-1,i]+w*(Pp[j,i]-Pp[j-1,i]); u[i]=xe/np.hypot(*xe)
    return te,lnS,u
def factors(Rr,f,Rp,Lp,Pp,T,X,V,c,Phi,Xo,Vo,use_grav=True,use_dop=True):
    te,lnS,u=emission(Rr,Rp,Lp,Pp,T)
    if np.all(np.isnan(te)): return None
    inner=np.array([rowpair(t,T,X,V)[0][0] if np.isfinite(t) else np.nan for t in te]); te=np.where(inner<Rr,te,np.nan)  # need an outward layer at R_rec then
    if np.all(np.isnan(te)): return None
    lnS=np.where(np.isfinite(te),lnS,np.nan)
    ve=np.array([np.interp(Rr,*rowpair(t,T,X,V)) if np.isfinite(t) else np.nan for t in te])
    ue=(1-f)*ve[:,None]*u/c; uo=np.array([(1-f)*Vo/c,0.0])
    if fw=='classical':
        vlos=np.sum((ue-uo)*nb,1); dop=1/(1+vlos)
        bmax=float(np.nanmax(np.hypot(*(ue-uo).T)))
    else:
        be=np.hypot(*ue.T); assert np.nanmax(be)<1 and abs(uo[0])<1,'speed>=c in SR'
        ge=1/np.sqrt(1-be**2); go=1/np.sqrt(1-uo[0]**2)
        dop=go*(1+uo@nb.T)/(ge*(1+np.sum(ue*nb,1))); vlos=None
        bmax=float(np.nanmax(be))
    gr=1+(float(Phi(Rr)[0])-float(Phi(Xo)[0]))/c**2
    S=np.exp(lnS); Tobs=3000/S*(dop if use_dop else 1)*(gr if use_grav else 1)
    return dict(te=te,S=S,dop=dop,grav=gr,T=Tobs,ue=ue,uo=uo,ve=ve,bmax=bmax)
def rowpair(t,T,X,V):
    j=min(max(np.searchsorted(T,t),1),len(T)-1); w=(t-T[j-1])/(T[j]-T[j-1]); return X[j-1]*(1-w)+X[j]*w, V[j-1]*(1-w)+V[j]*w
def wmean(y): m=np.isfinite(y); return float(np.sum(wgt[m]*y[m])/np.sum(wgt[m]))
def frac(y): return float(np.sum(wgt[np.isfinite(y)])/np.sum(wgt))
def ring(y): return np.r_[y,y[-2:0:-1]]          # 0..359 from 0..180
def stats(Tv,step):
    full=ring(Tv); bb=np.radians(np.arange(360.)); k=(np.arange(360)%step==0); d=(full-wmean(Tv))/wmean(Tv)
    d=d[k]; b=bb[k]; A=np.c_[np.ones_like(b),np.cos(b),np.sin(b)]; cf,*_=np.linalg.lstsq(A,d,rcond=None)
    return dict(full_amp_1e5=float((d.max()-d.min())/2/1e-5),dipole_1e5=float(np.hypot(cf[1],cf[2])/1e-5),
                hot_l=float((np.degrees(np.arctan2(cf[2],cf[1]))+p['l0'])%360),leftover_rms_1e5=float((d-A@cf).std()/1e-5))
sw=lambda y: float((np.nanmax(y)-np.nanmin(y))/2/np.nanmean(y)/1e-5)
res=dict(framework=fw,params=p,D_Gpc=Dphys); curves={}
for mass in (1,0):
    T,X,V,G,c,Phi,Xo,Vo,meta=layers(mass); case='with_mass' if mass else 'no_mass'
    res[case]=dict(meta=meta,observer_distance_D=Xo,observer_speed_over_c=Vo/c,cases={})
    for f in FS:
        Rp,Lp,Pp,out=trace(f,T,X,V,G,c,Xo); d=dict(f=float(f),rays_left_layer_grid=int(out.sum()))
        Rgrid=np.geomspace(1e-9,0.999*Xo,300)
        def diag(Rr):
            fa=factors(Rr,f,Rp,Lp,Pp,T,X,V,c,Phi,Xo,Vo)
            if fa is None: return 'no_emission_seen',None
            if fa['grav']<=0: return 'gravity_factor_negative',fa
            if np.nanmin(fa['dop'])<=0: return 'doppler_factor_negative',fa
            return 'ok',fa
        D_=[diag(r) for r in Rgrid]; cats=[x[0] for x in D_]; d['R_rec_scan']={k:int(cats.count(k)) for k in set(cats)}
        fr=np.array([frac(x[1]['T']) if x[1] is not None else 0 for x in D_]); d['max_sky_fraction_with_emission']=float(fr.max())
        if 'gravity_factor_negative' in cats:
            fa=D_[cats.index('gravity_factor_negative')][1]; d['example_gravity_factor']=float(fa['grav']); d['example_Phi_over_c2']=[float(Phi(Rgrid[cats.index('gravity_factor_negative')])[0]/c**2),float(Phi(Xo)[0]/c**2)]
        if 'doppler_factor_negative' in cats:
            fa=D_[cats.index('doppler_factor_negative')][1]; d['example_doppler_range']=[float(np.nanmin(fa['dop'])),float(np.nanmax(fa['dop']))]; d['example_max_through_space_speed_c']=fa['bmax']
        variants=[('all_effects',True,True),('without_gravity',False,True),('stretch_only',False,False)]
        d['variants']={}
        for vname,ug,ud in variants:
            vals=[]; frs=[]
            for r_,x in zip(Rgrid,D_):
                fa=factors(r_,f,Rp,Lp,Pp,T,X,V,c,Phi,Xo,Vo,ug,ud) if x[1] is not None else None
                good=fa is not None and np.nanmin(fa['T'])>0 and np.isfinite(np.nanmin(fa['T']))
                vals.append(wmean(fa['T'])-2.725 if good else np.nan); frs.append(frac(fa['T']) if good else 0)
            vals=np.array(vals); frs=np.array(frs); full=frs>0.999999
            ok=np.isfinite(vals)&(full if full.any() else True)
            roots=[i for i in range(len(Rgrid)-1) if ok[i] and ok[i+1] and np.sign(vals[i])!=np.sign(vals[i+1])]
            vd=dict(full_sky=bool(full.any()))
            if not roots:
                vd['not_computable']=('no valid R_rec (see R_rec_scan)' if not ok.any() else 'mean T never 2.725 K; valid mean T range %.4g..%.4g K'%(np.nanmin(vals[ok])+2.725,np.nanmax(vals[ok])+2.725))
                d['variants'][vname]=vd; continue
            i=roots[0]; Rr=np.exp(brentq(lambda z: wmean(factors(np.exp(z),f,Rp,Lp,Pp,T,X,V,c,Phi,Xo,Vo,ug,ud)['T'])-2.725,np.log(Rgrid[i]),np.log(Rgrid[i+1]),xtol=1e-12))
            fa=factors(Rr,f,Rp,Lp,Pp,T,X,V,c,Phi,Xo,Vo,ug,ud); assert abs(wmean(fa['T'])-2.725)<1e-6
            tunit=meta['t0_Gyr'] if fw=='classical' else Dphys/core.C_GPC_GYR
            rel=fa['ue']-fa['uo']; vd['sky_fraction_used']=frac(fa['T'])
            vd.update(n_roots=len(roots),R_rec_Mpc=float(Rr*Dphys*1e3),
              stretch_range=[float(np.nanmin(fa['S'])),float(np.nanmax(fa['S']))],swing_1e5_stretch=sw(1/fa['S']),
              doppler_range=[float(np.nanmin(fa['dop'])),float(np.nanmax(fa['dop']))],swing_1e5_doppler=sw(fa['dop']),
              gravity_factor=float(fa['grav']),Phi_over_c2_emitter=float(Phi(Rr)[0]/c**2),Phi_over_c2_observer=float(Phi(Xo)[0]/c**2),
              through_space_speed_c_toward_CoE=float(np.hypot(*rel[180])),through_space_speed_c_away_from_CoE=float(np.hypot(*rel[0])),
              los_speed_c_toward_CoE=float(np.sum(rel[180]*nb[180])),los_speed_c_away=float(np.sum(rel[0]*nb[0])),
              emitter_through_space_speed_c_max=float(np.nanmax(np.hypot(*fa['ue'].T))),observer_through_space_speed_c=float(fa['uo'][0]),
              travel_delay_Gyr=[float((T[-1]-np.nanmax(fa['te']))*tunit),float((T[-1]-np.nanmin(fa['te']))*tunit)],
              **({'deg6':stats(fa['T'],6),'deg1':stats(fa['T'],1)} if vd['sky_fraction_used']>0.999999 else {'partial_sky_swing_1e5':sw(fa['T'])}))
            d['variants'][vname]=vd; curves[f'{case}_f{f}_{vname}']=ring(fa['T']).tolist()
        res[case]['cases'][str(f)]=d; print(case,f,json.dumps(d,default=float)[:900],flush=True)
json.dump(res,open(f'cmb_{TAG}.json','w'),indent=1,default=float); json.dump(curves,open(f'cmb_{TAG}_curves.json','w'))
