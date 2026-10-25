"""GR (LTB) CMB: null geodesics from the comoving observer back to the surface R(r,t)=R_rec. 1+z = (dt/dlam)_e/(dt/dlam)_o for comoving emitter/observer."""
import numpy as np,json,sys
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from multiprocessing import Pool
import core,dyn
F=json.load(open('migkas_gr.json'))['1deg']; p=F['params']; tau=p['tau']
MT,LD,Rb=core.profile(p['DR'],p['alpha']); Dphys=p['DR']*core.RV; tunit=Dphys/core.C_GPC_GYR
RT=np.geomspace(1e-5,1e4,6000)
B=np.arange(0,181,1.0); wgt=np.sin(np.radians(B)); wgt[[0,-1]]=np.sin(np.radians(0.5))/4
def setup(mass):
    M=(mass*p['K']/2)*np.interp(RT,core.RG,MT)
    E,ok=dyn.ltb_E(RT,M,tau) if mass else (RT**2/(2*tau**2),np.ones_like(RT,bool))
    if not ok.all():   # tiny-r shells: bisection precision; use E ~ r^2 law from the first good shell (near-empty centre, Milne-like)
        bad=~ok; j=np.argmax(ok); assert np.all(RT[bad]<RT[j]*1.0001) and RT[j]<1e-3, ('LTB E solve failed at',RT[bad][-3:])
        E[bad]=E[j]*(RT[bad]/RT[j])**2
    assert np.all(1+2*E>0)
    return CubicSpline(RT,E),CubicSpline(RT,M)
from scipy.interpolate import RectBivariateSpline
_SP={}
def spline(mass):
    if mass in _SP: return _SP[mass]
    Es,Ms=setup(mass); Ef=lambda r: Es(r); Mf=lambda r: np.maximum(Ms(r),0)
    xg=np.linspace(np.log(1e-4),np.log(1e4),900); sg=np.linspace(np.log(tau)-25,np.log(tau)+1e-9,700)
    Lg=np.array([np.log(dyn.ltb_R(np.exp(xg),np.exp(s_),Ef,Mf)[0]) for s_ in sg]).T
    _SP[mass]=(RectBivariateSpline(xg,sg,Lg,kx=5,ky=5,s=0),Es); return _SP[mass]
def geo(args):
    b,mass=args; SP,Es=spline(mass)
    def AR(r,t):
        x=np.log(r); s_=np.log(t); L=SP(x,s_)[0,0]; Lx=SP(x,s_,dx=1)[0,0]; Lxx=SP(x,s_,dx=2)[0,0]; Ls=SP(x,s_,dy=1)[0,0]; Lxs=SP(x,s_,dx=1,dy=1)[0,0]
        R=np.exp(L); Lr=Lx/r; Lrr=(Lxx-Lx)/r**2; Rp=R*Lr; Rpp=R*(Lrr+Lr**2); Rd=R*Ls/t; Rdp=(Rp*Ls+R*Lxs/r)/t
        E=float(Es(r)); Ep=float(Es(r,1)); q=np.sqrt(1+2*E)
        return R,Rd,Rp/q,Rdp/q,Rpp/q-Rp*Ep/q**3
    Ef=lambda r: Es(r)
    class _M: pass
    R0,_,A0,_,_=AR(1.0,tau); J=R0*np.sin(np.radians(b))
    # independent variable s=ln t (going back from tau); w=ln|dt/dlam|, u=dr/dt
    def rhs(s,y):
        t=np.exp(s); r,ph,w,u,Ir,It=y; R,Rd,A,Ad,Ar=AR(max(abs(r),1.01e-4),t); Rp=A*np.sqrt(1+2*float(Ef(abs(r))))
        php=-J*np.exp(-w)/R**2
        dIr=-A*Ad*u*u; dIt=-R*Rd*php*php
        du=(R*Rp*php**2-A*Ar*u*u-2*A*Ad*u)/A**2+u*(A*Ad*u*u+R*Rd*php**2)
        return [t*u,t*php,t*(dIr+dIt),t*du,t*dIr,t*dIt]
    def stop(s,y): return AR(max(abs(y[0]),1.01e-4),np.exp(s))[0]-2e-5
    stop.terminal=True
    def far(s,y): return abs(y[0])-5e3
    far.terminal=True
    y0=[1.0,0.0,0.0,-np.cos(np.radians(b))/A0,0,0]
    sol=solve_ivp(rhs,(np.log(tau),np.log(tau)-25),y0,events=[stop,far],rtol=1e-10,atol=1e-12,method='DOP853')
    Y=sol.y; T=np.exp(sol.t); Rpath=np.array([AR(max(abs(r),1.01e-4),t)[0] for t,r in zip(T,Y[0])])
    def null(t,y):
        R,Rd,A,Ad,Ar=AR(max(abs(y[0]),1.01e-4),t); php=-J*np.exp(-y[2])/R**2; return -1+(A*y[3])**2+(R*php)**2
    nl=np.array([null(t,y) for t,y in zip(T,Y.T)])
    return dict(b=b,t=T,r=Y[0],ph=Y[1],td=-np.exp(Y[2]),Ir=Y[4],It=Y[5],R=Rpath,null=nl,Rmin=float(Rpath.min()),status=int(sol.status),hit_far=bool(len(sol.t_events[1])>0))
def at(path,Rr):
    R=path['R']; j=np.argmax(R<Rr)
    if R[j]>=Rr: return None
    w=(R[j-1]-Rr)/(R[j-1]-R[j]); g=lambda k: path[k][j-1]+w*(path[k][j]-path[k][j-1])
    return dict(te=g('t'),z1=abs(g('td')),Ir=g('Ir'),It=g('It'),re=g('r'),ph=g('ph'))
def wmean(y): return float(np.sum(wgt*y)/np.sum(wgt))
def ring(y): return np.r_[y,y[-2:0:-1]]
def stats(Tv,step):
    full=ring(Tv); bb=np.radians(np.arange(360.)); k=(np.arange(360)%step==0); d=(full-wmean(Tv))/wmean(Tv); d=d[k]; b=bb[k]
    A=np.c_[np.ones_like(b),np.cos(b),np.sin(b)]; cf,*_=np.linalg.lstsq(A,d,rcond=None)
    return dict(full_amp_1e5=float((d.max()-d.min())/2/1e-5),dipole_1e5=float(np.hypot(cf[1],cf[2])/1e-5),
                hot_l=float((np.degrees(np.arctan2(cf[2],cf[1]))+p['l0'])%360),leftover_rms_1e5=float((d-A@cf).std()/1e-5))
if __name__=='__main__':
    res=dict(framework='GR (LTB)',params=p,D_Gpc=Dphys,time_unit_Gyr=tunit,t_obs_Gyr=tau*tunit); curves={}
    for mass in (1,0):
        case='with_mass' if mass else 'no_mass'
        with Pool(8) as P: paths=P.map(geo,[(b,mass) for b in B])
        Rmin=np.array([pa['Rmin'] for pa in paths])
        d=dict(smallest_R_reached_per_direction_D_units={str(int(bb)):float(x) for bb,x in zip(B[::30],Rmin[::30])})
        def sel(Rr): return np.array([at(pa,Rr) is not None for pa in paths])
        def meanT(Rr,mask):
            e=[at(pa,Rr) for pa,m in zip(paths,mask) if m]
            if len(e)==0 or any(x is None for x in e): return np.nan
            z=np.array([x['z1'] for x in e]); return float(np.sum(wgt[mask]*3000/z)/np.sum(wgt[mask]))-2.725
        Rg=np.geomspace(3e-5,0.999,200); full=[r for r in Rg if sel(r).all()]
        d['full_sky_possible']=bool(len(full)>0)
        if full: mask=np.ones(len(B),bool); Rg2=np.array(full)
        else:
            # partial sky: directions that see an emission event at R_rec (labelled); mask fixed at the scan value giving most coverage
            cov=np.array([np.sum(wgt[sel(r)])/np.sum(wgt) for r in Rg]); d['sky_fraction_with_emission_vs_Rrec']={f'{r*Dphys*1e3:.4g} Mpc':float(c) for r,c in zip(Rg[::20],cov[::20])}
            res[case]=d; print(case,'no full-sky emission; max sky fraction %.3f'%cov.max(),flush=True)
            mask=None
        if mask is None:
            # solve R_rec with mean over directions that see emission at each trial (partial sky, labelled)
            def mT(Rr):
                m=sel(Rr); return meanT(Rr,m) if m.any() else np.nan
            v=np.array([mT(r) for r in Rg]); Rg2=Rg
        else: v=np.array([meanT(r,mask) for r in Rg2])
        ok=np.isfinite(v); roots=[i for i in range(len(Rg2)-1) if ok[i] and ok[i+1] and np.sign(v[i])!=np.sign(v[i+1])]
        if not roots:
            d['not_computable']='mean T never 2.725 K; valid range %s'%([float(np.nanmin(v)+2.725),float(np.nanmax(v)+2.725)] if ok.any() else 'none'); res[case]=d; print(case,d); continue
        i=roots[0]
        g_=(lambda z: meanT(np.exp(z),mask)) if mask is not None else (lambda z: meanT(np.exp(z),sel(np.exp(z))))
        Rr=np.exp(brentq(g_,np.log(Rg2[i]),np.log(Rg2[i+1]),xtol=1e-12))
        m=sel(Rr) if mask is None else mask; d['sky_fraction_used']=float(np.sum(wgt[m])/np.sum(wgt)); d['directions_with_emission_deg']=[float(B[m].min()),float(B[m].max())]
        nul=max(float(np.abs(pa['null'][pa['R']>=Rr]).max()) for pa,mm in zip(paths,m) if mm); d['max_null_error_to_emission']=nul; assert nul<1e-3,nul
        e=[at(pa,Rr) if mm else None for pa,mm in zip(paths,m)]; z1=np.array([x['z1'] if x else np.nan for x in e]); Tv=3000/z1; assert abs(np.nansum(wgt*Tv)/np.sum(wgt[m])-2.725)<1e-6
        Ir=np.array([x['Ir'] if x else np.nan for x in e]); It=np.array([x['It'] if x else np.nan for x in e]); assert np.allclose((Ir+It)[m],np.log(z1[m]),rtol=1e-4,atol=1e-4)
        te=np.array([x['te'] if x else np.nan for x in e])
        d.update(n_roots=len(roots),R_rec_Mpc=float(Rr*Dphys*1e3),one_plus_z_range=[float(np.nanmin(z1)),float(np.nanmax(z1))],
          split_radial_part_ln=[float(np.nanmin(Ir)),float(np.nanmax(Ir))],split_sideways_part_ln=[float(np.nanmin(It)),float(np.nanmax(It))],
          swing_1e5_radial_part=float((np.nanmax(np.exp(-Ir))-np.nanmin(np.exp(-Ir)))/2/np.nanmean(np.exp(-Ir))/1e-5),
          swing_1e5_sideways_part=float((np.nanmax(np.exp(-It))-np.nanmin(np.exp(-It)))/2/np.nanmean(np.exp(-It))/1e-5),
          travel_delay_Gyr=[float((tau-np.nanmax(te))*tunit),float((tau-np.nanmin(te))*tunit)],**({'deg6':stats(Tv,6),'deg1':stats(Tv,1)} if m.all() else {'partial_sky_swing_1e5':float((np.nanmax(Tv)-np.nanmin(Tv))/2/2.725/1e-5)}))
        res[case]=d; curves[case]=ring(Tv).tolist(); print(case,d,flush=True)
    json.dump(res,open('cmb_gr.json','w'),indent=1,default=float); json.dump(curves,open('cmb_gr_curves.json','w'))
