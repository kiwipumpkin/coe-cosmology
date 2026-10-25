"""CMB = our past null cone ∩ recombined matter (per-element 3000 K, T=2.725 K (rho/rho_obs,today)^{1/3}), full cycloid LTB.
usage: cmb3.py <law> <fwd|back> <gtc|le> [dstep]"""
import os,sys,json,time,numpy as np,warnings;warnings.filterwarnings('ignore')
law,mode,vm=sys.argv[1:4]; DST=float(sys.argv[4]) if len(sys.argv)>4 else 6.
os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
else: os.environ.pop('VMODE',None)
sys.argv=['x',law,mode]
exec(open('fit_lc.py').read().split("ck=f'fwd")[0])
import grltb as G, ltb2 as T
from scipy.integrate import solve_ivp
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
tag=f'{law}_{mode}_{vm}'; x=np.array(json.loads(os.environ['XOVR'])); OT=os.environ['OUTTAG']
P,o=build(x)
if o is None: json.dump(dict(x=list(x),failed='build rejected (constraint/unphysical)'),open(f'gr2/scan/{OT}.json','w')); sys.exit()
io=o['io']; xo=np.log(F.R0[io])
if mode=='fwd': L=T.LTB(P,'fwd'); tn=o['t']*F.C
else: L=T.LTB(P,'back',r_ref_obs=F.R0[io],t_ref=0.0); tn=0.0
X=np.array([xo]); Ro=L.R(tn,X)[0]; Uo=L.Rdot(tn,X)[0]; tbo=L.tB(X)[0]; rho_o=L.rho_rel(tn,X)[0]
LNTH=3*np.log(3000/2.725)
TB_MIN=float(np.nanmin(L.tB(np.linspace(np.log(1e-3)+1,np.log(3e4)-1,3000))))+1e-9*abs(tn)
f=G.fields(o); tr=G.trace(f,np.arange(0,360,4.),nstep=320,lbmax_frac=0.999); af=G.migkas(tr); chi,l0=F.chi2_l0(af)
xi=np.linspace(np.log(1e-3)+1,np.log(3e4)-1,4000); inter=L.interior(xi); x_int=xi[inter].max() if inter.any() else -np.inf
print(tag,'R_o',Ro,'U_o',Uo,'interior (E<-1/2) labels up to r0=',np.exp(x_int),flush=True)
def lnrho(t,y):
    v=L.rho_rel(t,np.array([y[0]]))[0]; return np.log(v/rho_o) if np.isfinite(v) and v>0 else (1e3 if np.isfinite(L.R(t,np.array([y[0]]))[0]) else -1e3)   # R'<=0 (shell-crossing caustic, density -> inf) counts as hot; unborn -> not recombined
def trace(th):
    Lang=-Ro*np.sin(th); fr=T.rhs(L,Lang); y=np.array([xo,0.,0.,-np.cos(th)]); t=tn; t1=TB_MIN
    e1=lambda t,y: lnrho(t,y)-LNTH; e1.terminal=True; e1.direction=1
    e2=lambda t,y: y[0]-(L.xmin+1.0); e2.terminal=True; e2.direction=-1
    e3=lambda t,y: y[0]-x_int-1e-3 if np.isfinite(x_int) else 1.0; e3.terminal=True; e3.direction=-1
    nc=0
    for _ in range(4):
        try: s=solve_ivp(fr,(t,t1),y,method='LSODA',rtol=1e-10,atol=1e-12,events=(e1,e2,e3))
        except ValueError:   # event bracketing failed at an undefined-metric point: integrate without events and scan
            s=solve_ivp(fr,(t,t1),y,method='LSODA',rtol=1e-10,atol=1e-12,dense_output=True); s.status=0; s.t_events=[[],[],[]]; s.y_events=[[],[],[]]
            for ev_i,ev in enumerate((e1,e2,e3)):
                vals=np.array([ev(tt,s.y[:,k]) for k,tt in enumerate(s.t)]); sg=np.where((vals[:-1]<0)&(vals[1:]>=0) if ev_i==0 else (vals[:-1]>0)&(vals[1:]<=0))[0]
                if len(sg):
                    k=sg[0]; tt=np.linspace(s.t[k],s.t[k+1],200); vv=np.array([ev(q,s.sol(q)) for q in tt]); j=np.argmax(vv>=0) if ev_i==0 else np.argmax(vv<=0)
                    s.t_events[ev_i]=[tt[j]]; s.y_events[ev_i]=[s.sol(tt[j])]
            hits=[(s.t_events[i][0],i) for i in range(3) if len(s.t_events[i])]
            if hits:
                tl,il=max(hits)   # earliest-reached along the backward ray = largest t
                s.status=1; s.t_events=[[tl] if i==il else [] for i in range(3)]; s.y_events=[[s.y_events[i][0]] if i==il else [] for i in range(3)]
        if s.status==1 and len(s.t_events[1]) and not len(s.t_events[0]):
            t=s.t_events[1][0]; y=s.y_events[1][0].copy(); y[1]+=np.pi; y[3]=-y[3]; nc+=1; t+=(t1-t)*1e-12; continue
        break
    if s.status==1 and len(s.t_events[0]): return 'rec',s.t_events[0][0],s.y_events[0][0],nc
    if s.status==1 and len(s.t_events[2]): return 'trapped',s.t_events[2][0],s.y_events[2][0],nc
    fin=np.isfinite(s.y).all(0); te=s.t[fin][-1]; ye=s.y[:,fin][:,-1]
    xx=np.array([ye[0]]); why='unborn' if te<L.tB(xx)[0]+1e-6*(tn-tbo) else ('crunched/fell_into_M' if not np.isfinite(L.R(te,xx)[0]) else ('shell_crossing' if L.geo(te,xx)[3][0]<=0 else 'integrator_stop'))
    return why,te,ye,nc
rows=[]; t0=time.time(); lam=F.local(o)['Hbar']/H70
for ell in np.arange(0,360,DST):
    why,te,ye,nc=trace(np.radians(ell)); r=dict(ell=float(ell),status=why,ncross=nc,t_e_Gyr=float((te-tbo)/F.C*lam),flight_Gyr=float((tn-te)/F.C*lam),r0_e=float(np.exp(ye[0])))
    xx=np.array([ye[0]]); r['R_e_over_D']=float(L.R(te,xx)[0]/Ro); r['phi_e_deg']=float(np.degrees(ye[1])); r['U_e']=float(L.Rdot(te,xx)[0])
    if why!='rec':
        g_=L.geo(te,xx); r['diag']=dict(W2=float(1+2*L.E(xx)[0]),A=float(g_[3][0]),lnrho_rel=float(lnrho(te,ye)),tau_since_bang_frac=float((te-L.tB(xx)[0])/(tn-tbo)))
    if why=='rec': r['z']=float(np.exp(ye[2])-1); r['T_obs_K']=3000/(1+r['z'])
    rows.append(r); print(ell,why,{k:round(v,5) if isinstance(v,float) else v for k,v in r.items()},'%.0fs'%(time.time()-t0),flush=True)
out=dict(tag=tag,chi2=chi,sigma=sig(chi,len(nm)+1),l0=l0,Uo_c=float(Uo),D_Gpc=float(Ro*lam),age_obs_Gyr=float((tn-tbo)/F.C*lam),interior_r0_max=float(np.exp(x_int)) if np.isfinite(x_int) else 0.,
  status_counts={k:sum(r['status']==k for r in rows) for k in set(r['status'] for r in rows)},rows=rows)
ok=[r for r in rows if r['status']=='rec']; curves={}
def nodip(S):
    th=np.radians(np.arange(360.)); Xm=np.c_[np.ones(360),np.cos(th),np.sin(th)]; c=np.linalg.lstsq(Xm,np.log(S),rcond=None)[0]; return S/np.exp(c[1]*np.cos(th)+c[2]*np.sin(th))
if ok:
    out['T_obs_mean_K']=float(np.mean([r['T_obs_K'] for r in ok]))
if len(ok)==len(rows):
    el=np.array([r['ell'] for r in ok]); S=1/(1+np.array([r['z'] for r in ok])); S360=np.interp(np.arange(360.),np.append(el,360),np.append(S,S[0]))
    for t_,SS in (('raw',S360),('nodip',nodip(S360))):
        So=np.interp((np.arange(360.)-l0)%360,np.arange(361.),np.append(SS,SS[0])); rr=CC.resid(CC.G_FULL,So); v=[score(SS,p)[0] for p in np.arange(0,360,1.)]
        out[t_]=dict(cancel_l0=100*(1-CC.rms(rr)/S0),res_l0_1e5=CC.rms(rr)/1e-5,cancel_free=100*(1-min(v)/S0),modes={int(m['m']):float(m['A']) for m in CC.fourier(So/So.mean()-1,6)['modes']})
        curves[t_]=(So/So.mean()-1,rr)
out['x']=list(x); json.dump(out,open(f'gr2/scan/{OT}.json','w'),indent=1,default=float); sys.exit()
fig,ax=plt.subplots(3,1,figsize=(16,16)); ell=np.arange(360.)
ax[0].errorbar(F.IMIG,F.Y,yerr=0.5*(F.MLO+F.MHI),fmt='o',ms=6,c='k',label='Migkas 71 bins');ax[0].plot(ell,F.shifted(af,l0),lw=3,c='C3',label=f'GR model χ²={chi:.2f} σ={out["sigma"]:.1f} l0={l0:.0f}°');ax[0].set_ylabel('A/A_all',fontsize=18);ax[0].legend(fontsize=15)
gl=lambda r:(r['ell']+l0)%360
if curves:
    for (k,(d,rr)),c in zip(curves.items(),['C0','C2']): ax[1].plot(ell,d,lw=3,c=c,label=f'{k} dT/T'); ax[1].plot(ell,rr,lw=2,ls=':',c=c,label=f"{k} G·S resid ({out[k]['cancel_l0']:.1f}%)")
    ax[1].plot(ell,CC.G_FULL-1,c='k',lw=1.5,label='G−1')
elif ok:
    e_=np.array([gl(r) for r in ok]); i=np.argsort(e_); S=1/(1+np.array([r['z'] for r in ok]))[i]; ax[1].plot(e_[i],S/S.mean()-1,'o-',lw=3,label='dT/T raw (visible directions only)')
for r in rows:
    if r['status']!='rec': ax[1].axvspan(gl(r)-DST/2,gl(r)+DST/2,color={'trapped':'0.6','unborn':'C1','crunched/fell_into_M':'C4'}.get(r['status'],'C3'),alpha=0.3)
ax[1].set_ylabel('dT/T, residual',fontsize=18);ax[1].legend(fontsize=13,ncol=2)
ax[1].set_title('shaded: no recombined matter (grey trapped E<−½, orange unformed shell, purple crunched into M, red other)',fontsize=14)
e_=np.array([gl(r) for r in rows]); i=np.argsort(e_)
ax[2].plot(e_[i],np.array([r['t_e_Gyr'] for r in rows])[i],'o-',lw=3,label="emission/stop time after our shell's bang (Gyr)")
ax[2].plot(e_[i],np.array([r['flight_Gyr'] for r in rows])[i],'s-',lw=3,label='flight time (Gyr)')
ax[2].set_ylabel('Gyr',fontsize=18);ax[2].legend(fontsize=14);ax[2].set_xlabel('galactic longitude ℓ (deg)',fontsize=18)
for a in ax: a.set_xlim(0,360);a.tick_params(labelsize=14);a.grid()
fig.suptitle(f"3000 K per-element recombination {tag}: D={out['D_Gpc']:.3g} Gpc U_o/c={Uo:.3g}  {out['status_counts']}",fontsize=16)
fig.tight_layout();fig.savefig(f'gr2/page3000_{tag}.png',dpi=80);print('DONE',tag,out['status_counts'],flush=True)
