"""Common-recombination-event CMB with the adaptive LTB integrator. usage: cmb2.py <law> <fwd|back> <gtc|le> [dstep]"""
import os,sys,json,time,numpy as np,warnings;warnings.filterwarnings('ignore')
law,mode,vm=sys.argv[1:4]; DST=float(sys.argv[4]) if len(sys.argv)>4 else 4.
os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
else: os.environ.pop('VMODE',None)
sys.argv=['x',law,mode]
exec(open('fit_lc.py').read().split("ck=f'fwd")[0])
import grltb as G, ltb2 as T
from scipy.optimize import brentq
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
tag=f'{law}_{mode}_{vm}'; x=np.array(json.load(open(f'gr/ck_{tag}.json'))['x'])
P,o=build(x); io=o['io']; xo=np.log(F.R0[io])
if mode=='fwd': L=T.LTB(P,'fwd'); tn=o['t']*F.C
else: L=T.LTB(P,'back',r_ref_obs=F.R0[io],t_ref=0.0); tn=0.0
Ro=L.R(tn,np.array([xo]))[0]; Uo=L.Rdot(tn,np.array([xo]))[0]
print('sanity R_o',Ro,'vs',o['r'][io],' U_o/c',Uo,'vs',o['U'][io]/F.C,flush=True)
f=G.fields(o); tr=G.trace(f,np.arange(0,360,4.),nstep=320,lbmax_frac=0.999); af=G.migkas(tr); chi,l0=F.chi2_l0(af)
tb=L.tB(np.array([xo]))[0]
t_rec=brentq(lambda t: L.R(t,np.array([xo]))[0]-Ro/1101,tb+1e-12*abs(tn-tb),tn,xtol=1e-14*abs(tn-tb))
print('tB_obs',tb,'t_rec',t_rec,'tn',tn,flush=True)
def fwd_from(t_e,xe,phe,q,Lang,lkt=0.0): return T.integ(L,Lang,t_e,[xe,phe,lkt,q],tn,rtol=1e-10,atol=1e-12)
def miss(s):
    xf,phf=s.y[0,-1],s.y[1,-1]
    if not np.isfinite(xf) or abs(s.t[-1]-tn)>1e-9*max(1,abs(tn-t_rec)): return np.array([np.nan,np.nan])
    _,_,Rp,A,_,_=L.geo(tn,np.array([xo])); dphi=(phf+np.pi)%(2*np.pi)-np.pi
    return np.array([A[0]*np.exp(xo)*(xf-xo),Ro*dphi])/Ro
rows=[]; ells=np.arange(0,360,DST); t0=time.time()
for ell in ells:
    th=np.radians(ell); r=dict(ell=float(ell))
    sb=T.ray(L,tn,xo,th,t_rec,rtol=1e-11,atol=1e-13)
    ok=abs(sb.t[-1]-t_rec)<1e-9*abs(tn-t_rec) and np.all(np.isfinite(sb.y[:,-1]))
    r['reached_rec']=bool(ok); r['ncross']=int(sb.ncross)
    if not ok:
        fin=np.isfinite(sb.y).all(0); xe=sb.y[0][fin][-1]; ts=sb.t[fin][-1]
        r['stop_t_frac_of_flight']=float((tn-ts)/(tn-t_rec)); r['stop_label_r0']=float(np.exp(xe)); r['stop_label_unborn_at_trec']=bool(t_rec<L.tB(np.array([xe]))[0])
        rows.append(r); print(ell,'not reached',r,flush=True); continue
    xe,phe,lkt,q=sb.y[:,-1]; Lang=-Ro*np.sin(th)
    r.update(z_back=float(np.exp(lkt)-1),x_e=float(xe),phi_e=float(phe),R_e=float(L.R(t_rec,np.array([xe]))[0]),unborn=bool(t_rec<L.tB(np.array([xe]))[0]))
    sf=fwd_from(t_rec,xe,phe,q,Lang,lkt); m0=miss(sf); r['land_err_before']=float(np.hypot(*m0))
    pe,qe=phe,q; sb_=sf
    for it in range(2):
        if not np.all(np.isfinite(m0)) or np.hypot(*m0)<1e-7: break
        h=1e-6; J=np.zeros((2,2))
        J[:,0]=(miss(fwd_from(t_rec,xe,pe+h,qe,Lang,lkt))-m0)/h
        J[:,1]=(miss(fwd_from(t_rec,xe,pe,qe-h*np.sign(qe),Lang,lkt))-m0)/(-h*np.sign(qe))
        try: d=np.linalg.solve(J,-m0)
        except Exception: break
        if not np.all(np.isfinite(d)): break
        pe+=d[0]; qe=np.clip(qe+d[1],-1,1); s2=fwd_from(t_rec,xe,pe,qe,Lang,lkt); m2=miss(s2)
        if np.all(np.isfinite(m2)) and np.hypot(*m2)<np.hypot(*m0): sb_,m0=s2,m2
        else: break
    r['land_err_after']=float(np.hypot(*m0)) if np.all(np.isfinite(m0)) else None
    r['z_fwd']=float(np.exp(lkt-sb_.y[2,-1])-1) if np.isfinite(sb_.y[2,-1]) else None
    rows.append(r); print('%5.0f z_b %.6g z_f %s land %.2e->%s unborn %s %.0fs'%(ell,r['z_back'],r['z_fwd'],r['land_err_before'],r['land_err_after'],r['unborn'],time.time()-t0),flush=True)
L_=F.local(o); lam=L_['Hbar']/H70
ok=[r for r in rows if r.get('z_fwd') is not None]
out=dict(tag=tag,x=list(x),chi2=chi,sigma=sig(chi,len(nm)+1),l0=l0,Uo_c=float(Uo),D_Gpc=float(Ro*lam),flight_time_Gyr=float((tn-t_rec)/F.C*lam),t_rec_after_obs_bang_Gyr=float((t_rec-tb)/F.C*lam),age_obs_Gyr=float((tn-tb)/F.C*lam),
  n_dir=len(rows),n_reached=sum(r['reached_rec'] for r in rows),n_unborn=sum(bool(r.get('unborn')) for r in rows),max_land_err_after=max([r['land_err_after'] for r in ok if r['land_err_after'] is not None],default=None),rows=rows)
out['all_emitters_formed']=len(ok)==len(rows) and all(not r['unborn'] for r in ok)
def nodip(S,th):
    X=np.c_[np.ones(len(th)),np.cos(th),np.sin(th)]; c=np.linalg.lstsq(X,np.log(S),rcond=None)[0]; return S/np.exp(c[1]*np.cos(th)+c[2]*np.sin(th))
curves={}
if len(ok)==len(rows):
    el=np.array([r['ell'] for r in ok]); S=1/(1+np.array([r['z_fwd'] for r in ok]))
    S360=np.interp(np.arange(360.),np.append(el,360),np.append(S,S[0]))
    for t_,SS in (('raw',S360),('nodip',nodip(S360,np.radians(np.arange(360.))))):
        So=np.interp((np.arange(360.)-l0)%360,np.arange(361.),np.append(SS,SS[0])); rr=CC.resid(CC.G_FULL,So)
        v=[score(SS,p)[0] for p in np.arange(0,360,1.)]
        out[t_]=dict(cancel_l0=100*(1-CC.rms(rr)/S0),res_l0_1e5=CC.rms(rr)/1e-5,cancel_free=100*(1-min(v)/S0),
            modes={int(m['m']):float(m['A']) for m in CC.fourier(So/So.mean()-1,6)['modes']},mean_z=float(np.mean(1/SS-1)))
        curves[t_]=(So/So.mean()-1,rr)
json.dump(out,open(f'gr2/{tag}.json','w'),indent=1,default=float)
fig,ax=plt.subplots(3,1,figsize=(16,15)); ell=np.arange(360.)
ax[0].errorbar(F.IMIG,F.Y,yerr=0.5*(F.MLO+F.MHI),fmt='o',ms=6,c='k',label='Migkas 71 bins');ax[0].plot(ell,F.shifted(af,l0),lw=3,c='C3',label=f'GR model χ²={chi:.2f} σ={out["sigma"]:.1f} l0={l0:.0f}°');ax[0].set_ylabel('A/A_all',fontsize=18);ax[0].legend(fontsize=15)
for (k,(d,rr)),c in zip(curves.items(),['C0','C2']):
    ax[1].plot(ell,d,lw=3,c=c,label=f'{k} dT/T'); ax[1].plot(ell,rr,lw=2,ls=':',c=c,label=f"{k} G·S resid ({out[k]['cancel_l0']:.1f}%)")
if curves: ax[1].plot(ell,CC.G_FULL-1,c='k',lw=1.5,label='G−1')
else: ax[1].text(0.03,0.5,f"CMB incomplete: {out['n_reached']}/{out['n_dir']} directions reach E_rec; {out['n_unborn']} emitters unformed",transform=ax[1].transAxes,fontsize=18)
ax[1].set_ylabel('dT/T, residual',fontsize=18);ax[1].legend(fontsize=13,ncol=2)
oe=[r for r in rows if 'R_e' in r]
if oe:
    e_=(np.array([r['ell'] for r in oe])+l0)%360; i=np.argsort(e_)
    ax[2].plot(e_[i],np.array([r['R_e'] for r in oe])[i]/Ro,'o-',lw=3,label='emitter radius R_e/D at t_rec')
    ax[2].plot(e_[i],np.degrees(np.array([r['phi_e'] for r in oe]))[i]/180,'s-',lw=2,label='emitter polar angle φ_e/180°')
    ax2=ax[2].twinx(); ax2.semilogy(e_[i],np.array([r['z_back'] for r in oe])[i],'k--',lw=2); ax2.set_ylabel('z (dashed)',fontsize=18); ax2.tick_params(labelsize=14)
ax[2].set_ylabel('R_e/D, φ_e/180°',fontsize=18);ax[2].legend(fontsize=13,loc='upper left')
ax[2].set_title(f"flight time t_n−t_rec = {out['flight_time_Gyr']:.4g} Gyr (same in all directions: common event)",fontsize=16)
for a in ax: a.set_xlim(0,360);a.tick_params(labelsize=14);a.grid()
ax[2].set_xlabel('galactic longitude ℓ (deg)',fontsize=18)
fig.suptitle(f"adaptive-LTB {tag}: D={out['D_Gpc']:.3g} Gpc U_o/c={Uo:.3g}  E_rec reached {out['n_reached']}/{out['n_dir']}",fontsize=18)
fig.tight_layout();fig.savefig(f'gr2/page_{tag}.png',dpi=80);print('DONE',tag,flush=True)
