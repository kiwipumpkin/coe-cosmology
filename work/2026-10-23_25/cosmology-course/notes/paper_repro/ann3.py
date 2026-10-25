import numpy as np,json
exec(open('cmbann.py').read().split("meanT=lambda")[0])
d=json.load(open('cmbann.json')); a=d['annulus']; sp=json.load(open('split.json')); qi,qo=sp['q_in'],sp['q_out']; RR=sp['R_rec_Mpc']/1e3/DG; d['R_rec_Mpc']=sp['R_rec_Mpc']
tr=sp['trec_Myr']/1e3/t0; b,S,dop,gm,te=np.load('cmbann_keep.npy')
times=[tr,0.25,0.5,0.75,1.0]
V=lambda t,lo=None,hi=None: r_of(hi if hi else qo,t)**3-r_of(lo if lo else qi,t)**3; V0=V(tr)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,axs=plt.subplots(8,1,figsize=(18,120),gridspec_kw=dict(height_ratios=[1]*5+[0.55,0.55,0.55]))
L=31; ph=np.linspace(0,2*np.pi,20000); lam0=2.0/a['stretch']   # drawn wavelength now = 2 Gpc
rows=[]
for ax,t in zip(axs,times):
    ri,ro,D=r_of(qi,t)*DG,r_of(qo,t)*DG,Dt(t)*DG; s=(V(t)/V0)**(1/3); R=C*(t-tr)*DG; Tk=3000/s
    if t==1.0: assert abs(ri-sp['subshells_now_Gpc'][0])<1e-6 and abs(ro-sp['subshells_now_Gpc'][2])<1e-6 and abs(D-sp['subshells_now_Gpc'][1])<1e-6
    if t==tr: assert abs(s-1)<1e-12
    si=(V(t,qi,1.0)/V(tr,qi,1.0))**(1/3); so=(V(t,1.0,qo)/V(tr,1.0,qo))**(1/3)
    assert ri<D<ro
    ax.fill(np.r_[D*np.cos(ph),ri*np.cos(ph[::-1])],np.r_[D*np.sin(ph),ri*np.sin(ph[::-1])],color='royalblue',alpha=.35)
    ax.fill(np.r_[ro*np.cos(ph),D*np.cos(ph[::-1])],np.r_[ro*np.sin(ph),D*np.sin(ph[::-1])],color='orange',alpha=.5)
    rm=(ri+ro)/2; lam=lam0*s; ncyc=2*np.pi*rm/lam
    if lam>0.05 and ro-ri>0.05:
        w=rm+0.4*(ro-ri)*np.sin(ncyc*ph); ax.plot(w*np.cos(ph),w*np.sin(ph),color='darkred',lw=1)
        wl=f'drawn wavelength {lam:.2f} Gpc'
    else: wl=f'wavelength {lam*1e3:.3f} Mpc, annulus {(ro-ri)*1e3:.2f} Mpc thick: too small to draw'
    ax.plot(0,0,'*',color='r',ms=30); ax.annotate(f'CoE + R_rec = {d["R_rec_Mpc"]:.1f} Mpc (too small to see)',xy=(0,0),xytext=(-29,-8),fontsize=20,arrowprops=dict(arrowstyle='->',lw=2))
    ax.plot(D,0,'o',color='k',ms=16,zorder=5)
    if R>0.05: ax.plot(D+R*np.cos(ph),R*np.sin(ph),color='green',lw=3)
    if ro<1.5: ax.annotate('annulus + us here',xy=(D,0),xytext=(8,8),fontsize=20,arrowprops=dict(arrowstyle='->',lw=2))
    ax.text(10,-22,'blue: inner sub-shell D−R..D\norange: outer D..D+R',fontsize=19)
    ax.text(-29,29,f"t = {t*t0*1e3:.2f} Myr" if t==tr else f"t = {t*t0:.3f} Gyr",fontsize=26,va='top',fontweight='bold')
    ax.text(-29,25,f"D = {D:.4g} Gpc\nannulus {ri:.4g} – {ro:.4g} Gpc\nV/V_emit = {s**3:.3g}\nstretch = {s:.4g}  (inner {si:.5g}, outer {so:.5g})\nT (stretch only) = {Tk:.4g} K\nR = ct* = {R*1e3:.1f} Mpc (green)\n{wl}",fontsize=21,va='top',bbox=dict(fc='w',alpha=.85))
    ax.set_xlim(-L,L); ax.set_ylim(-L,L); ax.set_aspect('equal'); ax.tick_params(labelsize=20); ax.set_xlabel('Gpc',fontsize=22); ax.set_ylabel('Gpc',fontsize=22)
    rows.append(dict(t_Gyr=t*t0,D=D,ri=ri,ro=ro,Vratio=s**3,stretch=s,T=Tk,R_Gpc=R))
B=axs[-3]; bb,Tsp,Tsg=np.load('split_T.npy'); B.plot(bb,Tsg,lw=4,label=f'traced-ray emitter layers, mean {Tsg.mean():.3f} K'); Tw=3000/sp['stretch_whole']; B.axhline(Tw,lw=4,ls='--',color='tab:green',label=f'single shell = our whole annulus, stretch {sp["stretch_whole"]:.2f}: {Tw:.5f} K'); B.plot(bb,Tsp,lw=4,label=f'split inner/outer (cosβ interpolation), mean {Tsp.mean():.3f} K')
B.axhline(2.725,color='k',ls=':'); B.set_xlim(0,360); B.tick_params(labelsize=20); B.legend(fontsize=20)
B.set_xlabel('arrival direction β (deg, 0 = away from CoE)',fontsize=24); B.set_ylabel('T_obs now (K)',fontsize=24); B.set_title('T(direction) now',fontsize=26)
Zp=axs[-2]; Zp.plot(bb,(Tsp-2.725)*1e3,lw=4,color='tab:orange'); Zp.axhline(0,color='k',ls=':'); Zp.axhline((3000/sp['stretch_whole']-2.725)*1e3,lw=4,ls='--',color='tab:green',label='single shell (whole annulus)'); Zp.plot([],[],lw=4,color='tab:orange',label='split inner/outer'); Zp.legend(fontsize=20)
Zp.set_xlim(0,360); Zp.tick_params(labelsize=20); Zp.set_xlabel('arrival direction β (deg, 0 = away from CoE)',fontsize=24)
Zp.set_ylabel('T − 2.725 K  (mK)',fontsize=24); Zp.set_title(f'Split-shell T(direction), zoomed: range {(Tsp.min()-2.725)*1e3:+.3f} to {(Tsp.max()-2.725)*1e3:+.3f} mK',fontsize=24)
M=axs[-1]; Q=json.load(open('/workspace/cosmology-course/notes/catwise-density-levelset-fit.json'))['migkas_fits']['soft_fixed15.2']
src2=open('/workspace/cosmology-course/notes/catwise-density-levelset-migkas.py').read(); n2={'__name__':'x'}; exec(src2[:src2.index('def main')],n2)
ar=n2['_pc_args']('soft',Q['DR'],Q['LD'],Q['shape']); rq=n2['pc'].run('soft',ar,Q['mc'],0.0); af=rq['af']/rq['af'].mean()
mod=n2['shifted'](af,Q['l0']); Y=n2['MA']; im=n2['IMIG']
R2=1-np.sum((Y-mod[im])**2)/np.sum((Y-Y.mean())**2); assert abs(R2-0.7427756)<1e-5
l0=Q['l0']; th=(np.arange(360)-l0)%360; o=np.argsort(th); thd=(n2['LDEG']-l0)%360
HS=json.load(open('/workspace/cosmology-course/notes/paper_repro/splitA.json')); bth=np.radians(np.arange(360))
Hb=HS['H_in']*(1-np.cos(bth))/2+HS['H_out']*(1+np.cos(bth))/2; afs=af*Hb/Hb.mean(); afs/=afs.mean()
modS=n2['shifted'](afs,Q['l0']); W=n2['W']
c2=float(np.sum(W*(mod[im]-Y)**2)); c2s=float(np.sum(W*(modS[im]-Y)**2)); R2s=1-np.sum((Y-modS[im])**2)/np.sum((Y-Y.mean())**2)
json.dump(dict(R2_orig=float(R2),chi2_orig=c2,R2_split=float(R2s),chi2_split=c2s,max_abs_diff=float(np.max(abs(modS-mod)))),open('/workspace/cosmology-course/notes/paper_repro/splitA_fit.json','w'),indent=1)
M.plot(th[o],mod[o],lw=4,color='tab:orange',label=f'soft-core fit (repro): R² = {R2:.3f}, χ² = {c2:.2f}')
M.plot(th[o],modS[o],lw=3,ls='--',color='tab:blue',label=f'with split-shell expansion: R² = {R2s:.3f}, χ² = {c2s:.2f}')
M.errorbar(thd,Y,yerr=[n2['MLO'],n2['MHI']],fmt='o',color='k',ms=7,capsize=3,label='Migkas A/A_all (71 bins)')
M.set_xlim(0,360); M.tick_params(labelsize=20); M.legend(fontsize=20)
M.set_xlabel('arrival direction β (deg, 0 = away from CoE)   [β = ℓ − 15.2°]',fontsize=24); M.set_ylabel('A / A_all',fontsize=24)
M.set_title('Migkas A/A_all vs the same angle',fontsize=26)
sec=M.secondary_xaxis('top',functions=(lambda x:(x+l0),lambda x:(x-l0))); sec.set_xlabel('Galactic ℓ (deg)',fontsize=20); sec.tick_params(labelsize=17)
fig.tight_layout(); fig.savefig('/workspace/cosmology-course/notes/paper_repro/cmb_annulus3.png',dpi=70)
for r in rows: print({k:round(v,5) for k,v in r.items()})
