import numpy as np,json
exec(open('cmbann.py').read().split("meanT=lambda")[0])
d=json.load(open('cmbann.json')); a=d['annulus']; sp=json.load(open('split.json')); qi,qo=sp['q_in'],sp['q_out']; RR=sp['R_rec_Mpc']/1e3/DG; d['R_rec_Mpc']=sp['R_rec_Mpc']
tr=sp['trec_Myr']/1e3/t0; b,S,dop,gm,te=np.load('cmbann_keep.npy')
times=[tr,0.25,0.5,0.75,1.0]
V=lambda t,lo=None,hi=None: r_of(hi if hi else qo,t)**3-r_of(lo if lo else qi,t)**3; V0=V(tr)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,axs=plt.subplots(6,1,figsize=(18,96),gridspec_kw=dict(height_ratios=[1]*5+[0.55]))
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
B=axs[-1]; bb,Tsp,Tsg=np.load('split_T.npy'); B.plot(bb,Tsg,lw=4,label=f'single shell (ray emitter layers), mean {Tsg.mean():.3f} K'); B.plot(bb,Tsp,lw=4,label=f'split inner/outer (cosβ interpolation), mean {Tsp.mean():.3f} K')
B.axhline(2.725,color='k',ls=':'); B.set_xlim(0,360); B.tick_params(labelsize=20); B.legend(fontsize=20)
B.set_xlabel('arrival direction β (deg, 0 = away from CoE)',fontsize=24); B.set_ylabel('T_obs now (K)',fontsize=24); B.set_title('T(direction) now',fontsize=26)
fig.tight_layout(); fig.savefig('cmb_annulus2.png',dpi=70)
for r in rows: print({k:round(v,5) for k,v in r.items()})
