import os,sys,json,numpy as np
sys.path.insert(0,'/workspace/cosmology-course/notes/age_refit'); os.chdir('/workspace/cosmology-course/notes/age_refit')
b=json.load(open('gr2/v4/best/best_eta_fwd_gtc_1.json'))
os.environ['XOVR']=json.dumps(b['x']); sys.argv=['cmb4.py','eta','fwd','gtc','30']
exec(open('cmb4.py').read().split("t0=time.time(); ells=")[0])
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
plt.rcParams.update({'font.size':19,'axes.titlesize':24,'axes.labelsize':21})
K=b['K_local_today_K']; RAYS=json.load(open('storyboard/sb_rays.json'))
OUT='storyboard'; Cg=F.C
Gyr=lambda t:(t-tbo)/Cg*lam            # time since our shell's bang, Gyr
tcode=lambda g: tbo+g*Cg/lam
Hk=lambda t: L.Rdot(t,X)[0]/L.R(t,X)[0]*Cg/lam*977.79   # km/s/Mpc
xs=np.linspace(np.log(1e-3)+1,np.log(3e4)-1,1200)
rs=2*L.m(np.array([L.xmin]))[0]*lam    # Schwarzschild radius of M, Gpc
def Telem(t): return K*(L.rho_rel(t,xs)/rho_o)**(1/3)
def page(fig,cap,n,tag='COMPUTED'):
    fig.text(0.5,0.015,cap,ha='center',va='bottom',fontsize=17,wrap=True,bbox=dict(fc='#f4f4f4',ec='0.6'))
    fig.text(0.99,0.985,tag,ha='right',va='top',fontsize=16,color='white',bbox=dict(fc='darkgreen' if tag.startswith('COMP') else 'darkorange'))
    fig.text(0.01,0.985,f'Storyboard {n}/8 — η forward, U>c (eta_fwd_gtc_1)',ha='left',va='top',fontsize=16,color='0.3')
    fig.savefig(f'{OUT}/sb_{n:02d}.png',dpi=100); plt.close(fig)
PAR=(f"M = {b['M']:.2e} kg   η density law   p = {b['x'][2]:.3f}   lf = {b['x'][3]:.3f}   lL = {b['x'][4]:.3f}   s = {b['x'][5]:.3f}   u = {b['x'][1]:.3f}\n"
     f"Outputs: D_today = {b['D_Gpc']:.2f} Gpc, U_obs/c = {b['Uo_c']:.2f}, age = {b['age_obs_Gyr']:.2f} Gyr, Migkas χ² = {b['chi2']:.1f} ({b['sigma']:.1f}σ), K = {K:.2f} K")
def snap(n,g,title,cap,show_rays=True,lim=None):
    t=tcode(g); R=L.R(t,xs)*lam; T_=Telem(t); Ro_=L.R(t,X)[0]*lam
    fig=plt.figure(figsize=(16,10)); ax=fig.add_axes([0.03,0.21,0.48,0.66],projection='polar'); a2=fig.add_axes([0.66,0.2,0.32,0.7])
    ok=np.isfinite(R)&np.isfinite(T_)&(R>0); Ri=R[ok]; Ti=T_[ok]; o_=np.argsort(Ri)
    lo,hi=-3.5,2.0                     # display radius = log10(R/Gpc)
    rr=np.linspace(lo,hi,300); Tr=np.exp(np.interp(rr*np.log(10),np.log(Ri[o_]),np.log(Ti[o_]),left=np.nan,right=np.nan))
    th=np.linspace(0,2*np.pi,181); cm=plt.cm.inferno; nrm=matplotlib.colors.LogNorm(1e3,1e5)
    pc=ax.pcolormesh(th,rr-lo,np.tile(Tr[:-1],(180,1)).T,cmap=cm,norm=nrm,shading='flat')
    hot=Ti>=3000
    if hot.any() and (~hot).any():
        k=np.where(np.diff(hot[o_].astype(int))!=0)[0]
        for kk in k: ax.plot(th,np.full_like(th,np.log10(Ri[o_][kk])-lo),'c--',lw=3)
        ax.text(np.radians(60),np.log10(Ri[o_][k[0]])-lo,' 3000 K',color='cyan',fontsize=17)
    ax.plot(0,np.log10(Ro_)-lo,'*',ms=26,color='lime',mec='k'); ax.text(np.radians(-12),np.log10(Ro_)-lo,'our shell',color='lime',fontsize=17)
    if show_rays:
        for k_,v in RAYS.items():
            q=v['rec']
            if q and Gyr(q['t'])<=g:
                ax.plot(q['phi'],np.log10(q['R']*lam)-lo,'o',ms=11,color='cyan',mec='k')
                if float(k_)%40==0 or float(k_)>=144: ax.text(q['phi'],np.log10(q['R']*lam)-lo+0.25,f"{float(k_):g}°",color='cyan',fontsize=13)
    ax.set_rlim(0,hi-lo); ax.set_rticks(np.arange(0,hi-lo+0.01,1)[1:]); ax.set_yticklabels([f'1e{int(v+lo) if (v+lo)%1==0 else v+lo:g}' if False else f'{10**(v+lo):.3g}' for v in np.arange(0,hi-lo+0.01,1)[1:]],color='w',fontsize=13)
    ax.set_rlabel_position(100); ax.set_thetagrids([0,90,180,270],['0° (away)','90°','180° (CoE dir.)','270°'],fontsize=15)
    ax.set_title(title,pad=28)
    ax.text(np.radians(225),0.15,'M',color='w',fontsize=16)
    cb=fig.colorbar(pc,cax=fig.add_axes([0.54,0.2,0.012,0.7])); cb.set_label('T_elem (K)')
    a2.loglog(Ri[o_],Ti[o_],lw=3,color='orangered'); a2.axhline(3000,color='c',ls='--',lw=2); a2.axvline(Ro_,color='lime',lw=2)
    a2.set_xlabel('shell radius R (Gpc)'); a2.set_title(f'T_elem(R) at t = {g*1e3:.3g} Myr',fontsize=20); a2.grid(alpha=.3)
    a2.text(0.03,0.05,f"our shell: R = {Ro_:.3g} Gpc\nH = Ṙ/R = {Hk(t):.3g} km/s/Mpc\nU = Ṙ = {L.Rdot(t,X)[0]:.3g} c",transform=a2.transAxes,fontsize=17,bbox=dict(fc='w',alpha=.85))
    page(fig,cap+"\nLeft: polar view, radius axis = log10(R/Gpc) (log-compressed; the CoE/M sits at the centre, below 3e-4 Gpc).",n)
# ---- 1
g1=3e-5
snap(1,g1,'1. Early: shells around the CoE, all plasma',
 "Exact LTB dust + central mass M:  (dR/dτ)² = 2Gm(r)/R + 2E(r).  Each shell's radiation: T_elem = K·(ρ/ρ_obs,today)^{1/3}.  Hot (>3000 K) = plasma.\n"+PAR,show_rays=False)
# ---- 2-4 recombination sweep
rt={k:Gyr(v['rec']['t']) for k,v in RAYS.items() if v['rec']}
for n,(g,tt) in zip([2,3,4],[(rt['0']*1.0005,'2. Recombination reaches the 0° side first'),(rt['135']*1.001,'3. Front moves toward the CoE side'),(rt['148.5']*1.001,'4. Last: rays just outside the CoE shadow')]):
    snap(n,g,tt,"Cyan dots = points where the light we receive today in direction θ last scattered (T_elem = 3000 K), shown once that event has happened.\n"
      f"Emission times: 0°: {rt['0']*1e3:.2f} Myr, 80°: {rt['80']*1e3:.2f}, 120°: {rt['120']*1e3:.2f}, 144°: {rt['144']*1e3:.2f}, 148.5°: {rt['148.5']*1e3:.2f} Myr → 0° side recombines before CoE side.\n"
      f"The cyan 3000 K front sits INSIDE the emission points: outer shells (behind us) cool first, the front moves inward toward the CoE.  "
      f"R_critical (0° emission point) = {RAYS['0']['rec']['R']/Ro:.3f} D_today: behind us, outside our shell.  Global U ≫ c there (U_e ≈ {b['rows'][0]['U_e_c']:.0f} c).",lim=None)
# ---- 5 today + H history
fig=plt.figure(figsize=(16,10)); ax=fig.add_axes([0.07,0.22,0.4,0.68]); a2=fig.add_axes([0.57,0.22,0.4,0.68])
gs=np.logspace(-6,np.log10(b['age_obs_Gyr']),200); Hs=[Hk(tcode(g)) for g in gs]; Rs=[L.R(tcode(g),X)[0]*lam for g in gs]; Us=[L.Rdot(tcode(g),X)[0] for g in gs]
ax.loglog(gs,Hs,lw=3); ax.set_xlabel('t since our shell’s bang (Gyr)'); ax.set_ylabel('H of our shell (km/s/Mpc)'); ax.grid(alpha=.3); ax.set_title('H starts high and falls')
ax.plot(gs[-1],Hs[-1],'*',ms=22,color='g'); ax.text(gs[-1],Hs[-1]*3,f'today {Hs[-1]:.1f}',ha='right',fontsize=18)
a2.loglog(gs,Rs,lw=3,label='R of our shell (Gpc)'); a2.loglog(gs,Us,lw=3,label='U = Ṙ  (units of c)'); a2.axhline(1,color='k',ls=':')
a2.set_xlabel('t since our shell’s bang (Gyr)'); a2.legend(); a2.grid(alpha=.3); a2.set_title('Our shell: D(t) and U(t) are outputs')
page(fig,f"5. Today: D = {b['D_Gpc']:.2f} Gpc from the CoE, age {b['age_obs_Gyr']:.2f} Gyr, shell speed U_obs = {b['Uo_c']:.2f} c (>c, global, coordinate-dependent).\n"
 f"Integrated forward from small D;  H = Ṙ/R of our shell today = {Hs[-1]:.1f} km/s/Mpc (sky-mean H scaled to 70).  Observer moves with the local dust flow (comoving),\n"
 "so the local observer speed relative to the flow is 0 < c here; our kinetic CMB dipole is removed (m=1) before scoring.",5)
# ---- 6 light paths
fig=plt.figure(figsize=(16,10)); ax=fig.add_axes([0.06,0.2,0.44,0.72]); a2=fig.add_axes([0.56,0.2,0.42,0.72])
for A,lim,ttl in [(ax,1.25,'Past light rays today → last scattering (R/D_today, φ)'),(a2,0.16,'Zoom near the mass')]:
    for k,v in RAYS.items():
        R=np.array(v['R'])/Ro; ph=np.array(v['phi']); m=np.isfinite(R)
        if v['rec']: 
            tr=np.array(v['t']); m&=tr>=v['rec']['t']
        col='crimson' if not v['rec'] else ('darkorange' if float(k)>=144 else 'steelblue')
        A.plot(R[m]*np.cos(ph[m]),R[m]*np.sin(ph[m]),lw=2.2,color=col)
        if v['rec']:
            q=v['rec']; A.plot(q['R']/Ro*np.cos(q['phi']),q['R']/Ro*np.sin(q['phi']),'o',ms=9,color=col,mec='k')
            if (A is ax and float(k)<144 and float(k)%40==0) or (A is a2 and float(k)>=144): A.text(q['R']/Ro*np.cos(q['phi'])+lim*0.02,q['R']/Ro*np.sin(q['phi'])+lim*0.03,f"{float(k):g}°: {rt[k]*1e3:.1f} Myr",fontsize=13)
    A.plot(0,0,'ko',ms=10); A.plot(1,0,'g*',ms=22); A.set_aspect('equal'); A.set_xlim(-0.15 if A is ax else -lim,1.2 if A is ax else lim); A.set_ylim(-0.25 if A is ax else -lim*0.5,0.8 if A is ax else lim); A.set_title(ttl,fontsize=20); A.grid(alpha=.3)
ax.text(1,-0.08,'us',color='g',ha='center'); ax.text(0,-0.08,'CoE/M',ha='center')
a2.text(0.02,0.95,'red = shadow rays (150°,165°,180°): no 3000 K point,\n   end falling into M (integrator stops)\norange = edge rays 144–148.5°: emitter azimuth\n   swept to 102–149° (bent around M)',transform=a2.transAxes,va='top',fontsize=14,bbox=dict(fc='w',alpha=.9))
sh=json.load(open('gr2/v4/shadow/sh_eta_fwd_gtc_1.json'))
page(fig,f"6. Light paths: past null geodesics from us (exact LTB, LSODA), 1+z = k^t(emit)/k^t(obs).  Axes = areal radius R at each point's time / D_today.  Upper half-plane only (sky is mirror-symmetric about the axis).\n"
 f"CoE shadow: θ = {sh['edges'][0]:.1f}°–{sh['edges'][1]:.1f}° (width {sh['shadow_width_deg']:.1f}°), rays fall into M.  Near the edge light is bent ~180° round M (primary/secondary image ring),\n"
 f"observed hot: T = {max(g['T_K'] for g in sh['near_edge']):.1f} K at the edge vs 2.725 K mean.  Secondary (flipped) image not traced separately. r_s(M) = {rs*1e3:.3g} Mpc.",6)
# ---- 7 sky
sp=b['score_partial']; gl=np.array(sp['gal_l']); dT=np.array(sp['dTT']); o=np.argsort(gl); gl,dT=gl[o],dT[o]
th=np.radians(gl); Xp=np.c_[np.ones(len(gl)),np.cos(th),np.sin(th)]; c=np.linalg.lstsq(Xp,dT,rcond=None)[0]; dTn=dT-Xp[:,1:]@c[1:]
md=b['migkas_data']; G=np.array(b['migkas_model'])
fig=plt.figure(figsize=(16,10)); ax=fig.add_axes([0.07,0.22,0.9,0.68])
ax.errorbar(md['l'],np.array(md['y'])-1,md['err'],fmt='o',color='0.5',ms=6,label='Migkas A/A_all − 1 (data)')
ax.plot(np.arange(360),G-1,lw=3,color='k',label=f'model G−1 (χ²={b["chi2"]:.1f})')
cl=(180+b['l0'])%360
ax.plot(np.arange(360),-(G-1),lw=2,ls=':',color='k',label='−(G−1): perfect cancellation target')
gp=np.r_[gl,cl];dp=np.r_[dTn,np.nan];oo=np.argsort(gp)
ax.plot(gp[oo],dp[oo],'o-',lw=3,color='crimson',label='CMB dT/T, m=1 removed (6° grid)')
ax.axvspan(cl-31.1,cl+31.1,color='0.3',alpha=0.3,label='shadow (no last-scattering)')
ax.axvline(cl,color='0.3',ls='--'); ax.text(cl,0.55,'CoE',ha='center',fontsize=18)
ax.set_ylim(-0.3,0.7); ax.set_xlim(0,360); ax.set_xlabel('galactic longitude l (deg)'); ax.set_ylabel('fractional deviation'); ax.legend(fontsize=15,loc='upper left'); ax.grid(alpha=.3)
ax.set_title('7. Resulting sky: CMB vs Migkas')
page(fig,f"7. Model CoE at l = 180 + l0 = {cl:.0f}° (fitted l0 = {b['l0']:.1f}°).  Away from the CoE the CMB runs opposite to G−1 (cancellation direction), but at 6° resolution the edge ring\n"
 f"(up to +{dTn.max()*100:.0f}%) dominates: rms dT/T (m=1 removed) = {sp['model_dTT_rms_m1rem']*100:.1f}% vs target {sp['target_rms_m1rem']*100:.1f}%; partial-sky cancellation {sp['cancel']:.0f}% over {sp['n']}/60 directions.\n"
 "Shadow directions excluded from the score (a choice, not physics).  Real CMB: ~1e-5.",7)
# ---- 8 comparison M scan
rows=[]
for f in sorted(os.listdir('gr2/v4/shadow')):
    if f.startswith('shM_') and f.endswith('.json'):
        s=json.load(open('gr2/v4/shadow/'+f))
        if 'M' not in s or s.get('shadow_width_deg') is None: continue
        rows.append((s['M'],s['shadow_width_deg'],s['chi2'],s['Uo_c'],max([g['T_K'] for g in s['near_edge']] or [np.nan]) if s['near_edge'] else np.nan))
rows.append((b['M'],sh['shadow_width_deg'],b['chi2'],b['Uo_c'],max(g['T_K'] for g in sh['near_edge'])))
rows=sorted(r for r in rows if r[1] is not None)
fig=plt.figure(figsize=(16,10)); ax=fig.add_axes([0.08,0.27,0.38,0.62]); a2=fig.add_axes([0.58,0.27,0.38,0.62])
Ms=[r[0] for r in rows]
ax.semilogx(Ms,[r[1] for r in rows],'o-',lw=3,ms=10); ax.set_xlabel('M (kg)'); ax.set_ylabel('shadow width (deg)'); ax.grid(alpha=.3); ax.set_title('Shadow vs mass (Migkas refit each M)')
for r in rows: ax.text(r[0],r[1]+2,f"χ²{r[2]:.1f}",fontsize=13,ha='center')
a2.semilogx(Ms,[r[4] for r in rows],'s-',lw=3,ms=10,color='crimson'); a2.axhline(2.725,ls=':',color='k'); a2.set_xlabel('M (kg)'); a2.set_ylabel('hottest edge ray on 2° grid (K)'); a2.grid(alpha=.3); a2.set_title('Edge-ring temperature')
page(fig,"8. Comparison: M scanned (u,p,lf,lL,s refit to Migkas, K held at 8.46 K).  M = 2.5e53 kg: shadow 6.4°, χ² 14.1 vs 62.2°, χ² 16.9 at the baseline 2.5e54 kg.\n"
 "Ring T is from the 2° shadow-probe grid (not edge-resolved, peak likely higher).  The full CMB rerun for 2.5e53 kg (K re-solved, fine edge sampling) is not in these files yet.\n"
 "M ≤ 2.5e52 kg: no buildable U>c universe found; 7.9e52 refit failed (χ² 373).",8)
print('done')
