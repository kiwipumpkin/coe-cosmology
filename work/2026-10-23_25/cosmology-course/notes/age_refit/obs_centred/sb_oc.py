import pickle,json,sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
sys.path.insert(0,'..')
for tag in ('oc_m254','oc_m253'):
    d=pickle.load(open(f'{tag}_rays.pkl','rb')); J=json.load(open(f'{tag}.json')); RD=d['RD'];tn=d['tn'];lam=d['lam']; C=299792.458 if False else None
    S=J['sync_chris']; 
    if 'rows' not in S: continue
    import forward as F; tst=S['t_star_Myr']/1e3/lam*F.C   # code units
    # emission point (x,y) relative to observer at each emission time: use phi and R along ray, CoE at origin, observer at (Ro,0)
    import ltb2
    fr=np.linspace(0,1,7); fig,ax=plt.subplots(2,4,figsize=(20,10)); ax=ax.ravel()
    for k,f_ in enumerate(fr):
        tobs=tst+f_*(tn-tst)   # observation time; circles = points on observer's past cone at t* seen at tobs
        a=ax[k]
        # approximation for drawing: light-travel radius c*(tobs-t*) around us (code units c=1), shown in units of D
        r=(tobs-tst)/d['Ro']; th=np.linspace(0,2*np.pi,300)
        a.plot(r*np.cos(th),r*np.sin(th),'r-',label='c(t-t*) circle')
        for rr in np.linspace(0,r,5)[1:-1]: a.plot(rr*np.cos(th),rr*np.sin(th),':',c='orange',lw=.8)
        a.plot(0,0,'b*',ms=12); a.plot(-1,0,'kx',ms=10); a.text(-1,0.05,'CoE')
        a.set_title(f't-t* = {(tobs-tst)/F.C*lam*1e3:.3g} Myr'); a.set_aspect('equal'); L_=max(1.3,1.1*r); a.set_xlim(-L_,L_); a.set_ylim(-L_,L_)
    # last panel: actual GR emission points of today's CMB (sync model), observer-centred coords
    a=ax[7]; rows=[r for r in S['rows'] if 'T_obs_K' in r]
    x=[]; y=[]; c=[]
    for r in rows:
        R=r['R_e_over_D']; p=np.radians(r['phi_e_deg']); x.append(R*np.cos(p)-1); y.append(R*np.sin(p)); c.append(r['T_obs_K'])
    sc=a.scatter(np.array(x)+1-1,y,c=c,cmap='coolwarm'); plt.colorbar(sc,ax=a,label='T_obs K')
    a.plot(0,0,'b*',ms=12); a.plot(-1,0,'kx'); a.set_aspect('equal'); a.set_title('GR: where today\'s CMB was emitted (at t*)\n(us at origin, CoE at x=-1, units D)')
    fig.suptitle(f'{tag}: growing visible circle from our recombination crossing t* (Chris construction); t*={S["t_star_Myr"]:.3g} Myr')
    plt.tight_layout(); plt.savefig(f'sb_{tag}.png',dpi=80)
print('ok')
