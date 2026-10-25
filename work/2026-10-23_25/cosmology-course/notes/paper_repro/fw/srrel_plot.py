import sys,json,numpy as np
sys.argv=['x','srrel','0.99']; src=open('cmb_flat.py').read(); exec(src[:src.index("res=dict(framework")])
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,(a1,a2)=plt.subplots(2,1,figsize=(12,10))
import core
a1.errorbar(core.LDEG,core.Y,yerr=[core.MLO,core.MHI],fmt='o',color='k',ms=4,alpha=.6,label='Migkas cluster data')
for fn,st,lab in (('migkas_srrel_0.99.json','m-','SR, front 0.99c, observer speed free'),('migkas_srrel_0.99_near0.9.json','r--','SR, front 0.99c, observer near front (speed ≥ 0.9 of front)')):
    m=json.load(open(fn))['1deg']; a1.plot(np.arange(360),m['curve'],st,lw=2.5,label=f"{lab}: R²={m['R2']:.3f}, χ²={m['chi2']:.1f}, D={m['D_Gpc']:.0f} Gpc, our speed {m['beta_obs_now']:.3f}c")
a1.set_xlim(0,360); a1.set_ylabel('A / A_all'); a1.grid(alpha=.3); a1.legend(fontsize=8); a1.set_title('Special relativity, fast layers (front at 0.99 of light speed)')
info={}
L=(np.arange(360)+p['l0'])%360; o=np.argsort(L)
for mass,col in ((1,'r'),(0,'b')):
    T,X,V,G,c,Phi,Xo,Vo,meta=layers(mass)
    for f,ls in ((1.0,'-'),(0.0,':')):
        Rp,Lp,Pp,out=trace(f,T,X,V,G,c,Xo); best=None
        for Rr in np.geomspace(1e-9,0.999*Xo,300):
            fa=factors(Rr,f,Rp,Lp,Pp,T,X,V,c,Phi,Xo,Vo)
            if fa is None: continue
            m_=wmean(fa['T'])
            k_=(round(frac(fa['T']),3),-m_)
            if best is None or k_>best[3]: best=(m_,Rr,fa,k_)
        m_,Rr,fa,_=best; full=ring(fa['T']); d=(full-m_)/m_/1e-5
        info[f'mass{mass}_f{f}']=dict(lowest_mean_T_K=m_,R_rec_Mpc=Rr*Dphys*1e3,sky_fraction=frac(fa['T']),grav=float(fa['grav']),Phi_c2=[float(Phi(Rr)[0]),float(Phi(Xo)[0])],
            doppler=[float(np.nanmin(fa['dop'])),float(np.nanmax(fa['dop']))],stretch=[float(np.nanmin(fa['S'])),float(np.nanmax(fa['S']))],
            swing_1e5=float((np.nanmax(full)-np.nanmin(full))/2/m_/1e-5),obs_beta=float(Vo),
            los_toward=float(np.sum((fa['ue']-fa['uo'])[180]*nb[180])) if np.isfinite(fa['te'][180]) else None,los_away=float(np.sum((fa['ue']-fa['uo'])[0]*nb[0])) if np.isfinite(fa['te'][0]) else None,
            delay_Gyr=[float((T[-1]-np.nanmax(fa['te']))*Dphys/core.C_GPC_GYR),float((T[-1]-np.nanmin(fa['te']))*Dphys/core.C_GPC_GYR)])
        a2.plot(L[o],d[o],col+ls,lw=2,label=f"{'With' if mass else 'No'} mass, {'all motion is space growing (f=1)' if f==1 else 'all motion through space (f=0)'}: light seen over {100*frac(fa['T']):.0f}% of sky, mean {m_:.0f} K (cannot reach 2.725 K)")
a2.set_xlim(0,360); a2.set_xlabel('Galactic longitude (degrees)'); a2.set_ylabel('Temperature change, parts per 100,000'); a2.grid(alpha=.3); a2.legend(fontsize=8)
a2.set_title('CMB vs direction at the recombination radius giving the most sky coverage (gaps = no recombination light)')
fig.tight_layout(); fig.savefig('cmb_sr_rel.png',dpi=130); json.dump(info,open('cmb_srrel_plotinfo.json','w'),indent=1,default=float); print(json.dumps(info,indent=1,default=float))
