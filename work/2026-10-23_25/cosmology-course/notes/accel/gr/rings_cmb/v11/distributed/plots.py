import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
obs=np.load('../../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs))
psi=np.arange(0,181,1.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(0.5))/4;w=w/w.sum();Gy=42*3.26156
def st(fn):
    d=np.load(fn);z=d['o'][:,1];assert np.isfinite(z).all();T=1/z;dd=T/np.sum(w*T)-1
    X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],dd*sw,rcond=None)[0];r=dd-X@c
    t0=float(d['t0']);tc=float(d['tc']) if 'tc' in d else (float(d['te']) if 'te' in d else np.nan)
    return dict(d=dd,r=r,rms=float(np.sqrt(np.sum(w*dd*dd))),rng=float(np.ptp(dd)),dip=float(abs(c[1])),rrng=float(np.ptp(r)),z=(float(z.min()),float(z.max())),t0=t0*Gy,tc=tc*Gy)
R={'own_nomass':st('../z4_nomass.npz'),'sim_nomass':st('../z4_sim.npz'),'own_point':st('../z4_m02.npz'),'sim_point':st('../z4_sim_m02.npz'),
   'own_a':st('z4_dist_a_own.npz'),'sim_a':st('z4_dist_a_sim.npz'),'own_b':st('z4_dist_b_own.npz'),'sim_b':st('z4_dist_b_sim.npz')}
for k,s in R.items(): print(f"{k:11s} range {s['rng']:.6f} ({s['rng']/PR:.2f}x) RMS {s['rms']:.6f} dipole {s['dip']:.6f} leftover {s['rrng']:.3e} ({s['rrng']/PR:.3f}x) 1+z {s['z'][0]:.6f}-{s['z'][1]:.6f} today {s['t0']:.3f} crossing {s['tc']:.3f}")
for n in ('a','b'):
    for sc in ('own','sim'):
        a,b=R[f'{sc}_{n}'],R[f'{sc}_nomass'];print(f"{n} {sc}: mass vs no mass RMS {100*(a['rms']/b['rms']-1):+.2f}% range {100*(a['rng']/b['rng']-1):+.2f}%")
json.dump({k:{q:v for q,v in s.items() if q not in('d','r')} for k,s in R.items()},open('dist_stats.json','w'),indent=1)
A=np.r_[psi,360-psi[-2::-1]];m=lambda y:np.r_[y,y[-2::-1]]
def lab(s,name):
    if s['rng']<1e-5: return f"{name}: 1+z = 2.000 everywhere, dT/T = 0 (range {s['rng']/PR:.1e}× Planck, numerical noise)"
    return f"{name}: range {s['rng']/PR:.1f}× Planck, RMS {s['rms']:.4f}, dipole {s['dip']:.4f}, left after dipole {s['rrng']/PR:.2f}× Planck"
for n,suf,mdesc,today in (('b','','distributed mass, Migkas-fitted amount (4.36 × critical inside our radius)','15.8'),('a','_equalmass','distributed mass, same amount inside our radius as the 0.2-of-critical point mass','29.0')):
    C=[('own_nomass',f'own crossing, no mass (today 34.2 Gyr)','C0'),(f'own_{n}',f'own crossing, {mdesc} (today {today} Gyr)','C3'),
       ('sim_nomass','simultaneous at our crossing 17.12 Gyr, no mass (today 34.2 Gyr)','C2'),(f'sim_{n}',f"simultaneous at our crossing {R[f'sim_{n}']['tc']:.2f} Gyr, same distributed mass (today {today} Gyr)",'C1')]
    for key,out,yl,tt in (('d',f'../dTT_4c_4scenarios_distributed{suf}.png','dT/T (relative to full-sky mean)','4 c, D = 42 Gpc, recombination radius 0.5 D; mass = soft core+tail law (CatWISE 3% lock), ρ ∝ 1/(1+(r/L)^0.3), L = 464 D; T ∝ 1/(1+z)'),
                          ('r',f'../dTT_4c_4scenarios_distributed{suf}_residual.png','dT/T left after removing the dipole (cosine fit)','same four runs, dipole removed')):
        fig,ax=plt.subplots(figsize=(13,8))
        for k_,nm,col in C: ax.plot(A,m(R[k_][key]),color=col,lw=2.5,label=lab(R[k_],nm))
        ax.axhline(0,color='k',lw=.6);ax.grid(alpha=.3);ax.set_xlim(0,360);ax.set_xticks(range(0,361,30))
        ax.set_xlabel('angle on the sky from the anti-CoE direction (0°/360° = away from CoE, 180° = toward CoE; mirrored)');ax.set_ylabel(yl)
        ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=9,frameon=False);ax.set_title(tt,fontsize=9.5,pad=-14)
        fig.savefig(out,dpi=110,bbox_inches='tight');plt.close(fig)
# density plot
LD=18159.848409348455/39.10676196793444;al=0.30000009809958916;Rb=1.5*(1+1/(42/14.30048287260589))
x=np.linspace(0,2.3,2000);rho=1/(1+(x/LD)**al)*(x<=Rb);rho1=1/(1+(1/LD)**al)
fig,ax=plt.subplots(figsize=(10,6));ax.plot(x,rho/rho1,'k',lw=2.5,label=f'soft core+tail: ρ/ρ(us) = [1+(1/464)^0.3]/[1+(r/(464 D))^0.3]; ρ(centre)/ρ(us) = {1/rho1:.3f}; cut to zero at {Rb:.2f} D (as in the fit code)')
ax.axvline(1,color='C0',ls='--');ax.text(1.02,0.45,'us (D = 42 Gpc)',color='C0');ax.axvline(0.5,color='C3',ls='--');ax.text(0.52,0.25,'recombination radius 0.5 D',color='C3')
ax.set_xlim(0,2.3);ax.set_ylim(0,1.3);ax.grid(alpha=.3);ax.set_xlabel('distance from CoE ÷ D');ax.set_ylabel('density ÷ density at our radius')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=9,frameon=False);fig.savefig('../density_soft_law_4c.png',dpi=110,bbox_inches='tight')
d=14.30048287260589/42;print('rho(0)/rho(D)',1/rho1*1,'... centre/us',(1+(1/LD)**al),'delta at D=42 with d=14.3',((1/(1+((1-d)/LD)**al))-(1/(1+((1+d)/LD)**al)))*(1+(1/LD)**al),'Rb',Rb)
