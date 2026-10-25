import healpy as hp,numpy as np,json
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':16})
alm=hp.read_alm('/workspace/inv/smica_alm_l30.fits');L=hp.Alm.getlmax(len(alm));T0=2.7255e6
ns=64
def lmap(ls):
 f=np.zeros(L+1);f[list(ls)]=1;return hp.alm2map(hp.almxfl(alm,f),ns,verbose=False) if False else hp.alm2map(hp.almxfl(alm,f),ns)
maps={l:lmap([l]) for l in range(2,7)};tot=lmap(range(2,L+1))
cl=hp.alm2cl(alm);print('Cl2..4 uK2',cl[2:5], 'units check D2',2*3*cl[2]/(2*np.pi))
t=np.radians(np.linspace(0,360,721))
def circ(pole_l,pole_b):
 # great circle with given pole (galactic deg)
 p=hp.ang2vec(pole_l,pole_b,lonlat=True);u=np.cross(p,[0,0,1.0]);
 if np.linalg.norm(u)<1e-6:u=np.array([1.0,0,0])
 u/=np.linalg.norm(u);v=np.cross(p,u);vec=np.outer(np.cos(t),u)+np.outer(np.sin(t),v)
 th_,ph_=hp.vec2ang(vec);return lambda m:1e6*hp.get_interp_val(m,th_,ph_)
# search circle where l2 and l3 most anti-correlated
best=None
for pl in range(0,360,10):
 for pb in range(0,90,10):
  c=circ(pl,pb);a,b=c(maps[2]),c(maps[3]);rho=np.ptp(a+b)/(np.ptp(a)+np.ptp(b))
  if best is None or rho<best[0]:best=(rho,pl,pb,np.ptp(a),np.ptp(b),np.ptp(a+b))
print('cancel',best)
sc=1.0 # map units
u=np.degrees(t)
fig,ax=plt.subplots(2,1,figsize=(13,12),sharex=True)
c=circ(0,0) # meridian through Galactic poles (l=0/180)
for l,col in zip(range(2,7),['tab:red','tab:blue','tab:green','tab:orange','tab:purple']):ax[0].plot(u,c(maps[l]),color=col,lw=2,label=f'ℓ={l}')
ax[0].plot(u,c(tot),'k',lw=3,label='sum ℓ=2–30')
ax[0].set_title('A: meridian great circle through Galactic poles (l=0°/180°)')
c=circ(best[1],best[2]);a,b=c(maps[2]),c(maps[3])
ax[1].plot(u,a,'tab:red',lw=2,label='ℓ=2');ax[1].plot(u,b,'tab:blue',lw=2,label='ℓ=3');ax[1].plot(u,a+b,'k',lw=3,label='ℓ=2 + ℓ=3');ax[1].fill_between(u,-1e3,1e3,where=a*b<0,color='grey',alpha=0.15);ax[1].set_ylim(1.2*min(a.min(),b.min(),(a+b).min()),1.2*max(a.max(),b.max(),(a+b).max()))
ax[1].set_title(f'B: great circle with pole (l,b)=({best[1]}°,{best[2]}°): shaded = ℓ2, ℓ3 opposite sign (cancelling)')
for x in ax:x.axhline(0,color='grey');x.set_ylabel('ΔT [µK]');x.legend(loc='upper center',bbox_to_anchor=(0.5,-0.08 if x is ax[0] else -0.12),ncol=6,frameon=False,fontsize=14)
ax[1].set_xlabel('angle along great circle [deg]')
fig.suptitle('Planck SMICA (real sky, a_ℓm to ℓ=30): ΔT along great circles by multipole',fontsize=17)
fig.tight_layout();fig.savefig('cmb_dT_vs_angle_multipoles.png',dpi=130)
json.dump(dict(best=best,ptp_tot=float(np.ptp(circ(0,0)(tot))),cl=cl[2:7].tolist()),open('nums.json','w'),default=float)
