import numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,sys
sys.path.insert(0,'../../rotfit');import fit
for T,lab in (('D42','D = 42 Gpc (our model): surface at ≈12.0 Gpc from us'),('D14p3','surface rescaled to ≈14.2 Gpc (Chris ring scale)')):
    m=np.load(f'map_{T}.npy');ns=hp.npix2nside(len(m));fig=plt.figure(figsize=(16,13))
    hp.mollview(m/1e-5,fig=fig.number,sub=(2,1,1),title=f'Smooth model + seed-3 web, full sky, dipole removed by motion [dT/T, 1e-5]\n{lab}',cmap='RdBu_r',min=-1500,max=1500,unit='dT/T ×1e5')
    hp.graticule(30,alpha=.3)
    for (l,b),s,c in (((264.02,48.25),'Planck hot','r'),((84.02,-48.25),'Planck cold','b'),((180,-15),'CoE','m'),((0,15),'anti-CoE','g')):
        hp.projscatter(l,b,lonlat=True,marker='*',s=250,c=c,edgecolors='k');hp.projtext(l,b+6,s,lonlat=True,fontsize=12,weight='bold')
    ax=fig.add_subplot(2,1,2);L=np.arange(0,360,0.25)
    # ring: great circle through CoE and anti-CoE (galactic b-tilted plane containing both, through Planck hot side)
    a=fit.vec(0,15);h=fit.vec(264.02,48.25);e2=h-(h@a)*a;e2/=np.linalg.norm(e2)
    th=np.radians(L);v=np.cos(th)[:,None]*a+np.sin(th)[:,None]*e2
    sm=np.asarray(hp.remove_dipole(np.load(f'web_{T}.npy')*0));
    y=hp.get_interp_val(m,*hp.vec2ang(v))/1e-5
    ax.plot(L,y,'k',lw=1.4,label='smooth + web (after motion)')
    ax.axhline(0,c='gray',lw=.6);ax.set_xlabel('angle along ring from anti-CoE (deg)  [ring: plane of anti-CoE and Planck hot]',fontsize=13);ax.set_ylabel('residual dT/T ×1e5',fontsize=13)
    ax.axvline(180,c='m',ls='--',label='CoE');ax.axvline(0,c='g',ls='--',label='anti-CoE');ax.legend(fontsize=12);ax.set_xlim(0,360)
    ax.set_title('Residual around one great circle (Fig. 55 style)',fontsize=14)
    plt.savefig(f'fullsky_web_{T}.png',dpi=110,bbox_inches='tight');plt.close()
