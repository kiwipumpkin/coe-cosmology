import pickle,numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':13})
R=pickle.load(open('rings_gal.pkl','rb'));TH=np.arange(0,360,2.0)
for (k,G) in [(1.6,12.0),(0.8,2.0),(1.6,0.0),(0.8,0.0)]:
  for c in (20.0,50.0):
    fig,ax=plt.subplots(2,1,figsize=(12,10),sharex=True)
    for tr in (2.85,2.9,2.95):
        x=R[((k,G),c,tr)]['res'];T=x[:,0]
        if np.isnan(T).all():continue
        d=T/np.nanmean(T)-1;z=np.nan_to_num(d);m=np.fft.rfft(z)/len(z)
        ax[0].plot(TH,d,lw=2,label=f'local recomb. t={tr}: RMS={np.sqrt(np.nanmean(d**2)):.2e}, range={np.nanmax(d)-np.nanmin(d):.2e}, dip={2*abs(m[1]):.2e}, quad={2*abs(m[2]):.2e}')
        ax[1].plot(TH,x[:,1],lw=2)
    ax[0].set_ylabel('ΔT/T');ax[1].set_ylabel('emission time');ax[1].set_xlabel('sight angle (deg; 0 = anti-CoE, 180 = toward CoE)')
    ax[0].legend(loc='lower center',bbox_to_anchor=(.5,1.02),frameon=False,fontsize=10);ax[0].set_title(f'k={k}, GM={G}, c={c}, observed at t=3',pad=70)
    fig.tight_layout();fig.savefig(f'rings_k{k}_GM{G}_c{c:.0f}.png',dpi=100);plt.close()
