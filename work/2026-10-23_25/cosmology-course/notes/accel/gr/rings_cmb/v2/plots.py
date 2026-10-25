import pickle,numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':13})
R=pickle.load(open('v2_weak.pkl','rb'));TH=np.arange(0,360,2.0)
for (k,G) in [(1.6,12.0),(0.8,2.0),(1.6,0.0)]:
    fig,ax=plt.subplots(2,1,figsize=(12,10),sharex=True)
    for gr in [('pow',0.05),('pow',0.1),('pow',0.2),('exp',50.0)]:
        x=R[((k,G),gr)];T=x[:,0];d=T/np.nanmean(T)-1;z=np.nan_to_num(d);m=np.fft.rfft(z)/len(z)
        nm=f'T0∝(r0/3)^-{gr[1]}' if gr[0]=='pow' else f'T0∝exp(-(r0-3)/{gr[1]:.0f})'
        ax[0].plot(TH,d,lw=2,label=f'{nm}: RMS={np.sqrt(np.nanmean(d**2)):.2e}, range={np.nanmax(d)-np.nanmin(d):.2e}, dip={2*abs(m[1]):.2e}, quad={2*abs(m[2]):.2e}, emission spread={np.nanmax(x[:,1])-np.nanmin(x[:,1]):.3f}')
        ax[1].plot(TH,x[:,1],lw=2)
    ax[0].set_ylabel('ΔT/T');ax[1].set_ylabel('emission time (CoE frame)');ax[1].set_xlabel('sight angle, observer frame (deg; 0 = anti-CoE, 180 = toward CoE)')
    ax[0].legend(loc='lower center',bbox_to_anchor=(.5,1.02),frameon=False,fontsize=9.5);ax[0].set_title(f'k={k}, GM={G}, c=20, observed t=3, observer recombines t=2.95, proper SR',pad=80)
    fig.tight_layout();fig.savefig(f'v2_k{k}_GM{G}.png',dpi=100);plt.close()
