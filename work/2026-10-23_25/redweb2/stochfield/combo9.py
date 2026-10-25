import json,glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=json.load(open('fit_l30_bg.json'));Df=np.load('Df_2500.npy');L=np.arange(2,2501)
bg=np.zeros(len(L));B=np.array(d['Dfix']);bg[:len(B)]=B
def h(l,a,c,s1,s2):return a*np.where(l<c,np.exp(-.5*((l-c)/s1)**2),np.exp(-.5*((l-c)/s2)**2))
tot=bg+Df+h(L,4973.3,219.6,83.9,92.4)
th=np.loadtxt(glob.glob('/workspace/planck/*minimum*')[0]);P=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt')
MC=np.load('../realweb/mc/mc_mean_passing.npz')['D4.0'];mc=MC.mean(0);q16,q84=np.percentile(MC,[16,84],0);lm=np.arange(MC.shape[1])
fig,ax=plt.subplots(2,1,figsize=(13,10))
a=ax[0];Pl=P[P[:,0]<=30];a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=4,label='Planck 2018 TT',zorder=5)
a.plot(th[:,0],th[:,1],'k-',lw=1.2,label='ΛCDM');s=L<=30;cv=np.sqrt(2/(2*L+1))*Df
a.fill_between(L[s],(tot-cv)[s],(tot+cv)[s],color='C3',alpha=.2,label='±1σ cosmic variance (field)')
a.plot(L[s],tot[s],'C3-o',ms=3,label='CoE bg + field + hump ensemble mean')
a.plot(L[s],bg[s],'r--',lw=1,label='deterministic CoE background')
a.set(xlim=(1.5,30.5),ylim=(0,2600),xlabel='ℓ',ylabel='D_ℓ [μK²]',title='Low multipoles, linear');a.legend(loc='upper right',fontsize=8)
a=ax[1];a.errorbar(P[:,0],P[:,1],yerr=[P[:,2],P[:,3]],fmt='.',color='0.65',ms=2,lw=.3,label='Planck 2018 TT',zorder=1)
a.plot(th[:,0],th[:,1],'k-',lw=1.2,label='ΛCDM');a.fill_between(L,tot-cv,tot+cv,color='C3',alpha=.2,label='±1σ cosmic variance');a.plot(L,tot,'C3-',lw=1.6,label='CoE: bg + field + ℓ220 hump')
a.set(xscale='log',xlim=(2,2500),ylim=(0,6500),xlabel='ℓ',ylabel='D_ℓ [μK²]',title='Full spectrum, log ℓ');a.legend(loc='upper left',fontsize=8)
fig.tight_layout();fig.savefig('combo9.png',dpi=110)
