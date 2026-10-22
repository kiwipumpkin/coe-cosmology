import json,glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from fwd import Fwd
d=json.load(open('fit_l30_bg.json'));p=json.load(open('chosen_l30_bg.json'))['p']
F=Fwd(R=4.0,lmax=400);Df=F.Dl(p);L=F.l;bg=np.zeros(len(L));B=np.array(d['Dfix']);bg[:len(B)]=B
th=np.loadtxt(glob.glob('/workspace/planck/*minimum*')[0]);Lt=th[:,1][np.searchsorted(th[:,0],L)]
P=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt')
def h(l,a,c,s1,s2):return a*np.where(l<c,np.exp(-.5*((l-c)/s1)**2),np.exp(-.5*((l-c)/s2)**2))
m=(L>=60)&(L<=330);po,_=curve_fit(h,L[m],(Lt-bg-Df)[m],p0=[5000,220,60,50])
print('hump a,c,s_left,s_right',np.round(po,1))
tot=bg+Df+h(L,*po);z=np.load('../realweb/mc/mc_mean_passing.npz');mc=z['D4.0'].mean(0)
fig,ax=plt.subplots(2,1,figsize=(12,10))
a=ax[0];a.errorbar(P[:,0],P[:,1],yerr=[P[:,2],P[:,3]],fmt='.',color='0.65',ms=2,lw=.4,label='Planck 2018 TT',zorder=1)
a.plot(th[:,0],th[:,1],'k-',lw=1.2,label='ΛCDM');a.plot(L,tot,'C3-',lw=1.8,label=f'CoE: bg + field + hump (ℓc={po[1]:.0f}, σL={po[2]:.0f}, σR={po[3]:.0f})')
a.set(xlim=(2,400),ylim=(0,6500),ylabel='D_ℓ [μK²]',title='Combined CoE spectrum, hump widened to ΛCDM near ℓ≈200');a.legend(loc='lower center',fontsize=8)
a=ax[1];Pl=P[P[:,0]<=30];a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=4,label='Planck 2018 TT',zorder=5)
a.plot(th[:,0],th[:,1],'k-',lw=1.2,label='ΛCDM')
s=L<=30;cv=np.sqrt(2/(2*L+1))*Df
a.fill_between(L[s],(tot-cv)[s],(tot+cv)[s],color='C3',alpha=.2,label='±1σ cosmic variance (field)')
a.plot(L[s],tot[s],'C3-o',ms=3,label='CoE bg + field ensemble mean')
a.plot(np.arange(len(mc))[2:],mc[2:],'C0-s',ms=3,label='CoE web, passing-universe mean R=4 Gpc')
a.plot(L[s],bg[s],'r--',lw=1,label='deterministic CoE background')
a.set(xlim=(1.5,30.5),ylim=(0,2600),xlabel='ℓ',ylabel='D_ℓ [μK²]',title='Low multipoles, linear');a.legend(loc='upper right',fontsize=8)
fig.tight_layout();fig.savefig('combo2.png',dpi=110)
