import json,glob,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from fwd import Fwd
d=json.load(open('fit_l30_bg.json'));p=json.load(open('chosen_l30_bg.json'))['p']
LM=2500;F=Fwd(R=4.0,lmax=LM,kmax=1000,nk=600);Df=F.Dl(p);L=F.l;np.save('Df_2500.npy',Df)
bg=np.zeros(len(L));B=np.array(d['Dfix']);bg[:len(B)]=B
def h(l,a,c,s1,s2):return a*np.where(l<c,np.exp(-.5*((l-c)/s1)**2),np.exp(-.5*((l-c)/s2)**2))
tot=bg+Df+h(L,4973.3,219.6,83.9,92.4)
th=np.loadtxt(glob.glob('/workspace/planck/*minimum*')[0]);P=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt')
fig,ax=plt.subplots(2,1,figsize=(13,10))
for a,sc in zip(ax,['linear','log']):
    a.errorbar(P[:,0],P[:,1],yerr=[P[:,2],P[:,3]],fmt='.',color='0.65',ms=2,lw=.3,label='Planck 2018 TT',zorder=1)
    a.plot(th[:,0],th[:,1],'k-',lw=1.2,label='ΛCDM');a.plot(L,tot,'C3-',lw=1.6,label='CoE: bg + field + ℓ220 hump')
    a.set(xscale=sc,xlim=(2,LM),ylim=(0,6500),ylabel='D_ℓ [μK²]',title=f'Full spectrum, {sc} ℓ');a.legend(loc='upper right',fontsize=8)
ax[1].set_xlabel('ℓ');fig.tight_layout();fig.savefig('combo3.png',dpi=110)
print([round(float(tot[l-2]),0) for l in (540,800,1100,2000)])
