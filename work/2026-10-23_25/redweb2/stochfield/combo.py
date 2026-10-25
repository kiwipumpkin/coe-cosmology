import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,sys
sys.path.insert(0,'/workspace/rings220/hp')
from fwd import Fwd
d=json.load(open('fit_l30_bg.json'));c=json.load(open('chosen_l30_bg.json'));p=c['p']
F=Fwd(R=4.0,lmax=400);Df=F.Dl(p);L=F.l
bg=np.zeros(len(L));B=np.array(d['Dfix']);bg[:len(B)]=B
# hump: Gaussian measured from map spectra (sigma 32 and 8), peak 5000
H32=5000*np.exp(-.5*((L-217)/31)**2);H8=5000*np.exp(-.5*((L-220)/8)**2)
P=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt')
th=np.loadtxt([f for f in __import__('glob').glob('/workspace/planck/*minimum*')][0])
fig,ax=plt.subplots(figsize=(12,5.5))
ax.errorbar(P[:,0],P[:,1],yerr=[P[:,2],P[:,3]],fmt='.',color='0.65',ms=2,lw=.4,label='Planck 2018 TT',zorder=1)
ax.plot(th[:,0],th[:,1],'k-',lw=1.2,label='ΛCDM')
ax.plot(L,(bg+Df+H32),'C3-',lw=1.8,label='CoE: background + stochastic field + ℓ220 hump (σℓ=32)')
ax.plot(L,(bg+Df+H8),'C1--',lw=1.2,label='same, hump σℓ=8')
ax.set(xscale='log',xlim=(2,400),ylim=(0,6500),xlabel='ℓ',ylabel='D_ℓ [μK²]',title='Combined CoE spectrum: low-ℓ model + characteristic hump')
ax.legend(loc='lower center',fontsize=8);fig.tight_layout();fig.savefig('combo.png',dpi=110)
print('Dfix shape',B.shape,'field D at l=30,100,200,400:',[round(float(Df[l-2]),1) for l in (30,100,200,399)],'hump32 at l=30:',H32[28],'l=100:',H32[98])
