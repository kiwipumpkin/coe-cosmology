import json,numpy as np,warnings;warnings.filterwarnings('ignore')
from fit import make,L as LL,Y,S
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':16})
from fit import MLO,MHI
fig,ax=plt.subplots(figsize=(13,8));o=np.argsort(LL);ax.errorbar(LL,Y,yerr=[MLO,MHI],fmt='o',color='k',ms=4,label='Migkas A/A_all (71 bins)')
nm={'soft':'soft core+tail','exp':'exponential','eta':'η-edge'}
for k,c in zip(['soft','exp','eta'],['C0','C2','C3']):
    d=json.load(open(f'fit_{k}.json'));m=make(k)[0](np.array(list(d['best'].values())))
    ax.plot(LL[o],m[o],c,lw=2.5,label=f"{nm[k]} (3% lock): χ²={d['chi2']:.2f}, R²={d['R2']:.3f}")
ax.set_xlabel('Galactic longitude (deg)');ax.set_ylabel('A / A_all');ax.legend(loc='lower center',bbox_to_anchor=(.5,1),frameon=False);fig.tight_layout();fig.savefig('density_fit_migkas.png',dpi=120)
