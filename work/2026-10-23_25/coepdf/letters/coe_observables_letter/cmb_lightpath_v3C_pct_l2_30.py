# v3C: share of l=2..30 power per l (monopole+dipole removed map) vs Planck TT.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from cmb_lightpath_v3_common import *
m,m1,cl,cl1,kind,mean=stats('C');l,pm,pp=pct(cl1)
fig,ax=plt.subplots(figsize=(10,4.5));ax.bar(l-.2,pm,.4,color='C3',label='v3C');ax.bar(l+.2,pp,.4,color='k',label='Planck TT')
ax.set_xlabel('l');ax.set_ylabel('% of sum_2^30 (2l+1)C_l');ax.legend();plt.savefig('cmb_lightpath_v3C_pct_l2_30.png',dpi=110,bbox_inches='tight')
