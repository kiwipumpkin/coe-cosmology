# v3B: Mollweide of total dT/T (monopole removed) and total with monopole+dipole removed [uK]. Data via cmb_lightpath_v3_common.py.
import numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from cmb_lightpath_v3_common import *
m,m1,cl,cl1,kind,mean=stats('B');fig=plt.figure(figsize=(14,4.6))
hp.mollview(m,fig=fig.number,sub=(1,2,1),title='v3B total dT/T, monopole removed',cmap='coolwarm')
hp.projscatter(180,-15,lonlat=True,marker='*',s=120,c='purple');hp.projscatter(260,60,lonlat=True,marker='x',s=80,c='k')
hp.mollview(m1*T0uK,fig=fig.number,sub=(1,2,2),title='v3B total, monopole+dipole removed [uK]',cmap='coolwarm')
hp.projscatter(180,-15,lonlat=True,marker='*',s=120,c='purple');hp.projscatter(260,60,lonlat=True,marker='x',s=80,c='k')
plt.savefig('cmb_lightpath_v3B_maps.png',dpi=110,bbox_inches='tight')
