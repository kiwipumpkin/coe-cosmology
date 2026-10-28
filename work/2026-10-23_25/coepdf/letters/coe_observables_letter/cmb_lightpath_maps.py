# Mollweide maps (Galactic, l increases left) of the light-path dT/T prediction, R_see run, from cmb_lightpath_maps.npz (cmb_lightpath_compute.py).
import numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
m=np.load('cmb_lightpath_maps.npz');T=2.7255e6
fig=plt.figure(figsize=(14,8.4))
P=[('geo','A. Geometric term (kinematic + gravitational at emission), dT/T',1,None),
   ('path','B. Along-path term (ISW/RS): identically 0 - model potential is static',1,(-1,1)),
   ('tot','C. Total minus monopole, dT/T',1,None),
   ('tot_nodip','D. Total, monopole+dipole removed, $\\mu$K',T,None)]
for i,(k,t,s,r) in enumerate(P):
    x=m[k]*s
    if k in('geo','tot'):x=x-x.mean()
    kw=dict(min=r[0],max=r[1]) if r else {}
    hp.mollview(x,fig=fig.number,sub=(2,2,i+1),title=t,cmap='coolwarm',**kw,notext=False)
    hp.projscatter(180,-15,lonlat=True,marker='*',s=120,c='purple');hp.projscatter(260,60,lonlat=True,marker='x',s=80,c='k')
    hp.graticule(dpar=30,dmer=60,alpha=.3)
fig.text(.5,.01,'R_see = 0.851 Gpc emission sphere; star = CoE (180,-15); x = axis of evil (260,60). Monopole of A: -828 (dT/T, linear).',ha='center')
plt.savefig('cmb_lightpath_maps.png',dpi=110,bbox_inches='tight')
