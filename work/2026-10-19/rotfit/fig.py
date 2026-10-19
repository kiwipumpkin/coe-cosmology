import numpy as np,healpy as hp,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from fit import sky,vec,lb,T1,V,T0,nside
S=json.load(open('rotfit_stats.json'));F=S['fit'];vm=vec(*F['motion_lb'])
steps=[('Panel 1: model frame, CoE axis at the pole\n(CoE at b = −90°, anti-CoE at b = +90°)',(0,-90)),
 ('Panel 2: after tilt by 75° (latitude move)\nCoE now at (0°, −15°)',(0,-15)),
 ('Panel 3: after yaw by 180° (longitude move)\nCoE now at (180°, −15°) as in the paper',(180,-15)),
 ('Panel 4: after roll about the CoE axis (fitted roll: none, the sky is\naxisymmetric so roll changes nothing) + 370 km/s motion added',(180,-15))]
MK=[((264.02,48.25),'r','^','Planck dipole hot (264.0°, +48.3°)'),((84.02,-48.25),'b','v','Planck dipole cold (84.0°, −48.3°)'),((257.73,57.4),'m','D','Planck quadrupole axis (257.7°, +57.4°)'),((77.73,-57.4),'m','D',None)]
fig=plt.figure(figsize=(20,12))
for i,(t,(cl,cb)) in enumerate(steps):
    a=vec((cl+180)%360,-cb);m=sky(a,vm if i==3 else None)*1e3
    hp.mollview(m,sub=(2,2,i+1),fig=fig.number,title=t,cmap='RdBu_r',unit='dT/T (×10⁻³)',min=-130,max=130)
    hp.graticule(dpar=30,dmer=60,alpha=.3)
    hp.projscatter(cl,cb,lonlat=True,c='k',marker='*',s=400,edgecolors='w');hp.projscatter((cl+180)%360,-cb,lonlat=True,c='k',marker='X',s=200,edgecolors='w')
    if i>=1:
        for (l,b),c,s,_ in MK: hp.projscatter(l,b,lonlat=True,c=c,marker=s,s=150,edgecolors='w')
    if i==3: hp.projscatter(*F['motion_lb'],lonlat=True,c='lime',marker='o',s=220,edgecolors='k')
h=[plt.Line2D([],[],c='k',marker='*',ls='',ms=16),plt.Line2D([],[],c='k',marker='X',ls='',ms=12)]+[plt.Line2D([],[],c=c,marker=s,ls='',ms=12) for _,c,s,n in MK[:3]]+[plt.Line2D([],[],c='lime',marker='o',ls='',ms=13,markeredgecolor='k')]
[a.title.set_fontsize(15) for a in fig.axes]
fig.subplots_adjust(hspace=.35,top=.86)
fig.legend(h,['CoE (hot end of model sky)','anti-CoE','Planck dipole hot end','Planck dipole cold end','Planck quadrupole axis (both ends)',f"fitted motion direction ({F['motion_lb'][0]:.1f}°, {F['motion_lb'][1]:+.1f}°)"],loc='upper center',ncol=3,fontsize=14,frameon=False,bbox_to_anchor=(.5,1.0))
fig.text(.5,.0,'Own crossing, distributed mass, normalisation b (Migkas-fitted), speed 4c. Galactic Mollweide, l = 0 at centre, l increasing to the left.',ha='center',fontsize=13)
fig.savefig('rotation_panels.png',dpi=80,bbox_inches='tight');plt.close(fig)
# summary
a=vec(0,15);m=sky(a,vm)*1e3;fig=plt.figure(figsize=(16,10))
hp.mollview(m,fig=fig.number,cmap='RdBu_r',unit='dT/T (×10⁻³)',title=f"After fit: model dipole hot end ({F['dipole_hot_lb'][0]:.1f}°, {F['dipole_hot_lb'][1]:+.1f}°), {F['dipole_offset_deg']:.1f}° from Planck hot; amplitude {F['dipole_amp_over_Planck']:.1f}× Planck\nmodel quadrupole axis ({F['quad_axis_lb'][0]:.1f}°, {F['quad_axis_lb'][1]:+.1f}°), {F['quad_offset_deg']:.1f}° from Planck's; chance of doing as well at random p = {S['MC']['p']:.2f} ({S['MC']['sigma_two_sided']:.1f}σ)",min=-130,max=130)
hp.graticule(dpar=30,dmer=60,alpha=.3)
hp.projscatter(180,-15,lonlat=True,c='k',marker='*',s=400,edgecolors='w');hp.projscatter(0,15,lonlat=True,c='k',marker='X',s=200,edgecolors='w')
for (l,b),c,s,_ in MK: hp.projscatter(l,b,lonlat=True,c=c,marker=s,s=170,edgecolors='w')
hp.projscatter(*F['motion_lb'],lonlat=True,c='lime',marker='o',s=220,edgecolors='k')
q=F['quad_axis_lb'];hp.projscatter([q[0],(q[0]+180)%360],[q[1],-q[1]],lonlat=True,c='orange',marker='s',s=170,edgecolors='k')
[a.title.set_fontsize(15) for a in fig.axes]
fig.legend(h+[plt.Line2D([],[],c='orange',marker='s',ls='',ms=12,markeredgecolor='k')],['CoE (180°, −15°)','anti-CoE (0°, +15°)','Planck dipole hot','Planck dipole cold','Planck quadrupole axis','fitted motion direction','model quadrupole axis'],loc='upper center',ncol=4,fontsize=13,frameon=False,bbox_to_anchor=(.5,1.12))
fig.savefig('summary_fit_sky.png',dpi=80,bbox_inches='tight')
