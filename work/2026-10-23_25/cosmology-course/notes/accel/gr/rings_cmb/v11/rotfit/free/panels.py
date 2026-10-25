import sys;sys.path.insert(0,'..');import numpy as np,healpy as hp,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from fit import sky,vec
F=json.load(open('free_fit.json'))
MK=[((264.02,48.25),'r','^'),((84.02,-48.25),'b','v'),((257.73,57.4),'m','D'),((77.73,-57.4),'m','D')]
for key,tag,fn in (('A_full_dipole','fit A (full model dipole matched to Planck dipole)','rotation_panels_A.png'),('B_axial_dipole_removed','fit B (axial dipole removed; residual matched to Planck dipole + quadrupole axis)','rotation_panels_B.png')):
    f=F[key];cl,cb=f['CoE_lb'];ml,mb=f['motion_lb']
    steps=[('Panel 1: model frame, CoE axis at the pole',(0,-90),None),(f'Panel 2: after tilt by {90+cb:.1f}°\nCoE now at (0°, {cb:+.1f}°)',(0,cb),None),
           (f'Panel 3: after yaw by {cl:.1f}°\nCoE now at ({cl:.1f}°, {cb:+.1f}°)',(cl,cb),None),(f'Panel 4: roll 0° (no effect, sky is axisymmetric)\n+ 370 km/s motion toward ({ml:.1f}°, {mb:+.1f}°)',(cl,cb),vec(ml,mb))]
    fig=plt.figure(figsize=(20,12))
    for i,(t,(l,b),v) in enumerate(steps):
        m=sky(vec((l+180)%360,-b),v)*1e3
        hp.mollview(m,sub=(2,2,i+1),fig=fig.number,title=t,cmap='RdBu_r',unit='dT/T (×10⁻³)',min=-130,max=130);hp.graticule(dpar=30,dmer=60,alpha=.3)
        hp.projscatter(l,b,lonlat=True,c='k',marker='*',s=400,edgecolors='w');hp.projscatter((l+180)%360,-b,lonlat=True,c='k',marker='X',s=200,edgecolors='w')
        if i: [hp.projscatter(*p,lonlat=True,c=c,marker=s,s=150,edgecolors='w') for p,c,s in MK]
        if v is not None: hp.projscatter(ml,mb,lonlat=True,c='lime',marker='o',s=220,edgecolors='k')
    [a.title.set_fontsize(15) for a in fig.axes];fig.subplots_adjust(hspace=.35,top=.84)
    h=[plt.Line2D([],[],c='k',marker='*',ls='',ms=16),plt.Line2D([],[],c='k',marker='X',ls='',ms=12)]+[plt.Line2D([],[],c=c,marker=s,ls='',ms=12) for _,c,s in MK[:3]]+[plt.Line2D([],[],c='lime',marker='o',ls='',ms=13,markeredgecolor='k')]
    fig.legend(h,['CoE (hot end of model sky)','anti-CoE','Planck dipole hot','Planck dipole cold','Planck quadrupole axis (both ends)','fitted motion'],loc='upper center',ncol=3,fontsize=14,frameon=False,bbox_to_anchor=(.5,.99))
    fig.text(.5,1.0,f"{tag}: quadrupole offset {f['quad_offset_deg']:.1f}°, dipole {f['dipole_amp_over_Planck']:.2f}× Planck",ha='center',fontsize=15)
    fig.savefig(fn,dpi=75,bbox_inches='tight');plt.close(fig)
