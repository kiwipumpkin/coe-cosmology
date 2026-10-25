import numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':13})
img=plt.imread('/workspace/migkas_src/Ho-map-75-FINAL-ANISOTROPY-HAMMER.png')
X=lambda lon:-np.radians(((np.asarray(lon)+180)%360)-180)
m=hp.read_map('/workspace/cosmology/cmb/commander_I_n128_uK.fits')
lo=np.linspace(-179.75,179.75,720);la=np.linspace(-90,90,361);LO,LA=np.meshgrid(lo,la)
cmb=m[hp.ang2pix(128,(LO%360),LA,lonlat=True)]
fig=plt.figure(figsize=(16,5.6),facecolor='white')
a1=fig.add_axes([0.0,0.13,0.5,0.78]);a1.imshow(img);a1.axis('off')
a1.set_title('LOW z: Migkas et al. 2020, A&A 636, A15, Fig. 21 (original, unmodified)\nH₀ from L_X–T, θ=75° cones, 313 clusters · Galactic, ℓ=0 centre, ℓ increases left',fontsize=12)
a2=fig.add_axes([0.52,0.2,0.46,0.68],projection='hammer')
pm=a2.pcolormesh(-np.radians(LO),np.radians(LA),cmb,cmap='RdBu_r',vmin=-300,vmax=300,shading='auto',rasterized=True)
for (l0,l1,ec,fc,al) in ((260,335,'steelblue','lightsteelblue',.45),(15,80,'goldenrod','gold',.35)):
  bb=np.radians(np.linspace(-90,90,181))
  for e in (l0,l1):a2.plot(X(np.full(181,e)),bb,color=ec,lw=2)
  LL,BB=np.meshgrid(np.linspace(l0,l1,131),bb)
  a2.pcolormesh(X(LL),BB,np.ones_like(LL),cmap=matplotlib.colors.ListedColormap([fc]),alpha=al,shading='auto',zorder=3)
a2.set_xticks(np.radians([-150,-90,0,90,150]));a2.set_xticklabels(['150°','90°','0°','270°','210°']);a2.grid(True,alpha=.4)
a2.set_title('HIGH z (z≈1100): Planck 2018 Commander CMB (Hammer, same orientation)',fontsize=12,pad=14)
c2=fig.add_axes([0.62,0.12,0.26,0.025]);fig.colorbar(pm,cax=c2,orientation='horizontal',label='ΔT [μK]')
h_=[plt.Rectangle((0,0),1,1,fc='gold',alpha=.4,ec='goldenrod',lw=2,label='ℓ≈50° deviation, ℓ=15–80° (left map: ℓ 15–80° lies just left of 0°, toward the 90° label)'),
    plt.Rectangle((0,0),1,1,fc='lightsteelblue',alpha=.5,ec='steelblue',lw=2,label='ℓ≈305° deviation, ℓ=260–335° (left map: from just right of 0° to past 270°, incl. Migkas low-H₀ region ~(280°,−20°))')]
fig.legend(handles=h_,loc='lower left',bbox_to_anchor=(0.01,0.0),ncol=1,frameon=False,fontsize=10.5)
fig.savefig('sky_migkas_orig_cmb.png',dpi=300,facecolor='white')
