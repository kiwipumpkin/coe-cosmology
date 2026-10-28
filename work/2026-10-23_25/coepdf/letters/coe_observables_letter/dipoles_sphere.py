# 3D Galactic sphere: each direction (filled) and its antipode (open), joined through the centre by a line of the same colour.
import os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__))
exec(open(H+'/dipoles_sky_map.py').read().split('def xy')[0].split('H=os.path')[1].split('\n',1)[1])  # loads P
def v(l,b):l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
fig=plt.figure(figsize=(16,8))
for k,(el,az) in enumerate([(0,-10),(0,170)]):
  ax=fig.add_subplot(1,2,k+1,projection='3d')
  u,w=np.mgrid[0:2*np.pi:60j,0:np.pi:30j]
  ax.plot_surface(np.cos(u)*np.sin(w),np.sin(u)*np.sin(w),np.cos(w),color='lightblue',alpha=.08,lw=0)
  t=np.linspace(0,2*np.pi,200)
  for B in (-60,-30,0,30,60):
    r=np.cos(np.radians(B));ax.plot(r*np.cos(t),r*np.sin(t),np.sin(np.radians(B))+0*t,'k',lw=1.2 if B==0 else .5,alpha=.6 if B==0 else .3)
  for L in range(0,360,30):
    bb=np.radians(np.linspace(-90,90,100));ax.plot(np.cos(bb)*np.cos(np.radians(L)),np.cos(bb)*np.sin(np.radians(L)),np.sin(bb),'k',lw=.5,alpha=.25)
  for L,lab in [(0,'l=0'),(90,'l=90'),(180,'l=180'),(270,'l=270')]:p=1.12*v(L,0);ax.text(*p,lab,fontsize=9,color='gray')
  ax.text(0,0,1.12,'NGP (b=+90)',fontsize=9,color='gray');ax.text(0,0,-1.18,'SGP (b=-90)',fontsize=9,color='gray');ax.scatter(0,0,0,color='k',s=30)
  for n,l,b,c,mk,anti in P:
    a=v(l,b);ax.scatter(*a,color=c,marker=mk,s=150 if mk!='*' else 300,edgecolors='k',depthshade=False,label=f'{n} ({l:g}°, {b:+g}°)' if k==0 else None)
    if anti:
      z=-a;ax.plot([a[0],z[0]],[a[1],z[1]],[a[2],z[2]],'-',color=c,lw=2.2)
      ax.scatter(*z,facecolors='none',edgecolors=c,marker=mk,s=130 if mk!='*' else 260,linewidths=2,depthshade=False)
    else:ax.plot([0,a[0]],[0,a[1]],[0,a[2]],':',color=c,lw=1.5)
  ax.set_box_aspect((1,1,1));ax.set_xlim(-1,1);ax.set_ylim(-1,1);ax.set_zlim(-1,1);ax.view_init(el,az);ax.set_axis_off()
  ax.set_title(f'edge-on view from l = {az%360}° (b = 0)',fontsize=11)
fig.legend(loc='lower center',ncol=4,fontsize=9,frameon=False)
fig.suptitle('Dipole/axis directions on the Galactic unit sphere. Filled: published; open: antipode; solid line: pair through centre; dotted: no antipode plotted',fontsize=11)
fig.savefig(H+'/dipoles_sphere.png',dpi=130,facecolor='white',bbox_inches='tight')
