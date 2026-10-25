import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
from matplotlib.lines import Line2D
ZA=[22,99,262,310];ZB=[14,104,266,321]
BANDS=[(10,85,'darkorange','gold',.40),(260,335,'navy','deepskyblue',.30)]
wrap=lambda l:((np.asarray(l,float)+180)%360)-180
def hammer(lam,b):
  d=np.sqrt(1+np.cos(b)*np.cos(lam/2));return np.cos(b)*np.sin(lam/2)/d, np.sin(b)/d   # normalised to [-1,1]
def lam_eff(l):  # Hammer longitude whose equator x matches the image's linear tick labels (0,90,270,edge=180)
  t=wrap(l)/180.;return 0 if t==0 else brentq(lambda L:hammer(L,0)[0]-t,-np.pi,np.pi)
def mig_xy(l,b):  # Migkas image pixels: ellipse x 170-1722, y 18-888
  l=np.atleast_1d(l);b=np.atleast_1d(b);L=np.array([lam_eff(v) for v in l])
  X,Y=hammer(L,np.radians(b));return 946-X*776, 453-Y*435
def moll(lam,b):
  th=np.array(b,float)
  for _ in range(50):th=th-(2*th+np.sin(2*th)-np.pi*np.sin(b))/(2+2*np.cos(2*th)+1e-12)
  return lam/np.pi*np.cos(th), np.sin(th)
def pl_xy(l,b):  # ESA image pixels: ellipse x 27-4630, y 26-2332 ; l increases left
  X,Y=moll(np.radians(wrap(l)),np.radians(b));return 2328.5-X*2301.5, 1179-Y*1153
def loop1():
  c=np.radians([329,17.5]);r=np.radians(58);t=np.linspace(0,2*np.pi,721)
  bl=np.arcsin(np.sin(c[1])*np.cos(r)+np.cos(c[1])*np.sin(r)*np.cos(t))
  dl=np.arctan2(np.sin(t)*np.sin(r)*np.cos(c[1]),np.cos(r)-np.sin(c[1])*np.sin(bl))
  return np.degrees(c[0]+dl)%360,np.degrees(bl)
def overlay(ax,f,lw):
  bb=np.linspace(-89.9,89.9,181)
  for L0,col in ((ZA,'red'),(ZB,'blue')):
    for z in L0:
      xa,ya=f(np.full_like(bb,z),bb);ax.plot(xa,ya,color=col,lw=lw*1.4,zorder=4)
  l,b=loop1();x,y=f(l,b);br=np.where(np.abs(np.diff(x))>200)[0]+1
  for s in np.split(np.arange(len(x)),br):ax.plot(x[s],y[s],'k--',lw=lw*.7,zorder=5)
  return f
fig=plt.figure(figsize=(16,5.4),facecolor='white')
a1=fig.add_axes([0.0,0.14,0.5,0.74]);a1.imshow(plt.imread('/workspace/migkas_src/Ho-map-75-FINAL-ANISOTROPY-HAMMER.png'));a1.axis('off')
overlay(a1,mig_xy,1.6);x,y=mig_xy(329,30);a1.text(x[0],y[0],'Loop I / NPS',ha='center',fontsize=10,fontweight='bold',zorder=6)
a1.set_title('LOW z: Migkas et al. 2020, A&A 636, A15, Fig. 21 — original H₀ map (θ=75° cones)\nzero-crossing meridians + Loop I overlaid; edges approximate (±few°) on this panel',fontsize=11.5)
a2=fig.add_axes([0.505,0.14,0.49,0.74]);a2.imshow(plt.imread('/workspace/planck_img/p.jpg'));a2.axis('off')
overlay(a2,pl_xy,1.6);x,y=pl_xy(329,30);a2.text(x,y,'Loop I / NPS',ha='center',fontsize=10,fontweight='bold',color='k',zorder=6)
for l in (0,90,270):x,y=pl_xy(l,0);a2.text(x,y+2600*0+1260,'' ,fontsize=1)
a2.set_title('HIGH z: ESA/Planck Collaboration, Planck 2018 CMB (official image)\nGalactic Mollweide, ℓ=0 centre, ℓ increases left',fontsize=11.5)
for l in (90,0,270):
  x,y=pl_xy(l,-90);a2.text(pl_xy(l,0)[0],2400,f'{l}°',ha='center',va='top',fontsize=10)
h_=[Line2D([],[],color='red',lw=2.2,label='A residual zero crossings: ℓ = 22°, 99°, 262°, 310°'),
    Line2D([],[],color='blue',lw=2.2,label='B residual zero crossings: ℓ = 14°, 104°, 266°, 321°'),
    Line2D([],[],color='k',ls='--',lw=1.4,label='Loop I / NPS: small circle centre (329°, +17.5°), r = 58° (Berkhuijsen+1971)')]
fig.legend(handles=h_,loc='lower center',ncol=3,frameon=False,fontsize=10.5)
fig.savefig('sky_migkas_planck_zeros.png',dpi=300,facecolor='white')
