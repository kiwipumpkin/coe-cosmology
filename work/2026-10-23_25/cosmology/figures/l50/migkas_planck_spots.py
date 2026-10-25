exec(open('/tmp/spots.py').read());OV=[]
MM=__import__('json').load(open('/tmp/M25.json'))
import json,sys;RR=json.load(open('/tmp/R.json'));X=sys.argv[1]
import matplotlib;matplotlib.rcParams['hatch.linewidth']=2.2
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
  G=np.interp(np.arange(0,360.5,1),np.r_[MM['l'][-1]-360,MM['l'],MM['l'][0]+360],np.r_[MM[X][-1],MM[X],MM[X][0]]);gm=np.abs(G).max()
  for q in range(360):
    if q in (179,):pass
    u0,u1=q+.02,q+.98;xa,ya=f(np.full_like(bb,u0),bb);xb,yb=f(np.full_like(bb,u1),bb)
    if np.abs(xa-xb).max()>300:continue
    ax.fill(np.r_[xa,xb[::-1]],np.r_[ya,yb[::-1]],color='red' if G[q]>0 else 'blue',alpha=.6*abs(G[q])/gm,lw=0,zorder=3)
  for L0,col in (((ZA,'red'),) if X=='A' else ((ZB,'blue'),)):
    for z in L0:
      xa,ya=f(np.full_like(bb,z),bb);ax.plot(xa,ya,color=col,lw=lw*1.4,zorder=4) if False else ax.plot(xa,ya,color='k',lw=lw*1.2,zorder=4)
  for e in (179.95,180.05):
    xe,ye=f(np.full_like(bb,e),bb);ax.plot(xe,ye,color='0.3',lw=1.2,zorder=4)
  xq,yq=f(np.arange(0.,360,1),np.zeros(360));ax.plot(xq,yq,color='0.6',lw=.8,ls=':',zorder=4)
  if ax in OV:
    for i,(t,sl,sb,c) in enumerate(SPOTS[X]):
      xs,ys=f(np.array([float(sl)]),np.array([float(sb)]));xs,ys=np.ravel(xs)[0],np.ravel(ys)[0]
      ax.plot(xs,ys,'o',ms=24,mfc=CC[c],mec='k',mew=2,zorder=7);ax.text(xs,ys,str(i+1),ha='center',va='center',color='white',fontsize=13,fontweight='bold',zorder=8);ax.annotate(t,(xs,ys),xytext=(-16,12),textcoords='offset points',ha='right',fontsize=10.5,fontweight='bold',color=CC[c],zorder=9,bbox=dict(fc='white',ec=CC[c],alpha=.85,boxstyle='round,pad=.2'))
  l,b=loop1();x,y=f(l,b);br=np.where(np.abs(np.diff(x))>200)[0]+1
  for s in np.split(np.arange(len(x)),br):ax.plot(x[s],y[s],'k--',lw=lw*.7,zorder=5)
  return f
fig=plt.figure(figsize=(16,10.2),facecolor='white')
t1=fig.add_axes([0.0,0.55,0.5,0.39]);t1.imshow(plt.imread('/workspace/migkas_src/Ho-map-75-FINAL-ANISOTROPY-HAMMER.png'));t1.axis('off')
t1.set_title('LOW z: Migkas et al. 2020, A&A 636, A15, Fig. 21 — original H₀ map (θ=75° cones)',fontsize=11.5)
t2=fig.add_axes([0.505,0.55,0.49,0.39]);t2.imshow(plt.imread('/workspace/planck_img/p_white.png'));t2.axis('off')
t2.set_title('HIGH z: ESA/Planck Collaboration, Planck 2018 CMB (official image)',fontsize=11.5)
fig.text(.5,.975,'Original maps (no overlay)',ha='center',fontsize=14,fontweight='bold');fig.text(.5,.505,f'With overlay: {X} residual high/low regions, zero crossings, Loop I',ha='center',fontsize=14,fontweight='bold')
a1=fig.add_axes([0.0,0.08,0.5,0.39]);a1.imshow(np.ones_like(plt.imread('/workspace/migkas_src/Ho-map-75-FINAL-ANISOTROPY-HAMMER.png')[...,:3]));a1.axis('off')
OV.append(a1);overlay(a1,mig_xy,1.6);x,y=mig_xy(329,30);a1.text(x[0],y[0],'Loop I / NPS',ha='center',fontsize=10,fontweight='bold',zorder=6)
a1.set_title('',fontsize=1) if True else a1.set_title('LOW z: Migkas et al. 2020, A&A 636, A15, Fig. 21 — original H₀ map (θ=75° cones)\nzero-crossing meridians + Loop I overlaid; edges approximate (±few°) on this panel',fontsize=11.5)
a2=fig.add_axes([0.505,0.08,0.49,0.39]);a2.imshow(np.full(plt.imread('/workspace/planck_img/p.jpg').shape,255,np.uint8));a2.axis('off')
OV.append(a2);overlay(a2,pl_xy,1.6);x,y=pl_xy(329,30);a2.text(x,y,'Loop I / NPS',ha='center',fontsize=10,fontweight='bold',color='k',zorder=6)
for l in (0,90,270):x,y=pl_xy(l,0);a2.text(x,y+2600*0+1260,'' ,fontsize=1)
a2.set_title('',fontsize=1) if True else a2.set_title('HIGH z: ESA/Planck Collaboration, Planck 2018 CMB (official image)\nGalactic Mollweide, ℓ=0 centre, ℓ increases left',fontsize=11.5)
for l in (90,0,270):
  x,y=pl_xy(l,-90);a2.text(pl_xy(l,0)[0],2400,f'{l}°',ha='center',va='top',fontsize=10)
h_=[plt.Rectangle((0,0),1,1,fc='red',alpha=.3,label=f'{X} high (shade ∝ mean residual, ±25°)'),plt.Rectangle((0,0),1,1,fc='blue',alpha=.3,label=f'{X} low (shade ∝ mean residual, ±25°)'),Line2D([],[],color='k',lw=2,label='zero crossings: '+', '.join(f'{z}°' for z in (ZA if X=='A' else ZB))),
    Line2D([],[],color='k',ls='--',lw=1.4,label='Loop I / NPS: small circle centre (329°, +17.5°), r = 58° (Berkhuijsen+1971)')]
h_+= [Line2D([],[],ls='',marker='o',mfc=CC[c],mec='k',ms=12,label=f'{i+1} {t} ({sl}°, {sb:+g}°)') for i,(t,sl,sb,c) in enumerate(SPOTS[X])];fig.legend(handles=h_,loc='lower center',ncol=4,frameon=False,fontsize=10.5)
fig.savefig('sky_spots_'+X+'.png',dpi=300,facecolor='white')
