import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
exec(open('/tmp/spots.py').read());MM=json.load(open('/tmp/M25.json'))
Z={'A':[22,99,262,310],'B':[14,104,266,321]}
wrap=lambda l:((np.asarray(l,float)+180)%360)-180
def moll(l,b):
  lam=np.radians(wrap(l));b=np.radians(np.asarray(b,float));th=b.copy()
  for _ in range(60):th=th-(2*th+np.sin(2*th)-np.pi*np.sin(b))/(2+2*np.cos(2*th)+1e-12)
  return -2*np.sqrt(2)/np.pi*lam*np.cos(th), np.sqrt(2)*np.sin(th)   # l increases left
def loop1():
  c=np.radians([329,17.5]);r=np.radians(58);t=np.linspace(0,2*np.pi,721)
  bl=np.arcsin(np.sin(c[1])*np.cos(r)+np.cos(c[1])*np.sin(r)*np.cos(t))
  dl=np.arctan2(np.sin(t)*np.sin(r)*np.cos(c[1]),np.cos(r)-np.sin(c[1])*np.sin(bl))
  return np.degrees(c[0]+dl)%360,np.degrees(bl)
bb=np.linspace(-89.9,89.9,181)
fig,axs=plt.subplots(2,1,figsize=(15,15),facecolor='white')
for ax,X in zip(axs,'AB'):
  G=np.interp(np.arange(0,360.5,1),np.r_[MM['l'][-1]-360,MM['l'],MM['l'][0]+360],np.r_[MM[X][-1],MM[X],MM[X][0]]);gm=np.abs(G).max()
  for q in range(360):
    if q==180 or q==179:
      pass
    xa,ya=moll(np.full_like(bb,q+.02),bb);xb,yb=moll(np.full_like(bb,q+.98),bb)
    if np.abs(xa-xb).max()>1:continue
    ax.fill(np.r_[xa,xb[::-1]],np.r_[ya,yb[::-1]],color='red' if G[q]>0 else 'blue',alpha=.6*abs(G[q])/gm,lw=0)
  t=np.linspace(0,2*np.pi,400);ax.plot(2*np.sqrt(2)*np.cos(t),np.sqrt(2)*np.sin(t),'k',lw=1.5)
  for lg in (0,90,270):
    x,y=moll(np.full_like(bb,lg),bb);ax.plot(x,y,color='0.75',lw=.8,ls=':')
    ax.text(moll(lg,0)[0],-1.55,f'{lg}°',ha='center',fontsize=12)
  x,y=moll(np.arange(0.,360,1),np.zeros(360));ax.plot(np.sort(x),np.zeros(360),color='0.75',lw=.8,ls=':')
  for z in Z[X]:x,y=moll(np.full_like(bb,z),bb);ax.plot(x,y,'k',lw=2)
  l,b=loop1();x,y=moll(l,b);br=np.where(np.abs(np.diff(x))>1)[0]+1
  for s in np.split(np.arange(len(x)),br):ax.plot(x[s],y[s],'k--',lw=1.3)
  h=[]
  for i,(nm,sl,sb,c) in enumerate(SPOTS[X]):
    xs,ys=moll(sl,sb);xs,ys=float(xs),float(ys)
    ax.plot(xs,ys,'o',ms=28,mfc=CC[c],mec='k',mew=2,zorder=7);ax.text(xs,ys,str(i+1),ha='center',va='center',color='white',fontsize=15,fontweight='bold',zorder=8)
    ax.annotate(nm,(xs,ys),xytext=(-20,14) if xs>0 else (20,14),textcoords='offset points',ha='right' if xs>0 else 'left',fontsize=12.5,fontweight='bold',color=CC[c],zorder=9,bbox=dict(fc='white',ec=CC[c],boxstyle='round,pad=.25'))
    h.append(Line2D([],[],ls='',marker='o',mfc=CC[c],mec='k',ms=13,label=f'{i+1} {nm} ({sl}°, {sb:+g}°)'))
  h+=[plt.Rectangle((0,0),1,1,fc='red',alpha=.45,label=f'{X} residual high'),plt.Rectangle((0,0),1,1,fc='blue',alpha=.45,label=f'{X} residual low'),
      Line2D([],[],color='k',lw=2,label='zero crossings: '+', '.join(f'{z}°' for z in Z[X])),Line2D([],[],color='k',ls='--',label='Loop I / NPS')]
  ax.legend(handles=h,loc='upper center',bbox_to_anchor=(.5,-.02),ncol=3,fontsize=11.5,frameon=False)
  ax.set_title(f'{X} residual after n=1–3 (shade ∝ ±25° mean)  —  '+('distance/velocity spots' if X=='A' else 'foreground/sample spots'),fontsize=15)
  ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-2.9,2.9);ax.set_ylim(-1.65,1.5)
fig.suptitle('Galactic Mollweide, ℓ = 0 at centre, ℓ increases to the left',fontsize=13,y=.995)
fig.tight_layout();fig.savefig('overlay_stack_AB.png',dpi=130,facecolor='white')
