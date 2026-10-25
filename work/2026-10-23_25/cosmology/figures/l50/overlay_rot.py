import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
exec(open('/tmp/spots.py').read());MM=json.load(open('/tmp/M25.json'))
Z={'A':[22,99,262,310],'B':[14,104,266,321]}
TILT,ROLL=15.0,80.0
def vec(l,b):l,b=np.radians(l),np.radians(b);return np.stack([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)],-1)
def rot(ax,a):
  ax=ax/np.linalg.norm(ax);a=np.radians(a);K=np.array([[0,-ax[2],ax[1]],[ax[2],0,-ax[0]],[-ax[1],ax[0],0]])
  return np.eye(3)+np.sin(a)*K+(1-np.cos(a))*K@K
R1=rot(np.array([0,1.,0]),TILT)           # CoE (180,0) -> b=+15
coe=R1@vec(180,0);R=rot(coe,ROLL)@R1      # +80 right-handed about CoE axis
def lm(l,b):m=vec(l,b)@R;return np.degrees(np.arctan2(m[...,1],m[...,0]))%360   # model-frame l
def moll(l,b):
  lam=np.radians(((np.asarray(l,float)+180)%360)-180);b=np.radians(np.asarray(b,float));th=b.copy()
  for _ in range(60):th=th-(2*th+np.sin(2*th)-np.pi*np.sin(b))/(2+2*np.cos(2*th)+1e-12)
  return -2*np.sqrt(2)/np.pi*lam*np.cos(th),np.sqrt(2)*np.sin(th)
# inverse mollweide grid
X,Y=np.meshgrid(np.linspace(-2*np.sqrt(2),2*np.sqrt(2),1400),np.linspace(-np.sqrt(2),np.sqrt(2),700))
ok=(X/(2*np.sqrt(2)))**2+(Y/np.sqrt(2))**2<=1
th=np.arcsin(np.clip(Y/np.sqrt(2),-1,1));B=np.degrees(np.arcsin(np.clip((2*th+np.sin(2*th))/np.pi,-1,1)))
L=np.degrees(-np.pi*X/(2*np.sqrt(2)*np.cos(th)+1e-12))%360
LM=lm(L,B);rows=[]
fig,axs=plt.subplots(2,1,figsize=(15,15),facecolor='white')
for ax,K in zip(axs,'AB'):
  v=np.interp(LM,MM['l'],MM[K],period=360);gm=np.abs(v[ok]).max()
  img=np.ones(X.shape+(4,));img[...,3]=0
  a=.6*np.abs(v)/gm;red=v>0
  img[...,0]=np.where(red,1,0);img[...,1]=0;img[...,2]=np.where(red,0,1);img[...,3]=np.where(ok,a,0)
  ax.imshow(img,extent=(X.min(),X.max(),Y.min(),Y.max()),origin='lower',interpolation='bilinear')
  for z in Z[K]:
    d=(LM-z+180)%360-180;ax.contour(X,Y,np.where(ok,d,np.nan),[0],colors='k',linewidths=2)
  t=np.linspace(0,2*np.pi,400);ax.plot(2*np.sqrt(2)*np.cos(t),np.sqrt(2)*np.sin(t),'k',lw=1.5)
  for lg in (0,90,270):ax.text(moll(lg,0)[0],-1.55,f'{lg}°',ha='center',fontsize=12)
  for nm_,(cl,cb),mk in [('CoE',(180,0),'X')]:
    c=R1@vec(cl,cb);c2=np.degrees([np.arctan2(c[1],c[0])%(2*np.pi),np.arcsin(c[2])]);xs,ys=moll(*c2);ax.plot(xs,ys,'kX',ms=16)
    ax.annotate('CoE (tilted +15°)',(float(xs),float(ys)),xytext=(10,8),textcoords='offset points',fontsize=12)
  h=[];sc=0
  for i,(nm,sl,sb,c) in enumerate(SPOTS[K]):
    val=float(np.interp(lm(sl,sb),MM['l'],MM[K],period=360));okk=(val>0)==(c=='hot');sc+=okk
    rows.append((K,nm,sl,sb,round(float(lm(sl,sb)),1),round(val,2),okk))
    xs,ys=moll(sl,sb);xs,ys=float(xs),float(ys)
    ax.plot(xs,ys,'o',ms=28,mfc=CC[c],mec='k',mew=2,zorder=7);ax.text(xs,ys,str(i+1),ha='center',va='center',color='white',fontsize=15,fontweight='bold',zorder=8)
    ax.annotate(f"{nm} {val:+.1f}% {'✓' if okk else '✗'}",(xs,ys),xytext=(-20,14) if xs>0 else (20,14),textcoords='offset points',ha='right' if xs>0 else 'left',fontsize=12.5,fontweight='bold',color=CC[c],zorder=9,bbox=dict(fc='white',ec=CC[c],boxstyle='round,pad=.25'))
    h.append(Line2D([],[],ls='',marker='o',mfc=CC[c],mec='k',ms=13,label=f'{i+1} {nm} ({sl}°, {sb:+g}°)'))
  h+=[plt.Rectangle((0,0),1,1,fc='red',alpha=.45,label=f'{K} residual high'),plt.Rectangle((0,0),1,1,fc='blue',alpha=.45,label=f'{K} residual low'),Line2D([],[],color='k',lw=2,label='zero crossings (rotated)')]
  ax.legend(handles=h,loc='upper center',bbox_to_anchor=(.5,-.02),ncol=3,fontsize=11.5,frameon=False)
  ax.set_title(f'{K} residual: CoE tilted up 15°, then rotated 80° about the CoE axis   —   {sc}/{len(SPOTS[K])} match',fontsize=15)
  ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-2.9,2.9);ax.set_ylim(-1.65,1.5)
fig.suptitle('Galactic Mollweide, ℓ = 0 at centre, ℓ increases to the left; spots at published (ℓ, b)',fontsize=13,y=.995)
fig.tight_layout();fig.savefig('overlay_rot_AB.png',dpi=110,facecolor='white')
for r in rows:print(r)
