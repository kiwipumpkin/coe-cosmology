src=open('/workspace/cosmology/figures/l50/boot_spots.py').read().split("M=json.load")[0]
exec(src)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
exec(open('/tmp/spots.py').read())
Z={'A':[22,99,262,310],'B':[14,104,266,321]}
fig,axs=plt.subplots(2,1,figsize=(18,11),facecolor='white')
cols={15:'#3b0f4f',25:'#2c6e8f',35:'#3cbc63'}
for ax,K in zip(axs,'AB'):
  r25=base[K][25];gm=np.abs(r25).max();ll=np.arange(0,360.5,.5)
  v=np.interp(ll,l,r25,period=360);f=np.clip((v/gm+1)/2,0,1)
  img=np.stack([.3*f+.7,np.full_like(f,.7),.3*(1-f)+.7],-1)[None]
  ax.imshow(img,extent=(0,360,-100,100),aspect='auto',zorder=0)
  h=[]
  for w in (15,25,35):
    y=base[K][w];o=np.argsort(l);ax.plot(np.r_[l[o],l[o][0]+360],np.r_[y[o],y[o][0]],color=cols[w],lw=2.6 if w==15 else 2,zorder=3)
    imx,imn=np.argmax(y),np.argmin(y)
    h.append(Line2D([],[],color=cols[w],lw=2.5,label=f'±{w}°  (max {y[imx]:+.1f}% at {l[imx]:.0f}°, min {y[imn]:+.1f}% at {l[imn]:.0f}°)'))
  ax.axhline(0,color='grey',lw=1)
  for z in Z[K]:ax.axvline(z,color='k',lw=1.5,ls='--',zorder=2)
  lo,hi=ax.get_ylim() if False else (min(base[K][15].min(),-1)*1.15,max(base[K][15].max(),1)*1.45)
  ax.set_ylim(lo,hi)
  MK=['o','s','D','^','v','p','h']
  for i,(nm,sl,sb,c) in enumerate(SPOTS[K]):
    lm_={'CMB dipole cold':97.2,'Bulk flow away':113.9,'CMB Cold Spot':156.2,'CMB dipole hot':277.2,'Bulk flow toward':293.9,'NVSS cold':58.5,'Hemi. asym. cold':29.9,'CatWISE cold':72.2,'NVSS hot':238.5,'Hemi. asym. hot':209.9,'CatWISE hot':252.2,'Loop I / NPS':335.1}[nm]
    ax.axvline(lm_,color=CC[c],lw=2,ls='-.',zorder=2)
    ax.plot(lm_,hi*.88,MK[i],ms=17,mfc=CC[c],mec='k',mew=1.5,zorder=5);ax.text(lm_,hi*.88,str(i+1),ha='center',va='center',color='white',fontsize=10,fontweight='bold',zorder=6)
    h.append(Line2D([],[],ls='',marker=MK[i],mfc=CC[c],mec='k',ms=11,label=f'{i+1} {nm}'))
  h+=[plt.Rectangle((0,0),1,1,fc=(1,.7,.7),label='residual high'),plt.Rectangle((0,0),1,1,fc=(.7,.7,1),label='residual low'),Line2D([],[],color='k',ls='--',label='zero crossings')]
  ax.legend(handles=h,loc='upper center',bbox_to_anchor=(.5,-.13),ncol=4,fontsize=11,frameon=False)
  ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.set_xlabel('model-frame angle ℓ_m [deg]',fontsize=13)
  ax.set_ylabel(f'{K} running mean residual\n% of {K} (after n=1–3)',fontsize=13)
fig.tight_layout();fig.savefig('/tmp/top_new.png',dpi=110,facecolor='white')
