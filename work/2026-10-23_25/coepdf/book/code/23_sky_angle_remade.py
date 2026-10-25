# REMADE (not original): A/B running-mean residual vs model-frame angle, spots at model-frame angle
# Based on /tmp/rn/recon_win_spotsAB.py (curves) + /tmp/r45c.py (lift 15, roll 45 rotation -> lm()).
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
MM=json.load(open('/tmp/M25.json'));RR=json.load(open('/tmp/R.json'))
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
r45=open('/tmp/r45c.py').read();exec(r45.split('def moll')[0].split("MM=json.load")[1].split('\n',1)[1])
exec(open('/tmp/spots.py').read())
plt.rcParams.update({'font.size':16});fig,axs=plt.subplots(2,1,figsize=(18,11),facecolor='white')
def rm(r,w):return np.array([np.mean(r[np.abs((l-c+180)%360-180)<=w]) for c in l])
W=[15,25,35];C=plt.cm.viridis(np.linspace(0,.7,3));MK=['o','s','D','^','v','p','h']
for a,(nm,Y,mb,k) in zip(axs,(('A',A0,AFl,0),('B',B0,BLl,1))):
  m=mb.copy()
  for n,p in steps2[:3]:m=m+T(n,p,th,k)
  r=100*(Y-m)/Y;a.axhline(0,color='grey',lw=1)
  G=np.interp(np.arange(0,360.5,1),np.r_[MM['l'][-1]-360,MM['l'],MM['l'][0]+360],np.r_[MM[nm][-1],MM[nm],MM[nm][0]]);gm=np.abs(G).max()
  for q in range(360):a.axvspan(q,q+1,color='red' if G[q]>0 else 'blue',alpha=.08+.3*abs(G[q])/gm,lw=0,zorder=0)
  h=[]
  for w,c in zip(W,C):
    M=rm(r,w);i,j=np.argmax(M),np.argmin(M)
    h.append(a.plot(l,M,'-',color=c,lw=2.2,label=f'±{w}°  (max {M[i]:+.1f}% at {l[i]:.0f}°, min {M[j]:+.1f}% at {l[j]:.0f}°)')[0])
  yt=a.get_ylim()[1]
  for i,(t,sl,sb,c) in enumerate(SPOTS[nm]):
    x=float(lm(sl,sb));a.axvline(x,color=CC[c],lw=1.8,ls='-.')
    a.plot(x,yt*.9,MK[i],ms=14,mfc=CC[c],mec='k',clip_on=False,zorder=5);a.text(x,yt*.9,str(i+1),ha='center',va='center',color='white',fontsize=8,fontweight='bold',zorder=6)
    h.append(Line2D([],[],ls='',marker=MK[i],mfc=CC[c],mec='k',ms=10,label=f'{i+1} {t}'))
  for z in Z[nm]:a.axvline(z,color='k',ls='--',lw=1.2)
  h+=[plt.Rectangle((0,0),1,1,fc='red',alpha=.3,label='residual high'),plt.Rectangle((0,0),1,1,fc='blue',alpha=.3,label='residual low'),Line2D([],[],color='k',ls='--',label='zero crossings')]
  a.set_xlim(0,360);a.set_xticks(range(0,361,30));a.set_xlabel('model-frame angle ℓ_m [deg]',fontsize=11)
  a.set_ylabel(f'{nm} running mean residual\n% of {nm} (after n=1–3)',fontsize=11)
  a.legend(handles=h,loc='upper center',bbox_to_anchor=(.5,-.1),ncol=3,fontsize=10,frameon=False)
fig.tight_layout();fig.subplots_adjust(hspace=.55);fig.savefig('sky_angle_remade.png',dpi=110,facecolor='white')
