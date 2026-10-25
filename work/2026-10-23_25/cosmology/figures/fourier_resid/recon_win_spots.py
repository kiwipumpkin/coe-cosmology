import json;MM=json.load(open('/tmp/M25.json'))
import json;RR=json.load(open('/tmp/R.json'))
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
plt.rcParams.update({'font.size':16});fig,axs=plt.subplots(2,1,figsize=(18,11),facecolor='white',sharex=True)
def rm(r,w):return np.array([np.mean(r[np.abs((l-c+180)%360-180)<=w]) for c in l])
W=[15,25,35];C=plt.cm.viridis(np.linspace(0,.7,3))
for a,(nm,Y,mb,k) in zip(axs,(('A',A0,AFl,0),('B',B0,BLl,1))):
  m=mb.copy()
  for n,p in steps2[:3]:m=m+T(n,p,th,k)
  r=100*(Y-m)/Y;a.axhline(0,color='grey',lw=1)
  G=np.interp(np.arange(0,360.5,1),np.r_[MM['l'][-1]-360,MM['l'],MM['l'][0]+360],np.r_[MM[nm][-1],MM[nm],MM[nm][0]]);gm=np.abs(G).max()
  for q in range(360):a.axvspan(q,q+1,color='red' if G[q]>0 else 'blue',alpha=.35*abs(G[q])/gm,lw=0,zorder=0)
  a.plot([],[],color='red',alpha=.25,lw=10,label='high (shade ∝ ±25° mean residual)');a.plot([],[],color='blue',alpha=.25,lw=10,label='low (shade ∝ ±25° mean residual)')
  for w,c in zip(W,C):
    M=rm(r,w);i,j=np.argmax(M),np.argmin(M);print(nm,w,'max',l[i],round(M[i],2),'min',l[j],round(M[j],2))
    a.plot(l,M,'-',color=c,lw=2.6,label=f'±{w}°  (max {M[i]:+.1f}% at {l[i]:.0f}°, min {M[j]:+.1f}% at {l[j]:.0f}°)')
    Z=[]
    for q in range(len(l)):
      q2=(q+1)%len(l);x1=l[q2] if q2 else l[0]+360
      if M[q]*M[q2]<0:Z.append((l[q]+(x1-l[q])*M[q]/(M[q]-M[q2]))%360)
    print(nm,w,'zeros',[round(z) for z in sorted(Z)]);[a.axvline(z,color=c,ls='--',lw=1.8,alpha=.9) for z in Z]
    a.plot(l[i],M[i],'^',color=c,ms=11);a.plot(l[j],M[j],'v',color=c,ms=11)
  a.set_ylabel(f'{nm} running mean residual\n% of {nm} (after n=1–3)')
  a.legend(loc='upper center',bbox_to_anchor=(.5,-.04 if k==0 else -.14),ncol=2,fontsize=12,frameon=False)
import csv
P=[('CMB dipole',264,'r'),('CMB dipole antipode',84,'b'),('CMB Cold Spot',209,'b'),('Axis of evil',260,'k'),('Hemispheric asym.',230,'r'),('CatWISE dipole',238,'r'),('NVSS dipole',196,'r'),('Migkas faintest',277,'b'),('Migkas brightest',32,'r'),('Migkas 842 aniso',303,'b'),('Migkas21 H0',280,'b'),('Bulk flow W09',287,'r'),('Loop I',329,'k')]
for a in axs:
  y0,y1=a.get_ylim()
  for i,(t,x,c) in enumerate(P):
    cc={'r':'darkred','b':'navy','k':'0.2'}[c];a.axvline(x,color=cc,lw=1.6,ls='-.',alpha=.8,zorder=2)
SP=sorted(P,key=lambda p:p[1]);key=[]
for i,(t,x,c) in enumerate(SP):
  cc={'r':'darkred','b':'navy','k':'0.2'}[c]
  for a in axs:
    yt=a.get_ylim()[1];a.annotate(str(i+1),(x,yt),xytext=(0,4+12*(i%2)),textcoords='offset points',ha='center',va='bottom',fontsize=10,fontweight='bold',color='white',bbox=dict(boxstyle='circle,pad=.2',fc=cc,ec='none'),annotation_clip=False)
  key.append(f'{i+1} {t} {x}°')
fig.text(.5,.005,'   '.join(key),ha='center',fontsize=10.5,wrap=True)
axs[1].set_xlim(0,360);axs[1].set_xticks(range(0,361,30));axs[1].set_xlabel('Galactic longitude ℓ [deg]')
fig.suptitle('Signed residual after n=1–3, with notable hot (dark red) and cold (navy) spot longitudes',y=1.0)
fig.tight_layout(rect=(0,.04,1,.95));fig.subplots_adjust(hspace=.32);fig.savefig('recon_win_spots.png',dpi=110,facecolor='white')
