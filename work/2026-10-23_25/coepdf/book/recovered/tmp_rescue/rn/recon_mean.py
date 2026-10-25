exec(open('recon_n3.py').read().split("cols=['purple'")[0])
plt.rcParams.update({'font.size':16});fig,axs=plt.subplots(2,1,figsize=(18,11),facecolor='white',sharex=True)
def rm(r):return np.array([np.mean(r[np.abs((l-c+180)%360-180)<=15]) for c in l])
def segs(mask):
  out=[];st=None
  for j in range(len(l)):
    if mask[j] and st is None:st=j
    if (not mask[j] or j==len(l)-1) and st is not None:out.append((l[st]-2.5,l[j if mask[j] else j-1]+2.5));st=None
  return out
for a,(nm,Y,mb,k) in zip(axs,(('A',A0,AFl,0),('B',B0,BLl,1))):
  m=mb.copy();r0=100*(Y-m)/Y
  for n,p in steps2[:3]:m=m+T(n,p,th,k)
  r3=100*(Y-m)/Y;M0=rm(r0);M3=rm(r3);s=np.std(M3)
  for x0,x1 in segs(M3>s):a.axvspan(x0,x1,color='orangered',alpha=.18,zorder=0)
  for x0,x1 in segs(M3<-s):a.axvspan(x0,x1,color='royalblue',alpha=.18,zorder=0)
  i,j=np.argmax(M3),np.argmin(M3);print(nm,'max',l[i],round(M3[i],2),'min',l[j],round(M3[j],2),'s',round(s,2),'hi',segs(M3>s),'lo',segs(M3<-s))
  i0,j0=np.argmax(M0),np.argmin(M0);print(nm,'before max',l[i0],round(M0[i0],2),'min',l[j0],round(M0[j0],2))
  a.axhline(0,color='grey',lw=1)
  a.plot(l,M0,'o-',color='k',lw=2.5,ms=4,label='before: CoE model only')
  a.plot(l,M3,'o-',color='tab:green',lw=2.5,ms=4,label='after: + n=1–3 harmonics')
  a.plot([],[],color='orangered',alpha=.3,lw=10,label='data above model (after > +1σ)')
  a.plot([],[],color='royalblue',alpha=.3,lw=10,label='data below model (after < −1σ)')
  a.set_ylabel(f'{nm} running mean residual (±15°)\n% of {nm}')
  a.legend(loc='upper center',bbox_to_anchor=(.5,-.04 if k==0 else -.14),ncol=4,fontsize=12,frameon=False)
axs[1].set_xlim(0,360);axs[1].set_xticks(range(0,361,30));axs[1].set_xlabel('Galactic longitude ℓ [deg]')
fig.suptitle('Signed residual (data − model), averaged over ±15°, before and after harmonics')
fig.tight_layout();fig.subplots_adjust(hspace=.28);fig.savefig('recon_mean.png',dpi=110,facecolor='white')
