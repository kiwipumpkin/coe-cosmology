exec(open('recon_n3.py').read().split("cols=['purple'")[0])
plt.rcParams.update({'font.size':16});fig,axs=plt.subplots(2,1,figsize=(18,11),facecolor='white',sharex=True)
def run(r):return np.array([np.sqrt(np.mean(r[np.abs((l-c+180)%360-180)<=15]**2)) for c in l])
def segs(mask):
  out=[];st=None
  for j in range(len(l)):
    if mask[j] and st is None:st=j
    if (not mask[j] or j==len(l)-1) and st is not None:out.append((l[st]-2.5,l[j if mask[j] else j-1]+2.5));st=None
  return out
for a,(nm,Y,mb,k) in zip(axs,(('A',A0,AFl,0),('B',B0,BLl,1))):
  S=Y;m=mb.copy();r0=Y-m
  for n,p in steps2[:3]:m=m+T(n,p,th,k)
  r3=Y-m;R0=run(100*r0/Y);R3=run(100*r3/Y);o0=np.sqrt(np.mean((100*r0/Y)**2));o3=np.sqrt(np.mean((100*r3/Y)**2))
  for x0,x1 in segs(R3>1.25*o3):a.axvspan(x0,x1,color='orangered',alpha=.18,zorder=0);print(nm,'HIGH',x0,x1)
  for x0,x1 in segs(R3<0.75*o3):a.axvspan(x0,x1,color='seagreen',alpha=.18,zorder=0);print(nm,'LOW',x0,x1)
  a.plot(l,R0,'o-',color='k',lw=2.5,ms=4,label=f'before: CoE model only (overall {o0:.1f}%)')
  a.plot(l,R3,'o-',color='tab:green',lw=2.5,ms=4,label=f'after: + n=1–3 harmonics (overall {o3:.1f}%)')
  a.axhline(o3,color='tab:green',ls=':',lw=1.5)
  a.plot([],[],color='orangered',alpha=.3,lw=10,label='high error (after > 1.25× overall)')
  a.plot([],[],color='seagreen',alpha=.3,lw=10,label='low error (after < 0.75× overall)')
  a.set_ylabel(f'{nm} running RMS error (±15°)\n% of {nm} itself');a.set_ylim(0,None)
  a.legend(loc='upper center',bbox_to_anchor=(.5,-.04 if k==0 else -.14),ncol=4,fontsize=12,frameon=False)
  print(nm,round(o0,1),round(o3,1))
axs[1].set_xlim(0,360);axs[1].set_xticks(range(0,361,30));axs[1].set_xlabel('Galactic longitude ℓ [deg]')
fig.suptitle('Error before and after the shared harmonics (reconstruction, 313 clusters)')
fig.tight_layout();fig.subplots_adjust(hspace=.28);fig.savefig('recon_err_pctdata.png',dpi=110,facecolor='white')
