exec(open('recon_n3.py').read().split("cols=['purple'")[0])
cols=['k','tab:blue','tab:red','tab:green']
plt.rcParams.update({'font.size':16});fig,axs=plt.subplots(2,1,figsize=(18,11),facecolor='white',sharex=True)
for a,(nm,Y,mb,k) in zip(axs,(('A',A0,AFl,0),('B',B0,BLl,1))):
  S=np.sqrt(np.mean((Y-1)**2))
  a.axhline(0,color='grey',lw=1)
  a.plot(l,100*(Y-1)/S,color='0.6',lw=6,alpha=.4,label='original signal (data − 1)')
  m=mb.copy()
  for i in range(4):
    if i>0:n,p=steps2[i-1];m=m+T(n,p,th,k)
    r=Y-m;pc=100*np.sqrt(np.mean(r**2))/S;print(nm,i,round(pc,1))
    a.plot(l,100*r/S,'o-' if i in(0,3) else '-',color=cols[i],lw=2.8 if i in(0,3) else 1.6,ms=4,alpha=1 if i in(0,3) else .7,label=('CoE model' if i==0 else f'+ n=1–{i}')+f'  (RMS = {pc:.0f}% of signal)')
  run=np.array([np.sqrt(np.mean(r[np.abs((l-c+180)%360-180)<=15]**2)) for c in l]);thr=1.25*np.sqrt(np.mean(r**2));hi=run>thr
  segs=[];st=None
  for j in range(len(l)):
    if hi[j] and st is None:st=j
    if (not hi[j] or j==len(l)-1) and st is not None:segs.append((l[st]-2.5,l[j if hi[j] else j-1]+2.5));st=None
  for x0,x1 in segs:a.axvspan(x0,x1,color='orange',alpha=.2,zorder=0);print(nm,'high-error',x0,x1)
  a.plot([],[],color='orange',alpha=.3,lw=10,label='highest error: running RMS (±15°) > 1.25× overall, after n=1–3')
  a.set_ylabel(f'{nm} residual, % of signal RMS\n(signal RMS = {S:.3f})');a.set_ylim(-250,250);a.legend(loc='upper right',ncol=3,fontsize=12,frameon=False)
axs[1].set_xlim(0,360);axs[1].set_xticks(range(0,361,30));axs[1].set_xlabel('Galactic longitude ℓ [deg]')
fig.suptitle('Residual as % of the original signal RMS (reconstruction, 313 clusters)')
fig.tight_layout();fig.savefig('recon_pct.png',dpi=110,facecolor='white')
