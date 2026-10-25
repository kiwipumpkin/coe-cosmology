exec(open('recon_n3.py').read().split("cols=['purple'")[0])
run=lambda r:np.array([np.sqrt(np.mean(r[np.abs((l-c+180)%360-180)<=15]**2)) for c in l])
cols=['k','tab:blue','tab:red','tab:green']
plt.rcParams.update({'font.size':16});fig,axs=plt.subplots(2,1,figsize=(18,11),facecolor='white',sharex=True)
for a,(nm,Y,mb,k) in zip(axs,(('A',A0,AFl,0),('B',B0,BLl,1))):
  a.axvspan(15,80,color='gold',alpha=.18);a.axvspan(260,335,color='lightsteelblue',alpha=.35);a.axhline(0,color='grey',lw=1)
  m=mb.copy()
  for i in range(4):
    if i>0:n,p=steps2[i-1];m=m+T(n,p,th,k)
    r=Y-m;lab=('CoE model' if i==0 else f'+ n=1–{i}')+f'  (RMS {np.sqrt(np.mean(r**2)):.3f})'
    a.plot(l,r,'o-' if i in (0,3) else '-',color=cols[i],lw=2.8 if i in(0,3) else 1.6,ms=4,alpha=1 if i in(0,3) else .7,label=lab)
    if i==3:a.plot(l,run(r),color=cols[i],ls=':',lw=1.5,label='running RMS ±15° (n=1–3)')
  a.set_ylabel(f'{nm} residual (data − model)');a.set_ylim(-0.15,0.15);a.legend(loc='upper right',ncol=3,fontsize=12,frameon=False)
axs[1].set_xlim(0,360);axs[1].set_xticks(range(0,361,30));axs[1].set_xlabel('Galactic longitude ℓ [deg]')
fig.suptitle('Residuals vs angle, reconstruction (313 clusters): CoE model → + shared n=1, 2, 3')
fig.tight_layout();fig.savefig('recon_n3_residuals.png',dpi=110,facecolor='white')
