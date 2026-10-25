import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('cum2.py').read().split("rA=Araw")[0])
d=json.load(open('../l50/curves.json'));g0=np.array(d['grid'])
P=lambda k,x:periodic_interp(g0,np.array(d[k]),x)
A0,B0,AFl,BLl=P('A0',l),P('B0',l),P('AF',l),P('BL',l);AFg,BLg=P('AF',g),P('BL',g)
rA=A0-AFl;rB=B0-BLl;steps2=[]
for n in range(1,4):
  h=np.column_stack([np.cos(n*th),np.sin(n*th)]);X=np.vstack([np.column_stack([o,Z,h]),np.column_stack([Z,o,h])]);s=np.r_[sA,sB]
  p=np.linalg.lstsq(X/s[:,None],np.r_[rA,rB]/s,rcond=None)[0];steps2.append((n,p));rA=rA-X[:71]@p;rB=rB-X[71:]@p
T=lambda n,p,x,k:p[k]+p[2]*np.cos(n*x)+p[3]*np.sin(n*x)
cols=['purple','tab:blue','tab:red','tab:green']
plt.rcParams.update({'font.size':16});fig,axs=plt.subplots(2,4,figsize=(26,11),facecolor='white',sharex=True,sharey='row')
for row,(nm,Y,sY,mb,mg,k) in enumerate((('A',A0,sA,AFl,AFg,0),('B',B0,sB,BLl,BLg,1))):
  m,mgg=mb.copy(),mg.copy()
  for i in range(4):
    if i>0:n,p=steps2[i-1];pm=mgg.copy();m=m+T(n,p,th,k);mgg=mgg+T(n,p,gt,k)
    a=axs[row,i];a.axvspan(15,80,color='gold',alpha=.18);a.axvspan(260,335,color='lightsteelblue',alpha=.35)
    a.errorbar(l,Y,yerr=sY,fmt='o',color='k',ms=4,alpha=.5,lw=1)
    if i>0:a.plot(g,pm,color='grey',ls='--',lw=2,alpha=.6)
    a.plot(g,mgg,color=cols[i],lw=3.5)
    r=Y-m;a.text(.03,.04,f'χ²={np.sum(r**2/sY**2):.2f}  r²={np.corrcoef(Y,m)[0,1]**2:.3f}  RMS={np.sqrt(np.mean(r**2)):.3f}',transform=a.transAxes,fontsize=14,color='0.2')
    if row==0:a.set_title('CoE model' if i==0 else f'+ shared n={i}',fontweight='bold')
    a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])
  axs[row,0].set_ylabel('A/A_all' if nm=='A' else 'B/B_all')
for a in axs[1]:a.set_xlabel('ℓ [deg]')
from matplotlib.lines import Line2D
h=[Line2D([],[],color='k',marker='o',ls='',alpha=.5,label='reconstruction, 313 clusters (±σ Fig. 6)'),Line2D([],[],color='purple',lw=3,label='CoE model (A: AF_smooth, B: sheet BL)')]+[Line2D([],[],color=cols[j],lw=3,label=f'+ shared n={j}') for j in (1,2,3)]+[Line2D([],[],color='grey',ls='--',lw=2,label='previous'),plt.Rectangle((0,0),1,1,color='gold',alpha=.3,label='ℓ 15–80°'),plt.Rectangle((0,0),1,1,color='lightsteelblue',alpha=.5,label='ℓ 260–335°')]
fig.legend(handles=h,loc='lower center',ncol=4,frameon=False,fontsize=15)
fig.tight_layout(rect=(0,.09,1,1));fig.savefig('recon_n3_build.png',dpi=95,facecolor='white')
for n,p in steps2:
  ph=np.degrees(np.arctan2(p[3],p[2]))%360;print(n,round(np.hypot(p[2],p[3]),4),sorted(round((ph/n+360*j/n)%360,1) for j in range(n)))
