import numpy as np,csv
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':17})
from matplotlib.lines import Line2D
rows=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
l=np.array([float(r['l_deg']) for r in rows]);th=np.radians(l)
d=np.genfromtxt('residuals.csv',delimiter=',',names=True);Araw,AF,sA=d['A_migkas'],d['A_axisym_AF_smooth'],d['sigma']
Bv=np.array([float(r['B_over_Ball']) for r in rows]);sB=np.array([(float(r['B_err_low'])+float(r['B_err_high']))/2 for r in rows])
wB=1/sB**2;Bm=np.sum(wB*Bv)/np.sum(wB)
rA=Araw-AF;rB=Bv-Bm;chi=lambda r,s:np.sum(r**2/s**2)
out=[];o=np.ones_like(th);Z=0*o
for n in range(1,5):
 h=np.column_stack([np.cos(n*th),np.sin(n*th)])
 X=np.vstack([np.column_stack([o,Z,h]),np.column_stack([Z,o,-h])]);s=np.r_[sA,sB]
 p=np.linalg.lstsq(X/s[:,None],np.r_[rA,rB]/s,rcond=None)[0]
 C=np.hypot(p[2],p[3]);pk=np.degrees(np.arctan2(p[3],p[2]))%360/n
 mA=X[:71]@p;mB=X[71:]@p
 out.append(dict(n=n,rA=rA.copy(),rB=rB.copy(),p=p,C=C,pk=pk,cA0=chi(rA,sA),cB0=chi(rB,sB),cA1=chi(rA-mA,sA),cB1=chi(rB-mB,sB)))
 rA=rA-mA;rB=rB-mB
print('row n C peakA(deg) chiA_before chiA_after chiB_before chiB_after comb_after')
for i,r in enumerate(out):print(i+2,r['n'],f"{r['C']:.4f} {r['pk']:.1f} {r['cA0']:.2f} {r['cA1']:.2f} {r['cB0']:.2f} {r['cB1']:.2f} {r['cA1']+r['cB1']:.2f}")
print('row1 chiA',chi(Araw-AF,sA),'chiB',chi(Bv-Bm,sB),'Bmean',Bm)
g=np.linspace(0,360,721);gt=np.radians(g);cols=['tab:blue','tab:red','tab:green','tab:orange']
gi=np.argsort(l)
fig,axs=plt.subplots(5,2,figsize=(20,27),sharex=True,sharey=False)
a=axs[0,0];a.errorbar(l,Araw,yerr=sA,fmt='o',color='k',ms=4,alpha=.6,capsize=2);a.plot(l[gi],AF[gi],color='purple',lw=3)
a.set_title(f'A/A_all data + axisymmetric fit (AF_smooth)\nχ²_A={chi(Araw-AF,sA):.2f}',fontsize=16)
a=axs[0,1];a.errorbar(l,Bv,yerr=sB,fmt='o',color='k',ms=4,alpha=.6,capsize=2);a.axhline(Bm,color='purple',lw=3,ls=':')
a.set_title(f'B/B_all data + weighted mean {Bm:.4f} (no B model exists)\nχ²_B={chi(Bv-Bm,sB):.2f}',fontsize=16)
for i,r in enumerate(out):
 n=r['n'];p=r['p'];hc=p[2]*np.cos(n*gt)+p[3]*np.sin(n*gt)
 for j,(res,s,c0,sg,b,af,nm) in enumerate(((r['rA'],sA,p[0],1,r['cA0'],r['cA1'],'A'),(r['rB'],sB,p[1],-1,r['cB0'],r['cB1'],'B'))):
  ax=axs[i+1,j];ax.errorbar(l,res,yerr=s,fmt='o',color='k',ms=4,alpha=.6,capsize=2);ax.axhline(0,color='grey',lw=1)
  ax.plot(g,c0+sg*hc,color=cols[i],lw=3)
  pk=r['pk'] if sg==1 else (r['pk']+180/n)%(360/n)
  ax.set_title(f"{nm} resid. + mutual n={n}{'' if sg==1 else ' (×−1)'}: C={r['C']:.3f}, peak {pk:.0f}°\nχ²_{nm} {b:.2f} → {af:.2f} after removal",fontsize=16)
for ax in axs[:,0]:ax.set_ylabel('A/A_all' if ax is axs[0,0] else 'residual')
axs[0,1].set_ylabel('B/B_all')
for ax in axs[-1]:ax.set_xlabel('Galactic longitude ℓ [deg]')
axs[0,0].set_xlim(0,360)
for ax in axs[1:].ravel():ax.set_ylim(-0.24,0.24)
for ax in axs[0]:ax.set_ylim(0.78,1.22)
h=[Line2D([],[],color='k',marker='o',ls='',label='data / residual ±σ'),Line2D([],[],color='purple',lw=3,label='A: axisymmetric AF_smooth'),Line2D([],[],color='purple',lw=3,ls=':',label='B: weighted mean (no B model exists)')]+[Line2D([],[],color=c,lw=3,label=f'mutual n={k+1} (B = −shared)') for k,c in enumerate(cols)]
fig.legend(handles=h,loc='lower center',ncol=3,fontsize=18,frameon=False)
fig.suptitle('Sequential mutual harmonics: each row fits ONE shared harmonic jointly to current A & B residuals\n(min χ²_A+χ²_B, 1/σ² weights, separate offsets, B = −shared); next row subtracts it',fontsize=19)
fig.tight_layout(rect=(0,0.06,1,0.96));fig.savefig('joint_sequential_mutual_harmonics.png',dpi=90)
