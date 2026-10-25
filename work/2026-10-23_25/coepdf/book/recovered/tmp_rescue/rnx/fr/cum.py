import numpy as np,csv
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':18})
from matplotlib.lines import Line2D
from make_A_vs_angle_axisymmetric_before_rotation import periodic_interp
rows=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
l=np.array([float(r['l_deg']) for r in rows]);th=np.radians(l)
d=np.genfromtxt('residuals.csv',delimiter=',',names=True);Araw,AF,sA=d['A_migkas'],d['A_axisym_AF_smooth'],d['sigma']
Bv=np.array([float(r['B_over_Ball']) for r in rows]);sB=np.array([(float(r['B_err_low'])+float(r['B_err_high']))/2 for r in rows])
gm=np.genfromtxt('sheet_BL_copy2_density_w_rotat.csv',delimiter=',',names=True)
g=np.linspace(0,360,721);gt=np.radians(g)
BG=periodic_interp(gm['ell_deg'],gm['BL_modelB_sheet'],l);BGg=periodic_interp(gm['ell_deg'],gm['BL_modelB_sheet'],g)
AFg=periodic_interp(l,AF,g)
SG=1.;chi=lambda r,s:np.sum(r**2/s**2);o=np.ones_like(th);Z=0*o
rA=Araw-AF;rB=Bv-BG;steps=[]
for n in range(1,5):
 h=np.column_stack([np.cos(n*th),np.sin(n*th)])
 X=np.vstack([np.column_stack([o,Z,h]),np.column_stack([Z,o,SG*h])]);s=np.r_[sA,sB]
 p=np.linalg.lstsq(X/s[:,None],np.r_[rA,rB]/s,rcond=None)[0]
 steps.append(dict(n=n,p=p,C=np.hypot(p[2],p[3]),pk=np.degrees(np.arctan2(p[3],p[2]))%360/n))
 rA=rA-X[:71]@p;rB=rB-X[71:]@p
def term(st,x,k,sg):n=st['n'];p=st['p'];return p[k]+sg*(p[2]*np.cos(n*x)+p[3]*np.sin(n*x))
cols=['purple','tab:blue','tab:red','tab:green','tab:orange']
for nm,Y,sY,base,baseg,k,sg in (('A',Araw,sA,AF,AFg,0,1.),('B',Bv,sB,BG,BGg,1,SG)):
 lab='A/A_all' if nm=='A' else 'B/B_all'
 cum=[(base.copy(),baseg.copy())]
 for st in steps:cum.append((cum[-1][0]+term(st,th,k,sg),cum[-1][1]+term(st,gt,k,sg)))
 fig,axs=plt.subplots(5,2,figsize=(22,30),facecolor='white');print(nm)
 for i in range(5):
  m,mg=cum[i];res=Y-m;c2=chi(res,sY);r2=np.corrcoef(Y,m)[0,1]**2;print(f' row{i+1} chi2={c2:.2f} r2={r2:.3f}')
  a=axs[i,0];a.errorbar(l,Y,yerr=sY,fmt='o',color='k',ms=5,alpha=.6,capsize=2)
  if i>0:a.plot(g,cum[i-1][1],color=cols[i-1],lw=2.2,ls='--',alpha=.55)
  a.plot(g,mg,color=cols[i],lw=3.2)
  hd='base model' if i==0 else 'base'+''.join(f' + n{j}' for j in range(1,i+1))
  a.set_title(f'{lab}: {hd}\nχ²_{nm}={c2:.2f},  Pearson r²={r2:.3f}',fontsize=18)
  a=axs[i,1];a.errorbar(l,res,yerr=sY,fmt='o',color='k',ms=5,alpha=.6,capsize=2);a.axhline(0,color='grey',lw=1.5)
  if i<4:
   st=steps[i];a.plot(g,term(st,gt,k,sg),color=cols[i+1],lw=3.2)
   pk=st['pk'] if sg>0 else (st['pk']+180/st['n'])%(360/st['n'])
   a.set_title(f'residual after row {i+1}; next: n={st["n"]}  C={st["C"]:.3f}, peak {pk:.0f}°\nχ²_{nm} {c2:.2f} → {chi(Y-cum[i+1][0],sY):.2f} after removal',fontsize=18)
  else:a.set_title(f'final residual (base + n1..n4)\nχ²_{nm}={c2:.2f}',fontsize=18)
  a.set_ylabel('residual');axs[i,0].set_ylabel(lab)
  for b_ in axs[i]:b_.set_xlim(0,360);b_.set_facecolor('white')
 for b_ in axs[:,0]:b_.set_ylim(0.78,1.22)
 for b_ in axs[:,1]:b_.set_ylim(-0.25,0.25)
 for b_ in axs[-1]:b_.set_xlabel('Galactic longitude ℓ [deg]')
 bn='AF_smooth (axisymmetric)' if nm=='A' else 'sheet BL: (1−BE)·AN'
 h=[Line2D([],[],color='k',marker='o',ls='',label='data / residual ±σ (Migkas Fig. 6)'),Line2D([],[],color='purple',lw=3,label=f'base: {bn}')]+[Line2D([],[],color=cols[j],lw=3,label=f'+ mutual n={j} (B=+shared)') for j in range(1,5)]+[Line2D([],[],color='grey',lw=2.2,ls='--',label='previous row model (left) ')]
 fig.legend(handles=h,loc='lower center',ncol=3,fontsize=18,frameon=False)
 fig.suptitle(f'{lab}: cumulative sequential mutual harmonics (same fit as joint_sequential_mutual_harmonics.png)\nadjacent 90° slices overlap → bins correlated, χ² overstates significance',fontsize=20)
 fig.tight_layout(rect=(0,0.05,1,0.965));fig.savefig(f'{nm}_cumulative_harmonics.png',dpi=75,facecolor='white')
