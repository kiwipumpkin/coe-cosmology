import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
N='/workspace/cosmology-course/notes/'
D,d,rhoD=58.2,14.3,0.3
lab={'soft':r'$\rho_0/(1+(r/L)^{0.3})$, R²=0.743, 5.1σ',
     'exp':r'$\rho_0\,e^{-r/L}$, R²=0.72, 5.3σ',
     'pl':r'$\rho_0(1+r/L)^{-1.2}$, R²=0.72, 5.0σ',
     'eta':r'$\rho_0\,e^{-(r/L)^{0.05}}$, R²=0.745, 5.0σ, ℓ₀=19.9°'}
df=pd.read_csv(N+'catwise-density-levelset.csv'); e=pd.read_csv(N+'catwise-eta-density.csv')
ys=df['rho_over_rhocrit_exp_L922.6'].values; r=df.r_Gpc.values
slope=np.gradient(ys,r)[np.argmin(abs(r-D))]
fig,ax=plt.subplots(figsize=(12,7.8))
ax.plot(r,df['rho_over_rhocrit_soft_L5110_alpha0.3_beta1'],'-',lw=3,label=lab['soft'])
ax.plot(r,ys,'--',lw=3,label=lab['exp'])
ax.plot(r,df['rho_over_rhocrit_pl_L1049_n1.2'],':',lw=3,label=lab['pl'])
ax.plot(e.r_Gpc,e['rho_over_rhocrit_eta0.0500009_d14.3005_L0.849'],'-',color='red',lw=2.5,alpha=.8,label=lab['eta'])
rt=np.linspace(D-d-8,D+d+8,50); ax.plot(rt,rhoD+slope*(rt-D),color='k',lw=2,label=f'Tangent at D, slope = {slope:.2e} ρ_crit/Gpc')
for x in (D-d,D+d): ax.axvline(x,color='purple',ls='--',lw=1.5)
ax.axvline(D,color='k',ls='-',lw=1.2)
ax.plot([D-d,D,D+d],[rhoD-slope*d,rhoD,rhoD+slope*d],'o',color='purple',ms=9)
ax.text(D,0.47,'Earth, D = 58.2',ha='center',fontsize=14,bbox=dict(fc='w',ec='none'))
ax.text(D-d,0.44,'D − d = 43.9',ha='center',fontsize=14,color='purple',bbox=dict(fc='w',ec='none'))
ax.text(D+d,0.44,'D + d = 72.5',ha='center',fontsize=14,color='purple',bbox=dict(fc='w',ec='none'))
ax.text(D-d-1,0.375,'[ρ(D−d) − ρ(D+d)] / ρ(D) = 3.1%',ha='right',fontsize=14,bbox=dict(fc='w',ec='gray'))
ax.set_xlim(0,80); ax.set_ylim(0.28,0.5)
ax.set_xlabel('R from CoE (Gpc)',fontsize=18); ax.set_ylabel('ρ / ρ_crit',fontsize=18)
ax.tick_params(labelsize=15); ax.grid(alpha=.3)
ax.legend(fontsize=13,loc='upper center',bbox_to_anchor=(0.5,-0.12),ncol=2)
ax.set_title('Density vs R, four laws under the 3% CatWISE lock',fontsize=18)
fig.tight_layout(); fig.savefig(N+'catwise-density-four-laws.png',dpi=130)

src=open(N+'catwise-density-levelset-migkas.py').read()
ns={'__name__':'x'}; exec(src[:src.index('Y = MA')],ns)
MA,MLO,MHI=ns['MA'],ns['MLO'],ns['MHI']
Ld=np.loadtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',skiprows=1,usecols=0)
fig,ax=plt.subplots(figsize=(12,7.5))
ax.errorbar(Ld,MA,yerr=[MLO,MHI],fmt='o',color='k',ms=4,alpha=.6,label='Migkas bins')
for f,ls in [('soft','-'),('exp','--'),('pl',':')]:
    q=pd.read_csv(N+f'catwise-density-levelset-migkas-{f}_fixed15.2.csv'); ax.plot(q.ell_deg,q.A_over_Aall,ls,lw=3,label=lab[f])
q=pd.read_csv(N+'catwise-eta-migkas-fit-d14.3.csv'); ax.plot(q.ell_deg,q.A_over_Aall,'-',color='red',lw=2.5,alpha=.8,label=lab['eta'])
ax.set_xlabel('Galactic longitude ℓ (deg)',fontsize=18); ax.set_ylabel('A / A_all',fontsize=18)
ax.set_xlim(0,360); ax.tick_params(labelsize=15); ax.grid(alpha=.3)
ax.legend(fontsize=14,loc='upper center',bbox_to_anchor=(0.5,-0.13),ncol=2)
ax.set_title('Migkas A/A_all, four laws (ℓ₀=15.2° unless marked)',fontsize=18)
fig.tight_layout(); fig.savefig(N+'catwise-AoverAall-four-laws.png',dpi=130)

fig=plt.figure(figsize=(12,8)); fig.text(0.04,0.93,'Density functions, all locked to the 3% CatWISE gradient',fontsize=22,weight='bold')
rows=[('Soft core + tail',r'$\rho(r)=\dfrac{\rho_0}{1+(r/L)^{\alpha}},\qquad \alpha=0.3$'),
('Exponential',r'$\rho(r)=\rho_0\,e^{-r/L}$'),
('Power law',r'$\rho(r)=\rho_0\left(1+\dfrac{r}{L}\right)^{-n},\qquad n=1.2$'),
('Stretched exponential',r'$\rho(r)=\rho_0\,e^{-(r/L)^{\eta}},\qquad \eta=0.05$'),
('Lock shared by all four',r'$\dfrac{\rho(D-d)-\rho(D+d)}{\rho(D)}=0.031,\qquad \rho(D)=0.3\,\rho_{\rm crit}$')]
y=0.80
for n_,f in rows:
    fig.text(0.05,y,n_,fontsize=19,weight='bold'); fig.text(0.07,y-0.075,f,fontsize=22); y-=0.165
fig.savefig(N+'catwise-density-functions.png',dpi=130)
print(slope)
