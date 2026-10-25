import json,numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':17})
from matplotlib.lines import Line2D
o=json.load(open('zshell_results.json'));sh=['z1','z2','z3']
zm=np.array([o[s]['zmed'] for s in sh]);dd=np.array([o[s]['d_Mpc'] for s in sh])
from scipy.integrate import quad
dc=lambda z:quad(lambda x:1/np.sqrt(.3*(1+x)**3+.7),0,z)[0]*299792.458/70
zz=np.linspace(0.02,0.18,80);dz=np.array([dc(q) for q in zz]);LS=258.5
fig,axs=plt.subplots(2,4,figsize=(26,13))
col={'A':'tab:blue','B':'tab:red'};off={'A':-0.003,'B':0.003}
for i in range(4):
  n=i+1;P=360/n;ax=axs[0,i];ax2=axs[1,i]
  for lab in 'AB':
    ph=np.array([o[s][lab][i]['peak'] for s in sh]);pe=np.array([o[s][lab][i]['peak_err'] for s in sh])
    C=np.array([o[s][lab][i]['C'] for s in sh]);Ce=np.array([o[s][lab][i]['C_err'] for s in sh])
    # unwrap relative to shell 1
    ph=ph[0]+((ph-ph[0]+P/2)%P-P/2)
    ax.errorbar(zm+off[lab],ph,yerr=pe,fmt='o',ms=11,color=col[lab],capsize=5,lw=2.5)
    ax.plot(zz,np.full_like(zz,ph[0]),color=col[lab],ls='-',lw=1.8,alpha=.6)
    # fixed-Mpc: offset from anchor LS shrinks as d1/d
    a=LS+((ph[0]-LS+P/2)%P-P/2);fm=LS+(ph[0]-a+ (a-LS))*0+ (ph[0]-LS - ((ph[0]-LS+P/2)%P-P/2)) 
    delta=((ph[0]-LS+P/2)%P-P/2);ax.plot(zz,ph[0]-delta+delta*dd[0]/dz,color=col[lab],ls='--',lw=2.2)
    ax2.errorbar(zm+off[lab],C,yerr=Ce,fmt='s',ms=11,color=col[lab],capsize=5,lw=2.5)
  ax.set_title(f'n={n}: peak longitude ℓ_n vs z\n(period {P:.0f}°)');ax.set_xlabel('median z of shell');ax.set_xlim(0.02,0.17)
  ax2.set_title(f'n={n}: amplitude C_n vs z');ax2.set_xlabel('median z of shell');ax2.set_xlim(0.02,0.17);ax2.set_ylim(0,0.55)
axs[0,0].set_ylabel('peak ℓ_n [deg] (unwrapped to shell 1)');axs[1,0].set_ylabel('C_n (residual amplitude)')
h=[Line2D([],[],color='tab:blue',marker='o',ls='',ms=11,label='A residual (shell A/A_all − AF_smooth)'),Line2D([],[],color='tab:red',marker='o',ls='',ms=11,label='B residual (shell B/B_all − sheet BL)'),
   Line2D([],[],color='k',lw=1.8,alpha=.6,label='fixed on sky: same ℓ_n at all z (anchored at shell 1)'),Line2D([],[],color='k',ls='--',lw=2.2,label='fixed Mpc: offset from ℓ*=258.5° (CoE-axis ℓ) shrinks ∝ d₁/d'),
   Line2D([],[],color='k',marker='s',ls='',ms=11,label='bulk flow expectation: C_n decays with z')]
fig.legend(handles=h,loc='lower center',ncol=2,fontsize=17,frameon=False)
fig.suptitle('Migkas+2020 313 clusters (VizieR J/A+A/636/A15) split in 3 z-shells (N≈104 each); Δℓ=90° slices, OLS L_X–T; 500 bootstraps',fontsize=19)
fig.tight_layout(rect=(0,0.13,1,0.96));fig.savefig('zshell_harmonic_phase_vs_z.png',dpi=80)
