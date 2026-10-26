# REGENERATED for rev7: CatWISE-locked density laws at the chosen background D=34.96 Gpc (rev5/6: D=58.2).
# Each law keeps its rev5 shape parameters (n, alpha, beta, eta) and its scale L is re-solved so that
# [rho(D-d)-rho(D+d)]/rho(D) = 0.031 with d = 14.3005 Gpc; rho(D)=0.3 rho_crit. eta-edge law: eta, L from the D=34.96 Migkas refit (p1).
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.optimize import brentq
D,d,dl=34.96,14.3005,0.031
P1=json.load(open('/workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json'))['34.96']
laws={'exp':(lambda r,L:np.exp(-r/L),'exp: $e^{-r/L}$'),'pl':(lambda r,L:(1+r/L)**-1.2,'power law $(1+r/L)^{-1.2}$'),
 'cored':(lambda r,L:1/(1+(r/L)**2),'cored $1/(1+(r/L)^2)$'),'gauss':(lambda r,L:np.exp(-0.5*(r/L)**2),'Gaussian $e^{-(r/L)^2/2}$'),
 'soft':(lambda r,L:(1+(r/L)**0.3)**-1.0,r'soft $(1+(r/L)^{0.3})^{-1}$'),'eta':(lambda r,L:np.exp(-(r/L)**0.25),r'stretched $e^{-(r/L)^{0.25}}$')}
r=np.linspace(0,80,1601);fig,ax=plt.subplots(1,2,figsize=(13,5.2));out={}
for k,(f,lab) in laws.items():
    g=lambda lL:(f(D-d,np.exp(lL))-f(D+d,np.exp(lL)))/f(D,np.exp(lL))-dl
    G=np.linspace(-3,15,2000);v=np.array([g(x) for x in G]);i=np.where(np.sign(v[:-1])!=np.sign(v[1:]))[0]
    lL=brentq(g,G[i[-1]],G[i[-1]+1]);L=np.exp(lL);rho=0.3*f(r,L)/f(D,L);out[k]=dict(L=L,rho0=float(rho[0]))
    for a in ax:a.plot(r,rho,lw=2,label=f'{lab}, L={L:.4g} Gpc')
e,Le=P1['eta'],P1['L'];re=0.3*np.exp((D/Le)**e-(r/Le)**e)
for a in ax:a.plot(r,re,'k--',lw=1.5,label=rf'$\eta$-edge (Migkas refit): $\eta$={e:.3f}, L={Le:.3f} Gpc')
for a in ax:a.axvline(D,color='grey',lw=1);a.axvspan(D-d,D+d,color='grey',alpha=.12,label=r'CatWISE lock window $D\pm d$');a.set_xlabel('r from CoE [Gpc]');a.set_ylabel(r'$\rho/\rho_{\rm crit}$')
ax[0].set_xlim(0,80);ax[0].set_ylim(0,1.0);ax[1].set_xlim(D-1.3*d,D+1.3*d);ax[1].set_ylim(0.28,0.32)
ax[0].set_title('Density laws at D=34.96 Gpc, $\\rho(D)=0.3\\rho_{\\rm crit}$');ax[1].set_title('Zoom: each law has $\\delta=0.031$ across $D\\pm d$')
h,l=ax[0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=2,frameon=False,fontsize=8.5);fig.tight_layout(rect=(0,0.22,1,1));fig.savefig('d1_migkas_laws.png',dpi=140,facecolor='white')
json.dump(out,open('d1_migkas_laws.json','w'),indent=1);print(out)
