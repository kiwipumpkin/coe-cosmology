# PROVENANCE COPY of /workspace/coepdf/dens_plots.py (generated rev6 Fig 15 web ingredients; not rerun for rev7)
import numpy as np,pandas as pd,json,sys,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.integrate import quad
plt.rcParams.update({'font.size':11})
def legbelow(fig,axs,n=3):
    h,l=[],[]
    for a in axs:
        hh,ll=a.get_legend_handles_labels();h+=hh;l+=ll
    d=dict(zip(l,h));fig.legend(d.values(),d.keys(),loc='lower center',ncol=n,frameon=False,fontsize=9)
N='/workspace/cosmology-course/notes/'
# 1 Migkas density laws
df=pd.read_csv(N+'catwise-density-levelset.csv');e=pd.read_csv(N+'catwise-eta-density.csv')
J=json.load(open(N+'catwise-density-levelset-fit.json'))['density_params']
lab={'exp':r'exp: $e^{-r/L}$, L=923','pl':r'power law: $(1+r/L)^{-n}$, L=1049, n=1.20','cored':r'cored: $1/(1+(r/L)^2)$, L=322',
 'gauss':r'Gaussian: $e^{-(r/L)^2/2}$, L=232','soft':r'soft core+tail: $(1+(r/L)^{\alpha})^{-\beta}$, L=5110, $\alpha$=0.30, $\beta$=1',
 'eta':r'stretched exp: $e^{-(r/L)^{\eta}}$, L=15246, $\eta$=0.25'}
r=df.r_Gpc.values;fig,ax=plt.subplots(1,2,figsize=(13,5.6))
for c in df.columns[1:]:
    k=c.split('_')[3]
    for a in ax:a.plot(r,df[c],lw=2,label=lab[k])
ce=[c for c in e.columns if 'eta0.05' in c][0]
for a in ax:a.plot(e.r_Gpc,e[ce],'k--',lw=1.5,label=r'$\eta$-edge: $\eta$=0.050, L=0.849 Gpc')
D,d=58.2,14.3005
for a in ax:
    a.axvline(D,color='grey',lw=1);a.axvspan(D-d,D+d,color='grey',alpha=.12,label=r'CatWISE lock window $D\pm d$')
    a.set_xlabel('r from CoE [Gpc]');a.set_ylabel(r'$\rho/\rho_{\rm crit}$')
ax[0].set_xlim(0,max(r));ax[1].set_xlim(D-1.3*d,D+1.3*d);ax[1].set_ylim(0.28,0.32)
ax[0].set_title('Density laws, all normalised to $\\rho(D)=0.3\\rho_{\\rm crit}$');ax[1].set_title(r'Zoom: every law has $[\rho(D-d)-\rho(D+d)]/\rho(D)=0.031$')
legbelow(fig,ax,2);fig.tight_layout(rect=(0,0.2,1,1));fig.savefig('fig/d1_migkas_laws.png',dpi=150)
# 2 web
sys.path.insert(0,'/tmp/coe/work/2026-10-20_21/redweb2/realweb');import cosmo as C
fig,ax=plt.subplots(1,3,figsize=(14,4.8))
M=np.geomspace(1e14,2e15,41);n,s=C.ps_dndlnM(M);ax[0].loglog(M,n,'b',lw=2,label='Press–Schechter $dn/d\\ln M$ (halos)')
ax[0].set_xlabel(r'M [M$_\odot$]');ax[0].set_ylabel(r'Mpc$^{-3}$')
Rv=np.geomspace(10,50,41);ax[1].loglog(Rv,C.svdw_dndlnR(Rv),'r',lw=2,label='Sheth–van de Weygaert $dn/d\\ln R$ (voids)');ax[1].set_xlabel('void radius R [Mpc]')
ax[2].stem([0],[-1],linefmt='r-',markerfmt='ro',basefmt=' ',label=r'void as painted: point deficit $-M_d$ at centre, $M_d=0.85\bar\rho_m(4\pi/3)R_v^3$')
ax[2].stem([1.25],[1],linefmt='b-',markerfmt='bo',basefmt=' ',label=r'compensating wall: 8 random points, $M_d/8$ each, at $1.25R_v$')
ax[2].axhline(0,color='k',lw=.5);ax[2].set_xlim(-0.2,1.6);ax[2].set_ylim(-1.3,1.3);ax[2].set_xlabel(r'$r/R_v$');ax[2].set_ylabel(r'mass / $M_d$')
legbelow(fig,ax,2);fig.tight_layout(rect=(0,0.16,1,1));fig.savefig('fig/d2_web.png',dpi=150)
# 2b modulation shp
CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937)
shp=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
# 3 stochastic
p=json.load(open('/tmp/coe/work/2026-10-21_22/stochfield/acoustic_half_fixed.json'))['r2free']['p']
k=np.linspace(0.5,600,3000)  # 1/Gpc
base=10**p['lgA']*k**p['n']*np.exp(-(k/p['kD'])**2)
ac=1+p['Aw']*np.cos(k*p['rs']+p['phi']);dm=p['B']*np.exp(-k/p['kB'])*np.cos(k*p['r2']+p['phi2'])
fig,ax=plt.subplots(1,2,figsize=(13,5))
km=k/1000
ax[0].plot(km,base,'grey',lw=1.5,label=r'$P_0k^ne^{-(k/k_D)^2}$');ax[0].plot(km,base*ac,'b',lw=1.2,label='+ acoustic term')
ax[0].plot(km,base*(ac+dm),'r',lw=2,label=r'full $P_\delta$ (+ damped $r_2$ term)')
ax[0].set_xlabel(r'k [Mpc$^{-1}$]');ax[0].set_ylabel(r'$P_\delta$ [code units, Gpc$^3$]');ax[0].set_yscale('log');ax[0].set_title('Fitted source spectrum ($r_2$ free)')
ch=np.linspace(2.0,5.2,400);ax[1].plot(ch,np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2),'r',lw=2,label=r'$W(\chi)=e^{-(\chi-R_c)^2/2w^2}$, $R_c$=3.572, $w$=0.336 Gpc')
ax[1].set_xlabel(r'$\chi$ [Gpc]');ax[1].set_ylabel('W');ax[1].set_title('Radial window')
legbelow(fig,ax,2);fig.tight_layout(rect=(0,0.12,1,1));fig.savefig('fig/d3_stoch.png',dpi=150)
# 4 radiation+matter
ae=1/3400;f=lambda a:a/np.sqrt(ae+a);T=quad(f,0,1)[0];sc=13.8e9/T
a=np.geomspace(1e-6,1,300);t=np.array([quad(f,0,x)[0]*sc for x in a])
fig,ax=plt.subplots(1,2,figsize=(13,5))
ax[0].loglog(a,ae/a**4,'orange',lw=2,label=r'radiation $\rho_r/\rho_{m,0}=a_{\rm eq}a^{-4}$');ax[0].loglog(a,a**-3.,'g',lw=2,label=r'matter $\rho_m/\rho_{m,0}=a^{-3}$')
ax[0].axvline(ae,color='k',ls='--',label=r'equality $a_{\rm eq}=1/3400$ (41 kyr)');ax[0].axvline(1/1090,color='m',ls=':',label='recombination z=1089 (297 kyr)')
ax[0].set_xlabel('a');ax[0].set_ylabel(r'$\rho/\rho_{m,0}$')
ax[1].loglog(t/1e3,a,'k',lw=2,label=r'$a(t)$ from $H\propto\sqrt{a_{\rm eq}a^{-4}+a^{-3}}$, $t_0$=13.8 Gyr')
ax[1].loglog(t/1e3,a[-1]*(t/t[-1])**(2/3),'g:',label=r'$t^{2/3}$');ax[1].loglog(t/1e3,ae*(t/40.8e3)**0.5,'orange',ls=':',label=r'$t^{1/2}$')
ax[1].set_ylim(1e-6,1.5);ax[1].set_xlabel('t [kyr]');ax[1].set_ylabel('a')
legbelow(fig,ax,2);fig.tight_layout(rect=(0,0.18,1,1));fig.savefig('fig/d4_radmat.png',dpi=150)
# 5 modulation
fig,ax=plt.subplots(figsize=(13,3.8));xx=np.linspace(0,3,400);ax.plot(xx*42,shp(xx),'b',lw=2,label=r'$m(r)=\exp\sum_k c_k(\min(r/42,R_b)-1)^{k+1}$, $c$=[−0.041, −0.035, 0.1353, 0], $R_b$=2.011')
ax.set_xlabel('distance from CoE axis point [Gpc]');ax.set_ylabel('m (before normalising)');legbelow(fig,[ax],1);fig.tight_layout(rect=(0,0.12,1,1));fig.savefig('fig/d5_mod.png',dpi=150)
print(p)
