# A/A_all data + locked CoE eta-edge fit, residuals, and Fourier decomposition of the residuals in l.
# Same data and model files as Fig04_amp_full.py.
import os,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__))
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high'];s=0.5*(el+eh)
c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1)
m=np.interp(l,c[:,0],c[:,1],period=360);r=y-m
chi2=float(np.sum(r**2/s**2));r2=float(1-np.sum(r**2)/np.sum((y-y.mean())**2))
# least-squares harmonic fit (points are 5 deg apart but not necessarily full-circle uniform)
M=8;L=np.radians(l);cols=[np.ones_like(L)]
for k in range(1,M+1):cols+=[np.cos(k*L),np.sin(k*L)]
X=np.stack(cols,1);coef=np.linalg.lstsq(X,r,rcond=None)[0]
a=coef[1::2];b=coef[2::2];amp=np.hypot(a,b);ph=np.degrees(np.arctan2(b,a))%360  # r ~ amp*cos(k(l-ph/k))
pw=amp**2/2;frac=pw/np.sum(pw);peak=(ph/np.arange(1,M+1))%(360/np.arange(1,M+1))
fig,ax=plt.subplots(3,1,figsize=(11,9.5),gridspec_kw=dict(height_ratios=[3,1.4,2]))
ax[0].errorbar(l,y,[el,eh],fmt='ko',ms=3.5,capsize=2,label='Migkas et al. (2020) $A/A_{\\rm all}$, 71 cones')
ax[0].plot(c[:,0],c[:,1],'r',lw=2.2,label=f'locked CoE $\\eta$-edge fit, $\\chi^2$={chi2:.3f}, $r^2$={r2:.3f}')
ax[0].set_ylabel('$A/A_{\\rm all}$');ax[0].legend(fontsize=9);ax[0].grid(alpha=.3);ax[0].set_xlim(0,360)
ax[1].errorbar(l,r,[el,eh],fmt='ko',ms=3,capsize=2);ax[1].axhline(0,color='r')
lf=np.linspace(0,360,721);ax[1].plot(lf,np.stack([np.ones_like(lf)]+sum([[np.cos(k*np.radians(lf)),np.sin(k*np.radians(lf))] for k in range(1,M+1)],[]),1)@coef,'g',lw=1.3,label=f'sum of m=1..{M}')
ax[1].set_ylabel('data − fit');ax[1].set_xlabel('Galactic longitude [deg]');ax[1].set_xlim(0,360);ax[1].grid(alpha=.3);ax[1].legend(fontsize=8)
k=np.arange(1,M+1);ax[2].bar(k,amp,color='steelblue')
for i in range(M):ax[2].text(k[i],amp[i],f'{100*frac[i]:.1f}%\npeak {peak[i]:.0f}°',ha='center',va='bottom',fontsize=8.5)
ax[2].set_xlabel('harmonic m (in Galactic longitude)');ax[2].set_ylabel('amplitude');ax[2].set_xticks(k);ax[2].set_ylim(0,amp.max()*1.45);ax[2].grid(alpha=.3,axis='y')
ax[2].set_title('Fourier decomposition of residuals: amplitude, fractional power, first-peak longitude',fontsize=10)
fig.tight_layout();fig.savefig(H+'/A_fit_residual_fourier.png',dpi=140,facecolor='white')
out=dict(chi2=chi2,r2=r2,offset=float(coef[0]),modes=[dict(m=int(k[i]),amp=float(amp[i]),phase_deg=float(ph[i]),first_peak_l=float(peak[i]),frac_power=float(frac[i])) for i in range(M)])
json.dump(out,open(H+'/A_fit_residual_fourier.json','w'),indent=1);print(json.dumps(out,indent=0))
