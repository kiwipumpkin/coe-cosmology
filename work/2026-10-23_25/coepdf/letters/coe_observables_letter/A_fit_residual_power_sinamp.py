# Power spectrum (amp^2/2 per m, m=1..35) of the Fig 11 residuals (data - locked eta-edge fit), with and without the sheet sin-amp transform (Amp=-0.06, phi=0.35 rad).
import os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__))
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall']
c0=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1)
AMP,PHI=-0.06,0.35;L=np.radians(l)
def spec(c,M):
  r=y-np.interp(l,c[:,0],c[:,1],period=360);X=np.stack([np.ones_like(L)]+sum([[np.cos(k*L),np.sin(k*L)] for k in range(1,M+1)],[]),1)
  co=np.linalg.lstsq(X,r,rcond=None)[0];return np.hypot(co[1::2],co[2::2])**2/2
c1=c0.copy();c1[:,1]*=1-AMP*np.sin(np.radians(c1[:,0])+PHI)
fig,ax=plt.subplots(1,2,figsize=(13,5))
for a,M in zip(ax,[8,35]):
  m=np.arange(1,M+1);p0,p1=spec(c0,M),spec(c1,M)
  a.semilogy(m,p0,'o-',color='gray',ms=4,label=f'Fig 11 original (total {p0.sum():.2e})')
  a.semilogy(m,p1,'s-',color='r',ms=4,label=f'with sin-amp (total {p1.sum():.2e})')
  a.set_title(f'Residual power spectrum, m=1..{M}');a.set_xlabel('harmonic m (Galactic longitude)');a.set_ylabel('power = amp$^2$/2');a.grid(alpha=.3,which='both');a.legend(fontsize=9)
  if M==8:a.set_xticks(m)
  print(M,'orig',np.round(p0[:8],5),'sinamp',np.round(p1[:8],5))
fig.tight_layout();fig.savefig(H+'/A_fit_residual_power_sinamp.png',dpi=140,facecolor='white')
