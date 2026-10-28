# Residuals between Migkas connect-the-dots Fourier power share (m>=2, dipole removed) and Planck 2018 TT (2l+1)C_l share (l>=2), matched m=l, ranges 2..35 and 2..2000.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];o=np.argsort(l);l,y=l[o],y[o]
N=2**16;g=np.arange(N)*360/N;A2=4*np.abs(np.fft.rfft(np.interp(g,l,y,period=360))/N)**2
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0].astype(int);Pl=(2*L+1)*2*np.pi*p[:,1]/(L*(L+1))
fig,ax=plt.subplots(2,2,figsize=(13,8))
for j,M in enumerate([35,2000]):
  m=np.arange(2,M+1);a=A2[2:M+1];a=a/a.sum();k=L<=M;b=Pl[k]/Pl[k].sum();r=a-b
  ax[0,j].plot(m,r,'b.-' if M<=35 else 'b-',lw=1);ax[0,j].axhline(0,color='k',lw=.7);ax[0,j].set(title=f'Residual: Migkas share - Planck share, m=$\\ell$=2..{M}',xlabel='m = $\\ell$',ylabel='residual (fraction of power)');ax[0,j].grid(alpha=.3)
  ax[1,j].semilogy(m,a/b,'b.-' if M<=35 else 'b-',lw=1);ax[1,j].axhline(1,color='k',lw=.7);ax[1,j].set(title='Ratio Migkas / Planck',xlabel='m = $\\ell$',ylabel='ratio');ax[1,j].grid(alpha=.3,which='both')
  if M>35:[x.set_xscale('log') for x in ax[:,j]]
  print(M,'sum|r|',np.abs(r).sum(),'rms r',np.sqrt((r**2).mean()));
  if M==35:
    for mm,rr,aa,bb in zip(m[:8],r[:8],a[:8],b[:8]):print(mm,round(aa,4),round(bb,4),round(rr,4))
fig.text(.5,.005,'Dipole removed; each spectrum normalised to 1 over its range. m azimuthal ring harmonic, $\\ell$ spherical multipole: not the same quantity.',ha='center',fontsize=9)
plt.tight_layout(rect=(0,.03,1,1));plt.savefig('A_lines_vs_planck_residuals.png',dpi=140)
