# Power spectrum (amp^2 per m) of the connect-the-dots A/A_all Fourier sums truncated at m=2, 8, 35, 2000.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];o=np.argsort(l);l,y=l[o],y[o]
N=2**16;g=np.arange(N)*360/N;f=np.interp(g,l,y,period=360);C=np.fft.rfft(f)/N;A2=4*np.abs(C)**2  # amp^2 of cos/sin term at m
fig,axs=plt.subplots(2,2,figsize=(12,8))
for a,(M,cl) in zip(axs.flat,[(2,'orange'),(8,'green'),(35,'purple'),(2000,'b')]):
  m=np.arange(1,M+1);p=A2[1:M+1]
  a.loglog(m,p,'o-' if M<=35 else '-',color=cl,ms=4,lw=1);a.plot(1,p[0],'r*',ms=10,label='m=1 dipole')
  a.set_title(f'Fourier sum m=1..{M}');a.set_xlabel('m (azimuthal)');a.set_ylabel('amplitude$^2$ (dimensionless)');a.grid(alpha=.3,which='both');a.set_xlim(.8,max(M*1.3,3));a.legend(fontsize=8)
plt.tight_layout();plt.savefig('A_lines_harm_spectra.png',dpi=140)
print(A2[1:9])
