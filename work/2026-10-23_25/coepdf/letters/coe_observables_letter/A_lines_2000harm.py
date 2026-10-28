# Migkas A/A_all connect-the-dots curve rebuilt from its Fourier series truncated at m=2000 (and m=2,8,35 for comparison).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high'];o=np.argsort(l);l,y,el,eh=l[o],y[o],el[o],eh[o]
N=2**16;g=np.arange(N)*360/N;f=np.interp(g,l,y,period=360);C=np.fft.rfft(f)
def rec(M):c=C.copy();c[M+1:]=0;return np.fft.irfft(c,N)
fig,ax=plt.subplots(figsize=(11,5));ax.errorbar(l,y,[el,eh],fmt='ko',ms=3,capsize=2,label='Migkas et al. (2020) $A/A_{\\rm all}$')
for M,cl,lw in [(2,'orange',1),(8,'green',1),(35,'purple',1),(2000,'b',1.5)]:
  r=rec(M);ax.plot(g,r,color=cl,lw=lw,label=f'Fourier sum m=0..{M} (max |diff| from lines {np.abs(r-f).max():.2g})')
ax.set(xlabel='Galactic longitude l [deg]',ylabel='$A/A_{\\rm all}$',xlim=(0,360));ax.set_xticks(range(0,361,60));ax.grid(alpha=.3);ax.legend(fontsize=8)
plt.tight_layout();plt.savefig('A_lines_2000harm.png',dpi=140)
