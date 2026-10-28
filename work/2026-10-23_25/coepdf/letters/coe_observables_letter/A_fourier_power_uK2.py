# Planck 2018 TT D_l in uK^2 (l=1..35, dipole from A=3362.08 uK) vs Migkas A/A_all weighted Fourier M=35 power per m (dimensionless, twin axis).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];s=(d['A_err_low']+d['A_err_high'])/2
M=35;t=np.radians(l);X=[np.ones_like(t)]
for m in range(1,M+1):X+=[np.cos(m*t),np.sin(m*t)]
X=np.array(X).T;w=1/s;c=np.linalg.lstsq(X*w[:,None],y*w,rcond=None)[0]
amp2=c[1::2]**2+c[2::2]**2;ms=np.arange(1,M+1)
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0];Dl=p[:,1];k=L<=35
Ad=3362.08;C1=4*np.pi/9*Ad**2;D1=1*2*C1/(2*np.pi)  # T=A cos(theta): a10=A*sqrt(4pi/3), C1=a10^2/3
fig,axs=plt.subplots(1,2,figsize=(13,5))
for a,dip in zip(axs,[True,False]):
  LL=np.r_[1,L[k]] if dip else L[k];DD=np.r_[D1,Dl[k]] if dip else Dl[k]
  a.semilogy(LL,DD,'s-',color='darkorange',label='Planck 2018 TT $D_\\ell$ [$\\mu$K$^2$]'+(' (dipole $A$=3362 $\\mu$K)' if dip else ''))
  a.set_xlabel('m or $\\ell$');a.set_ylabel('$D_\\ell=\\ell(\\ell+1)C_\\ell/2\\pi$  [$\\mu$K$^2$]');a.grid(alpha=.3,which='both')
  b=a.twinx();mm=ms if dip else ms[1:];aa=amp2 if dip else amp2[1:]
  b.semilogy(mm,aa,'o-',color='royalblue',label='Migkas Fourier M=35, amp$^2$ (dimensionless)');b.set_ylabel('Migkas amp$^2$ (dimensionless)')
  a.set_title('Dipole kept' if dip else 'Dipole removed');h1,l1=a.get_legend_handles_labels();h2,l2=b.get_legend_handles_labels();a.legend(h1+h2,l1+l2,fontsize=8,loc='lower left')
fig.text(.5,.01,'m = azimuthal ring harmonic; $\\ell$ = spherical multipole. Not the same quantity. Migkas axis is dimensionless; only Planck is in $\\mu$K$^2$.',ha='center',fontsize=9)
plt.tight_layout(rect=(0,.04,1,1));plt.savefig('A_fourier_power_uK2.png',dpi=140)
print(D1,Dl[:3])
