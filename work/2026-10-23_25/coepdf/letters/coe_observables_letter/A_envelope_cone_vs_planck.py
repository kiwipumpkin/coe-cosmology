# Smooth (unbinned) construction of Migkas et al. 2020 (A&A 636, A15) Fig.6 A/A_all vs Galactic longitude.
# Envelopes: upper y+eh, central y, lower y-el; bin max/min taken as the smooth system's max/min.
# KIND='cone'
#  step: each bin flat across its own 5 deg spacing (piecewise constant, full 360 deg), 0.01 deg grid.
#  cone: Migkas Sect. 5.2: "we consider regions centered at l with a width of Delta l = 90o ... the whole Galactic
#        latitude range b in (-90,+90) is covered by every region"; Fig.6 caption: "Every region covers a sky area of
#        Delta l = 90 and Delta b = 180". (Fig.6 is NOT a cone scan; cones of 45-90 deg are used for the 2D maps.)
#        So each bin's value is spread over its true 90 deg footprint; smooth curve = overlap-weighted mean of the
#        18 bins covering each longitude (equivalent to the step curve convolved with a 90 deg boxcar).
# FFT to m=2500, dipole m=1 removed. Top-hat steps produce sinc ringing at high m; the 90 deg overlap smooths.
import numpy as np,matplotlib,json;matplotlib.use('Agg');import matplotlib.pyplot as plt
KIND='cone'
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
N=36000;g=np.arange(N)*0.01
def smooth(v):
  W=5.0 if KIND=='step' else 90.0
  num=np.zeros(N);den=np.zeros(N)
  for li,vi in zip(l,v):
    dd=(g-li+180)%360-180;w=(np.abs(dd)<W/2)|(np.isclose(dd,-W/2));num+=w*vi;den+=w
  f=np.where(den>0,num/np.maximum(den,1),np.nan)
  if np.isnan(f).any():f[np.isnan(f)]=np.interp(g[np.isnan(f)],g[~np.isnan(f)],f[~np.isnan(f)],period=360)
  return f
MM=2500;m=np.arange(1,MM+1)
def spec(v):c=np.fft.rfft(smooth(v))/N;return 4*np.abs(c[1:MM+1])**2
S={'upper envelope (max)':(spec(y+eh),'crimson'),'central values':(spec(y),'k'),'lower envelope (min)':(spec(y-el),'royalblue')}
P=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');ll=P[:,0];C=P[:,1]*2*np.pi/(ll*(ll+1))
k=m>=2;fig,ax=plt.subplots(3,1,figsize=(11,13))
for a,xs in zip(ax[:2],['log','linear']):
  a.plot(ll,P[:,1],color='darkorange',lw=1,label='Planck 2018 TT $D_\\ell$');a.set_xscale(xs);a.grid(alpha=.3)
  a.set_ylabel('$D_\\ell$ [$\\mu$K$^2$]',color='darkorange');a.set_xlabel('m or $\\ell$ (%s)'%xs);t=a.twinx();mx=0
  for n,(s,c) in S.items():D=m[k]*(m[k]+1)*s[k]/(2*np.pi);mx=max(mx,D.max());t.plot(m[k],D,'-',color=c,lw=.7,label='Migkas '+n)
  t.set_ylim(0,mx*1.1);t.set_ylabel('Migkas m(m+1)amp$^2$/2$\\pi$');h1,l1=a.get_legend_handles_labels();h2,l2=t.get_legend_handles_labels();a.legend(h1+h2,l1+l2,fontsize=8,loc='upper right')
ax[0].set_title('Migkas A/A_all as smooth %s curve (%s), Fourier to m=2500, dipole removed, vs Planck'%(KIND,'5$^\\circ$ top-hats' if KIND=='step' else '90$^\\circ$ overlapping regions'))
ax[2].loglog(ll,C/C.sum(),color='darkorange',lw=1,label='Planck $C_\\ell$')
for n,(s,c) in S.items():ax[2].loglog(m[k],s[k]/s[k].sum(),'-',color=c,lw=.7,label='Migkas '+n)
ax[2].set(xlabel='m or $\\ell$',ylabel='fractional power over m,$\\ell\\geq2$');ax[2].grid(alpha=.3,which='both');ax[2].legend(fontsize=8)
fig.text(.5,.004,'m = azimuthal harmonic in Galactic longitude along a ring; $\\ell$ = spherical multipole. Not the same quantity.',ha='center',fontsize=8)
plt.tight_layout(rect=[0,.015,1,1]);plt.savefig('A_envelope_%s_vs_planck.png'%KIND,dpi=140)
out={}
for n,(s,c) in S.items():
  p=s[k]/s[k].sum();ms=m[k]
  out[n]={'m2_8_pct':list(np.round(100*p[:7],2)),'frac_m_ge_10':float(p[ms>=10].sum()),'frac_m_ge_36':float(p[ms>=36].sum()),'frac_m_ge_100':float(p[ms>=100].sum()),
   'p_at_m':{int(q):float(p[ms==q][0]) for q in [10,36,72,100,500,1000,2500]},'Dmax_m':int(ms[np.argmax(ms*(ms+1)*s[k])])}
  print(KIND,n,out[n])
json.dump(out,open('A_envelope_%s_vs_planck.json'%KIND,'w'),indent=1)
