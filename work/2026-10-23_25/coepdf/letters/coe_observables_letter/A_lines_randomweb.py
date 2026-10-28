# Random "cosmic web" test: random clusters (+, compact), masses (+, broader) and voids (-, large) on the sphere,
# averaged over Migkas-style 90-deg longitude strips (all latitudes, cos b weight) at the 71 bin centres,
# added to the connect-the-dots A/A_all, then FFT power (dipole removed) vs Planck 2018 TT. 30 realisations.
# Chosen (not fitted) toy settings: 60 clusters sigma 3deg amp +0.5, 30 masses sigma 8deg amp +0.3, 20 voids sigma 15deg amp -0.3 (fractional A).
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
rng=np.random.default_rng(1)
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];o=np.argsort(l);l,y=l[o],y[o]
lg=np.radians(np.arange(0.5,360,1));bg=np.radians(np.arange(-89.5,90,1));LL,BB=np.meshgrid(lg,bg);W=np.cos(BB)
xyz=np.stack([np.cos(BB)*np.cos(LL),np.cos(BB)*np.sin(LL),np.sin(BB)])
def blobs(n,sig,amp):
  m=np.zeros_like(LL)
  for _ in range(n):
    v=rng.normal(size=3);v/=np.linalg.norm(v);ang=np.arccos(np.clip(np.tensordot(v,xyz,1),-1,1))
    m+=amp*np.exp(-0.5*(ang/np.radians(sig))**2)
  return m
def strip(m):
  out=[]
  for c in np.radians(l):
    dl=np.angle(np.exp(1j*(LL-c)));k=np.abs(dl)<np.radians(45);out.append((m*W)[k].sum()/W[k].sum())
  return np.array(out)
N=2**16;g=np.arange(N)*360/N
def spec(v):f=np.interp(g,l,v,period=360);P=4*np.abs(np.fft.rfft(f)/N)**2;return P
M=2000;ms=np.arange(2,M+1)
P0=spec(y)[2:M+1];P0/=P0.sum()
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0];Pl=(2*L+1)*2*np.pi*p[:,1]/(L*(L+1));k=L<=M;Pl=Pl[k]/Pl[k].sum()
R=30;Ps=[];Pw=[];rms=[];ex=None
for r in range(R):
  mp=blobs(60,3,.5)+blobs(30,8,.3)+blobs(20,15,-.3);s=strip(mp);s-=s.mean();rms.append(s.std())
  yy=y+s;P=spec(yy)[2:M+1];Ps.append(P/P.sum());q=spec(s+1)[2:M+1];Pw.append(q/q.sum())
  if r==0:ex=(mp,yy)
Ps=np.array(Ps);Pw=np.array(Pw)
fig,ax=plt.subplots(1,2,figsize=(14,5.5))
for a,xl in zip(ax,[(2,35),(2,2000)]):
  a.loglog(L[k],Pl,'r-',lw=1,label='Planck 2018 TT, $\\ell\\geq2$')
  a.loglog(ms,P0,'k-',lw=1.2,label='Migkas lines (no web)')
  a.fill_between(ms,np.percentile(Ps,16,0),np.percentile(Ps,84,0),color='royalblue',alpha=.3,lw=0)
  a.loglog(ms,np.median(Ps,0),'b-',lw=1,label=f'Migkas + random web, median & 68% of {R}')
  a.loglog(ms,np.median(Pw,0),'g--',lw=1,label='random web alone (strip-averaged), median')
  a.set_xlim(*xl);a.set_xlabel('m or $\\ell$');a.set_ylabel('fraction of power (m or $\\ell\\geq2$)');a.grid(alpha=.3,which='both')
ax[0].set_ylim(1e-6,1);ax[0].legend(fontsize=8);ax[0].set_title('m, $\\ell$ = 2..35');ax[1].set_title('m, $\\ell$ = 2..2000')
fig.text(.5,.01,f'Dipole removed. Toy web: 60 clusters (3$^\\circ$, +0.5), 30 masses (8$^\\circ$, +0.3), 20 voids (15$^\\circ$, -0.3); 90$^\\circ$ strip average gives rms {np.mean(rms):.3f} in A/A_all. m azimuthal, $\\ell$ spherical: not the same quantity.',ha='center',fontsize=8)
plt.tight_layout(rect=(0,.04,1,1));plt.savefig('A_lines_randomweb_spectrum.png',dpi=140)
mp,yy=ex;fig,ax=plt.subplots(1,2,figsize=(14,4.5))
im=ax[0].imshow(mp,origin='lower',extent=[0,360,-90,90],cmap='RdBu_r',aspect='auto');plt.colorbar(im,ax=ax[0]);ax[0].set(title='One random web realisation (clusters, masses, voids)',xlabel='l [deg]',ylabel='b [deg]')
ax[1].plot(l,y,'k.-',label='Migkas lines');ax[1].plot(l,yy,'b.-',label='+ random web (strip-averaged)');ax[1].set(xlabel='l [deg]',ylabel='A/A_all',xlim=(0,360));ax[1].legend(fontsize=8);ax[1].grid(alpha=.3)
plt.tight_layout();plt.savefig('A_lines_randomweb_curve.png',dpi=140)
print('rms',np.mean(rms));print('share m2-8 no web',P0[:7].sum(),'with web median',np.median(Ps[:,:7].sum(1)),'web alone',np.median(Pw[:,:7].sum(1)))
