import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from build import shp,PC,GR
from cosmo import sigR,h
O='/workspace/redweb2/realweb/'
b={g['name']:np.load(O+f"bodies_g0_{g['name']}.npz") for g in GR}
Hp=np.vstack([b[n]['Hp'] for n in b]);Hm=np.concatenate([b[n]['Hm'] for n in b]);Vp=np.vstack([b[n]['Vp'] for n in b])
rng=np.random.default_rng(0);n=3000000;u=rng.random(n);r=(u*15**3)**(1/3);v=rng.normal(size=(n,3));v/=np.linalg.norm(v,axis=1)[:,None];dU=np.linalg.norm(r[:,None]*v-PC,axis=1)
E=np.linspace(28,56,15);c=0.5*(E[1:]+E[:-1]);Vs=np.histogram(dU,E)[0];ref=shp(c/42);ref/=np.average(ref,weights=Vs)
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
for lab,P_,w in (('halo mass',Hp,Hm),('halo number',Hp,None),('void number',Vp,None)):
  h_=np.histogram(np.linalg.norm(P_-PC,axis=1),E,weights=w)[0]/Vs;h_/=np.average(h_,weights=Vs);ax[0].plot(c,h_,'o-',ms=3,label=lab+' (group 0 = 1% sample)')
ax[0].plot(c,ref,'k-',lw=2,label='red ρ(r) shape');ax[0].set_xlabel('distance from CoE (Gpc)');ax[0].set_ylabel('shell density / mean');ax[0].legend(fontsize=8);ax[0].set_title('Shell-averaged body density vs red ρ(r) (web ball r<15 Gpc)')
# homogeneity: halo counts in cubes of side ~ (4/3 pi)^(1/3) * 100/h Mpc (volume = 100/h sphere), group0 fine region, shot noise subtracted
Rh=100/h;L=(4/3*np.pi)**(1/3)*Rh/1000;P0=b['fine']['Hp'];m0=b['fine']['Hm'];sel=np.all(np.abs(P0)<3.4,1);P0,m0=P0[sel],m0[sel]
nb=int(6.8//L);ed=np.linspace(-nb*L/2,nb*L/2,nb+1);Hn=np.histogramdd(P0,[ed]*3)[0];Hw=np.histogramdd(P0,[ed]*3,weights=m0)[0]
mn=Hn.mean();var_n=Hn.var()/mn**2-1/mn;w2=(m0**2).sum()/Hw.sum()**2*Hn.size;var_w=Hw.var()/Hw.mean()**2-w2
print('cube side Mpc',L*1000,'cells',Hn.size,'mean halos/cell (1%)',mn,'sigma_n(shot-subtracted)',np.sqrt(max(var_n,0)),'sigma_massweighted',np.sqrt(max(var_w,0)),'linear matter sigma(100/h)',sigR(Rh),'raw sigma_n',Hn.std()/mn)
ax[1].hist((Hn/mn).ravel(),40);ax[1].set_title(f'halo counts in {L*1000:.0f} Mpc cubes (1% sample): raw σ={Hn.std()/mn:.2f}, shot-noise-subtracted σ={np.sqrt(max(var_n,0)):.3f}; linear matter σ(100/h)={sigR(Rh):.3f}',fontsize=9)
plt.tight_layout();plt.savefig(O+'density_check.png',dpi=110);plt.close()
# map slice
fig,ax=plt.subplots(1,2,figsize=(16,7.5));d=np.load(O+'delta_fine_slice.npy');g=GR[0];N=g['N'];x=(np.arange(N)+0.5)*g['L']/N-g['L']/2
for a,lim in zip(ax,(5,1.5)):
  im=a.pcolormesh(x,x,d.T,cmap='RdBu_r',vmin=-0.3,vmax=0.3,shading='auto');s=lambda P_:(np.abs(P_[:,2])<0.026)&(np.abs(P_[:,0])<lim)&(np.abs(P_[:,1])<lim)
  sh=s(b['fine']['Hp']);a.scatter(b['fine']['Hp'][sh,0],b['fine']['Hp'][sh,1],s=2+20*(b['fine']['Hm'][sh]/2e15),c='k',label='halos (group 0)')
  sv=s(b['fine']['Vp']);a.scatter(b['fine']['Vp'][sv,0],b['fine']['Vp'][sv,1],s=(b['fine']['Vr'][sv]/10)**2*3,facecolors='none',edgecolors='g',label='voids (group 0)')
  t=np.linspace(0,2*np.pi,300)
  for R,cc in ((0.5,'m'),(0.9,'orange'),(2,'c'),(4,'purple')):a.plot(R*np.cos(t),R*np.sin(t),'--',color=cc,lw=1.2)
  a.annotate('',xy=(0.9*lim,0),xytext=(0,0),arrowprops=dict(arrowstyle='->',color='darkgreen',lw=2));a.text(0.5*lim,0.05*lim,'to CoE',color='darkgreen')
  a.plot(0,0,'k+',ms=12,mew=2);a.set_xlim(-lim,lim);a.set_ylim(-lim,lim);a.set_aspect('equal');a.set_xlabel('galactic x (Gpc)');a.set_ylabel('galactic y (Gpc)')
ax[0].legend(loc='lower left');plt.colorbar(im,ax=ax[1],label='linear δ (smoothed at 52 Mpc cell), z=0 slice')
ax[0].set_title('Realistic web: GRF δ (z≈0 slice) + group-0 bodies in |z|<26 Mpc; dashed R=0.5/0.9/2/4 Gpc');ax[1].set_title('zoom ±1.5 Gpc')
plt.tight_layout();plt.savefig(O+'realweb_map.png',dpi=100)
