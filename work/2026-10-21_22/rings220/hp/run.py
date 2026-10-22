import numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.special import sph_harm_y
ell=np.arange(1,351);lmax=350;nside=256
def unit(a):return a/np.linalg.norm(a,axis=-1,keepdims=True)
def build(width,seed):
    rg=np.random.default_rng(seed);target=np.exp(-.5*((ell-220)/width)**2)
    alm=np.zeros(hp.Alm.getsize(lmax),complex)
    for i,l in enumerate(ell):
        if target[i]<1e-9:continue
        axes=unit(rg.normal(size=(5,3)));w=rg.choice([-1.,1.],size=5)
        th=np.arccos(np.clip(axes[:,2],-1,1));ph=np.arctan2(axes[:,1],axes[:,0])
        from scipy.special import eval_legendre
        S=max(float(w@eval_legendre(int(l),np.clip(axes@axes.T,-1,1))@w),1e-12)
        A=np.sqrt(target[i]*(2*l+1)/(l*(l+1)*S))
        m=np.arange(l+1)
        for j in range(5):
            Y=sph_harm_y(l,m,th[j],ph[j])
            alm[hp.Alm.getidx(lmax,l,m)]+=A*w[j]*4*np.pi/(2*l+1)*np.conj(Y)
    cl=hp.alm2cl(alm);L=np.arange(lmax+1);D=L*(L+1)*cl/2/np.pi
    s=5000/D.max();alm*=np.sqrt(s);D*=s
    return target*5000,D,hp.alm2map(alm,nside)*1e0
P=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt')
fig,ax=plt.subplots(2,1,figsize=(11,8))
for k,(wd,sd) in enumerate([(32,32),(8,8)]):
    t,D,m=build(wd,sd);L=np.arange(351)
    ax[k].errorbar(P[:,0],P[:,1],yerr=[P[:,2],P[:,3]],fmt='.',color='0.6',ms=2,lw=.4,label='Planck 2018 TT',zorder=1)
    ax[k].plot(L[1:],D[1:],'C3-',lw=1.4,label='Measured map D_ℓ (scaled, peak 5000 μK²)')
    ax[k].plot(ell,t,'k--',lw=1,label='Prescribed target')
    ax[k].set(xlim=(2,600),ylim=(0,6500),ylabel='D_ℓ [μK²]',title=f'Peak ℓ=220, σℓ={wd}');ax[k].legend(loc='lower center',fontsize=8)
    np.save(f'map_s{wd}.npy',m)
    plt.figure(10+k);hp.mollview(m,title=f'Source surface, σℓ={wd} (scaled, μK-equivalent)',cmap='RdBu_r',unit='μK',min=-np.percentile(abs(m),99),max=np.percentile(abs(m),99));plt.savefig(f'surf_s{wd}.png',dpi=120)
    print(wd,'rms',m.std(),'peak l',L[D.argmax()])
ax[1].set_xlabel('ℓ');fig.tight_layout();fig.savefig('spec_5000.png',dpi=110)
