import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
import ff_prof as F
U=json.load(open('fullfree_results.json'))['universes'];Pl=np.loadtxt('../../../reproduce_cmb/planckTT.txt')
pick={'grid_C3_0.1_f1000':'mig','sim101':'mig360'}
for r in U:
    if r['name'] not in pick:continue
    k=pick[r['name']];m=r[k];th=np.array(m['theta']);tag='fullfree_%s_%s'%(r['name'],'l0m10_100' if k=='mig' else 'l0_0_360')
    plt.figure(figsize=(8,4.5));plt.errorbar(F.L,F.Y,F.S,fmt='o',ms=3,label='Migkas A/A_all');ll=F.L
    plt.plot(ll,F.model(th,r['ratio']),'r-',label='fit')
    plt.xlabel('Galactic longitude l (deg)');plt.ylabel('A/A_all');plt.title('%s  ratio %.4f  chi2 %.4f  R2 %.4f  l0 %.2f'%(r['name'],r['ratio'],m['chi2'],m['R2'],th[1]),fontsize=9);plt.legend();plt.tight_layout();plt.savefig(tag+'_migkas.png',dpi=130);plt.close()
    DR=r['DRV'];Rb=1.5*(1+1/DR);x=np.linspace(0,2.4,600)
    if r['fam']=='slc':lr=-0.041037*np.log(np.maximum(x,1e-6))-0.035*(x-1)**2
    else:lr=sum(c*(np.minimum(x,Rb)-1)**(i+1) for i,c in enumerate(r['C']))
    rho=np.where(x<=Rb,np.exp(lr),0)
    plt.figure(figsize=(7,4.5));plt.plot(x,rho);plt.axvline(1,ls=':',c='k');plt.axvline(Rb,ls='--',c='gray')
    plt.xlabel('R / D  (CoE at 0, observer at 1, mass edge Rb=%.3f)'%Rb);plt.ylabel('rho / rho(observer)');plt.title('%s  C=%s  frac=%g  D/R=%g'%(r['name'],r['C'],r['frac'],DR),fontsize=9);plt.tight_layout();plt.savefig(tag+'_density.png',dpi=130);plt.close()
    l=np.arange(2,31);D=np.array(r['Dl'])[2:31]
    plt.figure(figsize=(8,4.5));plt.errorbar(Pl[:,0],Pl[:,1],[Pl[:,2],Pl[:,3]],fmt='ko',ms=3,label='Planck 2018')
    plt.semilogy(l,np.maximum(D,1e-6),'r.-',label='model');plt.yscale('log');plt.xlabel('l');plt.ylabel('D_l (uK^2)')
    plt.title('%s  CMB chi2(l2-30) %.4g'%(r['name'],r['cmb']),fontsize=9);plt.legend();plt.tight_layout();plt.savefig(tag+'_Dl.png',dpi=130);plt.close()
    print(tag)
