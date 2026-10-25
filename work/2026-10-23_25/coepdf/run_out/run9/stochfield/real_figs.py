# One random sky per model (fixed seed 42): D_real = D_bg + D_field * chi2_{2l+1}/(2l+1); background is deterministic (no scatter).
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
import acoustic_fit as A
from acoustic_model import Dfield,Dfix,LMAX
from fwd import Fwd
def real(Df,l,seed=42):return Df*np.random.default_rng(seed).chisquare(2*l+1)/(2*l+1)
J=json.load(open('acoustic_fit.json'));Lt=A.Lt;Pf=A.Pf;Pb=A.Pb;ll=np.arange(LMAX+1)
for nm in ('free_W','fixed_shell'):
    S=J[nm];E=S['best'];p=E['p'];Df=Dfield(p)*np.array(S['corr']);Dm=Df+Dfix;Dr=Dfix+real(Df,np.maximum(ll,1))
    st={k:A.stats(Dr,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500)).items()};S['realization_seed42']=st
    fig,ax=plt.subplots(2,1,figsize=(12,10))
    for a,(l0,l1),key in zip(ax,((2,30),(2,LMAX)),('l2_30','l2_2500')):
        m=(ll>=l0)&(ll<=l1);a.plot(ll[m],Lt[m],'k-',lw=1.8,label=f'ΛCDM  χ²(Planck err)={J["lcdm"][key]["chi2_planck"]:.1f}, Pearson r²={J["lcdm"][key]["r2"]:.3f}')
        cv=np.sqrt(2/(2*np.maximum(ll,1)+1))*Df;a.fill_between(ll[m],(Dm-cv)[m],(Dm+cv)[m],color='r',alpha=.2,label='±1σ cosmic variance (field part)')
        a.plot(ll[m],Dm[m],'r-',lw=.7,alpha=.35,label='CoE ensemble mean (faint)')
        a.plot(ll[m],Dr[m],'r-',lw=1.0 if l1>30 else 1.6,marker='' if l1>30 else 'o',ms=3,label=f'CoE one realization (seed 42, {nm})  χ²(Planck err)={st[key]["chi2_planck"]:.1f}, Pearson r²={st[key]["r2"]:.3f}')
        mm=(Pf[:,0]>=l0)&(Pf[:,0]<=l1)
        if l1==30:a.errorbar(Pf[mm,0],Pf[mm,1],yerr=[Pf[mm,2],Pf[mm,3]],fmt='o',color='grey',ms=4,label='Planck 2018');a.set_ylim(0,2500)
        else:a.plot(Pf[mm,0],Pf[mm,1],'.',color='grey',ms=1.5,alpha=.5,label='Planck 2018 (unbinned)');a.set_xscale('log');a.set_ylim(0,7000)
        a.set_xlabel('ℓ');a.set_ylabel('D_ℓ (μK²)');a.legend(loc='lower center' if l1==30 else 'upper left',fontsize=7,framealpha=.85)
    ax[0].set_title(f'acoustic-feature stochastic field, R=4 Gpc, {nm}: one realization');fig.tight_layout();fig.savefig(f'acoustic_cl_{nm}_real.png',dpi=100);plt.close()
    print(nm,st)
json.dump(J,open('acoustic_fit.json','w'),indent=1)
# stage A (l<=30 / l<=100) panels with a realization
for tag in ('l30_bg','l30_nobg','l100_bg','l100_nobg'):
    R=json.load(open(f'fit_{tag}.json'));C=json.load(open(f'chosen_{tag}.json'));P_=np.array(R['planck']);l=P_[:,0].astype(int);F=Fwd(lmax=int(R['LMAX']))
    Df=F.Dl(C['p']);Dfx=np.array(R['Dfix']);Dm=Df+Dfx;Dr=Dfx+real(Df,l)
    e=np.where(Dr>P_[:,1],P_[:,3],P_[:,2]);c2=float(np.sum((Dr-P_[:,1])**2/e**2));r2=float(np.corrcoef(Dr,P_[:,1])[0,1]**2)
    fig,ax=plt.subplots(figsize=(11,5));ax.plot(l,R['lcdmD'],'k-',lw=2,label=f'ΛCDM χ²(Planck err)={R["lcdm"]["chi2_planck"]:.1f}, Pearson r²={R["lcdm"]["r2"]:.3f}')
    cv=np.sqrt(2/(2*l+1))*Df;ax.fill_between(l,Dm-cv,Dm+cv,color='r',alpha=.2,label='±1σ cosmic variance (field part)')
    ax.plot(l,Dm,'C0-',lw=.7,alpha=.35,label='ensemble mean (faint)');ax.plot(l,Dr,'C0-o',ms=3,lw=1.3,label=f'one realization (seed 42): χ²(Planck err)={c2:.1f}, Pearson r²={r2:.3f}')
    ax.errorbar(l,P_[:,1],yerr=[P_[:,2],P_[:,3]],fmt='o',color='grey',ms=3,label='Planck 2018');ax.set_ylim(0,2500 if R['LMAX']<=30 else 3300)
    ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)');ax.set_title(f'stochastic field, R=4 Gpc, one realization [{tag}]');ax.legend(loc='lower center',fontsize=7,framealpha=.85);fig.tight_layout();fig.savefig(f'a_cl_{tag}_real.png',dpi=100);plt.close()
    print(tag,'real chi2',round(c2,1),'r2',round(r2,3))
