import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.signal import find_peaks,savgol_filter
import acoustic_fit as A
from acoustic_model import Dfield,Dfix,EX,LMAX
J=json.load(open('acoustic_fit.json'));Lt=A.Lt;Pf=A.Pf;Pb=A.Pb;ll=np.arange(LMAX+1)
def peaks(D,lo=100):
    s=savgol_filter(D[lo:],61,3);pk,_=find_peaks(s,prominence=20);return (pk+lo)[:5].tolist()
pl=json.load(open('acoustic_fit.json')).get('peaks_planck_binned_spline',[]);J['peaks_lcdm']=peaks(Lt)
for nm in [k for k in ('free_W','fixed_shell') if k in J]:
    S=J[nm];C=np.array(S['corr']);E=min(S['Lcurve'],key=lambda e:e['chi2_fit']) ;p=E['p']
    Df=Dfield(p)*C;D=Df+Dfix;cv=np.sqrt(2/(2*ll+1))*Df;S['best']=E;S['peaks_model']=peaks(D)
    fig,ax=plt.subplots(2,1,figsize=(12,10))
    for a,(l0,l1) in zip(ax,((2,30),(2,LMAX))):
        m=(ll>=l0)&(ll<=l1);a.plot(ll[m],Lt[m],'k-',lw=1.8,label=f'ΛCDM  χ²(Planck err)={J["lcdm"]["l2_30" if l1==30 else "l2_2500"]["chi2_planck"]:.1f}, Pearson r²={J["lcdm"]["l2_30" if l1==30 else "l2_2500"]["r2"]:.3f}')
        a.fill_between(ll[m],(D-cv)[m],(D+cv)[m],color='r',alpha=.25,label='CoE stochastic field ±1σ cosmic variance')
        st=E['l2_30' if l1==30 else 'l2_2500'];a.plot(ll[m],D[m],'r-',lw=1.3,label=f'CoE stochastic field ({nm})  χ²(Planck err)={st["chi2_planck"]:.1f}, Pearson r²={st["r2"]:.3f}')
        mm=(Pf[:,0]>=l0)&(Pf[:,0]<=l1)
        if l1==30:a.errorbar(Pf[mm,0],Pf[mm,1],yerr=[Pf[mm,2],Pf[mm,3]],fmt='o',color='grey',ms=4,label='Planck 2018')
        else:a.plot(Pf[mm,0],Pf[mm,1],'.',color='grey',ms=1.5,alpha=.4,label='Planck 2018 (unbinned)');a.errorbar(Pb[:,0],Pb[:,1],yerr=Pb[:,2],fmt='o',color='dimgrey',ms=3,label='Planck binned');a.set_xscale('log');a.set_ylim(0,7000)
        if l1==30:a.set_ylim(0,2500)
        a.set_xlabel('ℓ');a.set_ylabel('D_ℓ (μK²)');a.legend(loc='lower center' if l1==30 else 'upper left',fontsize=7,framealpha=.85)
    ax[0].set_title(f'acoustic-feature stochastic field, R=4 Gpc, {nm}: peaks model {S["peaks_model"]} vs Planck {pl}')
    fig.tight_layout();fig.savefig(f'tmp_cl_{nm}.png',dpi=100);plt.close()
    fig,(a1,a2)=plt.subplots(1,2,figsize=(12,4.5));k=np.geomspace(0.1,1000,3000);a1.loglog(k,EX.Pd(p,k),'r');a1.set_xlabel('k (1/Gpc)');a1.set_ylabel('P_δ (Gpc³)')
    a1.set_title(f'P_δ: n={p["n"]:.2f}, A={p["Aw"]:.2f}, r_s={p["rs"]*1000:.1f} Mpc, φ={p["phi"]:.2f}, k_D={p["kD"]:.0f}/Gpc',fontsize=9)
    a2.plot(EX.ch,EX.W(p)[0],'r');a2.set_xlabel('χ (Gpc)');a2.set_ylabel('W');a2.set_title(f'window centre {p["chiw"]:.3f}, width {p["sw"]:.3f} Gpc (path ends at R=4)');fig.tight_layout();fig.savefig(f'tmp_pk_{nm}.png',dpi=100);plt.close()
    fig,a=plt.subplots(figsize=(12,4.5));Db=A.binned(D);Lb=A.binned(Lt)
    a.plot(Pb[:,0],(Pb[:,1]-Db)/Pb[:,2],'ro-',ms=3,label='Planck − CoE field');a.plot(Pb[:,0],(Pb[:,1]-Lb)/Pb[:,2],'ko-',ms=3,label='Planck − ΛCDM');a.axhline(0,color='grey');a.set_xlabel('ℓ (binned)');a.set_ylabel('residual / σ_Planck')
    np.save(f'resid_{nm}.npy',np.c_[Pb[:,0],(Pb[:,1]-Db)/Pb[:,2],(Pb[:,1]-Lb)/Pb[:,2]])
    a.legend(loc='lower center');a.set_title(f'residuals, binned ℓ≥30 ({nm})');fig.tight_layout();fig.savefig(f'acoustic_resid_{nm}_lin.png',dpi=100);plt.close()
    fig,a=plt.subplots(figsize=(6,5));X=[e['chi2_fit'] for e in S['Lcurve']];Y=[e['reg'] for e in S['Lcurve']];a.loglog(X,Y,'o-')
    for e in S['Lcurve']:a.annotate(f'λ={e["lam"]}',(e['chi2_fit'],e['reg']),fontsize=8)
    a.set_xlabel('χ² (fit measure)');a.set_ylabel('Reg');a.set_title(f'L-curve {nm}');fig.tight_layout();fig.savefig(f'tmp_lc_{nm}.png',dpi=100);plt.close()
json.dump(J,open('acoustic_fit.json','w'),indent=1);print('peaks planck',pl,'lcdm',J['peaks_lcdm'])
for nm in ('free_W','fixed_shell'):
    if nm in J:print(nm,J[nm]['peaks_model'],{k:J[nm]['best'][k] for k in ('p','chi2_fit','chi2_fit_cv','l2_30','l2_2500','l30_2500','sig_delta','lam')})
