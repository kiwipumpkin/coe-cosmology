# rev7 regeneration of rev6 Fig 17 (stage-A low-l fit). rev6 stage A was a regularised (lambda=100) GRF fit at R=4 Gpc. In rev7 the field
# is the unregularised multistart p8 fit at R_ls=0.18957 Gpc (joint_solve/phase2/p8_field_consistent.npz), so there is no lambda / L-curve;
# shown: ensemble mean Dfix+Df with ±1σ cosmic-variance band, the deterministic part Dfix, LCDM, Planck 2018, l=2-30. Stats from acoustic_fit.stats.
import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P='/workspace/coepdf/joint_solve/phase2/';z=np.load(P+'p8_field_consistent.npz');Df,Dfix=z['Df'],z['Dfix'];H=os.path.dirname(os.path.abspath(__file__))
os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.');import acoustic_fit as A
L=np.arange(len(Df));nu=2*np.maximum(L,1)+1;mean=Dfix+Df;cv=Df*np.sqrt(2/nu)
s=A.stats(mean,Df,2,30);sl=A.stats(A.Lt,0*A.Lt,2,30);k=slice(2,31)
fig,ax=plt.subplots(figsize=(10,4.8))
ax.fill_between(L[k],(mean-cv)[k],(mean+cv)[k],color='b',alpha=.2,label='±1σ cosmic variance')
ax.plot(L[k],mean[k],'b-',lw=2,label=f'CoE p8 ensemble mean: χ²_Planck+CV={s["chi2_cv"]:.1f}, χ²_Planck={s["chi2_planck"]:.1f}')
ax.plot(L[k],Dfix[k],'r--',label='deterministic part D_fix (background + bodies)');ax.plot(A.Lt[:31][k] if False else L[k],A.Lt[k],'k-',label=f'ΛCDM: χ²_Planck={sl["chi2_planck"]:.1f}')
pf=A.Pf[(A.Pf[:,0]>=2)&(A.Pf[:,0]<=30)];ax.errorbar(pf[:,0],pf[:,1],[pf[:,2],pf[:,3]],fmt='ko',ms=3,label='Planck 2018')
ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ [μK²]');ax.set_xlim(1.5,30.5);ax.legend(fontsize=8);ax.set_title('rev7 low-ℓ fit (R_ls=0.18957 Gpc, D=34.96 Gpc)',fontsize=10)
fig.tight_layout();fig.savefig(os.path.join(H,'R17_stageA.png'),dpi=140);json.dump(dict(coe=s,lcdm=sl),open(os.path.join(H,'R17_stageA.json'),'w'),default=float,indent=1);print(s,sl)
