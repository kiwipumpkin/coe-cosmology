# NEW (rev7): distribution of one-sky chi2 over 1000 realisations of the p8 field (Gaussian a_lm: D_l^sky = Dfix + Df*chi2_{2l+1}/(2l+1)),
# compared with LCDM (mean curve and 1000 LCDM realisations drawn the same way) and with seed 42. Bodies (compensation) not included: they change only l=2 of one sky. Planck errors as in Eq. (6).
import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P='/workspace/coepdf/joint_solve/phase2/';z=np.load(P+'p8_field_consistent.npz');Df,Dfix=z['Df'],z['Dfix']
os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.');import acoustic_fit as A
L=np.arange(len(Df));nu=2*np.maximum(L,1)+1;rg=dict(l2_30=(2,30),l2_2500=(2,2500))
c0={k:A.stats(A.Lt,0*A.Lt,*r)['chi2_planck'] for k,r in rg.items()}
rng=np.random.default_rng(2026);N=1000;res={k:[] for k in rg};resL={k:[] for k in rg}
for i in range(N):
    Dr=Dfix+Df*rng.chisquare(nu)/nu
    for k,r in rg.items():res[k].append(A.stats(Dr,Df,*r)['chi2_planck'])
    Dl=A.Lt*rng.chisquare(nu[:len(A.Lt)])/nu[:len(A.Lt)]
    for k,r in rg.items():resL[k].append(A.stats(Dl,0*Dl,*r)['chi2_planck'])
g=np.random.default_rng(42).chisquare(nu)/nu;s42={k:A.stats(Dfix+Df*g,Df,*r)['chi2_planck'] for k,r in rg.items()}
out={'N':N,'lcdm':c0,'seed42':s42,'lcdm_realisations':{k:dict(median=float(np.median(resL[k])),p16=float(np.percentile(resL[k],16)),p84=float(np.percentile(resL[k],84))) for k in rg}}
for k in rg:
    a=np.array(res[k]);out[k]=dict(median=float(np.median(a)),p16=float(np.percentile(a,16)),p84=float(np.percentile(a,84)),min=float(a.min()),
        frac_below_lcdm=float((a<c0[k]).mean()),frac_below_43=float((a<43).mean()) if k=='l2_30' else None,frac_below_100=float((a<100).mean()) if k=='l2_30' else None)
json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1)
fig,ax=plt.subplots(1,2,figsize=(10,3.6))
for a,k,t in zip(ax,rg,[r'$\ell=2$-30',r'$\ell=2$-2500']):
    v=np.array(res[k]);a.hist(v,bins=50,color='0.6',label='CoE p8 skies');a.hist(resL[k],bins=50,histtype='step',color='k',label=r'$\Lambda$CDM skies');a.axvline(c0[k],c='k',label=r'$\Lambda$CDM mean %.1f'%c0[k]);a.axvline(s42[k],c='r',ls='--',label='seed 42 %.1f'%s42[k])
    a.axvline(np.median(v),c='b',ls=':',label='median %.1f'%np.median(v));a.set_xlabel(r'one-sky $\chi^2$, '+t);a.legend(fontsize=8)
ax[0].set_ylabel('realisations (of 1000)');fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=150);print(json.dumps(out,indent=1))
