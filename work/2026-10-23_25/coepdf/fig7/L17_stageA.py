# PROVENANCE COPY of /workspace/redweb2/stochfield/figs.py (original generating script, rev6 figure; run in its original directory).
# Output figure copied here as L17_stageA.png. Not re-run for rev7: see Appendix B.
import sys,json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from fwd import Fwd
tag=sys.argv[1];LAMI=int(sys.argv[2]) if len(sys.argv)>2 else None
R=json.load(open(f'fit_{tag}.json'));Pf=np.array(R['planck']);l=Pf[:,0];Lt=np.array(R['lcdmD']);Dfix=np.array(R['Dfix']);LC=R['Lcurve']
# chosen: L-curve corner = max curvature in (log chi2, log reg); fallback lam=0 entry
x=np.log10([max(e['chi2_cv'],1e-9) for e in LC[1:]]);y=np.log10([max(e['reg'],1e-12) for e in LC[1:]])
if LAMI is None:
    dx,dy=np.gradient(x),np.gradient(y);ddx,ddy=np.gradient(dx),np.gradient(dy);k=(dx*ddy-dy*ddx)/np.maximum((dx**2+dy**2)**1.5,1e-12);LAMI=1+int(np.argmax(np.abs(k[1:-1])))+1
S=LC[LAMI];p=S['p'];F=Fwd(lmax=int(R['LMAX']));t=F.terms(p);Df=t['scalar']+t['mlens'];Dm=Df+Dfix;cv=np.sqrt(2/(2*l+1))*Df
lab=lambda c2,c2p,r2:f'χ²(Planck+CV)={c2:.1f}, χ²(Planck err)={c2p:.1f}, Pearson r²={r2:.3f}'
fig,ax=plt.subplots(figsize=(11,5));ax.plot(l,Lt,'k-',lw=2,label=f'ΛCDM  χ²(Planck err)={R["lcdm"]["chi2_planck"]:.1f}, r²={R["lcdm"]["r2"]:.3f}')
ax.fill_between(l,Dm-cv,Dm+cv,color='C0',alpha=.3,label='ensemble ±1σ cosmic variance');ax.plot(l,Dm,'C0-',lw=2,label='stochastic field ensemble mean: '+lab(S['chi2_cv'],S['chi2_planck'],S['r2']))
if R['BG']:ax.plot(l,Dfix,'r--',lw=1.2,label='deterministic CoE background (smooth red + RS_cc)')
ax.errorbar(l,Pf[:,1],yerr=[Pf[:,2],Pf[:,3]],fmt='o',color='k',ms=3,label='Planck 2018')
ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)');ax.set_ylim(0,max(2500,1.1*Pf[:,1].max()) if R['LMAX']>30 else 2500);ax.set_title(f'(a) stochastic field at R=4 Gpc, λ={S["lam"]:.3g}  [{tag}]');ax.legend(loc='lower center',fontsize=7,framealpha=.85)
fig.tight_layout();fig.savefig(f'a_cl_{tag}.png',dpi=100);plt.close()
# (b) P(k) and W
fig,(a1,a2)=plt.subplots(1,2,figsize=(12,4.5));kk=np.geomspace(1e-2,100,400)
for j,e in enumerate(LC[1::2]):a1.loglog(kk,F.Pd(e['p'],kk),color=plt.cm.viridis(j/6),lw=1,label=f'λ={e["lam"]:.2g}')
a1.loglog(kk,F.Pd(p,kk),'k',lw=2.5,label=f'chosen λ={S["lam"]:.2g}');a1.set_xlabel('k (1/Gpc)');a1.set_ylabel('P_δ(k) (Gpc³)');a1.legend(loc='lower center',fontsize=7,ncol=2)
a1.set_title(f'(b) P_δ: n={p["n"]:.2f}, k_c={p["kc"]:.3g}/Gpc, σ_δ(0.078)={S["reg_sig_delta"]:.3g}')
for j,e in enumerate(LC[1::2]):a2.plot(F.ch,F.W(e['p'])[0],color=plt.cm.viridis(j/6),lw=1)
a2.plot(F.ch,F.W(p)[0],'k',lw=2.5);a2b=a2.twinx();a2b.semilogy(F.ch,abs(F.KM),'C3:',label='|moving-lens kernel K_M|');a2b.legend(loc='lower center',fontsize=7)
a2.set_xlabel('distance from observer χ (Gpc)');a2.set_ylabel('W(χ)');a2.set_title(f'radial window: centre {p["chiw"]:.2f}, width {p["sw"]:.3g} Gpc (R=4)')
fig.tight_layout();fig.savefig(f'b_pk_w_{tag}.png',dpi=100);plt.close()
# (c) L-curve
fig,ax=plt.subplots(figsize=(6,5));X=[e['chi2_cv'] for e in LC[1:]];Y=[e['reg'] for e in LC[1:]];ax.loglog(X,Y,'o-')
for e in LC[1:]:ax.annotate(f'{e["lam"]:.2g}',(e['chi2_cv'],e['reg']),fontsize=7)
ax.loglog([S['chi2_cv']],[S['reg']],'r*',ms=14,label='chosen');ax.set_xlabel('χ² (Planck+CV)');ax.set_ylabel('regularizer');ax.set_title(f'(c) L-curve [{tag}]');ax.legend(loc='lower center')
fig.tight_layout();fig.savefig(f'c_lcurve_{tag}.png',dpi=100);plt.close()
json.dump(dict(chosen=S['lam'],p=p,chi2_cv=S['chi2_cv'],chi2_planck=S['chi2_planck'],r2=S['r2'],reg=S['reg'],sig_delta=S['reg_sig_delta']),open(f'chosen_{tag}.json','w'),indent=1)
print('chosen',S['lam'],p,S['chi2_cv'],S['chi2_planck'],S['r2'])
