# rev7 replacement of rev6 Fig 20 (passing domains). rev6 scanned the domain size R=4/8/14 Gpc; in rev7 R=R_ls is fixed by the
# background (p3_background_relations.json), so that scan has no meaning. Recomputed instead: 1000 one-sky realisations of the p8 field
# (same generator and seed as N5: D_l = Dfix + Df chi2_{2l+1}/(2l+1)), l=2-30, with the bands of the skies passing strict (<43) and loose (<100) cuts.
import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P='/workspace/coepdf/joint_solve/phase2/';z=np.load(P+'p8_field_consistent.npz');Df,Dfix=z['Df'],z['Dfix'];H=os.path.dirname(os.path.abspath(__file__))
os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.');import acoustic_fit as A
L=np.arange(len(Df));nu=2*np.maximum(L,1)+1;rng=np.random.default_rng(2026);S=[];C=[]
for i in range(1000):
    Dr=Dfix+Df*rng.chisquare(nu)/nu;S.append(Dr[:31]);C.append(A.stats(Dr,Df,2,30)['chi2_planck']);rng.chisquare(nu[:len(A.Lt)])  # keep N5's draw order (it interleaves LCDM draws)
S=np.array(S);C=np.array(C);k=np.arange(2,31);fig,ax=plt.subplots(1,2,figsize=(14,4.8));out={}
for lab,m,c in (('all 1000',C>-1,'0.7'),('loose χ²<100',C<100,'orange'),('strict χ²<43',C<43,'g')):
    ax[0].fill_between(k,np.percentile(S[m][:,2:],16,0),np.percentile(S[m][:,2:],84,0),color=c,alpha=.35,label=f'{lab}: {m.sum()} skies (16–84%)');out[lab]=int(m.sum())
ax[0].plot(k,A.Lt[2:31],'k-',label='ΛCDM');pf=A.Pf[(A.Pf[:,0]>=2)&(A.Pf[:,0]<=30)];ax[0].errorbar(pf[:,0],pf[:,1],[pf[:,2],pf[:,3]],fmt='ko',ms=3,label='Planck')
ax[0].set_xlabel('ℓ');ax[0].set_ylabel('D_ℓ [μK²]');ax[0].legend(fontsize=8);ax[0].set_title('one-sky bands, R fixed at R_ls=0.18957 Gpc',fontsize=10)
ax[1].hist(C,bins=60,color='0.6');ax[1].axvline(43,c='g',label='strict 43');ax[1].axvline(100,c='orange',label='loose 100');ax[1].set_xlabel('one-sky χ²_Planck, ℓ 2–30');ax[1].legend()
fig.tight_layout();fig.savefig(os.path.join(H,'R20_pass.png'),dpi=120);out['median']=float(np.median(C));json.dump(out,open(os.path.join(H,'R20_pass.json'),'w'));print(out)
