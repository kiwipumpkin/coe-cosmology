# Regenerated (rev7): TT spectrum of the current solve (p8 field, D=34.96, R=0.18957) vs Planck PR3 and LCDM.
# Input: p6_<tag>.npz from joint_solve/phase2/p6_tt_compensated.py (default tag p8).
import sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P='/workspace/coepdf/joint_solve/phase2/';tag=sys.argv[1] if len(sys.argv)>1 else 'p8'
z=np.load(P+f'p6_{tag}.npz');J=json.load(open(P+f'p6_{tag}.json'))
Pf=z['Pf'];l=Pf[:,0];d=Pf[:,1];e=0.5*(Pf[:,2]+Pf[:,3]);L=np.arange(len(z['D']))
fig,ax=plt.subplots(2,1,figsize=(10,6),sharex=True,gridspec_kw=dict(height_ratios=[2.2,1]))
ax[0].errorbar(l,d,e,fmt='.',ms=1.5,color='0.6',lw=.3,label='Planck PR3 TT',zorder=1)
ax[0].plot(L,z['Lt'],'k',lw=1,label='$\\Lambda$CDM best fit',zorder=5)
ax[0].plot(L,z['D'],'C3',lw=1,label='toy model, ensemble mean',zorder=4)
ax[0].plot(L,z['Drc'],'C0',lw=.6,alpha=.7,label='toy model, seed-42 compensated sky',zorder=3)
ax[0].set_xscale('log');ax[0].set_ylabel('$D_\\ell$ [$\\mu$K$^2$]');ax[0].legend(fontsize=8);ax[0].set_xlim(2,2500)
li=l.astype(int)
for k,c in (('D','C3'),('Lt','k')):ax[1].plot(l,(d-z[k][li])/e,'.',ms=1,color=c,alpha=.5)
ax[1].axhline(0,color='0.3',lw=.5);ax[1].set_ylabel('(data$-$model)/$\\sigma$');ax[1].set_xlabel('$\\ell$');ax[1].set_ylim(-6,6)
m=J['mean']['l2_2500']['chi2_planck'];s=J['seed42_comp']['l2_2500']['chi2_planck']
ax[0].set_title('$\\chi^2_{2-2500}$: mean %.0f, seed-42 compensated %.0f'%(m,s),fontsize=9)
plt.tight_layout();plt.savefig(__file__.replace('.py','.png'),dpi=150)
