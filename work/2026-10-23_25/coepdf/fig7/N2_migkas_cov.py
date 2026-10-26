# NEW (rev7): approximate covariance-aware chi2 for the 71 overlapping Migkas cones.
# Assumption (stated in paper): clusters roughly uniform in longitude, so the correlation of two
# 90-deg cones is their overlap fraction rho=max(0,1-|dl|/90) (circular). Not the bootstrap
# covariance from the 313-cluster catalogue (not on the box).
import json,numpy as np,importlib.util as iu,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.stats import chi2 as C2
sp=iu.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=iu.module_from_spec(sp);sp.loader.exec_module(cw)
P='/workspace/coepdf/joint_solve/phase2/p1_migkas_D34.96.json'
p1=json.load(open(P))['34.96'];m_coe=np.array(p1['curve'])[cw.IMIG]
y=cw.Y;s=cw.SIG;l=cw.LDEG;n=len(y)
dl=np.abs(l[:,None]-l[None,:]);dl=np.minimum(dl,360-dl)
rho=np.clip(1-dl/90.,0,None)
Cv=rho*np.outer(s,s)
# overlap matrix is near-singular (bins ~5 deg apart, cones 90 deg wide): keep eigenmodes with
# eigenvalue > 1% of the largest (pseudo-inverse); number kept is reported as n_modes.
ev,U=np.linalg.eigh(Cv);keep=ev>1e-2*ev.max();NMODES=int(keep.sum())
Ci=(U[:,keep]/ev[keep])@U[:,keep].T
def X2(r,cov=True):return float(r@Ci@r) if cov else float(np.sum(r**2/s**2))
def gls(cols):
    X=np.column_stack(cols);b=np.linalg.lstsq(X.T@Ci@X,X.T@Ci@y,rcond=None)[0];return X@b
L=np.radians(l);one=np.ones(n)
mods={'isotropy (k=0)':(one,0),'free dipole (k=3)':(gls([one,np.cos(L),np.sin(L)]),3),
 'dipole+quadrupole (k=5)':(gls([one,np.cos(L),np.sin(L),np.cos(2*L),np.sin(2*L)]),5),
 'CoE, D=34.96 (k=4: D,eta,m_c,l0)':(m_coe,4)}
# effective number of independent points: n_eff = n^2/sum(rho^2)
neff=float(n**2/np.sum(rho**2))
out={'assumption':'rho=max(0,1-|dl|/90) overlap model, pseudo-inverse eig>1e-2 max','n_modes':NMODES,'n_bins':n,'n_eff':neff,'models':{}}
for k,(m,p) in mods.items():
    r=y-m;c=X2(r);ci=X2(r,False)
    out['models'][k]=dict(k=p,chi2_indep=ci,chi2_cov=c,AIC=c+2*p,BIC=c+p*np.log(neff))
c0=out['models']['isotropy (k=0)']['chi2_cov']
for k,v in out['models'].items():
    if v['k']>0:
        d=c0-v['chi2_cov'];v['dchi2_vs_iso']=d;v['p_vs_iso']=float(C2.sf(d,v['k']))
json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1)
fig,ax=plt.subplots(1,2,figsize=(11,3.8),gridspec_kw=dict(width_ratios=[1.6,1]))
ax[0].errorbar(l,y,s,fmt='o',ms=3,color='k',lw=.7,label='Migkas $A/A_{\\rm all}$ (71 cones)')
for (k,(m,p)),c in zip(mods.items(),['grey','C0','C2','C3']):ax[0].plot(l,m,color=c,lw=1.5,label=k)
ax[0].set_xlabel('Galactic longitude [deg]');ax[0].set_ylabel('$A/A_{\\rm all}$');ax[0].legend(fontsize=7)
ks=list(mods);x=np.arange(len(ks))
ax[1].bar(x-.2,[out['models'][k]['chi2_indep'] for k in ks],.4,label='independent-error $\\chi^2$')
ax[1].bar(x+.2,[out['models'][k]['chi2_cov'] for k in ks],.4,label='overlap-covariance $\\chi^2$')
ax[1].set_xticks(x);ax[1].set_xticklabels(['iso','dip','dip+quad','CoE'],fontsize=8);ax[1].legend(fontsize=7)
ax[1].set_title('$n_{\\rm eff}=%.1f$'%neff,fontsize=9)
plt.tight_layout();plt.savefig(__file__.replace('.py','.png'),dpi=150)
print(json.dumps(out,indent=1))
