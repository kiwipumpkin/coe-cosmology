import json,sys,os,numpy as np
from scipy.stats import chi2 as C, norm
N='/workspace/cosmology-course/notes/'
src=open(N+'catwise-density-levelset-migkas.py').read()
ns={'__name__':'x'}; exec(src[:src.index('def main')],ns)
pc=ns['pc'];IMIG=ns['IMIG'];Y=ns['MA'];W=ns['W'];shifted=ns['shifted'];TH1=ns['TH1'];LDEG=ns['LDEG']
null=float(np.sum(W*(1-Y)**2))
F=json.load(open(N+'catwise-density-levelset-fit.json'))['migkas_fits']
E=json.load(open(N+'catwise-eta-migkas-fit.json'))['migkas_fits']['d14.3']; E['kind']='eta';E['shape']={'eta':E['eta']}
rows=[('Exponential, 3% lock, ℓ0=15.2',F['exp_fixed15.2'],2,0.716),
      ('Soft core 1/(1+(r/L)^0.3), lock, ℓ0=15.2',F['soft_fixed15.2'],3,0.743),
      ('η-law (η=0.25), lock, ℓ0=15.2',F['eta_fixed15.2'],3,None),
      ('Stretched exp η=0.05, lock, ℓ0 free',E,4,0.745)]
def r2(y,m): return 1-np.sum((y-m)**2)/np.sum((y-y.mean())**2)
out=[];curves={}
for name,p,nf,paper in rows:
    a=ns['_pc_args'](p['kind'],p['DR'],p['LD'],p.get('shape'))
    rr=pc.run(p['kind'],a,p['mc'],0.0); af=rr['af']/rr['af'].mean()
    for res in (6,1):
        if res==1: g=af
        else:
            t=TH1[::6]; g=np.interp(TH1,np.append(t,360),np.append(af[::6],af[0]))
        m=shifted(g,p['l0'])[IMIG]; c2=float(np.sum(W*(m-Y)**2)); d=null-c2
        s=float(norm.isf(C.sf(d,nf)))
        out.append(dict(model=name,res=res,R2=float(r2(Y,m)),chi2=c2,sigma=s,l0=p['l0'],D=p['D_Gpc'],
            stored_chi2=p['chi2'],stored_R2=p['R2_unweighted'],paper_R2=paper))
        if res==1: curves[name]=shifted(g,p['l0'])
json.dump(dict(null=null,rows=out),open('repro.json','w'),indent=1)
for o in out: print(o)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
MLO,MHI=ns['MLO'],ns['MHI']
fig,ax=plt.subplots(figsize=(14,8.5))
ax.errorbar(LDEG,Y,yerr=[MLO,MHI],fmt='o',color='k',ms=6,capsize=2,alpha=.7,label='Migkas A/A_all (71 bins)')
for (name,c),ls in zip(curves.items(),['--','-',':','-.']):
    o=[x for x in out if x['model']==name and x['res']==1][0]
    ax.plot(TH1,c,ls,lw=3.5,label=f"{name}: R²={o['R2']:.3f}, {o['sigma']:.1f}σ")
ax.axhline(1,color='gray',lw=1)
ax.set_xlim(0,360);ax.set_xlabel('Galactic longitude ℓ (deg)',fontsize=22);ax.set_ylabel('A / A_all',fontsize=22)
ax.tick_params(labelsize=18);ax.grid(alpha=.3);ax.legend(fontsize=14,loc='upper right')
ax.set_title('Paper geometry reproduced vs Migkas (scored at 1°)',fontsize=22)
fig.tight_layout();fig.savefig('paper_repro_migkas.png',dpi=130)
