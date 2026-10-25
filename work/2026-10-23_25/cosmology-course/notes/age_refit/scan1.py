"""Stage 1 of v4 scan: Migkas chi2 (71 bins, exact-GR light cone) on Sobol samples. usage: scan1.py law mode vm N seed"""
import os,sys,time,json,numpy as np,warnings;warnings.filterwarnings('ignore')
law,mode,vm,N,sd=sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),int(sys.argv[5]); os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
sys.argv=['x',law,mode]; exec(open('fit_lc.py').read().split("ck=f'fwd")[0]); import grltb as G
from scipy.stats import qmc
box=BX if mode=='back' else [B[n] for n in nm]
x0=np.array(json.load(open(f'gr/ck_{law}_{mode}_{vm}.json'))['x'])
# half the samples over the whole sampling box, half in a +-1.5 (log units; s: +-1) neighbourhood of the Migkas best fit
U=qmc.Sobol(len(nm),seed=sd).random(N)
lo=np.array([b[0] for b in box]); hi=np.array([b[1] for b in box]); w=np.array([1.0 if n=='s' else 1.5 for n in nm])
Xs=np.r_[[x0],lo+U[:N//2]*(hi-lo),x0-w+U[N//2:]*2*w]
res=[]
for x in Xs:
    try:
        P,o=build(x)
        if o is None: continue
        f=G.fields(o); tr=G.trace(f,np.arange(0,360,4.),nstep=320,lbmax_frac=0.999); af=G.migkas(tr)
        if not np.all(np.isfinite(af)): continue
        c,l0=F.chi2_l0(af); res.append(dict(x=list(x),chi2=float(c),l0=float(l0),M=float(P.get('M',np.nan))))
        if len(res)%5==0: json.dump(dict(names=nm,box=box,law=law,mode=mode,vm=vm,res=res,partial=True),open(f'gr2/v4/scan1_{law}_{mode}_{vm}.json','w'))
    except (Exception,MemoryError): pass
json.dump(dict(names=nm,box=box,law=law,mode=mode,vm=vm,res=res),open(f'gr2/v4/scan1_{law}_{mode}_{vm}.json','w'))
print('done',law,mode,vm,len(res))
