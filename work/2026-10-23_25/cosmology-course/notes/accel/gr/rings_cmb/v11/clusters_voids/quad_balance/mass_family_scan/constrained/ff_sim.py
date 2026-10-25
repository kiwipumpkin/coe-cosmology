import os,sys,json,subprocess,numpy as np,warnings;warnings.filterwarnings('ignore')
from ff_prof import fit
H=os.getcwd();sys.argv=[sys.argv[0]];import anP;os.chdir(H)
from multiprocessing import Pool
from scipy.stats import qmc
P=np.loadtxt('../../../reproduce_cmb/planckTT.txt');PL={int(r[0]):r[1:] for r in P}
def pchi(Dl,a=2,b=30):
    s=0.
    for l in range(a,b+1):
        d,lo,hi=PL[l];e=hi if Dl[l]>d else lo;s+=((Dl[l]-d)/e)**2
    return float(s)
def work(a):
    i,C,fr,dr=a;fn=H+'/ff_runs/s%03d.npz'%i
    try:
        if not os.path.exists(fn):
            e=dict(os.environ,DRV=str(dr),**{'C%d'%(k+1):str(C[k]) for k in range(4)})
            r=subprocess.run(['python3','famP.py','poly','0',str(fr),fn],env=e,capture_output=True,text=True,timeout=3000)
            if not os.path.exists(fn):return dict(i=i,C=C,frac=fr,DRV=dr,error=r.stderr[-300:])
        A=anP.an(fn);os.chdir(H);rat=float(A['ratio']);Dl=[float(x) for x in A['Dl'][:31]]
        m1=fit((rat,-10,100));m2=fit((rat,0,360))
        return dict(i=i,C=C,frac=fr,DRV=dr,ratio=rat,mig=m1,mig360=m2,Dl=Dl,cmb=pchi(Dl),cmb26=pchi(Dl,2,6))
    except Exception as ex:return dict(i=i,C=C,frac=fr,DRV=dr,error=repr(ex))
if __name__=='__main__':
    N=int(os.environ.get('N','48'));s=qmc.LatinHypercube(d=6,seed=7).random(N)
    lo=np.array([-0.3,-0.3,-0.3,-0.3,1,1.5]);hi=np.array([0.3,0.3,0.3,0.3,np.log10(50000),6]);X=lo+(hi-lo)*s
    jobs=[tuple(j) for j in json.loads(os.environ['JOBS'])] if 'JOBS' in os.environ else [(i,[round(v,5) for v in x[:4]],round(10**x[4],3),round(x[5],4)) for i,x in enumerate(X)]
    out='fullfree_sims.json';R=json.load(open(out)) if os.path.exists(out) else [];done={d['i'] for d in R}
    with Pool(3) as p:
        for d in p.imap_unordered(work,[j for j in jobs if j[0] not in done]):
            R.append(d);json.dump(R,open(out,'w'));print(d['i'],d.get('ratio'),d.get('mig',{}).get('chi2'),d.get('cmb'),d.get('error','')[:80],flush=True)
