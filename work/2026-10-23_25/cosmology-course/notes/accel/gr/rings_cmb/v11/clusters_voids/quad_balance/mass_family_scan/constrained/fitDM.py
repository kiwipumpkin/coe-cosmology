import os,sys,json,subprocess,numpy as np,warnings;warnings.filterwarnings('ignore')
H=os.getcwd();os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
exec(open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read().split('rng=np.random')[0]);os.chdir(H)
JA=sys.argv[1] if len(sys.argv)>1 else 'fitDM_jobs.json';sys.argv=[sys.argv[0]];import anP
from multiprocessing import Pool
P=np.loadtxt('../../../reproduce_cmb/planckTT.txt');PL={int(r[0]):r[1:] for r in P}
def mfit(rat,fix=None,n=40):
    rng=np.random.default_rng(12345);idx=[0,2,3,4,5];best=None
    f=(lambda p:chi(p,rat)) if fix is None else (lambda p:chi(np.insert(p,1,fix),rat))
    for k in range(n):
        th=LO+(HI-LO)*rng.random(6);th=th if fix is None else th[idx]
        r=minimize(f,th,method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    p=best.x
    for j in range(3):p=minimize(f,p,method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-9,fatol=1e-9)).x
    th=p if fix is None else np.insert(p,1,fix);m=model(th,rat)
    return float(chi(th,rat)),float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),th.tolist()
def pchi(Dl):
    s=0.
    for l in range(2,31):
        d,lo,hi=PL[l];e=hi if Dl[l]>d else lo;s+=((Dl[l]-d)/e)**2
    return float(s)
SH={'red':('poly',dict(C1='-0.041',C2='-0.035',C3='0.1353',C4='0')),'ref':('slc',dict(SL='0.041037'))}
def work(a):
    sh,dr,fr=a;fam,env=SH[sh];fn=H+'/fitDM_runs/%s_D%g_f%g.npz'%(sh,dr,fr)
    try:
        if not os.path.exists(fn):
            e=dict(os.environ,DRV=str(dr),**env);r=subprocess.run(['python3','famP.py',fam,'0.035' if fam=='slc' else '0',str(fr),fn],env=e,capture_output=True,text=True,timeout=3600)
            if not os.path.exists(fn):return dict(shape=sh,DRV=dr,frac=fr,error=r.stderr[-400:])
        A=anP.an(fn);c,R2,th=mfit(A['ratio']);c15,R15,th15=mfit(A['ratio'],15.6);Dl=A['Dl'][:31]
        return dict(shape=sh,DRV=dr,frac=fr,ratio=A['ratio'],a=A['a'],chi2=c,R2=R2,l0=th[1],theta=th,chi2_l015p6=c15,R2_l015p6=R15,theta_l015p6=th15,Dl=Dl,pschi2=pchi(Dl))
    except Exception as ex:return dict(shape=sh,DRV=dr,frac=fr,error=repr(ex))
if __name__=='__main__':
    out='fitDM_results.json';R=json.load(open(out)) if os.path.exists(out) else []
    done={(d['shape'],d['DRV'],d['frac']) for d in R}
    jobs=[tuple(j) for j in json.loads(open(JA).read()) if tuple(j) not in done]
    with Pool(int(os.environ.get('NP','7'))) as p:
        for d in p.imap_unordered(work,jobs):
            R.append(d);json.dump(R,open(out,'w'));print(d['shape'],d['DRV'],d['frac'],d.get('chi2'),d.get('error','')[:100],flush=True)
