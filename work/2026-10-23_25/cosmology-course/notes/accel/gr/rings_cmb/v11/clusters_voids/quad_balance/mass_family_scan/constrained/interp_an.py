import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.getcwd();RP=H+'/runs_poly/';REF=os.path.abspath('../runs_sl2/slc_0.041037_0.035.npz')
os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
exec(open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read().split('rng=np.random')[0]);os.chdir(H)
sys.argv=[sys.argv[0]];import anP
from multiprocessing import Pool
def mfit(rat,n=40):
    rng=np.random.default_rng(12345);best=None
    for k in range(n):
        r=minimize(chi,LO+(HI-LO)*rng.random(6),args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    th=best.x
    for j in range(3):th=minimize(chi,th,args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-9,fatol=1e-9)).x
    m=model(th,rat);return th,float(chi(th,rat)),float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2))
Ls=np.arange(0,361,1.0);Lsave=L.copy()
def curve(t,ra):
    global L;L=np.mod(Ls,360);L[-1]=360.0;v=model(t,ra);L=Lsave;return v
fmt=lambda v:('%g'%v)
C3=[0,0.02,0.04,0.06,0.08,0.1,0.12,0.1353];FR=[30,100,300,1000,3000,10000,20000]
runs={}
for c in C3:
    for f in FR:runs['grid_C3_%g_f%d'%(c,f)]=RP+'f%d_p_-0.041_-0.035_%s_0.npz'%(f,fmt(c))
path=[('path00_REF',REF)]
for k in range(1,10):
    s=k/10;path.append(('path%02d'%k,RP+'f%s_p_-0.041_-0.035_%s_0.npz'%(str(round(30*(20000/30)**s,3)),fmt(round(0.1353*s,6)))))
path.append(('path10_RED',RP+'f20000_p_-0.041_-0.035_0.1353_0.npz'))
for k,v in path:runs[k]=v
runs['poly_C3_0_f30']=RP+'p_-0.041_-0.035_0_0.npz'
for x in ['1.0','1.05','1.13','1.3']:
    for f in [20000,30]:runs['hyb_xs%s_f%d'%(x,f)]=RP+'hyb_xs%s_f%d.npz'%(x,f)
def work(kv):
    k,fn=kv
    if not os.path.exists(fn):return k,dict(run=fn,missing=True)
    try:
        A=anP.an(fn);th,c2,R2=mfit(A['ratio']);d=np.load(fn)
        A.update(frac=float(d['frac']),chi2=c2,R2=R2,l0=float(th[1]),l0_at_bound=bool(min(th[1]-LO[1],HI[1]-th[1])<0.5),l0_bounds=[float(LO[1]),float(HI[1])],theta=th.tolist(),
          shells_crossed=None);A['Dl']=A['Dl'][:31];A.pop('chi2_interp',None);return k,A
    except Exception as e:return k,dict(run=fn,error=repr(e))
if __name__=='__main__':
    with Pool(int(os.environ.get('NP','4'))) as p:R=dict(p.map(work,list(runs.items())))
    R['_meta']=dict(path='step k: C3=0.1353*k/10, frac=30*(20000/30)^(k/10), C1=-0.041,C2=-0.035; step0=REF slc 0.041037/0.035 frac30; step10=RED',
      hybrid='ln rho = (1-w) RED + w (REF + offset), w=0.5(1+tanh((x-xs)/0.01)), offset matches at xs',migkas='71 bins, 6 params, l0 bounds [-60,100], 40 NM starts seed 12345')
    json.dump(R,open('interp_results.json','w'),indent=1);print('saved')
