# Recombination when each comoving shell reaches fixed areal radius R_rec from the CoE. Back-trace photons from observer (today, t1).
import sys,os,pickle,warnings;warnings.filterwarnings('ignore')
G='/workspace/cosmology-course/notes/accel/gr';sys.path[:0]=[G+'/valley',G+'/ltb',G]
_h=os.path.dirname(os.path.abspath(__file__));os.chdir(G+'/valley');sys.path.insert(0,_h);import toyext as toy;os.chdir(_h)
import numpy as np;from ltbcore import photon
F=0.0178;FRAC=float(os.environ.get('RFRAC','0.3'))
def job(a):
    uD,f=a;Dt=uD;Mc=max(f,1e-12)*Dt**3/2;Mfun=lambda r:Mc+0*r;rm=0.02
    while True:
        try:L=toy.LTBM(Mfun,Dt,Dt,rm,rmax=max(150.,3*Dt),nr=900);break
        except ValueError: rm*=1.5
        if rm>0.5*Dt: return a,'no grid'
    Rrec=FRAC*Dt;T=L.t1;v0=0.00123;vo=L.frame(T,Dt,0.0,[v0,0.0]);out=[];ze=[]
    for th in np.radians(np.arange(360)):
        so=photon(L,Dt,0.0,th+np.pi,T,L.tg[0]*1.01)
        ts=np.linspace(T,so.t[-1],4000);y=so.sol(ts);Rv=L.Rs.ev(y[0],ts)-Rrec
        k=np.where(Rv<=0)[0]
        if len(k)==0: out.append(np.nan);ze.append(np.nan);continue
        k=k[0];from scipy.optimize import brentq
        te=brentq(lambda t:L.Rs.ev(so.sol(t)[0],t)-Rrec,ts[k],ts[k-1]);ye=so.sol(te)
        R,Rp,Rt,Gg,_,_=L.fields(T,Dt);vph=[np.cos(th+np.pi)/np.sqrt(Gg),np.sin(th+np.pi)/R]
        opz=L.kdu(T,Dt,vph,0.0,vo)/L.kdu(te,ye[0],ye[2:4],ye[4],[0.0,0.0])
        out.append(1/opz);ze.append(te)
    return a,dict(x=np.array(out),te=np.array(ze),Rrec=Rrec,rmin=rm,t1=T)
if __name__=='__main__':
    import multiprocessing as mp
    jobs=[(u,f) for u in (1.05,3.0,10.5) for f in (F,0.0)];out={}
    with mp.Pool(6) as P:
        for a,r in P.imap_unordered(job,jobs):
            out[a]=r;pickle.dump(out,open(f'recomb_{FRAC}.pkl','wb'))
            print(a,r if isinstance(r,str) else (r['rmin'],int(np.isnan(r['x']).sum()),np.nanstd(r['x']/np.nanmean(r['x']))),flush=True)
