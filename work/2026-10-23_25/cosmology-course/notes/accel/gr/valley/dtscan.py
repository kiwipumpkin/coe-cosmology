import sys,json,pickle,numpy as np,multiprocessing as mp
import toy, run_exactc as X
def job(dt):
    p=json.load(open('valley_points.json'))[2]; fac,Dt=toy.make(p); X.LTB=lambda core_frac=None,rmax=150.0,nr=900: fac(rmax=rmax,nr=nr)
    X.D=Dt; X.C0=np.array([Dt,0.0]); L=X.LTB(); x=X.run(None,0.00123,L.t1+dt); return dt,x
if __name__=='__main__':
    with mp.Pool(4) as P:
        out={}
        for dt,x in P.imap_unordered(job,[1.3,1.8,2.5,3.5]):
            out[dt]=x; pickle.dump(out,open('dtscan.pkl','wb')); n=np.where(np.isnan(x))[0]; print(dt,len(n),n[:3],n[-3:],flush=True)
