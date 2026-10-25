import sys,json,pickle,os,numpy as np,multiprocessing as mp
import szk, toy
RV=14.30048287260589
V={r['i']:r for r in json.load(open('../valley/valley_cmb.json'))}; P=json.load(open('../valley/valley_points.json'))
def job(a):
    i,C=a; fn=f'cmb_{i}_{C:+.2f}.pkl'
    if os.path.exists(fn): return a,'cached'
    p=P[i]; fac,Dt=toy.make(p); L=fac.func(*fac.args[:3],V[i]['rmin'])
    rin=(p['D']-p['kin']*RV)/RV; rout=(p['D']+p['kout']*RV)/RV; S=szk.SZ(L,C,rin,min(rout,L.rg[-1]))
    T=S.t1+1.3; ems,y=szk.ring(S,0.00123,T,Dt); x=szk.observe(S,y,ems,T)
    pickle.dump(dict(i=i,C=C,x=x,sig0=S.sig0,ratio=S.ratio),open(fn,'wb')); return a,int(np.isnan(x).sum())
if __name__=='__main__':
    jobs=[tuple((int(s.split(':')[0]),float(s.split(':')[1]))) for s in sys.argv[1:]]
    with mp.Pool(len(jobs)) as Pp:
        for a,r in Pp.imap_unordered(job,jobs): print(a,r,flush=True)
