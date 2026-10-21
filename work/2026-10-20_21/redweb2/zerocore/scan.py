import json,itertools,numpy as np,os,sys
from multiprocessing import Pool
import zc
def job(p):
    a,b,xm,fr,D=p;xp=xm/(a/b)**(1/b);o=zc.run(a,b,xp,fr,D);o['xm']=xm;o.pop('A',None) if False else None;return o
if __name__=='__main__':
    G=list(itertools.product([0.5,1,2,3],[1,2,4],[0.6,0.8,0.9,1.0,1.1,1.2,1.4],[30,1000,20000],[20,42,60]))
    out=[]
    with Pool(8) as P:
        for o in P.imap_unordered(job,G):
            out.append(o);json.dump(out,open('scan.json','w'))
            if len(out)%25==0:print(len(out),flush=True)
