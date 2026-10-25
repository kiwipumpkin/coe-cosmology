import json,itertools,numpy as np,os,sys
from multiprocessing import Pool
import zc
def job(p):
    a,b,xm,fr,D=p;xp=xm/(a/b)**(1/b);o=zc.run(a,b,xp,fr,D);o['xm']=xm;o.pop('A',None) if False else None;return o
if __name__=='__main__':
    G=list(itertools.product([0.25,0.5,1],[1,2,4],[0.6,0.8,1.0,1.2],[10,30,100],[0.1,0.2,0.3,0.5,1,2,5]))
    out=[]
    with Pool(5) as P:
        for o in P.imap_unordered(job,G):
            out.append(o);json.dump(out,open('scanB.json','w'))
            if len(out)%25==0:print(len(out),flush=True)
