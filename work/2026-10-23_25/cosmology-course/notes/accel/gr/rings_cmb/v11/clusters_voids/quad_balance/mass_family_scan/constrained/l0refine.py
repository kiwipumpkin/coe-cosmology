import os,json
from l0scan import work,R0
from multiprocessing import Pool
if __name__=='__main__':
    out='l0scan_results.json';D=json.load(open(out));T={}
    for d in D:T.setdefault(d['run'],{})[d['l0']]=d['chi2']
    jobs=[]
    for k,v in T.items():
        m=min(v,key=v.get)
        for dl in (-1.5,-1,-0.5,0.5,1,1.5):
            l=round(m+dl,1)
            if -10<=l<=100 and l not in v: jobs.append((k,l))
    print('jobs',len(jobs),flush=True)
    with Pool(7) as p:
        for d in p.imap_unordered(work,jobs):
            D.append(d);json.dump(D,open(out,'w'))
    print('done',len(D))
