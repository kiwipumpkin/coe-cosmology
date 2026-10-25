import os,sys,json,numpy as np
out={}
for LAW in ('exp','soft','eta','pl','exp_lock','soft_lock','eta_lock','pl_lock'):
    f=f'fwd/before/s12_{LAW}.json'
    if not os.path.exists(f): continue
    sys.argv=['x',LAW]; os.environ['RING']='low'; ns={}
    exec(open('unified2.py').read().split("if __name__=='__main__'")[0],ns)
    P=json.load(open(f))['params']; Q=dict(P); Q.pop('qc',None)
    e=ns['model'](dict(Q),False); e.pop('A',None)
    r={k:e[k] for k in ('cancel','rms','age','Hmean','vo','chi2')}
    o=json.load(open(f))['cmb_H70']; r['before_cancel']=o['cancel']; r['before_rms']=o['rms']; out[LAW]=r; print(LAW,r,flush=True)
json.dump(out,open('fwd/step2_redo.json','w'),indent=1,default=float)
