import sys,os,json,numpy as np
LAW=sys.argv[1]; sys.argv=['x',LAW]; os.environ['RING']='low'
src=open('unified2.py').read().split("if __name__=='__main__'")[0]; exec(src)
from scipy.optimize import differential_evolution, minimize
nm=['D','M','l0']+([] if LOCK else ['L'])+SH+['p']   # H fixed at 70 in step 1: Migkas sees only mc=2GM/(H^2D^3)
def up(z):
    P=dict(zip(nm,z)); P['l0']%=360
    for k in P:
        if k not in('l0','p'): P[k]=float(np.exp(P[k]))
    P['H']=70.0; P['qc']=0.3; return P
def c1(z):
    try: e=model(up(z),True)
    except Exception: e=None
    return 1e6 if e is None else e['chi2']
BX={'D':(np.log(1),np.log(2000)),'M':(80,140),'l0':(0,360),'L':(-6,8),'alpha':(-3,3),'eta':(-3,3),'n':(-3,6),'p':(-5,5)}
B=[BX[n] for n in nm]
x0=np.array([np.log(40.69),np.log(7.6e53),15.3]+([] if LOCK else [np.log(1.235)])+[0.0]*len(SH)+[1.0])
b=None
for st in range(2):
    z0=x0 if st==0 else x0+np.random.default_rng(st).normal(0,0.4,len(x0))*np.array([1 if n!='l0' else 15 for n in nm])
    r=minimize(c1,z0,method='Nelder-Mead',options=dict(maxiter=700,xatol=1e-7,fatol=1e-9,adaptive=True))
    print('start',st,r.fun,flush=True); json.dump(dict(start=st,fun=r.fun,x=list(r.x)),open(f's12_{LAW}_ckpt{st}.json','w'))
    if b is None or r.fun<b.fun: b=r
P=up(b.x); edges=[n for n,v,(lo,hi) in zip(nm,b.x,B) if n!='l0' and (v<lo+1e-3*(hi-lo) or v>hi-1e-3*(hi-lo))]
from scipy.stats import chi2 as C2, norm
mig=model(P,True); chi=mig['chi2']; k=len(nm); dchi=45.1366-chi; pval=C2.sf(max(dchi,0),k); sig=float(norm.isf(pval)) if dchi>0 else 0.0
out=dict(law=LAW,params=P,edges=edges,chi2=chi,R2=mig['R2'],sigma=sig,nfree=k,A=mig['A'])
# Step 2: same universe, no refit. H=70 (as in step 1) and H rescaled so age=13 at same mc (M scales as H^2)
for tag,H in (('H70',70.0),):
    Q=dict(P); Q.pop('qc')
    e=model(Q,False); out['cmb_'+tag]=e
age70=out['cmb_H70']['age']; H13=70*age70/13; Q=dict(P); Q.pop('qc'); Q['H']=H13; Q['M']=P['M']*(H13/70)**2
out['cmb_age13']=model(Q,False); out['cmb_age13']['H']=H13
for t in ('cmb_H70','cmb_age13'):
    e=out[t]; e.pop('A',None)
json.dump(out,open(f's12_{LAW}.json','w'),indent=1,default=float)
print('DONE',LAW,chi,mig['R2'],sig,edges,{k:out['cmb_H70'][k] for k in ('cancel','rms','age','Hmean','vo')},flush=True)
