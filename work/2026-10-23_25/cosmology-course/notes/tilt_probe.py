import numpy as np
from nobox_check import exact_ext, K0
FU=0.94; a=-0.06; ph=0.35
def parts(x,e,al,ka):
    L=1/x; S=(1-FU)*(6*L**3*np.exp(1/L)-6*L**3-6*L**2-3*L)+FU
    R=1/S; q=1-R; p=(1-FU)*x*R; k=al*ka
    return k,R,q,p
def variants(x,e,al,ka):
    k,R,q,p=parts(x,e,al,ka); H0=1+2*k-3*k*R; X=1.5*k*(2*q-p); Y=4.5*k*(4*q-p); C=6*k*q
    out={}
    for sg,name in ((1,'max1'),(-1,'max2')):
        d0=-a*np.cos(ph)*H0/C
        d1=-a*np.cos(ph)*(H0-sg*e*X)/(C-sg*e*Y)
        hh=(H0-sg*e*X); h2=-(C-sg*e*Y)
        d2=d1*(1+a*np.sin(ph)*(1-hh/h2))
        out[name]=[np.degrees(v) for v in (d0,d1,d2)]
    b=0.75*e; dp=-a*np.sin(ph)*(1-k)/C
    out['min1']=[np.degrees(b+dp), np.degrees(b+dp*(1+b/np.tan(ph))), np.degrees(b+dp*(1+b/np.tan(ph))*(1+a*np.cos(ph)*(1-(1-k)/C)))]
    out['min2']=[np.degrees(b+dp), np.degrees(b+dp*(1-b/np.tan(ph))), np.degrees(b+dp*(1-b/np.tan(ph))*(1+a*np.cos(ph)*(1-(1-k)/C)))]
    return out
def cd(u,v): return abs((u-v+180)%360-180)
def compare(x,e,al,ka,verbose=False):
    ex,_,_=exact_ext(x,e,al,ka,amp=a)
    v=variants(x,e,al,ka)
    pred={'max1':[d for d in v['max1']],'max2':[180-d for d in v['max2']],'min1':[90+d for d in v['min1']],'min2':[270-d for d in v['min2']]}
    mx=[l for t,l,A in ex if t=='max']; mn=[l for t,l,A in ex if t=='min']
    if len(mx)!=2 or len(mn)!=2: return None
    err={}
    for kname,arr in pred.items():
        pool=mx if 'max' in kname else mn
        err[kname]=[min(cd(l,pp) for l in pool) for pp in arr]
    if verbose: print(ex); print(pred)
    return err
if __name__=="__main__":
    compare(5.9,0.1,1.0,K0,True)
    SW=[("x",1.01,15,{}),("alpha",0.01,8,{}),("kappa",0.001,2,{}),("eps",0.005,0.999,{"alpha":1.3})]
    for key,lo,hi,over in SW:
        print("==",key)
        for v in np.linspace(lo,hi,25):
            kw=dict(x=5.9,e=0.1,al=1.0,ka=K0); m={'x':'x','alpha':'al','kappa':'ka','eps':'e'}
            for kk,vv in over.items(): kw[m[kk]]=vv
            kw[m[key]]=v
            r=compare(**kw)
            if r is None: print("  %.3g  extrema count != 2/2"%v); continue
            print("  %.3g  "%v+"  ".join("%s:%s"%(n,"/".join("%.2f"%t for t in r[n])) for n in r))
