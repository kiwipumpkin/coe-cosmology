import sys,json,warnings;warnings.filterwarnings('ignore');sys.path.insert(0,'.');from static import *
from scipy.optimize import brentq
FB=0.015426145753569655
psis=np.arange(0,180.001,0.05);w=np.sin(np.radians(psis));w[0]=w[-1]=np.sin(np.radians(0.025))/4;w/=w.sum()
out={}
for f in (0.0,FB):
    # threshold speeds: cap centre (s=0.5) and cap edge (s=sqrt(0.75)) reachable today iff t0(u) > s
    thr={}
    for nm,sd in (('centre',0.5),('edge',np.sqrt(0.75))):
        thr[nm]=brentq(lambda u:Model(u,f).t0-sd,0.1,10)
    rows=[]
    for u in (0.3,0.5,0.8,1.0,10**0.021927885054513485,1.1,1.2,1.5,1.9,2.0,2.5,4,5,10,30):
        m=Model(u,f);r=sky(m,psis);st=np.array([x['status'] for x in r])
        fr={k:float(w[st==k].sum()) for k in sorted(set(st))}
        opz=np.array([x.get('opz',np.nan) for x in r]);good=np.isfinite(opz)&(opz>0)
        rows.append(dict(u=u,t0=m.t0,tc=m.tc,frac=fr,frac_positive_redshift=float(w[good].sum()),
                         frac_emitted_before_our_crossing=float(w[[x.get('before_our_crossing',False) for x in r]].sum())))
        print(f'f={f:.4f} u={u:.4f} t0={m.t0:.4f} tc={m.tc:.4f} positive 1+z sky fraction {w[good].sum():.5f} ', {k:round(v,5) for k,v in fr.items()})
    out[str(f)]=dict(threshold_speeds=thr,rows=rows);print('thresholds',thr)
out['cap_fraction_geometric']=float((1-np.cos(np.radians(30)))/2);print('cap fraction',out['cap_fraction_geometric'])
json.dump(out,open('feasibility.json','w'),indent=1,default=float)
