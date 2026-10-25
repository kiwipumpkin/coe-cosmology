import json,numpy as np
from ff_prof import fit
from ff_sim import pchi
from multiprocessing import Pool
rows=json.load(open('/tmp/rows.json'));I=json.load(open('interp_results.json'))
def par(r):
    n=r['name']
    if r['src']=='sim':return r['C'],r['frac'],r['DRV'],'poly'
    if n.startswith('path'):
        k=int(n[4:6]);return ([-0.041,-0.035,0.1353*k/10,0] if k else 'slc SL=0.041037 p=0.035'),30*(20000/30)**(k/10),2.937,'poly' if k else 'slc'
    if n.startswith('grid'):
        c3=float(n.split('_')[2]);return [-0.041,-0.035,c3,0],float(n.split('_f')[-1]),2.937,'poly'
    if n.startswith('poly'):return [-0.041,-0.035,0,0],30,2.937,'poly'
    return n,I[n].get('frac'),2.937,'hyb'
def w(r):
    if 'mig' not in r:r['mig']=fit((r['ratio'],-10,100));r['mig360']=fit((r['ratio'],0,360))
    return r
if __name__=='__main__':
    for r in rows:r['C'],r['frac'],r['DRV'],r['fam']=par(r);r['cmb']=pchi(r['Dl']);r['cmb26']=pchi(r['Dl'],2,6)
    sel=[r for r in rows if r['src']=='sim' or 1.04<r['ratio']<1.1]
    with Pool(6) as P:sel=P.map(w,sel)
    R=json.load(open('fullfree_results.json'));R['universes']=sel;json.dump(R,open('fullfree_results.json','w'),indent=1)
    for key in ['mig','mig360']:
        print('==',key)
        for r in sorted(sel,key=lambda r:(round(r[key]['chi2'],4),r['cmb']))[:14]:
            t=r[key]['theta'];print(r['name'],r['fam'],r['C'],'%.4g'%r['frac'],r['DRV'],'rat %.4f'%r['ratio'],'chi2 %.5f R2 %.4f'%(r[key]['chi2'],r[key]['R2']),np.round(t,3).tolist(),'CMB %.4g  l2-6 %.4g'%(r['cmb'],r['cmb26']),'D2-6',['%.3g'%x for x in r['Dl'][2:7]])
