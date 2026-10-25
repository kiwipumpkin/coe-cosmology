# Exact projection everywhere: exact Fwd for l<=50; l>50 Limber x exact/Limber ratio on dense nodes (every ~40 l to 2500, log-interp).
# Ratio at l=50 matches Fwd exact to <0.2% -> continuous. Limber never reaches 1% agreement (ratio 0.96-0.99 at l 700-2500) so no Limber-only regime.
# Iterate: corr at current best -> Nelder-Mead polish -> repeat. Checkpoint fixed_ckpt_<mode>.json each step.
import sys,os,json,numpy as np
from multiprocessing import Pool
import half_fit as H, acoustic_model as M, acoustic_fit as A
from exact2 import node
NODES=np.r_[50,55,60,70,80,90,100,115,130,150,175,200:2501:40]
def corr(p):
    with Pool(8) as pl:r=np.array(pl.map(node,[(p,int(L),0.5) for L in NODES]))
    c=np.ones(M.LMAX+1);ll=np.arange(M.LMAX+1);s=ll>M.LX;c[s]=np.exp(np.interp(np.log(ll[s]),np.log(NODES),np.log(r)));return c,r
if __name__=='__main__':
    mode=sys.argv[1];npass=int(sys.argv[2]);ck=f'fixed_ckpt_{mode}.json'
    C=json.load(open(ck)) if os.path.exists(ck) else dict(passes=[])
    if not C['passes']:
        h=json.load(open(f'half_ckpt_{mode}.json'));xs=[s[1] for s in sorted(h['stage2'],key=lambda t:t[2])[:4]]
    else:xs=C['passes'][-1]['xs']
    while len(C['passes'])<npass:
        c,r=corr(H.P(xs[0],mode));H.CORR[:]=c;print('ratio',np.round(r,3).tolist(),flush=True)
        with Pool(4) as pl:R=sorted(pl.map(H.run,[(x,mode,1200) for x in xs]),key=lambda t:t[2])
        xs=[t[1] for t in R];C['passes'].append(dict(corr=c.tolist(),ratio=r.tolist(),xs=xs,f=[t[2] for t in R]));json.dump(C,open(ck,'w'));print('pass',len(C['passes']),R[0][2],flush=True)
    c,r=corr(H.P(xs[0],mode));H.CORR[:]=c;C['corr']=c.tolist();C['final_ratio']=r.tolist();C['best']=H.full(xs[0],mode)
    D=M.Dfield(C['best']['p'])*H.CORR+M.Dfix;j=np.abs(np.diff(D[2:]))/D[3:];ll=np.arange(3,M.LMAX+1)
    C['adjacent_step']=dict(max_30_100=float(j[(ll>=30)&(ll<=100)].max()),at_50_51=float(abs(D[51]-D[50])/D[50]),max_40_60=float(j[(ll>=40)&(ll<=60)].max()),median_30_100=float(np.median(j[(ll>=30)&(ll<=100)])))
    json.dump(C,open(ck,'w'),indent=1);print(C['best']['l2_2500'],C['best']['l2_30'],C['adjacent_step'])
