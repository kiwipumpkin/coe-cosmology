import json, multiprocessing as mp, solve as S
Ds=(14.7,15.0,17.0,500.0,1000.0,2000.0)
J=[(s,D,"median",False,0) for s in S.SL for D in Ds]
def job(a):
    try: r=S.solve(*a); r.pop("curve",None); return r
    except Exception as e: return dict(slice=a[0],D=a[1],err=repr(e))
if __name__=="__main__":
    with mp.Pool(7) as P: R=P.map(job,J,chunksize=1)
    json.dump(R,open("dscan_wide.json","w"),default=float)
    for r in R:
        if "err" in r: print(r["slice"],r["D"],"ERR",r["err"][:120]); continue
        print(r["slice"],"D %g resid %.1e mcc %.5f l0 %.1f amp %.4f R2 %.3f Mc %.3g"%(r["D"],r["resid"],r["mcc"],r["l0"],r["amp"],r["R2_point"],r["Mc_kg"]),flush=True)
