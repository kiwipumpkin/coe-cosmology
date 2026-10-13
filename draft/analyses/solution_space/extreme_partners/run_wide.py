import os,sys,json,time
os.environ["OUTDIR"]="raw_wide"
import run_ep_b as rb
rb.SETUPS[:]=[("overwide_3p+4C","overwide",3,"CCCC",10),("overwide_2p+4C4V","overwide",2,"CCCCVVVV",10)]
from multiprocessing import Pool
if __name__=="__main__":
    os.makedirs("raw_wide",exist_ok=True); tasks=[(nm,sg,n,k,i) for nm,sg,n,k,ns in rb.SETUPS for i in range(ns)]; res={}
    with Pool(int(sys.argv[1])) as pool:
        for r in pool.imap_unordered(rb.job,tasks):
            res.setdefault(r["setup"],[]).append(r); print(r["setup"],r["start"],round(r["left"],2),round(r["delta_us"],3),r["viol"],r["net"],flush=True)
            json.dump(res[r["setup"]],open(f"raw_wide/{r['setup']}.json","w"))
