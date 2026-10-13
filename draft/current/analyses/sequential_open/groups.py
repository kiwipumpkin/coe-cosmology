"""split each final layout's ring pattern into: bodies inside our sphere / in the CoE core / everywhere else"""
import sys,json,numpy as np
from common import *
def split(X,M,kernel):
    s=np.linalg.norm(X,axis=1); dc=np.linalg.norm(X-C_G,axis=1)
    G={"inside_us":s<RMIN,"core":(dc<A_IN_G)&(s>=RMIN),"rest":(dc>=A_IN_G)&(s>=RMIN)}
    W={g:(kernel(X[m],M[m]) if m.any() else np.zeros(360)) for g,m in G.items()}
    tot=sum(W.values()); out=dict(total=float(rms(R-tot)/FLOOR))
    for g in G: out[g]=dict(n=int(G[g].sum()),own_rms=float(rms(W[g])/FLOOR),without=float(rms(R-tot+W[g])/FLOOR),alone=float(rms(R-W[g])/FLOOR))
    return out
def deep_kernel(X,M):
    W=np.zeros(360)
    for i in range(0,len(X),3000): W+=M[i:i+3000]@seq_core.unit_pattern_fast(X[i:i+3000])
    return W
if __name__=="__main__":
    res={}
    for tag in sys.argv[1:]:
        f=np.load(f"final_{tag}.npz")
        if tag.startswith("shallow"):
            import shallow_kern as sk; X=f["X"]; res[tag]=split(X,np.zeros(len(X)),lambda X_,M_: sk.pat(X_))
        else: res[tag]=split(f["X"],f["M"],deep_kernel)
        print(tag,json.dumps(res[tag]))
    old=json.load(open("groups.json")) if os.path.exists("groups.json") else {}
    old.update(res); json.dump(old,open("groups.json","w"),indent=1)
