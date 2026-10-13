"""Demo null ensemble: random rule-following void+cluster networks through the existing CMB kernel.
usage: python demo_ensemble.py <Mmin> <n_real | start:end> [lcdm|hamaus|exact]  -> results/demo_<Mmin>_<mode>.json (appended after every realization)"""
import os,sys,json,time,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
ME=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,ME); sys.path.insert(0,os.path.join(ME,"..","joint_network"))
import vmn_gen as G
_c=os.getcwd(); import jn_core as J; os.chdir(_c)
def pattern(X,M,chunk=20000):
    out=np.zeros(360)
    for i in range(0,len(M),chunk): out+=M[i:i+chunk]@J.unit_field(X[i:i+chunk])
    return out-M.sum()*J.FLOC
if __name__=="__main__":
    Mmin=float(sys.argv[1]); a,b=(sys.argv[2].split(":") if ":" in sys.argv[2] else (0,sys.argv[2])); a,b=int(a),int(b); mode=sys.argv[3] if len(sys.argv)>3 else "lcdm"
    os.makedirs(os.path.join(ME,"results"),exist_ok=True)
    fn=os.path.join(ME,"results",f"demo_{Mmin:.0e}_{mode}_{a}-{b}.json"); res=json.load(open(fn)) if os.path.exists(fn) else []
    done={r["seed"] for r in res}
    for s in range(a,b):
        if s in done: continue
        t=time.time(); g=G.generate(Mmin=Mmin,seed=1000+s,comp_mode=mode); X,M=G.kernel_bodies(g); P=pattern(X,M)
        d=np.linalg.norm(X,axis=1); ins=d<J.RV
        Pin=pattern(X[ins],M[ins]); Pout=P-Pin
        r=dict(seed=s,Mmin=Mmin,mode=mode,n_bodies=int(len(M)),n_clusters=int(len(g["MC"])),n_voids=int(g["void_in"].sum()),
               left=float(J.rms(J.R-P)/J.FLOOR),pattern_rms=float(J.rms(P)/J.FLOOR),
               pattern_rms_inside=float(J.rms(Pin)/J.FLOOR),pattern_rms_outside=float(J.rms(Pout)/J.FLOOR),secs=time.time()-t)
        if s<3: np.save(os.path.join(ME,"results",f"pattern_{Mmin:.0e}_{mode}_{s}.npy"),P)
        res.append(r); json.dump(res,open(fn+".tmp","w")); os.replace(fn+".tmp",fn); print(r,flush=True)
