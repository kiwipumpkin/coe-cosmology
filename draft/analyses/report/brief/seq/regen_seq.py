"""Regenerate final layouts of the sequential-B runs used in Fig 7 (same code + seed; the original runs saved only trajectories).
usage: python3 regen_seq.py b <rule> <seed>   |   python3 regen_seq.py shallow cmb_aware <seed>
Writes final_<tag>.npz here (X, r, M, kind, w, left) and checks the final leftover against the saved json. Originals untouched."""
import sys,os,json,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); SEQ="/workspace/clump_scan/sequential"
mode,rule,seed=sys.argv[1],sys.argv[2],int(sys.argv[3])
if mode=="b":
    sys.path.insert(0,SEQ); os.chdir(SEQ)
    import scenario_b as sb
    g=sb.Grow(rule,seed,1.5).run(); X,M,k,left,w=sb.trajectory(g.order); Xa,ra,Ma,ka=g.pk.arrays()
    ref=json.load(open(f"{SEQ}/b_{rule}_Rv1.5_s{seed}.json"))["final_left"]
    np.savez(f"{HERE}/final_b_{rule}_s{seed}.npz",X=X,M=M,kind=k.astype(str),Xp=Xa,rp=ra,kp=np.array(ka).astype(str),w=w,left=left[-1],ref=ref)
    print("DONE",rule,seed,"final",left[-1],"saved",ref,flush=True)
else:
    sh=f"{SEQ}/shallow"; sys.path.insert(0,sh); os.chdir(sh)
    src=open(f"{sh}/run_shallow.py").read().split('if __name__')[0]
    src=src.replace('json.dump(out,open(os.path.join(HERE,f"grow_{rule}_s{seed}.json"),"w"))','out["_X"]=X; out["_w"]=C[-1]; return out')
    g={"__file__":f"{sh}/run_shallow.py"}; exec(src,g)
    out=g["grow"](rule,seed); ref=json.load(open(f"{sh}/grow_{rule}_s{seed}.json"))["final"]
    np.savez(f"{HERE}/final_shallow_{rule}_s{seed}.npz",X=out["_X"],Rv=1.5,w=out["_w"],left=out["final"],ref=ref)
    print("DONE shallow",rule,seed,"final",out["final"],"saved",ref,flush=True)
