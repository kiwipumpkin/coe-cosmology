"""Data for Fig. waterfall (b): real fitted layouts from ../../../solution_space (read-only), re-evaluated with ss_core.
Steps: no bodies (237x) -> best 4 masses (raw/CCCC.json) -> best 4 masses + 4 voids (raw/CCCCVVVV.json)
-> warm-start chain 9, 10, 11 bodies (extend_mixed.json; each adds one body to the previous layout and refits all)."""
import os,sys,json,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); SS=os.path.abspath(os.path.join(HERE,"..","..","..","solution_space"))
cwd=os.getcwd(); os.chdir(SS); sys.path.insert(0,SS)
from ss_core import Model,R,FLOOR,rms
from extend import encode
def best(f): d=json.load(open(f)); return min([r for r in d["runs"] if r["viol"]<=1e-3],key=lambda r:r["left"])
lay=[dict(kinds="",X=[],M=[],Rv=[],src="no bodies")]
for f in ["raw/C.json","raw/CCCC.json","raw/CCCCVVVV.json"]:
    b=best(f); lay.append(dict(kinds=b["runs"][0]["kinds"] if "runs" in b else os.path.basename(f)[:-5],X=b["X"],M=b["M"],Rv=b["Rv"],src=f))
for r in json.load(open("extend_mixed.json"))[1:]:
    lay.append(dict(kinds=r["kinds"],X=r["X"],M=r["M"],Rv=r["Rv"],src="extend_mixed.json N=%d"%r["N"]))
out=[]
for L in lay:
    k=L["kinds"]
    if not k: w=np.zeros(360)
    else: m=Model(k); w=m.pattern(encode(k,L["X"],L["M"],L["Rv"]))
    res=(R-w)/FLOOR
    out.append(dict(kinds=k,src=L["src"],X=L["X"],M=L["M"],Rv=L["Rv"],left=float(rms(R-w)/FLOOR),resid_x=res.tolist()))
    print(L["src"],k,round(out[-1]["left"],3))
os.chdir(cwd); json.dump(out,open(os.path.join(HERE,"placement_steps.json"),"w"))
