"""Collect curves + best example layouts for every point-object scenario from SAVED results, re-score residuals. -> collected.pkl"""
import os,sys,json,glob,pickle,numpy as np
os.environ["NCW"]="smooth"; os.environ["AO"]="smooth"
CS="/workspace/clump_scan"; SS=f"{CS}/solution_space"
for p in (SS,f"{CS}/inside_visible",f"{SS}/no_core_wall",f"{SS}/all_open",f"{SS}/extreme_partners"): sys.path.insert(0,p)
_c=os.getcwd()
from ss_core import R,FLOOR,TH,rms,amps,Model,C_G,A_IN_G,RV
from extend import encode
from inside_core import InModel
from ncw_core import NCW
from ao_core import AO
import ep_core
os.chdir(_c)
def load_runs(pattern,ok=lambda r:True):
    out=[]
    for f in glob.glob(pattern):
        d=json.load(open(f)); runs=d["runs"] if isinstance(d,dict) else d
        kinds=d.get("kinds") if isinstance(d,dict) else None
        for r in runs:
            r=dict(r); r.setdefault("kinds",kinds)
            if r.get("viol",0)<=1e-3 and ok(r): out.append(r)
    return out
def best_by_N(runs):
    b={}
    for r in runs:
        N=len(r["kinds"])
        if N not in b or r["left"]<b[N]["left"]: b[N]=r
    return b
def pick_N(curve):
    below=[n for n in sorted(curve) if curve[n]<1]; return below[0] if below else min(curve,key=curve.get)
def mode_info(v):
    a=amps(v)/FLOOR; rec=np.full(360,v.mean())
    for m in range(1,13): rec+=2*np.mean(v*np.cos(m*TH))*np.cos(m*TH)+2*np.mean(v*np.sin(m*TH))*np.sin(m*TH)
    return a,float(np.sqrt(2)*rms(v-rec)/FLOOR)
S={}   # scenario -> dict
def add(key,label,family,curve,ex=None,note="",xlabel="number of bodies N"):
    S[key]=dict(label=label,family=family,curve={int(k):float(v) for k,v in curve.items()},ex=ex,note=note,xlabel=xlabel)
def example(N,kinds,X,M,Rv,w,saved,extra=None,zoom=None):
    v=R-w; a,hi=mode_info(v); left=rms(v)/FLOOR
    X=np.array(X,float).reshape(-1,3)
    if zoom is None: zoom=15 if (len(X) and np.linalg.norm(X,axis=1).max()<RV) else 45
    e=dict(N=N,kinds=kinds,X=X,M=np.array(M,float),Rv=np.array(Rv,float),resid=v,amps=a,hi=hi,left=left,saved=saved,zoom=zoom)
    if extra: e.update(extra)
    print(f"   rescored {left:.4f} vs saved {saved:.4f}  {'OK' if abs(left-saved)<=1e-3*max(1,saved) else 'MISMATCH'}")
    return e
# ---------- outside (solution_space, hard wall) ----------
raw=load_runs(f"{SS}/raw/*.json")
for r in raw: pass
def ss_from_xyz(r):
    m=Model(r["kinds"]); p=encode(r["kinds"],r["X"],r["M"],r["Rv"]); return m.pattern(p)
ext_m=json.load(open(f"{SS}/extend_masses.json")); ext_x=json.load(open(f"{SS}/extend_mixed.json"))
cold=[r for f in ("cold_CCCCCCVVVVVV_polished.json","cold_CCCCCCVVVVV_polished.json") for r in json.load(open(f"{SS}/{f}")) if r["viol"]<=1e-3]
# masses only
mo=best_by_N([r for r in raw if set(r["kinds"])=={"C"}]); cm={n:r["left"] for n,r in mo.items()}
for s in ext_m[1:]: cm[s["N"]]=s["left"]
print("outside masses only"); s=ext_m[-1]
add("out_masses","Outside: masses only","outside",cm,example(s["N"],s["kinds"],s["X"],s["M"],s["Rv"],Model(s["kinds"]).pattern(np.array(s["p"])),s["left"]),
    "Hard wall at the CoE interior; 14.3–43 Gpc from us. Never reaches 1×; best at N=20 shown.")
# masses + voids
mx=best_by_N([r for r in raw if set(r["kinds"])=={"C","V"}]+[dict(r) for r in cold]); cx={n:r["left"] for n,r in mx.items()}
for s in ext_x[1:]: cx[s["N"]]=min(cx.get(s["N"],1e9),s["left"])
print("outside masses+voids"); s=ext_x[-1]
add("out_mixed","Outside: masses + deep voids (mirror layouts allowed)","outside",cx,
    example(s["N"],s["kinds"],s["X"],s["M"],s["Rv"],Model(s["kinds"]).pattern(np.array(s["p"])),s["left"]),"Hard wall; first below 1× at N=11 (layout came out mirror-symmetric).")
# voids only
vo=best_by_N([r for r in raw if set(r["kinds"])=={"V"}]); cv={n:r["left"] for n,r in vo.items()}
print("outside voids only"); Nb=min(cv,key=cv.get); r=vo[Nb]
add("out_voids","Outside: deep voids only","outside",cv,example(Nb,r["kinds"],r["X"],r["M"],r["Rv"],ss_from_xyz(r),r["left"]),"Voids capped at 4 Gpc; never beats the four clumps.")
# no mirror
asym=best_by_N(load_runs(f"{SS}/asym/raw/*.json",lambda r:r["mirror_matches"]==0 and r["on_plane"]==0)); ca={n:r["left"] for n,r in asym.items()}
print("outside no-mirror"); Nb=pick_N(ca); r=asym[Nb]
add("out_asym","Outside: masses + deep voids, no mirror layouts","outside",ca,
    example(Nb,r["kinds"],r["X"],r["M"],r["Rv"],Model(r["kinds"]).pattern(np.array(r["p"])),r["left"]),"No same-kind pair within 20° / 15% / 25% of each other's mirror image; nothing within 10° of the mirror plane.")
# giant partners
ep=[]
for f in glob.glob(f"{SS}/extreme_partners/raw*/*.json"):
    if "raw_wide" in f: continue
    ep+=[r for r in json.load(open(f)) if r["viol"]<=1e-3 and abs(r["net"])<=1e51]
cp={}
for r in ep:
    n=r["npart"]+len(r["kinds"]); cp[n]=min(cp.get(n,1e9),r["left"])
b=min(ep,key=lambda r:r["left"]); m=ep_core.EP(b["sign"],b["npart"],b["kinds"]); p=np.array(b["p"]); w=m.pattern(p)
print("giant partners")
e=example(b["npart"]+len(b["kinds"]),b["kinds"],b["Xo"],b["Mo"],b["Rvo"],w,b["left"],
    extra=dict(partners=dict(X=np.array(b["Xp"]),R=np.array(b["Rp"]),M=np.array(b["Mp"]),delta_us=b["delta_us"])),zoom=62)
gw=[r for r in ep if max(r["Rp"])>=8]; g=min(gw,key=lambda r:r["left"]); mg=ep_core.EP(g["sign"],g["npart"],g["kinds"])
print("giant partners (giant-void example)")
eg=example(g["npart"]+len(g["kinds"]),g["kinds"],g["Xo"],g["Mo"],g["Rvo"],mg.pattern(np.array(g["p"])),g["left"],
    extra=dict(partners=dict(X=np.array(g["Xp"]),R=np.array(g["Rp"]),M=np.array(g["Mp"]),delta_us=g["delta_us"])),zoom=62)
add("partners","Giant partners (our sphere over/underdense + 1–3 giant partners, ± ordinary bodies)","outside",cp,e,
    f"N counts partners + ordinary bodies. Best: {b['setup']} (δ_us={b['delta_us']:+.3f}).")
S["partners"]["ex2"]=eg; S["partners"]["ex2_note"]=f"Best layout that really uses a giant partner: {g['setup']} (δ_us={g['delta_us']:+.2f}, giant void R={max(g['Rp']):.1f} Gpc)."
# no wall
for var,lab in (("smooth","Outside, no wall: smooth-density limit"),("open","Outside, no wall: fully open")):
    B=best_by_N(load_runs(f"{SS}/no_core_wall/raw_{var}/*.json")); c={n:r["left"] for n,r in B.items()}
    print(lab); Nb=pick_N(c); r=B[Nb]
    add(f"nowall_{var}",lab,"nowall",c,example(Nb,r["kinds"],r["X"],r["M"],r["Rv"],NCW(r["kinds"]).pattern(np.array(r["p"])),r["left"]),
        "CoE-interior wall removed" + ("; mass ≤ 2e52·ρ(x)/ρ̄, supply volumes M/ρ may not overlap." if var=="smooth" else "; usual caps only."))
# inside
for tag,mhi,comp,key,lab,mixes in (("cap2e52",2e52,"local","in_masses","Inside our sphere: masses only",{"masses only"}),
                                    ("cap2e52",2e52,"local","in_mixed","Inside our sphere: masses + deep voids",{"masses + deep voids"}),
                                    ("cap2e52",2e52,"local","in_voids","Inside our sphere: deep voids only",{"deep voids only"}),
                                    ("cap1e53",1e53,"local","in_cap1e53","Inside, masses only, cap raised to 1e53 kg",{"masses only"}),
                                    ("nocomp",2e52,"none","in_nocomp","Inside, masses only, no compensation",{"masses only"}),
                                    ("surrcomp",2e52,"surround","in_surr","Inside, masses only, compensation taken from outside",{"masses only"})):
    res=json.load(open(f"{CS}/inside_visible/results_{tag}.json")); B={}
    for c_ in res.values():
        if c_["mix"] in mixes and (c_["N"] not in B or c_["best"]<B[c_["N"]]["best"]): B[c_["N"]]=c_
    c={n:v["best"] for n,v in B.items()}; Nb=pick_N(c); r=B[Nb]["best_run"]; print(lab)
    add(key,lab,"inside" if tag=="cap2e52" else "inside_var",c,
        example(Nb,r["kinds"],r["X"],r["M"],r["Rv"],InModel(r["kinds"],mhi,comp).pattern(np.array(r["p"])),r["left"],zoom=15),
        "Bodies 0.1–14.2 Gpc from us (never touching the ring); net mass removed smoothly over our sphere" if comp=="local" else
        ("no net-mass compensation" if comp=="none" else "net mass taken from the 14.3–30 Gpc surround"))
# all open
B=best_by_N(load_runs(f"{SS}/all_open/raw_smooth/*.json")); c={n:r["left"] for n,r in B.items()}; Nb=pick_N(c); r=B[Nb]; print("all open")
add("all_open","All open (anywhere 0.5–43 Gpc; smooth-density limit)","allopen",c,
    example(Nb,r["kinds"],r["X"],r["M"],r["Rv"],AO(r["kinds"]).pattern(np.array(r["p"])),r["left"]),"No wall, no visible-sphere exclusion; bodies may not straddle the ring.")
# ---------- sequential ----------
seq=f"{CS}/sequential"; lf=json.load(open(f"{seq}/local_fit.json")); bs=min(lf["runs"],key=lambda r:r["left"])
ca={i:t["left"] for i,t in enumerate(bs["traj"]) if i>0}
z=np.load(f"{seq}/local_frozen.npz"); print("seq A deep")
eA=example(48,"".join("C" if k>0 else "V" for k in z["M"]),z["X"],z["M"],np.where(z["M"]<0,z["r"],0),z["w"],bs["left"])
add("seqA_deep","Sequential A: near region fixed pair by pair (deep voids + compact clusters)","seq",ca,eA,
    "24 void+cluster pairs placed one at a time (14.3–30 Gpc), each chosen to lower the leftover most.",xlabel="bodies placed")
nf=json.load(open(f"{seq}/shallow/near_fit.json")); bb=nf["best"]; cs_={i:v for i,v in enumerate(bb["traj"]) if i>0}; print("seq A shallow")
Xs=np.array(bb["X"]); eAs=example(bb["n"],"V"*bb["n"],Xs,-np.ones(len(Xs)),np.array(bb["Rv"]),np.array(bb["w"]),bb["left"])
add("seqA_shallow","Sequential A: near region fixed, shallow (Clampitt & Jain) voids","seq",cs_,eAs,"Shallow compensated voids placed one at a time.",xlabel="bodies placed")
for rule,lab in (("cmb_aware","Sequential B: grow from empty, CMB-aware"),("foam_rim","Sequential B: grow from empty, foam rule"),("random","Sequential B: grow from empty, random (control)")):
    runs=[json.load(open(f)) for f in glob.glob(f"{seq}/b_{rule}_Rv1.5_s*.json")]; d=min(runs,key=lambda r:r["min_left"])
    add(f"seqB_{rule}",lab,"seqB",dict(zip(d["traj_obj"],d["traj_left"])),None,f"final {d['final_left']:.0f}×, lowest {d['min_left']:.1f}×, {d['n_obj']} bodies",xlabel="bodies placed")
    S[f"seqB_{rule}"]["final_amps"]=d["final_amps"]
g=[json.load(open(f)) for f in glob.glob(f"{seq}/shallow/grow_cmb_aware_s*.json")]; d=min(g,key=lambda r:r["min"])
add("seqB_shallow","Sequential B: grow from empty, CMB-aware, shallow voids","seqB",{5*i:v for i,v in enumerate(d["traj"]) if i>0},None,
    f"final {d['final']:.0f}×, lowest {d['min']:.0f}×, {d['n_void']} voids",xlabel="bodies placed")
pickle.dump(dict(S=S,R=R,FLOOR=FLOOR,TH=TH,C_G=C_G,A_IN_G=A_IN_G,RV=RV,start_amps=amps(R)/FLOOR,start_hi=mode_info(R)[1]),open("collected.pkl","wb"))
for k,v in S.items(): print(k,"| curve N:",min(v["curve"]),"-",max(v["curve"]),"best %.3g"%min(v["curve"].values()),"| example N",v["ex"]["N"] if v["ex"] else None)
