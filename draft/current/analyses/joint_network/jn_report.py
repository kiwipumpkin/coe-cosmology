"""Charts + results.json for the joint inside+outside network fits (reads raw/, null/; earlier runs read-only)."""
import os,sys,json,glob,numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from jn_core import *
CSd="/workspace/clump_scan"; SSd=f"{CSd}/solution_space"
plt.rcParams.update({"font.size":13,"axes.titlesize":14})
LVL=["free","cont","loose","web","strict"]
LAB={"free":"joint, no rules (core open, bodies anywhere)","cont":"joint + boundary continuity",
     "loose":"joint + continuity + web rules (loose)","web":"joint + continuity + web rules (strict), no density-trend cap",
     "strict":"joint + continuity + web rules + density-trend cap (strict)"}
COL={"free":"#1f77b4","cont":"#ff7f0e","loose":"#2ca02c","web":"#9467bd","strict":"#d62728"}
def score(x): lam=LEVELS[x["level"]]["lam"]; return x["left"]**2+(lam*x["soft_rms"])**2*len(x["p"])
def load(level,sym):
    out={}
    for f in glob.glob(f"raw/{level}_{sym}_N*.json"):
        r=[x for x in json.load(open(f)) if x["ok"]]
        if r:
            N=len(r[0]["kinds"]); met=[x for x in r if x["soft_max"]<=0.1]
            out[N]=dict(best_left=min(r,key=lambda x:x["left"]),best_obj=min(r,key=score),n=len(r),
                        met=min(met,key=lambda x:x["left"]) if met else None)
            out[N]["show"]=out[N]["met"] or out[N]["best_obj"]
    return dict(sorted(out.items()))
def bb(files):
    b={}
    for f in files:
        r=json.load(open(f)); r=r["runs"] if isinstance(r,dict) else r
        for x in r:
            if x.get("viol",0)>1e-3 or x.get("ok",True) is False: continue
            k=x.get("kinds") or x.get("name"); N=len(k)
            if N not in b or x["left"]<b[N]: b[N]=x["left"]
    return dict(sorted(b.items()))
def refs():
    ins=json.load(open(f"{CSd}/inside_visible/results_cap2e52.json")); bi={}
    for c in ins.values():
        x=c["best_run"]; N=len(x["kinds"])
        if x["viol"]<=1e-3 and (N not in bi or x["left"]<bi[N]): bi[N]=x["left"]
    R_={"inside our sphere only (earlier)":dict(sorted(bi.items())),
        "outside only, smooth limit (earlier)":bb(glob.glob(f"{SSd}/no_core_wall/raw_smooth/*.json")),
        "all open, smooth limit (earlier)":bb(glob.glob(f"{SSd}/all_open/raw_smooth/*.json"))}
    Rn={"outside only, no mirror (earlier)":bb(glob.glob(f"{SSd}/no_core_wall/raw_smooth_asym/*.json"))}
    return R_,Rn
def ok_rules(x,tol=0.05):
    """rules effectively met: every soft penalty <= tol (5 % past its tolerance)"""
    return x["soft_max"]<=tol
if __name__=="__main__":
    D={(lv,s):load(lv,s) for lv in LVL for s in "mn"}
    R_,Rn=refs()
    # ------------- best vs N -------------
    fig,axs=plt.subplots(2,1,figsize=(10,15))
    for ax,s,rf,ttl in ((axs[0],"m",R_,"(a) Mirror layouts allowed"),(axs[1],"n",Rn,"(b) No mirror layouts")):
        for (lab,c),st in zip(rf.items(),[":","--","-."]):
            N=[n for n in c if n>=4]; ax.plot(N,[c[n] for n in N],st,color="0.55",lw=1.8,label=lab)
        for lv in LVL:
            d=D[(lv,s)]
            if not d: continue
            N=list(d); y=[d[n]["show"]["left"] for n in N]; ax.plot(N,y,"-",color=COL[lv],lw=2.3,label=LAB[lv])
            for n,yy in zip(N,y):
                ax.plot(n,yy,"o",ms=8,color=COL[lv],mfc=COL[lv] if d[n]["met"] else "white",mew=2)
        ax.axhline(1,color="green",lw=2.2); ax.text(3.3,1.12,"1× target (CMB flat)",color="green",fontsize=12,va="bottom",ha="left")
        ax.axhline(237.3,color="0.6",lw=1); ax.text(33.2,250,"start 237×",color="0.4",fontsize=11,va="bottom",ha="right")
        ax.set_yscale("log"); ax.set_ylim(0.01,1500); ax.set_xlim(3,33.5); ax.set_xticks([4,8,12,16,20,24,32])
        ax.set_xlabel("number of bodies N (masses + voids, half each in the joint fits)"); ax.set_ylabel("best leftover (× floor, log scale)")
        ax.plot([],[],"o",color="0.3",ms=8,label="filled dot = all rules met (every penalty within 10 % of its tolerance)")
        ax.plot([],[],"o",color="0.3",mfc="white",mew=2,ms=8,label="open dot = no fit met every rule; best compromise shown")
        ax.set_title(ttl); ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.11),frameon=False,fontsize=11.5,ncol=1)
    fig.tight_layout(); fig.savefig("best_vs_N_logy.png",dpi=110); plt.close(fig)
    json.dump({f"{lv}_{s}":{N:dict(shown_left=v["show"]["left"],rules_met=v["met"] is not None,best_left=v["best_left"]["left"],best_left_soft_max=v["best_left"]["soft_max"],
              obj_left=v["best_obj"]["left"],obj_soft_max=v["best_obj"]["soft_max"],obj_breakdown=v["best_obj"]["breakdown"],
              n_ok=v["n"],props={k:vv for k,vv in v["show"]["props"].items() if k not in("links","link_len")}) for N,v in d.items()} for (lv,s),d in D.items()},
              open("results.json","w"),indent=1)
    print("refs",R_,Rn)
    for (lv,s),d in D.items():
        print(lv,s,{N:(round(v["show"]["left"],3),bool(v["met"]),round(v["best_left"]["left"],3),round(v["best_left"]["soft_max"],2),round(v["best_obj"]["left"],3),round(v["best_obj"]["soft_max"],2)) for N,v in d.items()})
