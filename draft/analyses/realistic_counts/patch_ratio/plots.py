import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from ps_core import *
R=json.load(open("results.json")); T=R["thresholds"]
plt.rcParams.update({"font.size":13})
def fmapB(B,p):
    df=B/p["Delta"]; f0=p["f0"]; return 100*(((f0+df)/(1-f0-df))/(f0/(1-f0))-1)
# ---------- chart 1: leftover vs allowed ratio shift
fig,ax=plt.subplots(figsize=(10,6.5))
cv=[("in05","extreme 50 Mpc","inside, 0.5 Gpc domains, patches +0.83 / −0.39","#c0392b","-"),
    ("in05","typical ±0.3, 50 Mpc","inside, 0.5 Gpc domains, patches ±0.3","#e67e22","--"),
    ("in20","extreme 50 Mpc","inside, 2 Gpc domains, patches +0.83 / −0.39","#2471a3","-"),
    ("both10","extreme 50 Mpc","inside + outside, no mirror, 1 Gpc, +0.83 / −0.39","#7d3c98","-")]
for c,pn,lab,col,ls in cv:
    if c not in T: continue
    p=pinfo(pn); B=[x[0] for x in T[c]["curve"]]; L=[max(x[1],1e-3) for x in T[c]["curve"]]
    if T[c]["B_1x"]:
        k=[i for i,b in enumerate(B) if b<T[c]["B_lower"]+1e-9]; B=[B[i] for i in k]+[T[c]["B_1x"]]; L=[L[i] for i in k]+[1.0]
    x=[fmapB(b,p) for b in B]; ax.plot(x,L,ls,color=col,lw=2.4,marker="o",ms=4,label=lab)
    if T[c]["B_1x"]: xs=fmapB(T[c]["B_1x"],p); ax.plot(xs,1,"o",ms=11,mfc="none",mec=col,mew=2.2)
ax.axhline(1,color="green",lw=2.5); ax.text(1.5,1.25,"1× = CMB flat",color="green",fontsize=13)
ax.set_yscale("log"); ax.set_xscale("log"); ax.set_ylim(5e-1,400); ax.set_xlim(1,150)
ax.set_xlabel("allowed shift of the over : under patch ratio from the cosmic mean (%)")
ax.set_ylabel("best CMB leftover (× floor)")
ax.set_title("Domains = fine mix of over- and under-dense patches:\nhow much must the mix ratio drift?")
ax.text(1.05,0.6,"open circle = first shift that reaches 1× (verified solution)",fontsize=11,color="0.35")
ax.legend(loc="upper right",frameon=False,fontsize=11.5); ax.grid(alpha=0.25,which="both")
fig.tight_layout(); fig.savefig("leftover_vs_ratio_shift.png",dpi=110); plt.close(fig)
# ---------- chart 2: required shift vs natural scatter
tab=R["table"]; fig,ax=plt.subplots(figsize=(10,6.5))
for pn,col,lab in [("extreme 30 Mpc","#922b21","patches 30 Mpc (+2.1 / −0.55)"),("extreme 50 Mpc","#c0392b","patches 50 Mpc (+0.83 / −0.39)"),
                   ("extreme 100 Mpc","#e59866","patches 100 Mpc (+0.24 / −0.25)"),("typical ±0.3, 50 Mpc","#2471a3","patches 50 Mpc, typical ±0.3")]:
    rr=sorted([t for t in tab if t["patch"]==pn and t["cfg"].startswith("in")],key=lambda t:t["H"])
    ax.plot([t["H"] for t in rr],[t["n_sigma_domain"] for t in rr],"-o",color=col,lw=2.4,ms=7,label=lab)
rr=sorted([t for t in tab if t["patch"]=="extreme 50 Mpc" and t["cfg"].startswith("in")],key=lambda t:t["H"])
ax.plot([t["H"] for t in rr],[t["n_sigma_LCDM"] for t in rr],"--s",color="0.3",lw=2,ms=7,label="vs ΛCDM σ at domain size")
ax.axhline(1,color="green",lw=1.5); ax.text(2.05,1.05,"1σ",color="green"); ax.axhline(3,color="green",lw=1.5,ls=":"); ax.text(2.05,3.1,"3σ",color="green")
ax.set_xscale("log",base=2); ax.set_yscale("log"); ax.set_xticks([0.5,1,2]); ax.set_xticklabels(["0.5","1","2"]); ax.set_xlim(0.45,2.4)
ax.set_xlabel("domain size (Gpc)"); ax.set_ylabel("required ratio drift per domain ÷ natural scatter (σ)")
ax.set_title("Is the needed drift per domain within ordinary randomness?\n(inside, mirror allowed; per domain, not yet counting coherence)")
ax.legend(loc="upper left",frameon=False,fontsize=11.5); ax.grid(alpha=0.25,which="both")
fig.tight_layout(); fig.savefig("required_vs_natural_scatter.png",dpi=110); plt.close(fig)
