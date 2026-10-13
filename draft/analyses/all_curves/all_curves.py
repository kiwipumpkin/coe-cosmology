"""All best-leftover-vs-N curves from saved results (no refits)."""
import json,glob,os,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
CS="/workspace/clump_scan"; SS=f"{CS}/solution_space"; IV=f"{CS}/inside_visible"
def inside(tag,mix):
    r=json.load(open(f"{IV}/results_{tag}.json")); d={}
    for c in r.values():
        if c["mix"]==mix: d[c["N"]]=min(d.get(c["N"],1e9),c["best"])
    return d
def outside():
    m,x,v={},{},{}
    for f in glob.glob(f"{SS}/raw/*.json"):
        d=json.load(open(f)); k=d["kinds"]
        if "S" in k: continue
        L=min(r["left"] for r in d["runs"] if r["viol"]<=1e-3); N=len(k); dd=m if set(k)=={"C"} else v if set(k)=={"V"} else x
        dd[N]=min(dd.get(N,1e9),L)
    for f,dd in (("extend_masses.json",m),("extend_mixed.json",x)):
        for s in json.load(open(f"{SS}/{f}")):
            if s["N"]>8: dd[s["N"]]=min(dd.get(s["N"],1e9),s["left"])
    for f in ("cold_CCCCCCVVVVVV_polished.json","cold_CCCCCCVVVVV_polished.json"):
        for r in json.load(open(f"{SS}/{f}")):
            if r["viol"]<=1e-3: N=len(r["kinds"]); x[N]=min(x.get(N,1e9),r["left"])
    return m,x,v
asym={int(N):r["left"] for N,r in json.load(open(f"{SS}/asym/asym_results.json"))["best"].items()}
ep=json.load(open(f"{SS}/extreme_partners/ep_results.json")); part={}
for nm,r in ep.items():
    np_=int(nm.split("_")[1][0]); extra=4 if "+4C" in nm else 0; extra+=4 if "4V" in nm else 0
    N=np_+extra; part[N]=min(part.get(N,1e9),r["left"])
C={"m":"#c0392b","x":"#7d3c98","v":"#2471a3","cap":"#e67e22","noc":"#f1948a","sur":"#6e2c00","asym":"#17a589","part":"#555555"}
IN=[("inside: masses only",inside("cap2e52","masses only"),C["m"],"o"),
    ("inside: masses + deep voids",inside("cap2e52","masses + deep voids"),C["x"],"o"),
    ("inside: deep voids only",inside("cap2e52","deep voids only"),C["v"],"o"),
    ("inside: masses, cap 1e53 kg",inside("cap1e53","masses only"),C["cap"],"^"),
    ("inside: masses, no compensation",inside("nocomp","masses only"),C["noc"],"v"),
    ("inside: masses, compensation from outside",inside("surrcomp","masses only"),C["sur"],"D")]
om,ox,ov=outside()
OUT=[("outside: masses only",om,C["m"],"s"),("outside: masses + deep voids",ox,C["x"],"s"),("outside: deep voids only",ov,C["v"],"s"),
     ("outside: masses + voids, no mirror layouts",asym,C["asym"],"P")]
def refs(ax,xt):
    for v,t,c in [(237.3,"start 237×","0.45"),(88.8,"four clumps 89×","0.45"),(1.0,"CMB flat target 1×","green")]:
        ax.axhline(v,color=c,ls="-" if c=="green" else ":",lw=2 if c=="green" else 1.2,zorder=0)
        if xt is not None: ax.text(xt,v*1.12,t,fontsize=11,color=c if c=="green" else "0.35",ha="right")
def setup(ax):
    ax.set_yscale("log"); ax.set_ylim(0.01,300); ax.set_xlim(0.5,20.5); ax.set_xticks([1,2,4,6,8,10,12,14,16,18,20])
    ax.set_yticks([0.01,0.1,1,10,100]); ax.set_yticklabels(["0.01","0.1","1","10","100"]); ax.set_xlabel("number of objects N")
def draw(ax,series,ls,lw=2.4):
    for lab,d,col,mk in series:
        N=sorted(d); ax.plot(N,[d[n] for n in N],ls,color=col,marker=mk,ms=6,lw=lw,label=lab)
def partners(ax):
    N=sorted(part); ax.plot(N,[part[n] for n in N],"*",color=C["part"],ms=15,ls="none",label="outside: giant-partner fits (best, total objects)")
# single chart
fig,ax=plt.subplots(figsize=(15,8))
draw(ax,IN,"-"); draw(ax,OUT,"--"); partners(ax); refs(ax,20.3); setup(ax)
ax.set_ylabel("best leftover (× 10⁻⁵ floor, log scale)"); ax.set_title("Best leftover vs number of objects — solid = inside our visible sphere, dashed = outside")
ax.legend(loc="center left",bbox_to_anchor=(1.01,0.5),frameon=False,fontsize=12)
fig.tight_layout(); fig.savefig("all_curves_logy.png",dpi=105); plt.close(fig)
# panels
fig,(a1,a2)=plt.subplots(1,2,figsize=(17,7.5),sharey=True)
draw(a1,IN,"-"); draw(a2,OUT,"--"); partners(a2)
for a in (a1,a2): setup(a)
refs(a1,None); refs(a2,20.3)
a1.set_title("Inside our visible sphere (14.3 Gpc)"); a2.set_title("Outside (14.3–43 Gpc, plus giant partners)")
a1.set_ylabel("best leftover (× 10⁻⁵ floor, log scale)")
a1.legend(loc="lower left",frameon=False,fontsize=11); a2.legend(loc="lower left",frameon=False,fontsize=11)
for v,t,c in [(237.3,"start 237×","0.35"),(88.8,"four clumps 89×","0.35"),(1.0,"CMB flat 1×","green")]: a1.text(20.3,v*1.12,t,fontsize=11,color=c,ha="right")
fig.tight_layout(); fig.savefig("all_curves_panels_logy.png",dpi=100); plt.close(fig)
# table
allc={l:d for l,d,_,_ in IN+OUT}; allc["outside: giant-partner fits"]=part
with open("curves.json","w") as f: json.dump({k:{str(n):v for n,v in sorted(d.items())} for k,d in allc.items()},f,indent=1)
for k,d in allc.items():
    first=min([n for n in d if d[n]<1],default=None); print(f"{k:45s} first<1x at N={first}  best={min(d.values()):.3f} (N={min(d,key=d.get)})  N4={d.get(4)}")
