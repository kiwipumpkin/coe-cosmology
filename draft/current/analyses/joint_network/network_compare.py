"""Fitted network vs random rule-obeying networks (null/*.json) + LCDM rarity of the fitted contrasts."""
import json,glob,numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from jn_core import *
from null_nets import best
plt.rcParams.update({"font.size":12.5})
ST=[("link_len_med","median link length (mean void diameters)"),("deg_mean_V","links per void"),("deg_mean_C","links per mass"),
    ("quad_axis_cos_CoE","|cos| main axis vs CoE direction (1 = aligned)"),("frac_inside","share of bodies inside our sphere"),
    ("soft_rms","rule penalty (rms; 0 = all rules met)"),("left","leftover (× floor)")]
def val(o,k):
    if k in("soft_rms","left"): return o[k]
    v=o["props"][k]; return np.nan if v is None else v
out={}
files=sorted(glob.glob("null/*.json"))
for f in files:
    key=f[5:-5]; lv,sym,N=key.split("_"); N=int(N[1:]); d=json.load(open(f)); fit=d["fit"]; nul=d["null"]
    b=best(lv,sym,N); ns=lcdm_nsigma(b["kinds"],b["X"],b["M"],b["Rv"],b["delta"])
    row=dict(fit_left=fit["left"],n_null=len(nul),null_ok=float(np.mean([o["ok"] for o in nul])),
             null_rules_met=float(np.mean([o["soft_max"]<=0.1 for o in nul])),lcdm_nsigma_median=float(np.median(ns["nsigma"])),
             lcdm_nsigma_max=float(np.max(ns["nsigma"])),lcdm_log10p_total=ns["log10p_total"],stats={})
    for k,lab in ST:
        x=np.array([val(o,k) for o in nul],float); x=x[np.isfinite(x)]; fv=val(dict(fit,**{"props":fit["props"]}),k) if k!="left" else fit["left"]
        if k=="soft_rms": fv=fit["soft_rms"]
        lo=float(np.mean(x<=fv)); row["stats"][k]=dict(fit=float(fv),null_median=float(np.median(x)),null_p05=float(np.percentile(x,5)),null_p95=float(np.percentile(x,95)),
                                                      frac_null_le_fit=lo,two_sided_p=float(min(1,2*min(lo,1-lo+np.mean(x==fv)))) )
    out[key]=row
    print(key,{k:(round(v["fit"],3),round(v["null_median"],3),round(v["frac_null_le_fit"],3)) for k,v in row["stats"].items()},
          "nsig med/max",round(row["lcdm_nsigma_median"]),round(row["lcdm_nsigma_max"]),"log10p",round(row["lcdm_log10p_total"]))
json.dump(out,open("network_compare.json","w"),indent=1)
# chart: one case per column? keep plain: stacked panels for the main case(s)
cases=[c for c in ("cont_m_N12","web_m_N24","strict_m_N24") if c in out]
show=["link_len_med","deg_mean_V","quad_axis_cos_CoE","frac_inside","left"]
fig,axs=plt.subplots(len(show),1,figsize=(10,3.0*len(show)+1.6))
col={"cont_m_N12":"#ff7f0e","web_m_N24":"#9467bd","strict_m_N24":"#d62728"}
lab={"cont_m_N12":"continuity only, N=12","web_m_N24":"continuity + strict web rules, N=24","strict_m_N24":"+ density-trend cap, N=24"}
for ax,k in zip(axs,show):
    for c in cases:
        nul=json.load(open(f"null/{c}.json"))["null"]; x=np.array([val(o,k) for o in nul],float); x=x[np.isfinite(x)]
        fv=out[c]["stats"][k]["fit"]
        if k=="left":
            bins=np.logspace(-2,3.5,45); ax.set_xscale("log")
        else: bins=30
        ax.hist(x,bins=bins,histtype="step",lw=2,color=col[c],label=f"random networks, {lab[c]}")
        ax.axvline(fv,color=col[c],lw=2.5,ls="--",label=f"fitted layout, {lab[c]}")
    if k=="left": ax.axvline(1,color="green",lw=2,label="1× target")
    if k in("frac_inside","quad_axis_cos_CoE"): ax.set_xlim(0,1)
    ax.set_xlabel(dict(ST)[k]); ax.set_ylabel("count")
h,l=axs[-1].get_legend_handles_labels()
fig.legend(h,l,loc="lower center",ncol=1,frameon=False,fontsize=11)
fig.suptitle("Fitted layouts (dashed) vs 100 random networks obeying the same rules (lines)",fontsize=14)
fig.tight_layout(rect=(0,0.15,1,0.97)); fig.savefig("network_vs_random.png",dpi=100); plt.close(fig)
