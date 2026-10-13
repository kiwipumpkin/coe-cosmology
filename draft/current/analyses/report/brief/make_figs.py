"""Paper-size figures for the brief, from ../collected.pkl and the styles in ../build_pdf.py (no new fits)."""
import os,sys,shutil,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("../build_pdf.py").read(); src=src[:src.index("with PdfPages(")]
src=src.replace("s=30+260*min(abs(M[i])/2e52,1.5)","s=SZ*(30+260*min(abs(M[i])/2e52,1.5))").replace('"k+",ms=14,mew=2','"k+",ms=8,mew=1.2').replace('fontsize=11,fontweight="bold"','fontsize=6.5,fontweight="bold"').replace("edgecolors=\"k\",lw=0.5","edgecolors=\"k\",lw=0.3")
SZ=0.3
os.chdir(".."); exec(src); os.chdir("brief")
plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,"savefig.dpi":300,"pdf.fonttype":42})
OUT="figures"
def save(fig,name): fig.savefig(f"{OUT}/{name}.pdf",bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.png",bbox_inches="tight",dpi=200); plt.close(fig)
# ---- Fig 1 fingerprint
th=TH; r=R/FLOOR; deg=np.degrees(th); m=np.arange(1,13)
a=np.array([2*np.mean(r*np.cos(k*th)) for k in m]); b=np.array([2*np.mean(r*np.sin(k*th)) for k in m]); amp=np.hypot(a,b)
fig=plt.figure(figsize=(7.1,2.35)); gs=fig.add_gridspec(1,3,width_ratios=[1.6,1.05,0.95],wspace=0.32)
ax=fig.add_subplot(gs[0]); ax.plot(deg,r,"k",lw=1.3); ax.axhline(0,color="green",lw=1); ax.set_xlim(0,360); ax.set_xticks(range(0,361,90))
ax.set_xlabel(r"angle around the ring $\ell$ (deg)"); ax.set_ylabel(r"leftover ($\times10^{-5}$)"); ax.set_title("(a) leftover around the ring (rms 237×)")
bx=fig.add_subplot(gs[1]); bx.bar(m,amp,color=["#c0392b" if k in (3,4,5) else "0.6" for k in m]); bx.set_xticks([1,3,5,7,9,11]); bx.set_ylim(0,235)
bx.set_xlabel("mode $m$"); bx.set_ylabel(r"amplitude ($\times10^{-5}$)"); bx.set_title("(b) the same, by mode")
px=fig.add_subplot(gs[2],projection="polar"); base=1.0; sc=0.6/np.abs(r).max(); t=np.append(th,th[0]); rr=np.append(base+sc*r,base+sc*r[0])
px.plot(np.linspace(0,2*np.pi,400),np.full(400,base),color="0.6",lw=0.8,ls="--")
px.fill_between(t,base,rr,where=rr>=base,color="#c0392b",alpha=0.5,lw=0); px.fill_between(t,base,rr,where=rr<base,color="#2471a3",alpha=0.5,lw=0)
px.plot(t,rr,"k",lw=1); px.set_ylim(0,1.65); px.set_yticks([]); px.set_xticks(np.radians(range(0,360,90))); px.tick_params(labelsize=6.5,pad=0)
px.plot(0,0,"k+",ms=7); px.set_title("(c) drawn around the sky",pad=10)
save(fig,"fig1_fingerprint")
comp=[a[i]*np.cos(m[i]*th)+b[i]*np.sin(m[i]*th) for i in range(12)]
fig,axs=plt.subplots(3,4,figsize=(7.1,3.9),sharex=True)
wave_panels(axs,deg,comp,amp,(3,4,5),fs=7.5,lw_big=1.3,lw=1.0,barw=0.4,ybar=0.04)
for axi in axs[-1]: axi.set_xlabel(r"$\ell$ (deg)"); axi.set_xticks([0,180,360])
for axi in axs[:,0]: axi.set_ylabel(r"$\times10^{-5}$")
fig.tight_layout(h_pad=0.6,w_pad=0.6); save(fig,"fig1b_waves")
# ---- curves
def curves(ax,keys,xmax=20,ymin=0.01,logx=False,leg=True,title=None,ncol=1,legloc="upper right"):
    for k in keys:
        c,ls,mk=STY[k]; N=sorted(S[k]["curve"]); y=[S[k]["curve"][n] for n in N]
        if ls=="none": ax.plot(N,y,mk,color=c,ms=8,ls="none",label=SHORT[k])
        else: ax.plot(N,y,ls=ls,marker=None if logx else mk,color=c,lw=1.4,ms=3,label=SHORT[k])
    ax.axhline(1,color="green",lw=1.5); ax.axhline(237.3,color="0.5",ls=":",lw=0.8)
    if logx: ax.set_xscale("log"); ax.set_xlim(1,7e4); ax.set_xlabel("bodies placed")
    else: ax.set_xlim(0.5,xmax+0.5); ax.set_xticks([t for t in [1,2,4,6,8,10,12,14,16,18,20] if t<=xmax]); ax.set_xlabel("number of added bodies $N$")
    ax.set_yscale("log"); ax.set_ylim(ymin,2000 if logx else 500); ax.set_ylabel(r"best leftover ($\times$ floor)")
    if title: ax.set_title(title)
    if leg: ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.17),frameon=False,ncol=1,handlelength=2.6,fontsize=7)
allk=["out_masses","out_mixed","out_voids","out_asym","partners","nowall_smooth","nowall_open","in_masses","in_mixed","in_voids","in_cap1e53","in_nocomp","in_surr","all_open"]
fig,ax=plt.subplots(figsize=(7.1,5.6)); curves(ax,allk,leg=False)
ax.text(19.9,1.12,"flat target 1×",color="green",ha="right",fontsize=7); ax.text(19.9,260,"start 237×",color="0.4",ha="right",fontsize=7)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.11),frameon=False,handlelength=2.6,ncol=3,fontsize=7); save(fig,"fig2_overview")
fig,axs=plt.subplots(2,2,figsize=(7.1,9.0))
curves(axs[0,0],["out_masses","out_mixed","out_voids","out_asym","partners"],title="(a) outside our sphere, CoE wall in place")
curves(axs[0,1],["out_mixed","nowall_smooth","nowall_open","all_open"],xmax=12,title="(b) removing the walls")
curves(axs[1,0],["in_masses","in_mixed","in_voids","all_open"],xmax=12,title="(c) inside our sphere, and all open")
curves(axs[1,1],["in_masses","in_cap1e53","in_nocomp","in_surr"],xmax=12,title="(d) inside: mass cap and compensation")
fig.tight_layout(h_pad=0.8); save(fig,"fig3_families")
fig,ax=plt.subplots(figsize=(3.4,3.6)); curves(ax,["seqA_deep","seqA_shallow","seqB_cmb_aware","seqB_foam_rim","seqB_random","seqB_shallow"],logx=True,ymin=10,leg=False)
ax.axhline(88.8,color="0.6",ls=":",lw=0.8); ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.2),frameon=False,ncol=1,fontsize=6.5); save(fig,"fig4_sequential")
# ---- examples
def ex_row(axs,key,e,lab):
    wall=S[key]["family"] in ("outside","seq")
    topview(axs[0],e,wall); axs[0].tick_params(labelsize=6); axs[0].set_xlabel("Gpc",fontsize=7); axs[0].set_ylabel("Gpc",fontsize=7)
    axs[0].set_title(f"{lab}. {NAME9[key]}: {e['left']:.2f}×\nN = {e['N']}: {mixstr(e['kinds'])}"+(f" + {len(e['partners']['R'])} partners" if 'partners' in e else ""),fontsize=7.5)
    a=axs[1]; mm=np.arange(1,13); a.bar(mm,D["start_amps"],color="none",edgecolor="0.5",lw=0.6,ls="--")
    a.bar(mm,np.maximum(e["amps"],1e-4),color="#7d3c98"); a.bar(13.8,max(e["hi"],1e-4),color="0.55",width=1.3)
    a.axhline(1,color="green",lw=1.2); a.set_yscale("log"); a.set_ylim(1e-3,400); a.set_xticks([1,3,5,7,9,11,13.8]); a.set_xticklabels(["1","3","5","7","9","11","13+"])
    a.set_xlabel("mode $m$"); a.set_ylabel(r"$\times$ floor"); a.set_title("leftover by mode (dashed = start)",fontsize=7.5)
    c=axs[2]; c.plot(np.degrees(TH),e["resid"]/FLOOR,color="#7d3c98",lw=0.9); c.axhspan(-1,1,fc=(0.0,0.5,0.0,0.25),ec=(0.0,0.5,0.0,0.9),lw=0.6,zorder=1.5)
    c.set_xlim(0,360); c.set_xticks(range(0,361,90)); c.set_xlabel(r"$\ell$ (deg)"); c.set_ylabel(r"$\times10^{-5}$"); c.set_title(f"residual after fit: {e['left']:.2f}×",fontsize=7.5)
NAME9=dict(seqA_deep="sequential, near us",partners="giant partners",out_asym="outside, no mirror",nowall_smooth="outside, no wall (smooth limit)",
           in_masses="inside, masses only",out_mixed="outside, masses + deep voids",in_cap1e53="inside, masses, cap 1e53 kg",all_open="all open",in_mixed="inside, masses + deep voids")
def ex_fig(rows,name,start):
    n=len(rows); fig,axs=plt.subplots(n,3,figsize=(7.1,2.3*n),gridspec_kw=dict(width_ratios=[0.9,1,1.25]))
    for i,(k,e) in enumerate(rows): ex_row(axs[i],k,e,str(start+i))
    fig.tight_layout(h_pad=1.0,w_pad=1.0); save(fig,name)
# Chris: one global sequence over the three example figures, largest leftover first (unrounded), rows 1-3 / 4-6 / 7-9
ALL=[(k,S[k]["ex"]) for k in ["out_mixed","out_asym","partners","nowall_smooth","all_open","in_mixed","in_masses","in_cap1e53","seqA_deep"]]
ALL.sort(key=lambda r:-r[1]["left"])
for j,(k,e) in enumerate(ALL): print(j+1,k,e["N"],e["kinds"],repr(e["left"]))
ex_fig(ALL[0:3],"fig5_examples_outside",1); ex_fig(ALL[3:6],"fig6_examples_open_inside",4); ex_fig(ALL[6:9],"fig7_examples_more",7)
# ---- core fully open (../../core_open, collected by core_open/collect_core_open.py): examples 10-13 + best-vs-N
import pickle as _pk
CO=_pk.load(open("core_open/core_open.pkl","rb"))
CONAME=dict(out="core open: masses + deep voids",asym="core open: no mirror",ep="core open: giant partners",ao="core open: all open")
for k_ in CO: S["co_"+k_]=dict(family="outside",ex=CO[k_]["ex"]); NAME9["co_"+k_]=CONAME[k_]   # 'outside' -> core edge drawn dashed
ROWS=sorted([("co_"+k_,CO[k_]["ex"]) for k_ in CO],key=lambda r:-r[1]["left"])
for j,(k_,e_) in enumerate(ROWS): print(10+j,k_,e_["N"],e_["kinds"],repr(e_["left"]))
ex_fig(ROWS,"fig_examples_core_open",10)
fig,ax=plt.subplots(figsize=(7.1,5.6))   # full page width (figure*)
for k_,old,lab in [("out","out_mixed","outside, masses + deep voids"),("asym","out_asym","outside, no mirror"),("ep","partners","giant partners"),("ao","all_open","all open")]:
    c_=STY[old][0]; N_=sorted(S[old]["curve"]); ax.plot(N_,[S[old]["curve"][n] for n in N_],ls=":",color=c_,lw=1.1)
    N_=sorted(CO[k_]["curve"]); ax.plot(N_,[CO[k_]["curve"][n] for n in N_],"-o",color=c_,lw=1.5,ms=3,label=lab+", core open")
ax.plot([],[],ls=":",color="0.4",lw=1.1,label="same scenario, earlier runs\n(CoE wall in place; all open:\nreach limited to 43 Gpc)")
ax.axhline(1,color="green",lw=1.5); ax.axhline(237.3,color="0.5",ls=":",lw=0.8)
ax.text(0.7,1.15,"flat target 1×",color="green",ha="left",fontsize=8.5); ax.text(12.4,260,"start 237×",color="0.4",ha="right",fontsize=8.5)
ax.set_yscale("log"); ax.set_ylim(0.03,500); ax.set_xlim(0.5,12.5); ax.set_xticks([1,2,4,6,8,10,12]); ax.set_xlabel("number of added bodies $N$"); ax.set_ylabel(r"best leftover ($\times$ floor)")
ax.set_title("CoE core fully open")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.11),frameon=False,ncol=3,handlelength=2.6,fontsize=8); save(fig,"fig_core_open_best_vs_N")
# ---- reuse existing figures
CS="/workspace/clump_scan"
for s,d in [("inside_visible/smooth_field/tradeoff_leftover_vs_rms_delta.png","fig8b_smooth_tradeoff.png"),
            ("realistic_counts/objects_needed_by_scenario.png","fig9a_realistic_counts.png"),("realistic_counts/mixed/leftover_vs_ratio_deviation.png","fig9b_mixed_ratio.png"),
            ("density_vs_scale/density_vs_scale_fits_radius_v2.png","fig10_density_vs_scale_v2.png"),("density_vs_scale/distribution/density_vs_size_heatmap.png","fig11a_heatmap.png"),
            ("density_vs_scale/distribution/density_hist_by_size.png","fig11b_hist.png")]:
    shutil.copy(f"{CS}/{s}",f"{OUT}/{d}")
import subprocess; subprocess.run([sys.executable,"fig_smooth.py"],check=True)
print("ok")
