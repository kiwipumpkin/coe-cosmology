import pickle,numpy as np,matplotlib,textwrap
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams.update({"font.size":12,"axes.titlesize":13})
D=pickle.load(open("collected.pkl","rb")); S=D["S"]; R=D["R"]; FLOOR=D["FLOOR"]; TH=D["TH"]; C_G=D["C_G"]; A_IN=D["A_IN_G"]; RV=D["RV"]
CS="/workspace/clump_scan"; W,H=11,8.5
STY={"out_masses":("#c0392b","--","s"),"out_mixed":("#7d3c98","--","s"),"out_voids":("#5dade2","--","s"),"out_asym":("#17a589","--","P"),
     "partners":("#444444","none","*"),"nowall_smooth":("#e67e22","--","D"),"nowall_open":("#f5b041",":","D"),
     "in_masses":("#c0392b","-","o"),"in_mixed":("#7d3c98","-","o"),"in_voids":("#5dade2","-","o"),"in_cap1e53":("#a04000","-","^"),
     "in_nocomp":("#f1948a","-","v"),"in_surr":("#7b241c","-","D"),"all_open":("#1a5276","-.","D"),
     "seqA_deep":("#c0392b","-",None),"seqA_shallow":("#e59866","-",None),"seqB_cmb_aware":("#7d3c98","-",None),"seqB_foam_rim":("#5dade2","-",None),
     "seqB_random":("#7f8c8d","-",None),"seqB_shallow":("#d7bde2","-",None)}
SHORT={"out_masses":"outside: masses only","out_mixed":"outside: masses + deep voids","out_voids":"outside: deep voids only","out_asym":"outside: no mirror layouts",
       "partners":"giant partners (+ ordinary bodies)","nowall_smooth":"outside, no wall: smooth-density limit","nowall_open":"outside, no wall: fully open",
       "in_masses":"inside: masses only","in_mixed":"inside: masses + deep voids","in_voids":"inside: deep voids only","in_cap1e53":"inside: masses, cap 1e53 kg",
       "in_nocomp":"inside: masses, no compensation","in_surr":"inside: masses, compensation from outside","all_open":"all open (smooth-density limit)",
       "seqA_deep":"A: near fix, deep voids + clusters","seqA_shallow":"A: near fix, shallow voids","seqB_cmb_aware":"B: growth, CMB-aware",
       "seqB_foam_rim":"B: growth, foam rule","seqB_random":"B: growth, random","seqB_shallow":"B: growth, CMB-aware, shallow voids"}
def refs(ax,x_txt,four=True):
    ax.axhline(1,color="green",lw=2.2); ax.text(x_txt,1.12,"CMB flat target 1×",color="green",fontsize=11,ha="right",va="bottom")
    ax.axhline(237.3,color="0.5",ls=":",lw=1.2); ax.text(x_txt,237.3*1.12,"start 237×",color="0.35",fontsize=11,ha="right",va="bottom")
    if four: ax.axhline(88.8,color="0.65",ls=":",lw=1); ax.text(x_txt,88.8*1.12,"four clumps 89×",color="0.45",fontsize=10,ha="right",va="bottom")
def curve_page(pdf,keys,title,note,xmax=20,ymin=0.01,logx=False):
    fig,ax=plt.subplots(figsize=(W,H))
    for k in keys:
        c,ls,mk=STY[k]; N=sorted(S[k]["curve"]); y=[S[k]["curve"][n] for n in N]
        if ls=="none": ax.plot(N,y,mk,color=c,ms=15,ls="none",label=SHORT[k])
        else: ax.plot(N,y,ls=ls,marker=mk if not logx else None,color=c,lw=2.4,ms=6,label=SHORT[k])
    if logx:
        ax.set_xscale("log"); ax.set_xlim(1,7e4); refs(ax,6.5e4)
    else:
        ax.set_xlim(0.5,xmax+0.5); ax.set_xticks([1,2,4,6,8,10,12,14,16,18,20][:sum(1 for t in [1,2,4,6,8,10,12,14,16,18,20] if t<=xmax)]); refs(ax,xmax+0.4)
    ax.set_yscale("log"); ax.set_ylim(ymin,2000 if logx else 400)
    yt=[t for t in [0.01,0.1,1,10,100,1000] if t>=ymin]; ax.set_yticks(yt); ax.set_yticklabels([("%g"%t) for t in yt])
    ax.set_xlabel("number of bodies placed (log scale)" if logx else "number of added bodies N"); ax.set_ylabel("best leftover (× 10⁻⁵ floor, log scale)")
    ax.set_title(title,fontsize=15)
    ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.1),ncol=2,frameon=False,fontsize=11)
    fig.text(0.5,0.015,note,ha="center",fontsize=10.5,color="0.3",wrap=True)
    fig.tight_layout(rect=(0,0.04,1,1)); pdf.savefig(fig); plt.close(fig)
def topview(ax,e,old_wall=True):
    lim=e["zoom"]; t=np.linspace(0,2*np.pi,400)
    ax.set_aspect("equal"); ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    if "partners" in e and abs(e["partners"]["delta_us"])>0.05:
        ax.fill(RV*np.cos(t),RV*np.sin(t),color="#c0392b" if e["partners"]["delta_us"]>0 else "#2471a3",alpha=0.10)
    ax.plot(RV*np.cos(t),RV*np.sin(t),"k",lw=1.3)
    if old_wall: ax.plot(C_G[0]+A_IN*np.cos(t),C_G[1]+A_IN*np.sin(t),color="0.65",ls="--",lw=1)
    k=np.array(list(e["kinds"])); X=e["X"]; M=e["M"]; Rv=e["Rv"]
    for i in range(len(k)):
        if k[i]=="C": ax.scatter(X[i,0],X[i,1],s=30+260*min(abs(M[i])/2e52,1.5),c="#c0392b",edgecolors="k",lw=0.5,zorder=3)
        else: ax.add_patch(plt.Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.4))
    if "partners" in e:
        P=e["partners"]
        for i in range(len(P["R"])):
            c="#c0392b" if P["M"][i]>0 else "#2471a3"
            ax.add_patch(plt.Circle(P["X"][i,:2],P["R"][i],facecolor=c,alpha=0.22,edgecolor=c,lw=2)); ax.text(P["X"][i,0],P["X"][i,1],f"P{i+1}",ha="center",va="center",fontsize=11,fontweight="bold")
    ax.plot(0,0,"k+",ms=14,mew=2); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")
def mixstr(k):
    c=k.count("C"); v=k.count("V"); s=[]
    if c: s.append(f"{c} mass"+("es" if c>1 else ""))
    if v: s.append(f"{v} void"+("s" if v>1 else ""))
    return " + ".join(s)
def example_page(pdf,key,e,note,title=None):
    wall=S[key]["family"] in ("outside","seq")
    fig=plt.figure(figsize=(W,H)); gs=fig.add_gridspec(2,2,width_ratios=[1.05,1],height_ratios=[1,1],left=0.06,right=0.98,top=0.86,bottom=0.08,wspace=0.22,hspace=0.4)
    ax=fig.add_subplot(gs[:,0]); topview(ax,e,wall)
    ax.set_title("Side view (us = +; black circle = our 14.3 Gpc sphere;\n"+("grey dashed = CoE-interior wall; " if wall else "")+"red = masses sized by mass;\nblue = voids at true radius; out-of-plane bodies projected)",fontsize=10.5)
    a=fig.add_subplot(gs[0,1]); m=np.arange(1,13)
    a.bar(m,np.maximum(e["amps"],1e-4),color="#7d3c98"); a.bar(13.6,max(e["hi"],1e-4),color="0.55",width=1.2)
    a.bar(m,D["start_amps"],color="none",edgecolor="0.4",lw=1,ls="--",label="start")
    a.axhline(1,color="green",lw=2); a.set_yscale("log"); a.set_ylim(1e-3,400); a.set_xticks(list(m)+[13.6]); a.set_xticklabels([str(i) for i in m]+["13+"],fontsize=10)
    a.set_title("Leftover per mode m (bars; dashed outline = start;\ngrey '13+' = all higher modes combined; green = 1×)",fontsize=10.5); a.set_ylabel("× floor",fontsize=10)
    deg=np.degrees(TH)
    c=fig.add_subplot(gs[1,1]); c.plot(deg,e["resid"]/FLOOR,color="#7d3c98",lw=1.8); c.axhline(1,color="green",ls="--",lw=1); c.axhline(-1,color="green",ls="--",lw=1); c.axhline(0,color="0.85",lw=1)
    c.set_title(f"Residual after the fit ({e['left']:.2f}×; green dashed = ±1)",fontsize=10.5); c.set_xlim(0,360); c.set_xticks(range(0,361,90)); c.set_xlabel("angle around the ring ℓ (deg)",fontsize=10); c.set_ylabel("10⁻⁵",fontsize=10)
    fig.suptitle(title or f"{S[key]['label']}\nExample: N = {e['N']} ({mixstr(e['kinds'])}{', '+str(len(e['partners']['R']))+' giant partners' if 'partners' in e else ''}), leftover {e['left']:.2f}×",fontsize=14)
    fig.text(0.5,0.015,note+f"   [re-scored {e['left']:.4f}× = saved {e['saved']:.4f}×]",ha="center",fontsize=9.5,color="0.3")
    pdf.savefig(fig); plt.close(fig)
def text_page(pdf,title,lines,fs=13):
    fig=plt.figure(figsize=(W,H)); fig.text(0.06,0.93,title,fontsize=20,weight="bold",va="top")
    y=0.82 if "\n" not in title else 0.78
    for ln in lines:
        for i,w in enumerate(textwrap.wrap(ln,115) or [""]):
            fig.text(0.07 if i==0 else 0.09,y,w,fontsize=fs); y-=0.038 if fs>=12 else 0.031
        y-=0.008 if fs>=12 else 0.004
    pdf.savefig(fig); plt.close(fig)
def image_page(pdf,title,imgs,caps,note=""):
    fig=plt.figure(figsize=(W,H)); n=len(imgs)
    for i,(p,cp) in enumerate(zip(imgs,caps)):
        ax=fig.add_axes([0.02+i*0.98/n,0.08,0.96/n-0.01,0.8]); ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(cp,fontsize=11)
    fig.suptitle(title,fontsize=15); fig.text(0.5,0.03,note,ha="center",fontsize=10.5,color="0.3",wrap=True)
    pdf.savefig(fig); plt.close(fig)
def wave_panels(axs,deg,comp,amp,big,fs=11,lw_big=2,lw=1.4,barw=0.32,ybar=0.86):
    ref=amp[2]
    for i,axi in enumerate(axs.flat):
        k=i+1; red=k in big; c="#c0392b" if red else "0.35"
        axi.plot(deg,comp[i],color=c,lw=lw_big if red else lw); axi.axhline(0,color="0.85",lw=0.8)
        A=amp[i]; axi.set_ylim(-2.1*A,2.1*A); axi.set_xlim(0,360); axi.set_xticks(range(0,361,90))
        axi.text(0.02,0.97,f"m = {k}  ({A:.0f}×)" if A>=10 else f"m = {k}  ({A:.1f}×)",transform=axi.transAxes,va="top",fontsize=fs,
                 weight="bold" if red else "normal",color=c)
        x0=0.98-barw; y0=ybar; h=0.08
        axi.add_patch(plt.Rectangle((x0,y0),barw,h,transform=axi.transAxes,fc="0.92",ec="0.6",lw=0.5,clip_on=False))
        axi.add_patch(plt.Rectangle((x0,y0),barw*max(A/ref,0.004),h,transform=axi.transAxes,fc=c,ec="none",clip_on=False))
        axi.tick_params(axis="y",labelsize=fs-3)
def fingerprint_pages(pdf):
    th=TH; r=R/FLOOR; deg=np.degrees(th); m=np.arange(1,13)
    a=np.array([2*np.mean(r*np.cos(k*th)) for k in m]); b=np.array([2*np.mean(r*np.sin(k*th)) for k in m])
    amp=np.hypot(a,b); comp=[a[i]*np.cos(m[i]*th)+b[i]*np.sin(m[i]*th) for i in range(12)]
    big=[3,4,5]; col=["#c0392b" if k in big else "0.6" for k in m]
    fig=plt.figure(figsize=(W,H)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.25],left=0.07,right=0.97,top=0.86,bottom=0.14,hspace=0.75,wspace=0.35)
    ax=fig.add_subplot(gs[0,:]); ax.plot(deg,r,color="k",lw=2); ax.axhline(0,color="0.8",lw=1)
    ax.axhspan(-1,1,color="green",alpha=0.6)
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle around the ring ℓ (deg)"); ax.set_ylabel("× 10⁻⁵")
    ax.set_title("1. The leftover pattern around the ring before any bodies are added (rms 237×)\nthe flat target (±1×) is a green band around 0, too thin to see at this scale",fontsize=12.5)
    bx=fig.add_subplot(gs[1,0]); bx.bar(m,amp,color=col); bx.set_ylim(0,235); bx.set_xticks(m); bx.set_xlabel("mode m (number of waves around the ring)"); bx.set_ylabel("amplitude (× 10⁻⁵)")
    for k in big: bx.text(k,amp[k-1]+5,f"{amp[k-1]:.0f}",ha="center",fontsize=10,color="#c0392b")
    bx.set_title("2. The same pattern split into its waves (modes):\nm = 3, 4, 5 (red) carry almost all of it",fontsize=12)
    px=fig.add_subplot(gs[1,1],projection="polar"); base=1.0; sc=0.6/np.abs(r).max(); rr=base+sc*r
    t=np.append(th,th[0]); rr2=np.append(rr,rr[0])
    px.plot(np.linspace(0,2*np.pi,400),np.full(400,base),color="0.6",lw=1,ls="--")
    px.fill_between(t,base,rr2,where=rr2>=base,color="#c0392b",alpha=0.5,lw=0); px.fill_between(t,base,rr2,where=rr2<base,color="#2471a3",alpha=0.5,lw=0)
    px.plot(t,rr2,"k",lw=1.5); px.set_ylim(0,1.65); px.set_yticks([]); px.set_xticks(np.radians(range(0,360,45))); px.tick_params(labelsize=9)
    px.annotate("CoE direction",xy=(np.radians(195.2),1.6),fontsize=9,ha="center",color="0.3")
    px.plot(0,0,"k+",ms=12,mew=2)
    px.set_title("3. The fingerprint drawn around the sky\n(dashed = flat; red = above, blue = below)",fontsize=12,pad=22)
    fig.suptitle("The fingerprint: the leftover CMB pattern every scenario tries to cancel",fontsize=16)
    fig.text(0.5,0.02,"This leftover pattern is the fingerprint. Every scenario in this report is an attempt to cancel it by placing masses and voids.\n"
             "1× = flat: the pattern squeezed down to the green band.",ha="center",fontsize=12)
    pdf.savefig(fig); plt.close(fig)
    fig,axs=plt.subplots(4,3,figsize=(W,H),sharex=True)
    wave_panels(axs,deg,comp,amp,big,fs=11,lw_big=2,lw=1.4)
    for axi in axs[-1]: axi.set_xlabel("ℓ (deg)")
    for axi in axs[:,0]: axi.set_ylabel("× 10⁻⁵",fontsize=10)
    fig.suptitle("The fingerprint, wave by wave (m = 1…12)",fontsize=15)
    fig.text(0.5,0.015,"Each panel has its own scale; the numbers give the size, and the bar at top right shows it relative to m = 3 (full bar = 200×).\n"
             "Adding these 12 waves back together gives the curve on the previous page. m = 3, 4, 5 (red) are the big ones.",ha="center",fontsize=11)
    fig.tight_layout(rect=(0,0.06,1,0.95)); pdf.savefig(fig); plt.close(fig)
    return amp
def first_below(k):
    c=S[k]["curve"]; b=[n for n in sorted(c) if c[n]<1]; return b[0] if b else None
with PdfPages("all_scenarios.pdf") as pdf:
    text_page(pdf,"CoE shell model: canceling the leftover CMB residual\nwith added bodies — all scenarios",[
        "All the scenarios we have run, in one place: leftover vs number of added bodies, an example universe for each (side view), and the residual that remains.",
        "How to read the charts:",
        "• Leftover = the rms of the CMB pattern that is still unexplained after the fit, in units of the 10⁻⁵ flatness floor. The green line at 1× is the target ('CMB flat').",
        "• Start = 237×: the leftover after special relativity + the central mass, before any bodies are added. The four clumps reached 89× (grey dotted line).",
        "• Bodies: compact masses (red dots, sized by mass) and deep voids (blue circles, δ = −0.7, drawn at their true radius). 'Inside' = inside our 14.3 Gpc visible sphere; 'outside' = 14.3–43 Gpc from us.",
        "• 'Mode m' = how many waves fit around the ring (m = 1 dipole … m = 12); '13+' collects all finer ripple.",
        "• Side views look face-on at the ring (galactic) plane; bodies above or below it are projected onto it.",
        "• Every example was re-scored from the saved layout with the saved model code; the re-scored leftover matches the saved value in every case (shown at the bottom of each page).",
        "• No new fits were run for this report. Above ~5 bodies, many different layouts fit about equally well (overfitting), so the example layouts are illustrations, not unique answers."],fs=12.5)
    amp=fingerprint_pages(pdf); print("amps",np.round(amp,1),"saved",np.round(D["start_amps"],1),"recon rms",np.sqrt(np.mean((R/FLOOR)**2)))
    allk=["out_masses","out_mixed","out_voids","out_asym","partners","nowall_smooth","nowall_open","in_masses","in_mixed","in_voids","in_cap1e53","in_nocomp","in_surr","all_open"]
    curve_page(pdf,allk,"Overview: best leftover vs number of added bodies, every point-body scenario",
        "Solid = inside our sphere, dashed = outside, dash-dot = all open, stars = giant partners. The next pages split these into smaller groups.")
    curve_page(pdf,["out_masses","out_mixed","out_voids","out_asym","partners"],"Outside our sphere (14.3–43 Gpc, hard wall at the CoE interior)",
        "Masses alone stall near 4×; masses + deep voids reach 1× at 11 bodies (12 with mirror layouts forbidden); voids alone stay near 85×. Stars: giant partners (N = partners + ordinary bodies).")
    curve_page(pdf,["out_mixed","nowall_smooth","nowall_open","all_open"],"Removing the walls",
        "Removing the CoE-interior wall gets to 1× one body sooner (N = 10). Opening our visible sphere too ('all open') gets there at N = 8 — the bodies move inside our sphere.",xmax=12)
    curve_page(pdf,["in_masses","in_mixed","in_voids","all_open","nowall_smooth"],"Inside our visible sphere, and everything open",
        "Inside our sphere: 7 masses reach 1× (0.62×), 4 masses + 4 voids reach 0.18×. Caveat: for bodies between us and the CMB the model leaves out ISW/lensing effects.")
    curve_page(pdf,["in_masses","in_cap1e53","in_nocomp","in_surr"],"Inside our sphere: variants (mass cap and net-mass compensation)",
        "Raising the mass cap to 1e53 kg reaches 1× with only 4 masses; dropping compensation changes little; taking the compensating mass from outside our sphere stalls near 6×.",xmax=12)
    curve_page(pdf,["seqA_deep","seqA_shallow","seqB_cmb_aware","seqB_foam_rim","seqB_random","seqB_shallow"],"Sequential placement: one body at a time, never overlapping",
        "A: the near region (14.3–30 Gpc) fixed pair by pair (best seed). B: growth from empty until the region jams (best seed per rule; coarse 1.5 Gpc voids). None reaches 1×.",logx=True,ymin=10)
    order=["out_masses","out_mixed","out_voids","out_asym","partners","nowall_smooth","nowall_open","in_masses","in_mixed","in_voids","in_cap1e53","in_nocomp","in_surr","all_open","seqA_deep","seqA_shallow"]
    for k in order:
        example_page(pdf,k,S[k]["ex"],S[k]["note"])
        if k=="partners": example_page(pdf,k,S[k]["ex2"],S[k]["ex2_note"],title=f"{S[k]['label']}\nSecond example: best layout with a giant partner, leftover {S[k]['ex2']['left']:.2f}×")
    # sequential B page
    fig=plt.figure(figsize=(W,H)); ax=fig.add_axes([0.0,0.08,0.6,0.82]); ax.imshow(plt.imread(f"{CS}/sequential/b_example_slices.png")); ax.axis("off")
    ax.set_title("Example grown universes (slices, saved chart)",fontsize=11)
    a=fig.add_axes([0.68,0.18,0.30,0.55]); m=np.arange(1,len(S["seqB_cmb_aware"]["final_amps"])+1)
    a.bar(m-0.2,S["seqB_cmb_aware"]["final_amps"],width=0.4,color="#7d3c98",label="CMB-aware, when full")
    a.bar(m+0.2,S["seqB_foam_rim"]["final_amps"],width=0.4,color="#5dade2",label="foam rule, when full")
    a.axhline(1,color="green",lw=2); a.set_yscale("log"); a.set_ylim(0.1,2000); a.set_xticks(m); a.set_xlabel("mode m"); a.set_ylabel("× floor"); a.legend(frameon=False,fontsize=10)
    a.set_title("Leftover per mode at the end of growth",fontsize=11)
    fig.suptitle("Sequential B: grow from empty, one body at a time (example)",fontsize=15)
    fig.text(0.5,0.03,"Growth layouts (tens of thousands of bodies) were not saved body-by-body, so no ring residual is re-scored here;\nthe per-mode leftover is the saved final value.",ha="center",fontsize=10,color="0.3")
    pdf.savefig(fig); plt.close(fig)
    # field type
    F=[(f"{CS}/inside_visible/smooth_field/map_minimal_field_2Gpc.png","Field-type scenario: smallest smooth density field that reaches 1× (2 Gpc smoothing)","Instead of individual bodies, a smooth density field inside our sphere is solved for directly (net mass zero). Red = slightly overdense, blue = slightly underdense; rms about 1.2%."),
       (f"{CS}/inside_visible/smooth_field/tradeoff_leftover_vs_rms_delta.png","Smooth field: leftover vs how strong the field has to be","Leftover falls from 237× to 1× as the field's rms density contrast goes from about 0.5% to 1.0–1.4% (depending on smoothing)."),
       (f"{CS}/realistic_counts/objects_needed_by_scenario.png","Realistic counts: how many extra clusters / voids a 1× solution corresponds to","Counts are 'extra objects beyond normal' (a rearrangement of the normal cosmic web), per scenario, compared with how many each region holds anyway."),
       (f"{CS}/realistic_counts/mixed/leftover_vs_ratio_deviation.png","Mixed cluster:void ratio: leftover vs how far the local ratio may differ from the cosmic mean","Inside our sphere (mirror allowed) about ±6.5% reaches 1×; inside + outside with no mirror needs somewhere between ±20% and ±60%.")]
    for p,t,n in F: image_page(pdf,t,[p],[""],n)
    lines=[]
    for k in order+["seqB_cmb_aware","seqB_foam_rim","seqB_random","seqB_shallow"]:
        c=S[k]["curve"]; fb=first_below(k); bN=min(c,key=c.get)
        lines.append(f"• {SHORT[k]}: "+(f"reaches 1× at N = {fb} ({c[fb]:.2f}×); best {c[bN]:.3g}× at N = {bN}." if fb else f"never reaches 1×; best {c[bN]:.3g}× at N = {bN}."))
    lines.append("• Field-type: a smooth field of ~1–1.4% rms reaches 1×; a mixed cluster:void ratio needs ~±6.5% (inside) or ±20–60% (no mirror). See the preceding four pages.")
    text_page(pdf,"Takeaways: fewest bodies to reach 1× (or best reached)",lines,fs=10.5)
print("ok")
