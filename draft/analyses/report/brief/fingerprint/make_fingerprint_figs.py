"""Brief figures 1-3 + Table I (fingerprint). Adapted from /workspace/cosmology-course/notes/cmb-modes-sheet.py (copied here,
original untouched). Numbers: cmb-cancel-story.json (copy) and the model S_model/resid from cmb-cancel.py (copy, read-only data)."""
import json,pickle,os,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Wedge
HERE=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(HERE,"..","figures")
plt.rcParams.update({"font.size":8,"axes.titlesize":8,"savefig.dpi":300,"pdf.fonttype":42})
d=json.load(open(os.path.join(HERE,"cmb-cancel-story.json"))); c=np.array(d["c_m"]); phi=d["params"][4]
cG=np.array(d["lnG_m"]); cS=np.array(d["lnInvS_m"])
src=open(os.path.join(HERE,"cmb-cancel_COPY.py")).read().split('if __name__')[0]; g={}; exec(src,g)
S=g["S_model"](*d["params"])[0]; R=g["resid"](g["G_FULL"],S)
M=12; A=np.abs(c); pk=[(phi+(0 if c[m-1]>0 else 180/m))%(360/m) for m in range(1,M+1)]; vf=(c**2/2)/np.mean(R**2)
# consistency with the clump_scan residual
D=pickle.load(open(os.path.join(HERE,"..","..","collected.pkl"),"rb"))
chk=dict(rms_sheet=float(np.sqrt(np.mean(R**2))),rms_clump=float(np.sqrt(np.mean(D["R"]**2))),maxabs=float(np.abs(R).max()),
         corr=float(np.corrcoef(R,D["R"])[0,1]),maxdiff=float(np.abs(R-D["R"]).max()),
         amps_sheet_x=list(np.round(A/1e-5,1)),amps_clump_x=list(np.round(D["start_amps"],1)),
         cdiff=list(np.round((cG-cS-c)/1e-5,3)),sum_vf_1_12=float(vf.sum()),vf35=float(vf[2:5].sum()))
json.dump(chk,open(os.path.join(HERE,"consistency.json"),"w"),indent=1); print(chk)
def save(fig,n): fig.savefig(f"{FIG}/{n}.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/{n}.png",bbox_inches="tight",dpi=200); plt.close(fig)
# ---------- Fig 1 geometry (side view = ring plane)
CG=np.array([-56.17,-15.22]); DCoE=np.hypot(*CG); RIN,ROUT,RV=43.9,72.5,14.3; u=CG/DCoE
fig,axs=plt.subplots(1,2,figsize=(7.1,3.7),gridspec_kw=dict(width_ratios=[1.25,1]))
ax=axs[0]
ax.add_patch(Wedge(CG,ROUT,0,360,width=ROUT-RIN,color="0.85",lw=0))
ax.add_patch(Circle(CG,RIN,fill=False,color="0.4",lw=0.8)); ax.add_patch(Circle(CG,ROUT,fill=False,color="0.4",lw=0.8))
ax.plot(*CG,"o",color="#c0392b",ms=5); ax.text(CG[0],CG[1]-6,"CoE",ha="center",color="#c0392b")
ax.add_patch(Circle((0,0),RV,fc="#d6eaf8",ec="#2471a3",lw=1.2)); ax.plot(0,0,"k+",ms=6,mew=1.2)
ax.plot([CG[0],0],[CG[1],0],"k--",lw=0.7); ax.text(*(CG*0.5+np.array([0,4])),"58.2 Gpc",ha="center",fontsize=7,rotation=np.degrees(np.arctan2(CG[1],CG[0]))-180)
ax.annotate("",xy=CG+RIN*np.array([0,1]),xytext=CG,arrowprops=dict(arrowstyle="->",lw=0.7,color="0.3")); ax.text(CG[0]+1.5,CG[1]+22,"43.9",fontsize=7,color="0.3")
ax.annotate("",xy=CG+ROUT*np.array([-np.sin(0.6),np.cos(0.6)]),xytext=CG,arrowprops=dict(arrowstyle="->",lw=0.7,color="0.3")); ax.text(CG[0]-27,CG[1]+48,"72.5",fontsize=7,color="0.3")
ax.text(CG[0],CG[1]-58,"shell (matter)",ha="center",fontsize=7.5,color="0.3"); ax.text(CG[0]-14,CG[1]-22,"CoE interior",ha="center",fontsize=7,color="0.4")
ax.text(0,-21,"our visible sphere",fontsize=7,color="#2471a3",ha="center",va="top")
ax.add_patch(plt.Rectangle((-17,-17),34,34,fill=False,ec="k",lw=0.6,ls=":"))
ax.set_xlim(-135,22); ax.set_ylim(-92,62); ax.set_aspect("equal"); ax.set_xlabel("x (Gpc)"); ax.set_ylabel("y (Gpc)")
ax.set_title("(a) side view, to scale: CoE, shell and our sphere")
ax=axs[1]; L=24
ax.add_patch(Circle(CG,ROUT,fc="0.85",ec="0.4",lw=0.8)); ax.add_patch(Circle(CG,RIN,fc="white",ec="0.4",lw=0.8))
ax.add_patch(Circle((0,0),RV,fill=False,ec="#2471a3",lw=2)); ax.plot(0,0,"k+",ms=9,mew=1.5); ax.text(1,-2.6,"us",fontsize=7.5)
ax.plot([u[0]*L*1.5,-u[0]*L*1.5],[u[1]*L*1.5,-u[1]*L*1.5],"k--",lw=0.7)
ax.text(-19.5,1.2,"to CoE\n$\\ell=195^\\circ$",fontsize=7,ha="center")
ax.text(-u[0]*20,-u[1]*20+1.5,"$\\ell=15^\\circ$",fontsize=7,ha="center")
for l in (0,90,180,270):
    t=np.radians(l); ax.plot([RV*np.cos(t)],[RV*np.sin(t)],"o",color="#2471a3",ms=2.5)
    ax.text(11.6*np.cos(t),11.6*np.sin(t),f"{l}°",ha="center",va="center",fontsize=6.5,color="#2471a3")
ax.annotate("",xy=(RV*np.cos(np.radians(120)),RV*np.sin(np.radians(120))),xytext=(0,0),arrowprops=dict(arrowstyle="->",lw=0.8,color="#2471a3"))
ax.text(-6.2,6.0,"14.3 Gpc",fontsize=7,color="#2471a3",ha="right")
ax.text(10.5,-15.5,"CMB ring",fontsize=7.5,color="#2471a3")
ax.text(-18.5,-17,"inner edge\nof shell",fontsize=6.5,color="0.35",ha="center"); ax.text(15.5,10.5,"outer edge\nof shell",fontsize=6.5,color="0.35",ha="center")
ax.set_xlim(-L,L); ax.set_ylim(-L,L); ax.set_aspect("equal"); ax.set_xlabel("x (Gpc)"); ax.set_ylabel("y (Gpc)")
ax.set_title("(b) zoom on our sphere (dotted box in a)")
fig.tight_layout(w_pad=1.5); save(fig,"fig1_geometry")
# ---------- Fig 2 source of the residual
m=np.arange(1,M+1); fig,axs=plt.subplots(1,2,figsize=(7.1,4.2),sharey=True)
ax=axs[0]; w=0.38
ax.bar(m-w/2,np.abs(cG),w,color="#2c3e50",label=r"gravity (shell): $|(\ln G)_m|$")
ax.bar(m+w/2,np.abs(cS),w,color="#e67e22",label=r"special relativity + central mass: $|(\ln 1/S)_m|$")
ax.set_yscale("log"); ax.set_ylim(1e-8,0.3); ax.set_xticks(m); ax.set_xlabel("mode $m$"); ax.set_ylabel("size of the mode")
ax.legend(frameon=False,loc="upper right",fontsize=6.8); ax.set_title("(a) the two large patterns, mode by mode")
ax=axs[1]; ax.bar(m,A,0.6,color=["#c0392b" if k in (3,4,5) else "0.55" for k in m])
ax.axhline(1e-5,color="green",lw=1.3); ax.text(12.4,0.06,"green line: flat floor $10^{-5}$",color="green",fontsize=7,ha="right")
ax.set_xticks(m); ax.set_xlabel("mode $m$"); ax.set_title(r"(b) what is left: $|c_m| = |(\ln G)_m-(\ln 1/S)_m|$")
fig.tight_layout(w_pad=1.0); save(fig,"fig2_source")
# ---------- Fig 3 polar modes (Chris's preferred illustration)
Lr=np.radians(np.arange(0,360.5,0.5)); th=Lr-np.radians(phi); cols=list(plt.cm.tab10.colors)+["k","goldenrod"]
def polar(ax,f,col,title,k=0.6,lw=1.4):
    rad=1.0+f/np.max(np.abs(f))*k
    ax.plot(Lr,rad,color=col,lw=lw); ax.plot(Lr,np.ones_like(Lr),color="gray",lw=0.5,ls=":")
    ax.plot([np.radians(phi)]*2,[0,1.7],color="k",lw=0.6,ls="--")
    ax.set_ylim(0,1.7); ax.set_yticklabels([]); ax.set_thetagrids([0,90,180,270]); ax.tick_params(labelsize=6,pad=-2); ax.set_title(title,fontsize=6.8,pad=9,linespacing=1.15)
fig=plt.figure(figsize=(7.1,9.3)); gs=fig.add_gridspec(5,4,height_ratios=[1,1,1,0.12,1.9],hspace=1.0,wspace=0.45)
for k in range(1,M+1):
    ax=fig.add_subplot(gs[(k-1)//4,(k-1)%4],projection="polar")
    polar(ax,c[k-1]*np.cos(k*th),cols[k-1],f"m = {k}  (period {360/k:.1f}°)\n$|c_{{{k}}}|$ = {A[k-1]:.2e} ({A[k-1]/1e-5:.3g}×)\npeak at {pk[k-1]:.1f}°")
ax=fig.add_subplot(gs[4,0:2],projection="polar"); f=sum(c[k-1]*np.cos(k*th) for k in range(1,M+1))
_bw=0.3/np.max(np.abs(f))*1e-5; ax.fill_between(Lr,1-_bw,1+_bw,color="green",alpha=0.35,lw=0); ax.plot(Lr,np.full_like(Lr,1.0),color="green",lw=0.8)
polar(ax,f,"k",f"Sum of m = 1–12 (radius = 1 + scaled residual)\nrms {np.sqrt(np.mean(R**2)):.3e} (= {np.sqrt(np.mean(R**2))/1e-5:.0f}×), max {np.abs(R).max():.2e}",k=0.3,lw=2)
ax=fig.add_subplot(gs[4,2:4],projection="polar")
for k in range(1,7): ax.plot(Lr,1+c[k-1]*np.cos(k*th)/A.max()*0.6,color=cols[k-1],lw=1.3,label=f"m={k}")
ax.plot([np.radians(phi)]*2,[0,1.7],color="k",lw=0.6,ls="--")
ax.set_ylim(0,1.7); ax.set_yticklabels([]); ax.set_thetagrids([0,90,180,270]); ax.tick_params(labelsize=6,pad=-2); ax.set_title("m = 1–6 overlaid, same scale",fontsize=7.3,pad=8)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.06),ncol=6,frameon=False,fontsize=6.5,handlelength=1.5,columnspacing=0.8)
save(fig,"fig3_modes_polar")
# ---------- Table (LaTeX booktabs)
def sci(x):
    e=int(np.floor(np.log10(abs(x)))); mnt=x/10**e; return f"${mnt:+.2f}\\times10^{{{e}}}$"
rows=[]
for k in range(1,M+1):
    org="fit matches (cancels)" if k<=2 else ("shell vs Doppler mismatch" if k<=6 else "Doppler tail only")
    rows.append(f"{k} & {360/k:.1f} & {sci(c[k-1])} & {A[k-1]/1e-5:.3g} & {pk[k-1]:.2f} & {100*vf[k-1]:.2f} & {sci(cG[k-1])} & {sci(cS[k-1])} & {org}\\\\")
open(os.path.join(FIG,"table_modes.tex"),"w").write("\n".join(rows)+"\n")
print("ok")
# ---------- Fig (waterfall): adapted from /workspace/cosmology-course/notes/cmb-waterfall.py (copy kept here).
# Steps recomputed from the model instead of hard-coded, and checked against the original's numbers.
G0=g["G_FULL"]; Lw=np.radians(np.arange(360.0)); thw=Lw-np.radians(phi)
def strip(Rv,ms): return Rv-sum(c[k-1]*np.cos(k*thw) for k in ms)
rr=[float(np.sqrt(np.mean((G0-1)**2))),float(np.sqrt(np.mean(R**2))),float(np.sqrt(np.mean(strip(R,range(1,6))**2))),
    float(np.sqrt(np.mean(strip(R,range(1,7))**2))),float(np.sqrt(np.mean(strip(R,range(1,13))**2)))]
orig=[5.751e-02,2.373e-03,6.113e-04,2.168e-04,1.833e-06]; print("waterfall recomputed",[f"{v:.4e}" for v in rr],"original",orig)
json.dump(dict(recomputed=rr,original=orig),open(os.path.join(HERE,"waterfall_check.json"),"w"),indent=1)
x=[v/1e-5 for v in rr]; lab=["shell fit\n(start)","+ SR\n+ mass","solve\nm = 1–5","solve\nm = 1–6","solve\nm = 1–12"]
colw=["k","C3","C1","C2","C0"]
from matplotlib.lines import Line2D
ST=json.load(open(os.path.join(HERE,"placement_steps.json")))   # from make_placement_steps.py (solution_space fits)
pick=[0,2,3,6]; ttl=["Step 1: no bodies yet","Step 2: 4 masses","Step 3: 4 masses + 4 voids","Step 4: 6 masses + 5 voids"]
fig=plt.figure(figsize=(7.1,5.7))
gs=fig.add_gridspec(3,4,height_ratios=[1.45,1.35,0.9],hspace=0.5,wspace=0.42,left=0.08,right=0.985,top=0.955,bottom=0.075)
a=fig.add_subplot(gs[0,0:2]); a.bar(range(5),x,color=colw,width=0.7); a.set_yscale("log"); a.set_ylim(0.05,8e4)
a.axhline(1,ls="--",color="purple",lw=1.2); a.text(4.45,1.3,"CMB level 1×",color="purple",fontsize=6.5,ha="right")
for i,v in enumerate(x): a.text(i,v*1.35,f"{v:,.0f}×" if v>=1 else f"{v:.2f}×",ha="center",fontsize=7)
for i in range(1,5): a.text(i-0.5,2.6e4,f"÷{x[i-1]/x[i]:.1f}",ha="center",fontsize=6.5,color="dimgray",fontweight="bold")
a.set_xticks(range(5)); a.set_xticklabels(lab,fontsize=6.3); a.set_ylabel("rms anisotropy\n(× $10^{-5}$, log)")
a.set_title("(a) anisotropy left after each step")
k_=fig.add_subplot(gs[0,2:4]); k_.axis("off")
k_.set_title("(b) placing masses and voids (side view, below)",loc="left")
hd=[Line2D([],[],marker="o",ls="",color="#c0392b",ms=6,label="mass (dot size ∝ mass,\n0.4–2×10$^{52}$ kg)"),
    Line2D([],[],marker="o",ls="",mfc="#aed6f1",mec="#2471a3",ms=9,label="void (drawn to size,\nradius 3.5–4 Gpc)"),
    Line2D([],[],marker="o",ls="",mfc="none",mec="k",mew=1.2,ms=10,label="body added in this step"),
    Line2D([],[],color="#2471a3",lw=2,label="CMB ring (edge of our sphere)"),
    plt.Rectangle((0,0),1,1,fc="#fadbd8",ec="none",label="CoE interior (no bodies)")]
k_.legend(handles=hd,loc="center left",bbox_to_anchor=(0.02,0.45),frameon=False,fontsize=6.8,labelspacing=1.1,handletextpad=0.6)
Ld=np.arange(360)
for j,(ix,tt) in enumerate(zip(pick,ttl)):
    s_=ST[ix]; ax=fig.add_subplot(gs[1,j]); W=36
    ax.add_patch(Circle(CG,RIN,fc="#fadbd8",ec="none",zorder=0))
    ax.add_patch(Circle((0,0),RV,fill=False,ec="#2471a3",lw=1.6,zorder=2)); ax.plot(0,0,"k+",ms=5,mew=1,zorder=3)
    ax.annotate("",xy=u*(W-3),xytext=u*(RV+2),arrowprops=dict(arrowstyle="->",lw=0.7,color="#c0392b"))
    if j==0:
        ax.text(*(u*24+np.array([0,-6])),"to CoE",fontsize=6,color="#c0392b",ha="center")
        ax.text(0,-RV-6.5,"our sphere\n14.3 Gpc",fontsize=6,color="#2471a3",ha="center",va="top")
    prev=ST[pick[j-1]]["kinds"] if j>0 else ""
    new=range(len(prev),len(s_["kinds"]))
    for n,(kk,X,M,Rv) in enumerate(zip(s_["kinds"],s_["X"],s_["M"],s_["Rv"])):
        if kk=="C": ax.plot(X[0],X[1],"o",color="#c0392b",ms=3+5*M/2e52,zorder=4)
        else: ax.add_patch(Circle(X[:2],Rv,fc="#aed6f1",ec="#2471a3",lw=0.8,zorder=3))
        if n in new: ax.add_patch(Circle(X[:2],(Rv if kk!="C" else 2.4)+1.3,fill=False,ec="k",lw=0.9,zorder=5))
    ax.set_xlim(-W,W); ax.set_ylim(-W,W); ax.set_aspect("equal"); ax.tick_params(labelsize=6)
    ax.set_xticks([-30,0,30]); ax.set_yticks([-30,0,30])
    if j==0: ax.set_ylabel("y (Gpc)",fontsize=7)
    ax.set_xlabel("x (Gpc)",fontsize=6.5,labelpad=1); ax.set_title(tt,fontsize=7.2)
    r_=fig.add_subplot(gs[2,j]); res=np.array(s_["resid_x"])
    r_.axhline(0,color="0.6",lw=0.5,zorder=1); r_.axhspan(-1,1,fc=(0.0,0.5,0.0,0.25),ec=(0.0,0.5,0.0,0.9),lw=0.6,zorder=1.5)
    r_.plot(Ld,res,color="k",lw=1.0); lim=max(1.6,1.15*np.abs(res).max()); r_.set_ylim(-lim,lim)
    r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=6)
    v=s_["left"]; r_.set_title(f"leftover {v:,.0f}× the CMB level" if v>=10 else f"leftover {v:.2g}× the CMB level",fontsize=6.8,color="#c0392b" if v>1 else "green")
    r_.set_xlabel("longitude ℓ (°) around the ring",fontsize=6.3,labelpad=1)
    if j==0: r_.set_ylabel("leftover\n(× $10^{-5}$)",fontsize=6.5)
fig.text(0.985,0.005,"Bottom row: leftover around the ring; each panel has its own vertical scale. Green band: ±1× (CMB level).",ha="right",fontsize=6.2,color="0.3")
save(fig,"fig_waterfall")
print("placement steps", [(ST[i]["kinds"],round(ST[i]["left"],3)) for i in pick])
