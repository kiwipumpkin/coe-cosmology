"""v3 one-panel-per-page structure images (Chris, 11 Oct 23:06). Physics = resonant_struct_pair.py run with GROWTH_P=3, BETA=1
(late cubic onset, growth amplitude proportional to rho/rho_unmod, bounded symmetric contrast |delta| <= 0.999, exact per-shell mean);
the computation is executed from that file unchanged (its page drawing is not run). One page per (mode, mass case): 8 images.
Under the wedge: a 0.15 in neutral zone strip (compressed = dark grey diagonal hatch where psi < 0, i.e. shells narrower than the mean;
expanded = light grey where psi > 0), then three equal-height graphs: standing-wave width, density / unmodulated, clumpiness.
No coloured shading behind any graph; thin grey dashed lines mark the zone boundaries (psi = 0) through the graphs.
Nothing in main.tex or any PDF is touched."""
import os
os.environ["GROWTH_P"]="3"; os.environ["BETA"]="1"
src=open("resonant_struct_pair.py").read().split("# ---------------- geometry (inches) ----------------")[0]
exec(src)
from matplotlib.patches import Rectangle
# ---------------- geometry ----------------
W=7.4; LM=0.95; AW=6.0; HW=AW*np.sin(np.radians(SPAN))*1.01; HZ=0.15; GZ=0.06; HG=0.70; GP=0.08
TOP=0.66; XL=0.36; LEG=0.55; CBH=0.82
H=TOP+HW+HZ+GZ+3*HG+2*GP+XL+LEG+CBH
TH=np.linspace(0,np.radians(SPAN),400)
def zones(psi):
    xs=np.linspace(0,1,20001)[1:-1]; sg=np.sign(psi(xs)); sg[sg==0]=1   # interior points only (psi = 0 at both ends)
    cut=[0.0]+[float(xs[j]) for j in range(1,len(xs)) if sg[j]!=sg[j-1]]+[1.0]
    out=[]
    for a,b in zip(cut[:-1],cut[1:]):
        out.append((R0/R+a*(R-R0)/R,R0/R+b*(R-R0)/R,psi(0.5*(a+b))<0))
    return out   # (x0, x1 in r/R, compressed?)
def page(key,title,lab,psi,tag):
    d=RES[key]; e=d["e"]; xe=e/R; S=d[tag]["S"]; rows=S["rows"]; rho=d[tag]["rho"]; w=d["w"]; cw=d["c"]
    fig=plt.figure(figsize=(W,H)); Ax=lambda y,h: fig.add_axes([LM/W,y/H,AW/W,h/H])
    mname="equal mass per shell" if tag=="eq" else "mass rising outward"
    fig.text(0.5,1-0.22/H,f"Late structure growth — mode {title.split(':')[0]}, {mname}",ha="center",va="center",fontsize=12,fontweight="bold")
    fig.text(0.5,1-0.46/H,f"{title.split(': ')[1]};  {lab},  A = {A:g};  rms = 0.8 (i/23)³ ρ/ρ$_{{unmod}}$,  |δ| ≤ 1"+("" if tag=="eq" else f";  mass ×{QM:.4f}/shell"),ha="center",va="center",fontsize=8.4,color="0.25")
    yw=H-TOP-HW; aw=Ax(yw,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
    img=CMAP(NORM(np.clip(np.nan_to_num(S["D"],nan=1.0),VMIN,VMAX))); img[np.isnan(S["D"])]=(1,1,1,0); img[S["core"]]=(0.86,0.86,0.86,1)
    aw.imshow(img,origin="lower",extent=[0,1,0,S["asp"]],interpolation="antialiased",zorder=1)
    for b in xe: aw.plot(b*np.cos(TH),b*np.sin(TH),color="white",lw=0.45,zorder=3)
    aw.text(0.0,HW/AW*0.99,f"{N} shells, inner radius {R0:g}, R = {R:g}\nfirst shells flat; structure only in outer shells\nouter-shell clumpiness {rows[-1]['sigma_measured']:.2f}"+(" (at the |δ| ≤ 1 limit)" if rows[-1].get("capped") else ""),fontsize=9.3,va="top",ha="left",zorder=4)
    # zone strip
    az=Ax(yw-HZ-GZ,HZ); az.set_xlim(0,1); az.set_ylim(0,1); az.set_xticks([]); az.set_yticks([])
    for s_ in az.spines.values(): s_.set_visible(False)
    az.add_patch(Rectangle((xe[0],0),1-xe[0],1,facecolor="none",edgecolor="0.5",lw=0.6,zorder=5))
    Z=zones(psi)
    for x0,x1,comp in Z:
        az.add_patch(Rectangle((x0,0),x1-x0,1,facecolor="white" if comp else "0.88",edgecolor="0.25" if comp else "none",hatch="////" if comp else None,lw=0))
    az.set_ylabel("zone",rotation=0,fontsize=8,labelpad=14,va="center")
    plt.rcParams["hatch.linewidth"]=0.6
    axs=[]; y=yw-HZ-GZ-GP
    for _ in range(3): y-=HG; axs.append(Ax(y,HG)); y-=GP
    aW,aD,aC=axs
    for k in range(N): aW.fill_between([xe[k],xe[k+1]],0.4,w[k]/cw,color=CMAP(NORM(rho[k])),alpha=0.55,lw=0)
    for k in range(N): aW.hlines(w[k]/cw,xe[k],xe[k+1],color="0.25",lw=0.9)
    rr_=np.linspace(R0,R,600); aW.plot(rr_/R,1+A*psi(X(rr_)),color="#6c3483",lw=1.3); aW.axhline(1,color="0.6",lw=0.6)
    aW.set_ylim(0.4,1.6); aW.set_yticks([0.5,1,1.5]); aW.set_ylabel("shell width\n1 + Aψ")
    dc=S["dcon"]
    for k in range(N): aD.hlines(dc[k],xe[k],xe[k+1],color=CMAP(NORM(rho[k])),lw=2.8,zorder=3)
    aD.axhline(1,color="0.5",lw=0.7); aD.set_yscale("log"); aD.set_ylim(*QYL); aD.set_yticks([0.5,1,2]); aD.set_yticklabels(["0.5","1","2"]); aD.yaxis.set_minor_formatter(mt.NullFormatter())
    aD.set_ylabel("density ÷\nunmodulated")
    sm=[r["sigma_measured"] for r in rows]; st=[r["sigma_target"] for r in rows]
    for k in range(N): aC.hlines(sm[k],xe[k],xe[k+1],color="#7b241c",lw=2.6,zorder=3)
    st=np.minimum(st,1.06); aC.step(np.r_[xe[:-1],xe[-1]],np.r_[st,st[-1]],where="post",color="0.35",lw=0.8,ls="--",zorder=2)
    aC.axhline(1,color="k",lw=0.6,ls=(0,(1,1))); aC.set_ylim(0,1.08); aC.set_yticks([0,0.5,1]); aC.set_ylabel("clumpiness\nrms(δ)")
    for a in axs:
        for b in xe: a.axvline(b,color="0.85",lw=0.5,zorder=0)
        a.set_xlim(0,1); a.grid(axis="y",alpha=0.25)
    for a in axs+[az]:
        for x0,x1,comp in Z[1:]: a.axvline(x0,color="0.45",lw=0.8,ls=(0,(4,3)),zorder=1)
    for a in axs[:-1]: plt.setp(a.get_xticklabels(),visible=False)
    aC.set_xlabel("normalised radius  r / R   (light grey lines = shell boundaries)",labelpad=2)
    hh=[Patch(facecolor="white",edgecolor="0.25",hatch="////",lw=0.5),Patch(facecolor="0.88",edgecolor="0.6",lw=0.5),Line2D([],[],color="0.45",lw=0.8,ls=(0,(4,3))),
        Line2D([],[],color="#6c3483",lw=1.3),Line2D([],[],color="0.3",lw=2.8),Line2D([],[],color="#7b241c",lw=2.6),Line2D([],[],color="0.35",lw=0.8,ls="--"),Line2D([],[],color="k",lw=0.6,ls=(0,(1,1)))]
    ll=["compressed (narrower shells)","expanded (wider shells)","zone boundary","standing wave (bars: shell widths)","shell density ÷ unmodulated (colour = density)",
        "measured clumpiness","growth rule target","|δ| = 1 limit"]
    fig.legend(hh,ll,loc="upper center",bbox_to_anchor=(0.5,(y+GP-XL-0.02)/H),ncol=3,frameon=False,fontsize=8.6,borderaxespad=0,columnspacing=1.3,handlelength=2.0,labelspacing=0.35)
    cax=fig.add_axes([LM/W,0.45/H,AW/W,0.15/H])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation="horizontal",extend="min")
    cb.set_label("density ÷ innermost shell of the unmodulated equal-mass set  (one log scale for all eight v3 pages)",fontsize=8.8,labelpad=2); cb.ax.tick_params(labelsize=8)
    out=f"resonant_sine_struct_v3_{key}_{tag}.png"; fig.savefig(out,dpi=250); plt.close(fig); return out
for key,title,lab,psi in SC:
    for tag in ("eq","mo"): print(page(key,title,lab,psi,tag))
print("page",W,round(H,2),"norm",VMIN,VMAX)
