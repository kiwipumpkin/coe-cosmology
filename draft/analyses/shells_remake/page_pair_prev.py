"""Two-panels-per-page versions (PAIR=A: (a) constant + (b) widening; PAIR=B: (c) narrowing + (d) comparison), 45 deg wedges, all graphs
under a wedge the same height. Derived from page_2x2.py. Original description of the 2x2 version follows.
One-page 2x2 versions of the four shell figure sets (Chris, 11 Oct): 22 deg wedges instead of 45 deg.
Physics, fits, colour scales and structure rules are IMPORTED unchanged from shells_remake.py (equal-mass fit A x^-n + C)
and shells_massout.py (mass-rising fit A x^-n + C x^p; structure builder sbuild). Only change: the wedge opening angle SPAN = 22 deg.
Structure fields are rebuilt on the 22 deg wedge (same seed 11, same nx = 1100) so each shell's clumpy field averages EXACTLY to its
shell density over the area actually drawn; the noise realisation therefore differs from the 45 deg images and the per-shell statistics
shift slightly (saved to page_2x2_22deg.json).
Layout: page 7.1 x 9.8 in (one portrait page of the two-column brief at full text width), quadrants a-c = constant / widening / narrowing,
d = the comparison plot. Each wedge axes has the same x-range (r/R = 0..1) and width as the graphs under it, equal aspect,
point at r = 0, flat edge on the top of the first graph, so each shell edge drops straight onto its grey boundary line.
Legends are below the graphs. Nothing in main.tex or any PDF is touched."""
import json,numpy as np,matplotlib,matplotlib.ticker as mt
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Patch
from matplotlib.lines import Line2D
import shells_massout as SM
import os
SPAN=float(os.environ.get("SPAN","45")); SM.SPAN=SPAN; TAG=f"{SPAN:g}deg"
plt.rcParams.update({"font.size":9.5,"axes.labelsize":9.5,"xtick.labelsize":8.5,"ytick.labelsize":8.5,"axes.titlesize":10.5})
PAIR=os.environ.get("PAIR","A"); SLOT=dict(constant=0,widening=1) if PAIR=="A" else dict(narrowing=0,d=1)
N=SM.N; CMAP=SM.CMAP; NORM_W=SM.NORM_W; NORM_S=SM.NORM_S; CASES=SM.CASES
COLS={"constant":"#555555","widening":"#1747a6","narrowing":"#b03050"}
LET=dict(constant="a",widening="b",narrowing="c")
def fit_eq(xe,rho):   # equal-mass fit of shells_remake.py: A x^-n + C
    from scipy.optimize import curve_fit
    xm=0.5*(xe[1:]+xe[:-1]); f=lambda x,lA,n,C: np.log(np.exp(lA)*x**(-n)+C)
    p,_=curve_fit(f,xm,np.log(rho),p0=[np.log(rho[0]*xm[0]**2),2.0,1e-3],bounds=([-30,0,0],[30,10,10]),maxfev=20000)
    A,n,C=np.exp(p[0]),p[1],p[2]
    return (lambda x: A*x**(-n)+C),dict(A=float(A),n=float(n),C=float(C))
def lab_eq(pf): return rf"smooth fit  $\rho={pf['A']:.3g}\,x^{{-{pf['n']:.2f}}}$"+(rf"$\,+\,{pf['C']:.3g}$" if pf['C']>1e-6 else "")
def lab_mo(pf): return rf"smooth fit  $\rho={pf['A']:.3g}\,x^{{-{pf['n']:.2f}}}+{pf['C']:.3g}\,x^{{{pf['p']:.2f}}}$"
# ---------- page geometry (inches) ----------
W=7.4; LMc=0.82; AW=2.55; GAPC=0.42; X0=[LMc,LMc+AW+GAPC+LMc]
HW=AW*np.sin(np.radians(SPAN))*1.01; GP=0.08; HG=1.1   # every graph under a wedge = wave-graph height (Chris, 11 Oct)
def page(title,sub,hs,legend_h):
    # row block: panel title 0.22, wedge HW, graphs hs (gaps GP), x label 0.30, per-panel legend legend_h, row gap 0.12
    blk=0.30+HW+sum(hs)+GP*(len(hs)-1)+0.40+legend_h+0.12
    top=0.45+(0.22 if sub else 0); bot=1.35   # shared legend + colour bar strip
    H=top+blk+bot; fig=plt.figure(figsize=(W,H))
    fig.text(0.5,1-0.2/H,title,ha="center",va="center",fontsize=12,fontweight="bold")
    if sub: fig.text(0.5,1-0.47/H,sub,ha="center",va="center",fontsize=9.5,color="0.25")
    A=lambda x,y,w,h: fig.add_axes([x/W,y/H,w/W,h/H])
    slots=[]
    for r in range(1):
        ytop=H-top-r*blk
        for c in range(2):
            y=ytop-0.30-HW; aw=A(X0[c],y,AW,HW); axs=[]
            for h in hs: y-=h; axs.append(A(X0[c],y,AW,h)); y-=GP
            y+=GP; slots.append(dict(aw=aw,axs=axs,ytop=ytop,ybot=y,leg_y=y-0.42,x0=X0[c]))
    return fig,H,slots,A
def wedge_axes(aw):
    aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
TH=np.linspace(0,np.radians(SPAN),300)
def band(r0,r1): return np.r_[np.c_[r1*np.cos(TH),r1*np.sin(TH)],np.c_[r0*np.cos(TH[::-1]),r0*np.sin(TH[::-1])]]
def common(axs,xe):
    for a in axs:
        for b in xe: a.axvline(b,color="0.8",lw=0.7,zorder=0)
        a.set_xlim(0,1); a.grid(axis="y",alpha=0.25)
    for a in axs[:-1]: plt.setp(a.get_xticklabels(),visible=False)
    axs[-1].set_xlabel("normalised radius  r / R",labelpad=2)
def info_txt(Rout,q,extra=""):
    return f"R = {Rout:.2f}"+("" if q==1 else f"\neach shell\n{q:.4f}× as wide")+extra
def panel_title(s,key,title): s["aw"].set_title(f"({LET[key]}) {title}",loc="left",pad=4,fontweight="bold")
def below_legend(fig,s,H,handles,labels,ncol=1,fs=8.8):
    fig.legend(handles,labels,loc="upper left",bbox_to_anchor=((s["x0"]-0.62)/W,s["leg_y"]/H),frameon=False,fontsize=fs,ncol=ncol,handlelength=2.0,borderaxespad=0,columnspacing=1.0,labelspacing=0.3)
def cbar(fig,H,norm,label):
    cax=fig.add_axes([X0[0]/W,0.50/H,(X0[1]+AW-X0[0])/W,0.15/H])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=CMAP),cax=cax,orientation="horizontal"); cb.set_label(label,fontsize=9.5,labelpad=2); cb.ax.tick_params(labelsize=8.5)
# ======================= waves pages =======================
def wave_page(massout,out):
    masses=SM.Mi if massout else np.ones(N)
    T="Shells with mass rising outward (outermost 3× innermost, same total)" if massout else "Equal-mass concentric shells"
    fig,H,S,A=page(T,None,[HG,HG,HG],0.22)
    res={}
    for key,title,Rout in CASES:
        q=SM.q_for(Rout); e=SM.edges(q); xe=e/e[-1]; rho=SM.dens(e,masses)
        g,pf=(SM.fitrho(xe,rho) if massout else fit_eq(xe,rho)); res[key]=dict(q=q,xe=xe,rho=rho,rhoeq=SM.dens(e,np.ones(N)),fit=pf)
        if key not in SLOT: continue
        s=S[SLOT[key]]
        aw=s["aw"]; wedge_axes(aw); panel_title(s,key,title.replace(" (geometric growth)","").replace(" (geometric shrink)",""))
        aw.add_patch(Polygon(np.r_[[[0,0]],np.c_[xe[0]*np.cos(TH),xe[0]*np.sin(TH)]],fc="#dcdcdc",ec="none"))
        for k in range(N): aw.add_patch(Polygon(band(xe[k],xe[k+1]),fc=CMAP(NORM_W(rho[k])),ec="white",lw=0.7))
        aw.text(0.0,HW/AW*0.98,info_txt(Rout,q),fontsize=9,va="top",ha="left")
        a1,a2,a3=s["axs"]; x=np.linspace(0,1,3000); wv=SM.waves(e,x); rf=np.where(x>=xe[0],g(np.maximum(x,1e-9)),np.nan); env=rf/g(xe[0])
        if massout: a1.step(np.r_[xe[:-1],xe[-1]],np.r_[res[key]["rhoeq"],res[key]["rhoeq"][-1]],where="post",color="0.55",lw=0.9,ls=":")
        for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM_W(rho[k])),lw=3)
        a1.plot(x,rf,"--",color="0.15",lw=1.0)
        a1.set_yscale("log"); a1.set_ylim(10**-2.3,1.8); a1.set_yticks([1e-2,1e-1,1]); a1.set_ylabel("density\n(÷ inner)" if not massout else "density\n(÷ eq. inner)",fontsize=9)
        a2.plot(x,wv,color="#2c3e50",lw=1.0); a2.set_ylabel("wave",fontsize=9); a2.set_ylim(-1.3,1.3); a2.set_yticks([-1,0,1])
        a3.plot(x,env,color="#b03050",lw=0.9,ls="--"); a3.plot(x,-env,color="#b03050",lw=0.9,ls="--"); a3.plot(x,env*wv,color="#2c3e50",lw=1.0)
        a3.set_ylabel("density fit\n× wave",fontsize=9); a3.set_ylim(-1.3,1.3); a3.set_yticks([-1,0,1])
        for a in (a2,a3): a.axhline(0,color="0.6",lw=0.6)
        common([a1,a2,a3],xe)
        below_legend(fig,s,H,[Line2D([],[],ls="--",color="0.15",lw=1.0)],[lab_mo(pf) if massout else lab_eq(pf)])
    # (d) comparison of the three density profiles, spanning the same block height
    if "d" in SLOT:
        s=S[SLOT["d"]]; [a.remove() for a in s["axs"]]; s["aw"].remove()
        ax=A(s["x0"],s["ybot"]+0.62,AW,s["ytop"]-0.30-s["ybot"]-0.30-0.62)
        ax.set_title("(d) The three density profiles",loc="left",pad=4,fontweight="bold")
        hs,ls=[],[]
        for key,title,Rout in CASES:
            r=res[key]; xe=r["xe"]
            h,=ax.step(np.r_[xe[:-1],xe[-1]],np.r_[r["rho"],r["rho"][-1]],where="post",color=COLS[key],lw=1.8); hs.append(h)
            ls.append(f"{key}, R = {Rout:.2f}"+(f", outer/inner {r['rho'][-1]/r['rho'][0]:.3f}" if not massout else ", mass rising"))
            if massout:
                h2,=ax.step(np.r_[xe[:-1],xe[-1]],np.r_[r["rhoeq"],r["rhoeq"][-1]],where="post",color=COLS[key],lw=0.9,ls="--",alpha=0.55,zorder=1)
        if massout: hs.append(Line2D([],[],color="0.4",lw=0.9,ls="--",alpha=0.55)); ls.append("faint dashed: equal-mass profiles")
        ax.set_yscale("log"); ax.set_xlim(0,1); ax.set_ylim(10**-2.3,1.8); ax.set_xlabel("normalised radius  r / R",labelpad=2)
        ax.set_ylabel("shell density ÷ innermost" + (" equal-mass shell" if massout else " shell")); ax.grid(alpha=0.3,which="both")
        fig.legend(hs,ls,loc="upper left",bbox_to_anchor=((s["x0"]-0.62)/W,(s["ybot"]+0.62-0.48)/H),frameon=False,fontsize=9,ncol=1,borderaxespad=0,labelspacing=0.3)
    # shared legend for the panels + colour bar
    hh=[Line2D([],[],color="#2c3e50",lw=1.0),Line2D([],[],color="#b03050",lw=0.9,ls="--")]; ll=["wave (one cycle per shell);  fit × wave","± density fit"]
    hh.append(Line2D([],[],color="0.8",lw=0.7)); ll.append("shell boundaries")
    if massout: hh.append(Line2D([],[],color="0.55",lw=0.9,ls=":")); ll.append("dotted: equal mass")
    fig.legend(hh,ll,loc="center",bbox_to_anchor=(0.5,1.08/H),ncol=2,frameon=False,fontsize=9,borderaxespad=0,columnspacing=1.5)
    cbar(fig,H,NORM_W,"shell density ÷ innermost"+(" equal-mass" if massout else "")+" shell  (one logarithmic colour scale for every shell figure)")
    fig.savefig(out,dpi=250); plt.close(fig); return H,{k:dict(fit=v["fit"]) for k,v in res.items()}
# ======================= structure pages =======================
def struct_page(massout,out):
    masses=SM.Mi if massout else np.ones(N)
    T="Structure in shells with mass rising outward" if massout else "Structure in equal-mass shells"
    sub="innermost shell smooth, clumpier outward; each shell's clumpy density averages exactly to its shell density"
    fig,H,S,A=page(T,sub,[HG,HG,HG],0.22)
    res={}
    for key,title,Rout in CASES:
        q=SM.q_for(Rout); e=SM.edges(q); xe=e/e[-1]; rho=SM.dens(e,masses)
        g,pf=(SM.fitrho(xe,rho) if massout else fit_eq(xe,rho))
        B=SM.sbuild(e,rho); rows=B["rows"]; Dv=B["D"][~np.isnan(B["D"])]
        res[key]=dict(q=q,xe=xe,fit=pf,rows=rows,frac_above_scale=float(np.mean(Dv>NORM_S.vmax)),frac_below_scale=float(np.mean(Dv<NORM_S.vmin)))
        if key not in SLOT: continue
        s=S[SLOT[key]]
        aw=s["aw"]; wedge_axes(aw); panel_title(s,key,title.replace(" (geometric growth)","").replace(" (geometric shrink)",""))
        img=CMAP(NORM_S(np.clip(np.nan_to_num(B["D"],nan=1.0),NORM_S.vmin,NORM_S.vmax))); img[np.isnan(B["D"])]=(1,1,1,0); img[B["core"]]=(0.86,0.86,0.86,1)
        aw.imshow(img,origin="lower",extent=[0,1,0,B["asp"]],interpolation="antialiased",zorder=1)
        for b in xe: aw.plot(b*np.cos(TH),b*np.sin(TH),color="white",lw=0.6,zorder=3)
        aw.text(0.0,HW/AW*0.98,info_txt(Rout,q),fontsize=9,va="top",ha="left",zorder=4)
        a1,a2,a3=s["axs"]
        for k in range(N):
            p=rows[k]["pct"]; c=CMAP(NORM_S(rho[k]))
            a1.fill_between([xe[k],xe[k+1]],p["p2.5"],p["p97.5"],color=c,alpha=0.18,lw=0,zorder=1)
            a1.fill_between([xe[k],xe[k+1]],p["p16"],p["p84"],color=c,alpha=0.45,lw=0,zorder=2)
            if k>0: a1.errorbar([0.5*(xe[k]+xe[k+1])],[rho[k]],yerr=[[rho[k]-p["min"]],[p["max"]-rho[k]]],fmt="none",ecolor="0.25",elinewidth=0.7,capsize=2.0,capthick=0.7,zorder=3)
        for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM_S(rho[k])),lw=3,zorder=4)
        xx=np.linspace(xe[0],1,600); a1.plot(xx,g(xx),"--",color="0.15",lw=0.9,zorder=5)
        xm=0.5*(xe[1:]+xe[:-1]); a1.plot(xm,[r["mean_over_rho"]*r["rho"] for r in rows],"o",ms=3,mfc="white",mec="k",mew=0.7,zorder=6)
        a1.set_yscale("log"); a1.set_ylim(10**-3.3,30); a1.set_yticks([1e-3,1e-1,10]); a1.set_ylabel("density\n(÷ inner)" if not massout else "density\n(÷ eq. inner)",fontsize=9)
        a1.yaxis.set_minor_formatter(mt.NullFormatter())
        sm=[r["sigma_measured"] for r in rows]; sp=[r["spacing_measured"] for r in rows]
        for k in range(N): a2.hlines(sm[k],xe[k],xe[k+1],color="#7b241c",lw=2.5)
        a2.set_ylim(0,1.6); a2.set_yticks([0,0.5,1,1.5]); a2.set_ylabel("rms\nclumpiness",fontsize=9)
        for k in range(1,N): a3.hlines(sp[k],xe[k],xe[k+1],color="#1f3a93",lw=2.5)
        a3.text(0.5*(xe[0]+xe[1]),0.3,"smooth",ha="center",va="center",fontsize=8.5,color="0.4")
        a3.set_yscale("log"); a3.set_ylim(0.2,3); a3.set_yticks([0.3,1,2]); a3.set_yticklabels(["0.3","1","2"]); a3.yaxis.set_minor_formatter(mt.NullFormatter())
        a3.set_ylabel("clump\nspacing",fontsize=9)
        common([a1,a2,a3],xe)
        below_legend(fig,s,H,[Line2D([],[],ls="--",color="0.15",lw=0.9)],[lab_mo(pf) if massout else lab_eq(pf)])
    if "d" in SLOT:
        s=S[SLOT["d"]]; [a.remove() for a in s["axs"]]; s["aw"].remove()
        yb=s["ybot"]+0.62; tot=s["ytop"]-0.30-yb-0.30; h=(tot-0.14)/2
        b0=A(s["x0"],yb+h+0.14,AW,h); b1=A(s["x0"],yb,AW,h)
        b0.set_title("(d) Clumpiness and spacing",loc="left",pad=4,fontweight="bold")
        hs,ls=[],[]
        for key,title,Rout in CASES:
            r=res[key]; xe=r["xe"]; sm=[q_["sigma_measured"] for q_ in r["rows"]]; sp=[q_["spacing_measured"] or np.nan for q_ in r["rows"]]
            hh,=b0.step(np.r_[xe[:-1],xe[-1]],np.r_[sm,sm[-1]],where="post",color=COLS[key],lw=1.8); hs.append(hh); ls.append(f"{key} (R = {Rout:.2f})")
            b1.step(np.r_[xe[:-1],xe[-1]],np.r_[sp,sp[-1]],where="post",color=COLS[key],lw=1.8)
        b0.set_ylabel("clumpiness\nrms(ρ/ρ_shell − 1)"); b0.set_ylim(0,1.6); b0.grid(alpha=0.3); b0.set_xlim(0,1); plt.setp(b0.get_xticklabels(),visible=False)
        b1.set_ylabel("typical clump spacing\n(radius units)"); b1.set_yscale("log"); b1.set_yticks([0.3,0.5,1,2]); b1.set_yticklabels(["0.3","0.5","1","2"]); b1.yaxis.set_minor_formatter(mt.NullFormatter())
        b1.set_ylim(0.3,2.3); b1.grid(alpha=0.3,which="both"); b1.set_xlim(0,1); b1.set_xlabel("normalised radius  r / R",labelpad=2)
        fig.legend(hs,ls,loc="upper left",bbox_to_anchor=((s["x0"]-0.62)/W,(yb-0.48)/H),frameon=False,fontsize=9,ncol=1,borderaxespad=0,labelspacing=0.3)
    hh=[Line2D([],[],color="0.3",lw=3),Line2D([],[],ls="none",marker="o",ms=3.5,mfc="white",mec="k",mew=0.7),Patch(color="0.45",alpha=0.45),Patch(color="0.45",alpha=0.18),Line2D([],[],color="0.25",lw=0.7,marker="_",ms=6)]
    ll=["shell density","mean of clumpy field","16–84%","2.5–97.5%","min–max","shell boundaries"]; hh.append(Line2D([],[],color="0.8",lw=0.7))
    fig.legend(hh,ll,loc="center",bbox_to_anchor=(0.5,1.08/H),ncol=3,frameon=False,fontsize=9,borderaxespad=0,columnspacing=1.2,handlelength=1.6)
    cbar(fig,H,NORM_S,"density ÷ innermost"+(" equal-mass" if massout else "")+" shell  (one logarithmic colour scale for every structure figure)")
    fig.savefig(out,dpi=250); plt.close(fig)
    return H,{k:dict(fit=v["fit"],frac_above_scale=v["frac_above_scale"],frac_below_scale=v["frac_below_scale"],rows=v["rows"]) for k,v in res.items()}
if __name__=="__main__":
    import sys; which=sys.argv[1:] or ["equal_waves","equal_structure","massout_waves","massout_structure"]
    info={}
    for w in which:
        mo=w.startswith("massout"); f=wave_page if w.endswith("waves") else struct_page
        H,inf=f(mo,f"page_pair_{w}_{PAIR}.png"); info[w]=dict(fig_in=[W,H],**inf); print(w,"size in",W,round(H,2))
    try: old=json.load(open(f"page_pair_{PAIR}.json"))
    except Exception: old={}
    old.update(info); json.dump(old,open(f"page_pair_{PAIR}.json","w"),indent=1)
