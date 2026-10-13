"""Structure-formation versions of the sine resonant pair pages (Chris, 11 Oct 22:25).
Geometry/mass exactly as resonant_pair.py MODE=sine V2=1 (24 shells, r_in = 3, R = 15, w_i = c(1 + A psi(x_mid)), A = 0.5,
psi = +/- sin(pi x), +/- sin(2 pi x); left equal mass, right mass x 3^(1/23) per shell, same total).
Structure rule = structure_growth.py / shells_massout.sbuild, adapted to 24 shells:
    sigma_i = rms(rho/rho_i - 1) = 0.8 * i/(N-1) * (rho_i / rho_ref,i)^(1/3),  i = 0..23  (shell 1 smooth)
    rho_ref,i = density of shell i in the UNMODULATED equal-mass 24-shell set (equal widths 0.5) -- the analogue of the old
    'same shell in the constant-spacing equal-mass case'.
    clump spacing lambda_i = 0.6 w_i; Gaussian random field with power in a log-normal band (width 0.35 in ln k) around 2 pi/lambda_i,
    same white-noise seed (11) and grid (nx = 1100 over R = 15) as before; lognormal transform, rescaled so each shell's mean over its
    drawn area is EXACTLY rho_i. Measured clumpiness = rms over the shell's pixels; spacing = 2 pi / k at the power-spectrum peak.
    The old code asserted a strictly monotonic sigma; here the formula is applied unchanged and NOT forced monotonic: in two equal-mass
    cases (n=2 +, n=2 -) the target dips by <= 0.0032 over 1-2 shells where the resonance makes shells much less dense than the reference.
v2 (Chris, 11 Oct 22:50: 'pattern formation starts too soon; flat after recombination, takes time to form'), run with GROWTH_P=3:
    sigma_i = 0.8 * (i/(N-1))^P * (rho_i / rho_unmod(r_mid,i))^BETA,  P = 3, BETA = 1   (P = 1 reproduces v1 / the brief's rule)
    rho_unmod(r) = continuum density of the UNMODULATED equal-mass 24-shell set at the same radius (= the 'density / unmodulated'
    graph of the v2 sine pages), so the causal chain is explicit: wave compresses/expands a shell -> its density relative to the
    unmodulated set rises/falls -> its structure grows faster/slower (growth amplitude proportional to the density contrast, BETA = 1;
    v1 used a cube root).
    FIELD (v2, Chris 22:56): bounded symmetric contrast instead of lognormal: rho = rho_shell (1 + delta), delta = tanh(a g) minus its shell
    mean, rescaled if needed so |delta| <= 0.999. Mean of delta over the shell = 0 exactly (excess = deficit), voids never below zero,
    peaks never above 2x the shell mean; a solved per shell so rms(delta) = target (target capped at 95% of the bounded maximum,
    flagged 'capped'). Denser shells therefore reach larger ABSOLUTE amplitude (bound = rho_shell).
    Shell index i = time since recombination (innermost = earliest, exactly flat). Everything else unchanged; per-shell mean of the
    final drawn field is checked against the shell density. Output names get the suffix _v2.
Nothing in main.tex or any PDF is touched."""
import os,json,numpy as np,matplotlib,matplotlib.ticker as mt
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm,LinearSegmentedColormap
from scipy.optimize import brentq
plt.rcParams.update({"font.size":9,"axes.labelsize":8.8,"xtick.labelsize":8,"ytick.labelsize":8})
N=24; R0=3.0; R=15.0; A=0.5; SPAN=45.0; M=1.0; QM=3.0**(1/(N-1)); S_OUT=0.8; LFAC=0.6; BW=0.35
P=float(os.environ.get("GROWTH_P","1")); BETA=float(os.environ.get("BETA","1")); SUF="" if P==1 else "_v2"
MI_MO=QM**np.arange(N); MI_MO=N*MI_MO/MI_MO.sum()
CMAP=LinearSegmentedColormap.from_list("shells_v2",["#06103d","#123a9a","#2c78cf","#6ea6df","#b9a3d6","#ee8f95","#e0584a","#a8141c","#5c0812"])  # = v2 sine pair map
SC=[("half_pos","n = 1 (+): half wavelength, widest in the middle",r"$\psi=\sin(\pi x)$",lambda x:np.sin(np.pi*x)),
    ("half_neg","n = 1 (−): half wavelength, narrowest in the middle",r"$\psi=-\sin(\pi x)$",lambda x:-np.sin(np.pi*x)),
    ("full_pos","n = 2 (+): full wavelength, widens first",r"$\psi=\sin(2\pi x)$",lambda x:np.sin(2*np.pi*x)),
    ("full_neg","n = 2 (−): full wavelength, narrows first",r"$\psi=-\sin(2\pi x)$",lambda x:-np.sin(2*np.pi*x))]
X=lambda r:(r-R0)/(R-R0)
def march(c,psi):
    e=[R0]
    for i in range(N):
        w=c
        for _ in range(200):
            wn=c*(1+A*psi(X(e[-1]+w/2)))
            if abs(wn-w)<1e-14: break
            w=wn
        e.append(e[-1]+w)
    return np.array(e)
def vol(e): return 4/3*np.pi*(e[1:]**3-e[:-1]**3)
EREF=np.linspace(R0,R,N+1); RHO_REF0=M/vol(EREF)[0]; RHOREF=M/vol(EREF)/RHO_REF0
DB=0.999   # |delta| bound (1 = void reaches zero density; positive side bounded equally)
def bounded_delta(gm,target):
    def make(a):
        u=np.tanh(a*gm); u=u-u.mean(); mx=np.abs(u).max()
        return u*(DB/mx) if mx>DB else u
    rms=lambda a: np.sqrt(np.mean(make(a)**2))
    rmax=rms(60.0); t=min(target,0.95*rmax)
    a=brentq(lambda a: rms(a)-t,1e-8,60.0,xtol=1e-12); return make(a),bool(target>0.95*rmax)
def sbuild(e,rho,nx=1100,seed=11):   # = shells_massout.sbuild with N=24 and no monotonicity assert
    w=np.diff(e); Rout=e[-1]; rm=0.5*(e[1:]+e[:-1]); DCON=rho/(M/(4*np.pi*rm**2*((R-R0)/N))/RHO_REF0)
    sig=S_OUT*np.arange(N)/(N-1)*(rho/RHOREF)**(1/3) if P==1 else S_OUT*(np.arange(N)/(N-1))**P*DCON**BETA; lam=LFAC*w
    ny=int(nx*np.sin(np.radians(SPAN)))+2; dx=Rout/nx; Xg=(np.arange(nx)+0.5)*dx; Yg=(np.arange(ny)+0.5)*dx
    XX,YY=np.meshgrid(Xg,Yg); RR=np.hypot(XX,YY); TT=np.degrees(np.arctan2(YY,XX))
    px,py=nx+512,ny+512; K=np.hypot(*np.meshgrid(2*np.pi*np.fft.fftfreq(px,dx),2*np.pi*np.fft.fftfreq(py,dx)))
    Wn=np.fft.fft2(np.random.default_rng(seed).normal(size=(py,px))); D=np.full((ny,nx),np.nan); rows=[]
    kb=np.geomspace(2*np.pi/(px*dx)*2,np.pi/dx,240); idx=np.digitize(K.ravel(),kb); cnt=np.maximum(np.bincount(idx,minlength=len(kb)+1)[1:len(kb)],1)
    for i in range(N):
        m=(RR>=e[i])&(RR<e[i+1])&(TT<=SPAN)
        if sig[i]==0:
            D[m]=rho[i]; rows.append(dict(shell=1,rho=float(rho[i]),sigma_target=0.,sigma_measured=0.,spacing_target=float(lam[i]),spacing_measured=None,pct={k_:float(rho[i]) for k_ in ("min","p2.5","p16","p84","p97.5","max")},mean_over_rho=1.0)); continue
        ki=2*np.pi/lam[i]; T=np.where(K>0,np.exp(-np.log(np.maximum(K,1e-12)/ki)**2/(4*BW**2))/np.maximum(K,1e-12),0.)
        G=np.real(np.fft.ifft2(Wn*T)); g=G[:ny,:nx]; g=(g-g[m].mean())/g[m].std()
        if P==1:
            s=np.sqrt(np.log(1+sig[i]**2)); L=np.exp(s*g-s*s/2); L=L/L[m].mean(); D[m]=rho[i]*L[m]; capped=False
        else:   # bounded, symmetric contrast: delta in [-DB, +DB], mean exactly 0 -> excess = deficit, no negative density, no tails
            dl,capped=bounded_delta(g[m],sig[i]); L=np.ones_like(g); L[m]=1+dl; D[m]=rho[i]*L[m]
        pw=np.bincount(idx,weights=(np.abs(np.fft.fft2(G))**2*K**2).ravel(),minlength=len(kb)+1)[1:len(kb)]/cnt
        kpk=np.sqrt(kb[:-1]*kb[1:])[np.argmax(pw)]; r_=L[m]-1
        rows.append(dict(shell=i+1,rho=float(rho[i]),sigma_target=float(sig[i]),sigma_measured=float(np.sqrt(np.mean(r_**2))),max_abs_delta=float(np.abs(r_).max()),capped=capped,spacing_target=float(lam[i]),
            spacing_measured=float(2*np.pi/kpk),mean_over_rho=float(L[m].mean()),
            pct=dict(zip(("min","p2.5","p16","p84","p97.5","max"),[float(v) for v in rho[i]*np.r_[L[m].min(),np.percentile(L[m],[2.5,16,84,97.5]),L[m].max()]]))))
    return dict(rows=rows,D=D,core=(RR<e[0])&(TT<=SPAN),asp=ny/nx,dcon=DCON)
RES={}
for key,title,lab,psi in SC:
    c=brentq(lambda c: march(c,psi)[-1]-R,0.05,2.0,xtol=1e-14); e=march(c,psi)
    RES[key]=dict(c=c,e=e,w=np.diff(e))
    for tag,mass in (("eq",np.full(N,M)),("mo",MI_MO)):
        rho=mass/vol(e)/RHO_REF0; RES[key][tag]=dict(rho=rho,mass=mass,S=sbuild(e,rho)); print(key,tag,"done",flush=True)
# ---------------- colour scale: v2 scale (0.033..1.35) widened only to hold the clump extremes ----------------
allD=np.concatenate([RES[k][t]["S"]["D"][~np.isnan(RES[k][t]["S"]["D"])] for k in RES for t in ("eq","mo")])
lo,hi=np.percentile(allD,[0.1,99.9]); print("pixel range",allD.min(),allD.max(),"0.1/99.9 pct",lo,hi)
if P==1: VMIN=min(0.033,10**(np.floor(np.log10(lo)*10)/10)); VMAX=max(1.35,hi)
else:
    lo,hi=np.percentile(allD,[0.2,99.8]); VMIN=10**(np.floor(np.log10(lo)*20)/20); VMAX=max(10**(np.ceil(np.log10(hi)*20)/20),1.03*max(RES[k][t]['rho'].max() for k in RES for t in ('eq','mo')))   # v2 (never below the densest smooth shell): fit the new data (0.2/99.8 pct, 0.05-dex rounding)
if P!=1:   # bounded field: voids reach ~0, so the log scale is set by the shell means and the 2x caps, voids below 0.01 shown as darkest
    VMIN=0.01; VMAX=1.35
NORM=LogNorm(VMIN,VMAX)
FA=float(np.mean(allD>VMAX)); FB=float(np.mean(allD<VMIN)); print("NORM",VMIN,VMAX,"frac above",FA,"below",FB)
allp=[]; 
for k in RES:
    for t in ("eq","mo"): allp+= [r["pct"]["min"] for r in RES[k][t]["S"]["rows"]]+[r["pct"]["max"] for r in RES[k][t]["S"]["rows"]]
YL=(min(allp)/1.4,max(allp)*1.4) if P==1 else (0.006,3.5)
allq=np.concatenate([RES[k][t]["S"]["dcon"] for k in RES for t in ("eq","mo")]); QYL=(allq.min()/1.12,allq.max()*1.12)
allsm=[r["sigma_measured"] for k in RES for t in ("eq","mo") for r in RES[k][t]["S"]["rows"]]; SYL=(0,max(allsm)*1.12) if P==1 else (0,1.08)
allsp=[r["spacing_measured"] for k in RES for t in ("eq","mo") for r in RES[k][t]["S"]["rows"] if r["spacing_measured"]]; PYL=(min(allsp)/1.25,max(allsp)*1.25)
# ---------------- geometry (inches) ----------------
W=7.4; LM=0.72; GAPC=0.42; AW=(W-LM-GAPC-0.2)/2; X0=[LM,LM+AW+GAPC]; HW=AW*np.sin(np.radians(SPAN))*1.01
HG=0.62; GP=0.07; NG=3 if P==1 else 5; TOP=0.72; CT=0.26; XL=0.34; LEG=0.62+(0.42 if P!=1 else 0); CBH=0.85
H=TOP+CT+HW+NG*HG+(NG-1)*GP+XL+LEG+CBH
TH=np.linspace(0,np.radians(SPAN),400)
def page(key,title,lab,psi):
    d=RES[key]; e=d["e"]; xe=e/R
    fig=plt.figure(figsize=(W,H)); Ax=lambda x,y,ww,h: fig.add_axes([x/W,y/H,ww/W,h/H])
    fig.text(0.5,1-0.22/H,(f"Late structure growth — sine mode {title.replace('half wavelength, ','').replace('full wavelength, ','')}" if P!=1 else f"Structure, sine mode {title}"),ha="center",va="center",fontsize=11.5 if P!=1 else 12,fontweight="bold")
    fig.text(0.5,1-0.48/H,f"equal mass (left) vs mass rising outward (right);  {N} shells, w ∝ 1 + Aψ,  {lab},  A = {A:g}",ha="center",va="center",fontsize=8.5,color="0.25")
    for col,(tag,ctitle) in enumerate([("eq","(a) equal mass per shell"),("mo",f"(b) mass rising outward (×{QM:.4f}/shell)")]):
        x0=X0[col]; S=d[tag]["S"]; rows=S["rows"]; rho=d[tag]["rho"]
        yw=H-TOP-CT-HW; aw=Ax(x0,yw,AW,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
        aw.set_title(ctitle,loc="left",fontsize=9.2,fontweight="bold",pad=3)
        img=CMAP(NORM(np.clip(np.nan_to_num(S["D"],nan=1.0),VMIN,VMAX))); img[np.isnan(S["D"])]=(1,1,1,0); img[S["core"]]=(0.86,0.86,0.86,1)
        aw.imshow(img,origin="lower",extent=[0,1,0,S["asp"]],interpolation="antialiased",zorder=1)
        for b in xe: aw.plot(b*np.cos(TH),b*np.sin(TH),color="white",lw=0.35,zorder=3)
        aw.text(0.0,HW/AW*0.99,f"outer-shell clumpiness {rows[-1]['sigma_measured']:.2f}\nclump spacing {min(r['spacing_measured'] for r in rows if r['spacing_measured']):.2f} – {max(r['spacing_measured'] or 0 for r in rows):.2f}",fontsize=8.3,va="top",ha="left",zorder=4)
        axs=[]; y=yw
        for _ in range(NG): y-=HG; axs.append(Ax(x0,y,AW,HG)); y-=GP
        if P==1: a1,a2,a3=axs
        else:
            aW,aD,a1,a2,a3=axs; w=d["w"]; cw=d["c"]
            for k in range(N): aW.fill_between([xe[k],xe[k+1]],0.4,w[k]/cw,color=CMAP(NORM(rho[k])),alpha=0.55,lw=0)
            for k in range(N): aW.hlines(w[k]/cw,xe[k],xe[k+1],color="0.25",lw=0.9)
            rr_=np.linspace(R0,R,600); aW.plot(rr_/R,1+A*psi(X(rr_)),color="#6c3483",lw=1.2); aW.axhline(1,color="0.6",lw=0.6); aW.set_ylim(0.4,1.6); aW.set_yticks([0.5,1,1.5])
            dc=S["dcon"]
            for k in range(N): aD.hlines(dc[k],xe[k],xe[k+1],color=CMAP(NORM(rho[k])),lw=2.4,zorder=3)
            if tag=="mo": de=d["eq"]["S"]["dcon"]; aD.step(np.r_[xe[:-1],xe[-1]],np.r_[de,de[-1]],where="post",color="0.45",lw=0.8,ls="--",alpha=0.6,zorder=2)
            aD.axhline(1,color="0.5",lw=0.7); aD.set_yscale("log"); aD.set_ylim(*QYL); aD.set_yticks([0.5,1,2]); aD.set_yticklabels(["0.5","1","2"]); aD.yaxis.set_minor_formatter(mt.NullFormatter())
            for a in axs:   # faint Doppler shading: blue = compressed (narrower than mean), red = expanded
                for k in range(N):
                    pk=psi(X(0.5*(e[k]+e[k+1]))); a.axvspan(xe[k],xe[k+1],color=("#2c78cf" if pk<0 else "#e0584a"),alpha=0.06+0.06*abs(pk),lw=0,zorder=0)
        for k in range(N):
            p=rows[k]["pct"]; cc=CMAP(NORM(rho[k]))
            a1.fill_between([xe[k],xe[k+1]],p["p2.5"],p["p97.5"],color=cc,alpha=0.20,lw=0,zorder=1)
            a1.fill_between([xe[k],xe[k+1]],p["p16"],p["p84"],color=cc,alpha=0.50,lw=0,zorder=2)
            if k>0: a1.errorbar([0.5*(xe[k]+xe[k+1])],[rho[k]],yerr=[[rho[k]-p["min"]],[p["max"]-rho[k]]],fmt="none",ecolor="0.25",elinewidth=0.5,capsize=1.2,capthick=0.5,zorder=3)
        for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM(rho[k])),lw=2.0,zorder=4)
        a1.set_yscale("log"); a1.set_ylim(*YL); a1.yaxis.set_major_locator(mt.LogLocator(numticks=5)); a1.yaxis.set_minor_formatter(mt.NullFormatter())
        sm=[r["sigma_measured"] for r in rows]; st=[r["sigma_target"] for r in rows]
        for k in range(N): a2.hlines(sm[k],xe[k],xe[k+1],color="#7b241c",lw=2.2)
        a2.step(np.r_[xe[:-1],xe[-1]],np.r_[st,st[-1]],where="post",color="0.35",lw=0.7,ls="--")
        if tag=="mo":
            se=[r["sigma_measured"] for r in d["eq"]["S"]["rows"]]; a2.step(np.r_[xe[:-1],xe[-1]],np.r_[se,se[-1]],where="post",color="#7b241c",lw=0.8,ls=":",alpha=0.6)
        if P!=1:
            md=[r.get("max_abs_delta",0) for r in rows]; a2.step(np.r_[xe[:-1],xe[-1]],np.r_[md,md[-1]],where="post",color="#c0392b",lw=0.9)
            a2.axhline(1,color="k",lw=0.6,ls=(0,(1,1)))
            for k in range(N): a1.hlines(2*rho[k],xe[k],xe[k+1],color="k",lw=0.6,ls=(0,(1,1)),zorder=5)
            a1.text(0.01,0.04,"floor ρ = 0\nbelow axis",transform=a1.transAxes,ha="left",va="bottom",fontsize=7,color="0.3")
        a2.set_ylim(*SYL)
        sp=[r["spacing_measured"] or np.nan for r in rows]
        for k in range(1,N): a3.hlines(sp[k],xe[k],xe[k+1],color="#1f3a93",lw=2.2)
        a3.step(np.r_[xe[:-1],xe[-1]],np.r_[LFAC*d["w"],LFAC*d["w"][-1]],where="post",color="0.35",lw=0.7,ls="--")
        a3.set_yscale("log"); a3.set_ylim(*PYL); a3.set_yticks([0.2,0.3,0.5]); a3.set_yticklabels(["0.2","0.3","0.5"]); a3.yaxis.set_minor_formatter(mt.NullFormatter())
        if col==0:
            a1.set_ylabel("density\n(÷ ref.)"); a2.set_ylabel("clumpiness\n(rms)"); a3.set_ylabel("clump\nspacing")
            if P!=1: aW.set_ylabel("width\n1 + Aψ"); aD.set_ylabel("ρ ÷ ρ$_{unmod}$")
        for a in axs:
            for bb in xe: a.axvline(bb,color="0.82",lw=0.5,zorder=0)
            a.set_xlim(0,1); a.grid(axis="y",alpha=0.25)
        for a in axs[:-1]: plt.setp(a.get_xticklabels(),visible=False)
        a3.set_xlabel("normalised radius  r / R",labelpad=2); ylow=y+GP
    hh=[Line2D([],[],color="0.3",lw=2.0),Patch(color="0.45",alpha=0.5),Patch(color="0.45",alpha=0.2),Line2D([],[],color="0.25",lw=0.5,marker="_",ms=5),
        Line2D([],[],color="#7b241c",lw=2.2),Line2D([],[],color="#7b241c",lw=0.8,ls=":",alpha=0.6),Line2D([],[],color="#1f3a93",lw=2.2),Line2D([],[],color="0.35",lw=0.7,ls="--"),Line2D([],[],color="0.82",lw=0.5)]
    ll=["shell density (= mean)","16–84% of pixels","2.5–97.5% of pixels","min–max in shell",
        "clumpiness rms(δ)","equal-mass rms (in b)","clump spacing (0.6 × width)",("rule: 0.8 (i/23)³ ρ/ρ$_{unmod}$; 0.6 w" if P!=1 else "rule target (dashed)"),"shell boundaries"]
    if P!=1:
        hh=[Line2D([],[],color="#6c3483",lw=1.2),Line2D([],[],color="0.5",lw=0.7),Line2D([],[],color="0.45",lw=0.8,ls="--",alpha=0.6),Patch(color="#2c78cf",alpha=0.15),Patch(color="#e0584a",alpha=0.15)]+hh
        hh+=[Line2D([],[],color="#c0392b",lw=0.9),Line2D([],[],color="k",lw=0.6,ls=(0,(1,1)))]
        ll=["wave 1 + Aψ (bars: widths)","1 = unmodulated density","equal-mass value (in b)","compressed (narrow shells)","expanded (wide shells)"]+ll+["max |δ| in shell","bounds: 2× mean; |δ| = 1"]
    fig.legend(hh,ll,loc="upper center",bbox_to_anchor=(0.5,(ylow-XL-0.02)/H),ncol=3,frameon=False,fontsize=7.8,borderaxespad=0,columnspacing=1.0,handlelength=1.8,labelspacing=0.3)
    cax=fig.add_axes([X0[0]/W,0.45/H,(X0[1]+AW-X0[0])/W,0.14/H])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation="horizontal",extend="both")
    cb.set_label("density ÷ innermost shell of the unmodulated equal-mass set  (one log scale for all four structure pages)",fontsize=8.5,labelpad=2); cb.ax.tick_params(labelsize=8)
    fig.savefig(f"resonant_sine_struct_pair_{key}{SUF}.png",dpi=250); plt.close(fig)
for s in SC: page(*s)
info=dict(norm=[VMIN,VMAX],frac_pixels_above=FA,frac_pixels_below=FB,page_in=[W,H],rho_ref="unmodulated equal-mass 24-shell set, same shell index",cases={})
for k in RES:
    info["cases"][k]={t:dict(rows=RES[k][t]["S"]["rows"]) for t in ("eq","mo")}
    for t in ("eq","mo"):
        rw=RES[k][t]["S"]["rows"]; print(k,t,"sigma meas",[round(rw[j]["sigma_measured"],4) for j in (0,1,2,3,5,7,11,17,23)],"maxd",[round(rw[j].get("max_abs_delta",0),3) for j in (0,3,7,11,17,23)],"capped",[r_["shell"] for r_ in rw if r_.get("capped")],"target last %.3f"%rw[-1]["sigma_target"],"max mean err %.1e"%max(abs(r["mean_over_rho"]-1) for r in rw))
# mass conservation check on the final drawn field
import numpy as _np
ny_,nx_=RES[SC[0][0]]["eq"]["S"]["D"].shape; dx_=R/nx_; XX_,YY_=_np.meshgrid((_np.arange(nx_)+0.5)*dx_,(_np.arange(ny_)+0.5)*dx_); RR_=_np.hypot(XX_,YY_); TT_=_np.degrees(_np.arctan2(YY_,XX_))
merr=0.0
for k in RES:
    e=RES[k]["e"]
    for t in ("eq","mo"):
        Dd=RES[k][t]["S"]["D"]; rho=RES[k][t]["rho"]
        for i in range(N):
            m=(RR_>=e[i])&(RR_<e[i+1])&(TT_<=SPAN); merr=max(merr,abs(Dd[m].mean()/rho[i]-1))
info["mass_conservation_max_rel_err"]=float(merr); info["growth_P"]=P; print("mass conservation max rel err",merr)
json.dump(info,open(f"resonant_sine_struct{SUF}.json","w"),indent=1); print("page",W,round(H,2))
