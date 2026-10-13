"""Structure growth in the three equal-mass shell cases (Chris, 11 Oct). Shells = shells_remake.py (12 shells, inner radius 3, first width 1;
constant / widening q=1.0970 / narrowing q=0.9537; R = 15 / 24 / 12.37; equal mass per shell).
MODEL (simple, stated):
 seed      : at recombination every shell is smooth with the same Gaussian random field, rms 1e-5. Statistics use a 1D RING per shell
             (the shell's mid-radius circle, 8192 points, same random seed); the wedge texture uses a 2D map per shell (same spectrum, same seed)
             for display only. Primordial spectrum: power per log k proportional to k^2 (blue) up to 4 k_J of the densest shell.
 background: all shells expand together (homologous, Einstein-de Sitter-like), a = 1 at recombination, a = 1100 today, t = t_rec a^1.5,
             t_rec = 0.38 Myr. Clump sizes are quoted at recombination size (multiply by a for later physical sizes).
 local density sets growth: Omega_i = rho_i / rho_ref, rho_ref = mean density of the whole constant-spacing system (total mass / total volume,
             r = 3..15), the same reference for all three cases (shell 1 is identical in all three).
             linear growth in an EdS background with local source 4 pi G rho_i = 1.5 Omega_i H^2:
                 delta'' + 2H delta' = 1.5 Omega_eff H^2 delta  ->  delta_k = delta_k,0 a^p,  p = Re[(-1/2 + sqrt(1/4 + 6 Omega_eff)) / 2]
             (p = 1 for Omega = 1, faster for denser shells, slower for thinner ones);
             pressure (Jeans) cut: Omega_eff(k) = Omega_i (1 - (k/k_J,i)^2), k_J,i = 2 pi / lambda_J,i, lambda_J,i = 0.15 (rho_i/rho_1)^-1/2
             (Jeans length at fixed sound speed: thinner shells only grow on larger scales -> structure spreads out; modes below the
             Jeans length oscillate and decay as a^-1/4).
 nonlinear : the time a shell goes nonlinear is when its linear rms reaches 1 (clumps form); beyond that the texture is capped at rms 1.5
             and drawn lognormal, rho/rho_bar = exp(g - s^2/2).
 outputs   : rms clumpiness today (linear theory; > 1 = nonlinear, clumps formed), t_nl, clump spacing = 2 pi / k_peak of the power per log k
             (growing modes only, k <= k_J; at the moment it goes nonlinear, or today if it never does).
             The wedge texture also shows the growing modes only (sub-Jeans seed ripples are left out of the picture), scaled to the shell's rms.
"""
import json,numpy as np,matplotlib,matplotlib.ticker
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,SymLogNorm,LinearSegmentedColormap
from matplotlib.patches import Polygon
from scipy.optimize import brentq,curve_fit
plt.rcParams.update({"font.size":12})
N=12; R0=3.0; W0=1.0; AF=1100.0; TREC=0.38; LJ0=0.15; DELTA0=1e-5
CASES=[("constant","Constant shell spacing",15.00),("widening","Widening shell spacing (geometric growth)",24.00),("narrowing","Narrowing shell spacing (geometric shrink)",12.37)]
def edges(q): return R0+np.r_[0,np.cumsum(W0*q**np.arange(N))]
def q_for(Rout): return 1.0 if abs(Rout-15)<1e-9 else brentq(lambda q: edges(q)[-1]-Rout,0.5,2.0)
def rho_of(e): return 1/(4/3*np.pi*(e[1:]**3-e[:-1]**3))
e15=edges(1.0); RHO1=rho_of(e15)[0]; RHOREF=N/(4/3*np.pi*(e15[-1]**3-e15[0]**3))
def fitrho(xe,rho):
    xm=0.5*(xe[1:]+xe[:-1]); f=lambda x,lA,n,C: np.log(np.exp(lA)*x**(-n)+C)
    p,_=curve_fit(f,xm,np.log(rho),p0=[np.log(rho[0]*xm[0]**2),2.0,1e-3],bounds=([-30,0,0],[30,10,10]),maxfev=20000)
    return (lambda x: np.exp(p[0])*x**(-p[1])+p[2]),dict(A=float(np.exp(p[0])),n=float(p[1]),C=float(p[2]))
def pexp(Om,k,kJ):
    Oe=Om*(1-(k/kJ)**2); return np.real((-0.5+np.sqrt(0.25+6*Oe+0j))/2)
KMAX=4*2*np.pi/LJ0
def ring(Om,rmid,kJ,seed=7,n=8192):
    L=2*np.pi*rmid; m=np.arange(1,n//2); k=m/rmid; ok=k<KMAX
    rng=np.random.default_rng(seed); z=(rng.normal(size=len(m))+1j*rng.normal(size=len(m)))/np.sqrt(2)
    amp=np.where(ok,np.sqrt(k),0.0)                       # power per mode ~ k  -> power per log k ~ k^2
    amp*=DELTA0/np.sqrt(np.sum(amp**2))                   # rms 1e-5 at recombination (sum over modes of |d_k|^2 = variance)
    p=pexp(Om,k,kJ)
    def sig(a): return float(np.sqrt(np.sum((amp*a**p)**2)))
    s_today=sig(AF)
    a_nl=float(np.exp(brentq(lambda la: np.log(sig(np.exp(la))),0,np.log(AF)))) if s_today>1 else np.inf
    a_ref=min(a_nl,AF); Pl=k*(amp*a_ref**p)**2               # power per log k
    g_=k<=kJ; kp=k[g_][np.argmax(Pl[g_])] if g_.any() else k[0]; return s_today,a_nl,2*np.pi/kp,amp,p
# ---------------- per case
RES={}
for key,title,Rout in CASES:
    q=q_for(Rout); e=edges(q); rho=rho_of(e); Om=rho/RHOREF; lam=LJ0*(rho/RHO1)**-0.5; kJ=2*np.pi/lam; rm=0.5*(e[1:]+e[:-1])
    rows=[]
    for i in range(N):
        s,anl,sp,_,_=ring(Om[i],rm[i],kJ[i]); rows.append(dict(shell=i+1,r_in=float(e[i]),r_out=float(e[i+1]),rho_rel=float(rho[i]/rho[0]),Omega=float(Om[i]),
            jeans_len=float(lam[i]),rms_today_linear=s,nonlinear=bool(s>1),a_nl=anl,t_nl_Myr=float(TREC*anl**1.5) if np.isfinite(anl) else None,
            spacing=float(sp),spacing_over_shell_width=float(sp/(e[i+1]-e[i])),p_growth_large_scale=float(pexp(Om[i],0.0,1.0))))
    RES[key]=dict(q=float(q),R=float(Rout),rows=rows)
json.dump(dict(model=__doc__,rho_ref_over_rho1=float(RHOREF/RHO1),cases=RES),open("structure_growth.json","w"),indent=1)
# ---------------- figures
DCM=plt.cm.RdBu_r; DN=SymLogNorm(linthresh=1e-5,vmin=-1,vmax=30,base=10)
CMAPD=LinearSegmentedColormap.from_list("shells",["#08154f","#1747a6","#4f86cf","#8a6db0","#b45a8e","#d45a6a","#9b0f1a"]); NORMD=LogNorm(10**-2.1,1.0)
def texture(key,e,Rout,rows,AW_px=1100,SPAN=45.0):
    nx=AW_px; ny=int(nx*np.sin(np.radians(SPAN)))+2; dx=Rout/nx
    X=(np.arange(nx)+0.5)*dx; Y=(np.arange(ny)+0.5)*dx; XX,YY=np.meshgrid(X,Y); RR=np.hypot(XX,YY); TT=np.degrees(np.arctan2(YY,XX))
    kx=2*np.pi*np.fft.fftfreq(nx,dx); ky=2*np.pi*np.fft.fftfreq(ny,dx); K=np.hypot(*np.meshgrid(kx,ky))
    rng=np.random.default_rng(11); W=np.fft.fft2(rng.normal(size=(ny,nx)))
    img=np.full((ny,nx,4),0.0); img[...,:3]=1
    for i,rw in enumerate(rows):
        m=(RR>=e[i])&(RR<e[i+1])&(TT<=SPAN)
        if not m.any(): continue
        p=pexp(rw["Omega"],K,2*np.pi/rw["jeans_len"]); a=min(rw["a_nl"],AF) if rw["nonlinear"] else AF
        T=np.where((K>0)&(K<=2*np.pi/rw["jeans_len"]),a**p,0.0)   # white 2D seed x growth (growing modes only): power per log k ~ k^2 a^2p, like the ring
        g=np.real(np.fft.ifft2(W*T)); g=g[m]; s_target=min(rw["rms_today_linear"],1.5); g=g/g.std()*s_target
        d=np.exp(g-s_target**2/2)-1
        cols=DCM(DN(d)); sub=img[m]; sub[:]=cols; img[m]=sub
    core=(RR<e[0])&(TT<=SPAN); img[core]=(0.86,0.86,0.86,1)
    return img,ny/nx
for key,title,Rout in CASES:
    D=RES[key]; q=D["q"]; e=edges(q); xe=e/e[-1]; rows=D["rows"]; rho=np.array([r["rho_rel"] for r in rows]); g,pf=fitrho(xe,rho)
    SPAN=45.0; Wf=10.0; LM=1.1; AW=8.6; HW=AW*np.sin(np.radians(SPAN))*1.01
    hcb,gcb,h3,gp,h2,h1=0.22,1.05,2.2,0.22,2.2,2.4; base=1.0; top=0.95; H=base+hcb+gcb+h3+gp+h2+gp+h1+HW+top
    fig=plt.figure(figsize=(Wf,H)); ax_=lambda y,h: fig.add_axes([LM/Wf,y/H,AW/Wf,h/H])
    y3=base+hcb+gcb; y2=y3+h3+gp; y1=y2+h2+gp; yw=y1+h1
    fig.suptitle(f"Structure grown in equal-mass shells — {title}",fontsize=15,fontweight="bold",y=1-0.3/H)
    aw=ax_(yw,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
    img,asp=texture(key,e,Rout,rows); aw.imshow(img,origin="lower",extent=[0,1,0,asp],interpolation="nearest",zorder=1)
    th=np.linspace(0,np.radians(SPAN),300)
    for b in xe: aw.plot(b*np.cos(th),b*np.sin(th),color="0.25",lw=0.8,zorder=3)
    aw.plot([0,np.cos(th[-1])],[0,np.sin(th[-1])],color="0.25",lw=0.8,zorder=3); aw.plot([0,1],[0,0],color="0.25",lw=0.8,zorder=3)
    nl=[r for r in rows if r["nonlinear"]]
    aw.text(0.0,HW/AW*0.97,"structure today (13.8 Gyr after recombination)\nseeded with the same 10⁻⁵ ripples in every shell\n"+
            (f"shells 1–{len(nl)} went nonlinear (clumps formed)" if nl else "no shell went nonlinear"),fontsize=11.5,va="top",ha="left",zorder=4)
    a1=ax_(y1,h1); a2=ax_(y2,h2); a3=ax_(y3,h3); cax=fig.add_axes([LM/Wf,base/H,AW/Wf,hcb/H])
    x=np.linspace(xe[0],1,800)
    for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAPD(NORMD(rho[k])),lw=4)
    a1.plot(x,g(x),"--",color="0.15",lw=1.3,label=rf"smooth fit $\rho={pf['A']:.3g}\,(r/R)^{{-{pf['n']:.2f}}}$"+(rf"$\,+\,{pf['C']:.3g}$" if pf['C']>1e-6 else ""))
    a1.set_yscale("log"); a1.set_ylim(10**-2.3,1.6); a1.set_ylabel("shell density\n(÷ innermost shell)"); a1.legend(loc="upper right",fontsize=10,frameon=False)
    s=np.array([r["rms_today_linear"] for r in rows]); sp=np.array([r["spacing"] for r in rows])
    a2.axhspan(1,1e9,color="#fdebd0",zorder=0); a2.axhline(1,color="#ca6f1e",lw=1)
    for k in range(N): a2.hlines(s[k],xe[k],xe[k+1],color="#7b241c",lw=3.5)
    for k,r in enumerate(rows):
        if r["nonlinear"]: a2.text(0.5*(xe[k]+xe[k+1]),min(s[k]*2.2,2e5),f"{r['t_nl_Myr']:.0f}",ha="center",fontsize=8.5,color="#7b241c")
    a2.set_yscale("log"); a2.set_ylim(1e-7,1e6); a2.set_ylabel("rms clumpiness today\n(linear theory)")
    a2.text(0.015,0.92,"shaded: nonlinear, clumps formed; numbers = Myr after recombination when it went nonlinear",transform=a2.transAxes,fontsize=9,color="#ca6f1e",va="top")
    for k in range(N): a3.hlines(sp[k],xe[k],xe[k+1],color="#1f3a93",lw=3.5)
    a3.set_yscale("log"); a3.set_ylim(0.3,4); a3.set_yticks([0.5,1,2,3]); a3.set_yticklabels(["0.5","1","2","3"]); a3.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a3.set_ylabel("typical clump spacing\n(radius units, at recombination)")
    for a_ in (a1,a2,a3):
        for b in xe: a_.axvline(b,color="0.8",lw=0.8,zorder=0)
        a_.set_xlim(0,1); a_.grid(axis="y",alpha=0.25)
    plt.setp(a1.get_xticklabels(),visible=False); plt.setp(a2.get_xticklabels(),visible=False)
    a3.set_xlabel("normalised radius  r / R   (grey lines = shell boundaries)")
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=DN,cmap=DCM),cax=cax,orientation="horizontal",ticks=[-1,-1e-2,-1e-4,0,1e-4,1e-2,1,10])
    cb.set_label("wedge: density contrast of the grown structure, δ = ρ/ρ̄ − 1\n(same symmetric-log colour scale in all three figures; red = clump, blue = gap)",fontsize=11)
    fig.savefig(f"structure_{key}.png",dpi=130); plt.close(fig)
# comparison
fig,ax=plt.subplots(figsize=(10,6.5)); cols={"constant":"#555555","widening":"#1747a6","narrowing":"#b03050"}
ax.axhspan(1,1e9,color="#fdebd0",zorder=0); ax.axhline(1,color="#ca6f1e",lw=1)
for key,title,Rout in CASES:
    D=RES[key]; e=edges(D["q"]); xe=e/e[-1]; s=[r["rms_today_linear"] for r in D["rows"]]
    nnl=sum(r["nonlinear"] for r in D["rows"])
    ax.step(np.r_[xe[:-1],xe[-1]],np.r_[s,s[-1]],where="post",color=cols[key],lw=2.4,label=f"{key} (R = {Rout:.2f}): {nnl} of 12 shells nonlinear by today")
ax.set_yscale("log"); ax.set_ylim(1e-7,1e5); ax.set_xlim(0,1); ax.grid(alpha=0.3); ax.set_xlabel("normalised radius  r / R"); ax.set_ylabel("rms clumpiness today (linear theory)")
ax.text(0.01,0.97,"shaded: nonlinear (clumps formed)",transform=ax.transAxes,va="top",fontsize=10,color="#ca6f1e")
ax.set_title("Structure growth in equal-mass shells: clumpiness vs radius",fontsize=14,fontweight="bold")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.13),frameon=False,fontsize=11)
fig.savefig("structure_compare.png",dpi=130,bbox_inches="tight"); plt.close(fig)
for key in RES:
    print(key); [print("  %2d rho %.4f Om %.3f p %.2f rms %.2e tnl %s spacing %.3f"%(r["shell"],r["rho_rel"],r["Omega"],r["p_growth_large_scale"],r["rms_today_linear"],r["t_nl_Myr"],r["spacing"])) for r in RES[key]["rows"]]
