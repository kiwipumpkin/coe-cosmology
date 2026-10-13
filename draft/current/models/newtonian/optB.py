import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("heat_combo.py").read().split("N=42")[0])          # calc, vesc, P, curve, ang
exec("th="+open("cone_smooth.py").read().split("\nth=",1)[1].split("names=[]")[0])
cases=[dict(m=1,Dv=58.2,sv=14.3,rst=59.364,col="#1f77b4",n="curve 2"),dict(m=2,Dv=40,sv=8,rst=41,col="#17becf",n="curve 3"),dict(m=4,Dv=40,sv=2,rst=40.05,col="#e6550d",n="curve 7")]
for q in cases:
    q["v"]=curve(q["m"],q["Dv"],q["sv"],q["rst"]); q["h"]=q["v"][:i180]; q["fr"],q["mn"]=stats(q["h"]); q["y"]=np.sqrt(1-q["Dv"]/q["rst"])
    q["h60"]=K60@q["h"]; q["v60"]=np.concatenate([q["h60"],q["h60"][::-1][1:len(q["v"])-i180+1]])
    q["lab"]="%s: M=%g×, D=%g, R=%g, stall=%.4g Gpc"%(q["n"],q["m"],q["Dv"],q["sv"],q["rst"])
    print(q["lab"],"frac %.1f mean %.4f y %.3f cone60 %.1f"%(q["fr"],q["mn"],q["y"],stats(q["h60"])[0]))
def base(ax):
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.axhline(0,color="k",lw=.8); ax.set_xlabel("angle θ (deg): 0° = away from CoE, 180° = toward CoE")
# 1: one panel per case
f,axs=plt.subplots(1,3,figsize=(16,4.6),sharey=True)
for ax,q in zip(axs,cases):
    v=q["v"]; ax.plot(ang,v,color=q["col"],lw=2.2); ax.fill_between(ang,v,0,where=v>0,color="#f6b8b8"); ax.fill_between(ang,v,0,where=v<0,color="#b8cdf6")
    ax.axhline(q["mn"],color="k",ls=":",lw=1.5); base(ax); ax.set_title(q["lab"]+"\n%.0f%% receding, sky mean %+.3g c (dotted)"%(q["fr"],q["mn"]),fontsize=9)
axs[0].set_ylabel("v_LOS / c  (+ = receding)"); f.suptitle("Option B1 — collapsing shells, LOS velocity vs angle (red = redshift, blue = blueshift window)",fontsize=11)
f.tight_layout(); f.savefig("optB1_vlos_angle.png",dpi=115); plt.close(f)
# 2 overlay
f,ax=plt.subplots(figsize=(9.5,5))
for q in cases:
    ax.plot(ang,q["v"],color=q["col"],lw=1.3,alpha=.6); ax.plot(ang,q["v60"],color=q["col"],lw=2.6,label="%s  — raw %.0f%%, 60° cone %.0f%% receding"%(q["n"],q["fr"],stats(q["h60"])[0]))
base(ax); ax.set_ylabel("v_LOS / c"); ax.set_title("Option B2 — raw (thin) vs 60° Migkas cone average (thick)",fontsize=11)
ax.legend(frameon=False,fontsize=9,loc="upper center",bbox_to_anchor=(.5,-.15)); f.savefig("optB2_raw_vs_cone60.png",dpi=115,bbox_inches="tight"); plt.close(f)
# 3 Mollweide
lon=np.radians(np.linspace(-180,180,361)); lat=np.radians(np.linspace(-90,90,181)); LO,LA=np.meshgrid(lon,lat)
theta=180-np.degrees(np.arccos(np.clip(np.cos(LA)*np.cos(LO),-1,1)))
f=plt.figure(figsize=(17,4.6)); a=max(abs(q["h"]).max() for q in cases)
for k,q in enumerate(cases):
    ax=f.add_subplot(1,3,k+1,projection="mollweide"); Vm=np.interp(theta,th,q["h"])
    pc=ax.pcolormesh(LO,LA,Vm,cmap="RdBu_r",vmin=-a,vmax=a,shading="auto"); ax.contour(LO,LA,Vm,[0],colors="k",linewidths=.7)
    ax.plot(0,0,"k*",ms=12); ax.text(0,np.radians(10),"CoE",ha="center",fontweight="bold"); ax.set_xticklabels([]); ax.set_yticklabels([]); ax.grid(alpha=.3)
    ax.set_title(q["lab"]+"\n%.0f%% receding"%q["fr"],fontsize=9)
cb=f.colorbar(pc,ax=f.axes,orientation="horizontal",fraction=.05,pad=.06,label="v_LOS / c  (red = receding, blue = approaching); CoE at centre, anti-CoE at edges")
cb.ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5)); cb.ax.tick_params(labelsize=9); f.suptitle("Option B3 — all-sky maps",fontsize=11); f.savefig("optB3_mollweide.png",dpi=115,bbox_inches="tight"); plt.close(f)
# 5 Hubble per case
f,axs=plt.subplots(1,3,figsize=(16,4.6)); cols=plt.cm.coolwarm(np.linspace(1,0,5))
for ax,q in zip(axs,cases):
    M,D,rst=q["m"],q["Dv"],q["rst"]; vf=lambda r:-np.sqrt(2*GM0*M/GPC*np.maximum(0,1/r-1/rst))/c_ if False else None
    def vfun(r): return -np.sqrt(2*GM0*M/GPC*np.maximum(0,1/r-1/rst))/2.998e8
    xo=np.array([D,0.]); vo=vfun(D)*xo/D; Rs=np.linspace(0.2,min(25,0.9*D),150)
    for an,cl in zip((0,45,90,135,180),cols):
        nn=np.array([np.cos(np.radians(an)),np.sin(np.radians(an))])
        vv=[((vfun(np.linalg.norm(xo+r*nn))*(xo+r*nn)/np.linalg.norm(xo+r*nn)-vo)@nn) for r in Rs]; ax.plot(Rs,vv,color=cl,lw=2,label="θ = %d°"%an)
    ax.axhline(0,color="k",lw=.8); ax.axvline(q["sv"],color="0.6",ls=":"); ax.set_xlabel("galaxy distance R (Gpc)"); ax.set_title(q["lab"],fontsize=9)
axs[0].set_ylabel("v_LOS / c  (+ = receding)"); axs[0].legend(frameon=False,fontsize=8); f.suptitle("Option B5 — Hubble diagram by direction (dotted = R used in other plots)",fontsize=11)
f.tight_layout(); f.savefig("optB5_hubble.png",dpi=115); plt.close(f)
# 6 region
X=np.linspace(0.02,0.9,50); Y=np.linspace(0.02,0.99,50); F=np.zeros((50,50)); Sm=F.copy()
for i_,y_ in enumerate(Y):
    for j_,x_ in enumerate(X): F[i_,j_],Sm[i_,j_]=stats(curve(1,58.2,x_*58.2,58.2/(1-y_**2))[:i180])
f,ax=plt.subplots(figsize=(7.5,5.5)); ok=((F>50)&(Sm>0)).astype(float)
ax.contourf(X,Y,ok,[.5,1.5],colors=["#f6b8b8"]); ax.contour(X,Y,ok,[.5],colors="#c0392b",linewidths=1.5)
for q in cases: ax.plot(q["sv"]/q["Dv"],q["y"],"*",color=q["col"],mec="k",ms=17,label="%s (R/D=%.2f, |v|/v_esc=%.2f)"%(q["n"],q["sv"]/q["Dv"],q["y"]))
ax.text(.45,.08,"looks expanding\n(>50% redshifted AND sky-mean redshift)",fontsize=9.5,color="#7b1f1f"); ax.text(.45,.8,"blueshift dominated",fontsize=10,color="#2e5fa8")
ax.set_xlabel("R / D  (galaxy distance / our distance to CoE)"); ax.set_ylabel("|v| / v_esc  (our infall speed / escape speed)"); ax.legend(frameon=False,fontsize=8.5,loc="upper left")
ax.set_title("Option B6 — where a collapsing universe looks expanding",fontsize=10); f.savefig("optB6_region.png",dpi=115,bbox_inches="tight"); plt.close(f)
