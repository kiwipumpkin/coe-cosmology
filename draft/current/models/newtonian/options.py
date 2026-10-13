import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("heat_combo.py").read().split("N=42")[0])          # calc, vesc, P, curve, ang
exec("th="+open("cone_smooth.py").read().split("\nth=",1)[1].split("names=[]")[0])
from matplotlib.patches import Circle
M,D=1.0,58.2; RD=0.25; R=RD*D; rst=1.05*D; ve=vesc(M,D); yv=np.sqrt(1-D/rst)
v=curve(M,D,R,rst); h=v[:i180]; fr,mn=stats(h); print("case y=%.3f frac %.1f mean %.4f vus %.4f"%(yv,fr,mn,yv*ve))
lab="collapsing: stall 1.05·D, R/D = 0.25, D = 58.2 Gpc, M = 7.5×10⁵² kg (|v|/v_esc = %.2f)"%yv
def base(ax):
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg): 0° = away from CoE, 180° = toward CoE"); ax.axhline(0,color="k",lw=.8)
# 1
f,ax=plt.subplots(figsize=(9,4.8)); ax.plot(ang,v,"k",lw=2.2,label="v_LOS / c")
ax.fill_between(ang,v,0,where=v>0,color="#f6b8b8",label="receding (redshift) %.0f%% of sky"%fr); ax.fill_between(ang,v,0,where=v<0,color="#b8cdf6",label="blueshift window")
ax.axhline(mn,color="k",ls=":",lw=1.5,label="sky mean %+.3g c"%mn); base(ax); ax.set_ylabel("v_LOS / c  (+ = receding)")
ax.set_title("Option 1 — LOS velocity vs angle\n"+lab,fontsize=10); ax.legend(frameon=False,fontsize=9,loc="upper center",bbox_to_anchor=(.5,-.15),ncol=2)
f.savefig("option1_vlos_angle.png",dpi=120,bbox_inches="tight"); plt.close(f)
# 2
v60=full2=np.concatenate([K60@h,(K60@h)[::-1][1:len(v)-i180+1]])
f,ax=plt.subplots(figsize=(9,4.8)); ax.plot(ang,v,"0.5",lw=1.8,label="raw (%.0f%% receding)"%fr); ax.plot(ang,v60,"#c0392b",lw=2.4,label="60° Migkas cone average (%.0f%% receding)"%stats(K60@h)[0])
base(ax); ax.set_ylabel("v_LOS / c"); ax.set_title("Option 2 — raw vs 60° cone-averaged\n"+lab,fontsize=10); ax.legend(frameon=False,fontsize=9,loc="upper center",bbox_to_anchor=(.5,-.15),ncol=2)
f.savefig("option2_raw_vs_cone60.png",dpi=120,bbox_inches="tight"); plt.close(f)
# 3 Mollweide, CoE at centre
lon=np.radians(np.linspace(-180,180,361)); lat=np.radians(np.linspace(-90,90,181)); LO,LA=np.meshgrid(lon,lat)
sep=np.degrees(np.arccos(np.clip(np.cos(LA)*np.cos(LO),-1,1)))   # angle from CoE
theta=180-sep; Vm=np.interp(theta,th,h)
f=plt.figure(figsize=(9,5)); ax=f.add_subplot(111,projection="mollweide"); a=np.abs(Vm).max()
pc=ax.pcolormesh(LO,LA,Vm,cmap="RdBu_r",vmin=-a,vmax=a,shading="auto"); ax.contour(LO,LA,Vm,[0],colors="k",linewidths=.8)
ax.plot(0,0,"k*",ms=14); ax.text(0,np.radians(8),"CoE",ha="center",fontsize=10,fontweight="bold"); ax.text(np.radians(178),np.radians(-8),"anti-CoE",ha="right",fontsize=9)
ax.set_xticklabels([]); ax.grid(alpha=.3); plt.colorbar(pc,orientation="horizontal",pad=.05,fraction=.05,label="v_LOS / c  (red = receding, blue = approaching)")
ax.set_title("Option 3 — all-sky map (CoE at centre)\n"+lab,fontsize=10); f.savefig("option3_mollweide.png",dpi=120,bbox_inches="tight"); plt.close(f)
# 4 schematic
f,ax=plt.subplots(figsize=(9,6)); vf=lambda r:-np.sqrt(2*GM0*M/GPC*np.maximum(0,1/r-1/rst))/c
xo=np.array([D,0.]); vo=vf(D)*xo/D; S=60
ax.add_patch(Circle((0,0),D,fill=False,ls="--",color="0.7")); ax.plot(0,0,"k*",ms=16); ax.text(0,-4,"CoE",ha="center",va="top",fontweight="bold")
ax.plot(*xo,"ko",ms=9); ax.text(D+1.5,-3,"us",fontweight="bold")
ax.arrow(*xo,*(vo*S),color="0.55",width=.35,length_includes_head=True)
for an in range(0,360,45):
    nn=np.array([np.cos(np.radians(an)),np.sin(np.radians(an))]); xp=xo+R*nn; rp=np.linalg.norm(xp); vg=vf(rp)*xp/rp; vr=vg-vo; los=vr@nn
    ax.plot(*xp,"o",color="0.3",ms=6); (ax.arrow(*xp,*(vg*S),color="0.55",width=.3,length_includes_head=True) if np.linalg.norm(vg)>1e-6 else ax.text(xp[0],xp[1]-2.2,"at rest",fontsize=7,color="0.4",ha="center"))
    ax.arrow(*xp,*(vr*S),color="#c0392b" if los>0 else "#2e5fa8",width=.35,length_includes_head=True,zorder=5)
    ax.plot([D,xp[0]],[0,xp[1]],color="0.85",lw=.8,zorder=0); ax.text(*(xp+np.array([0,3.2])),"%d°"%an,ha="center",va="center",fontsize=8)
ax.plot([],[],color="0.55",lw=3,label="absolute infall velocity (all toward CoE)"); ax.plot([],[],color="#c0392b",lw=3,label="relative to us, receding"); ax.plot([],[],color="#2e5fa8",lw=3,label="relative to us, approaching")
ax.set_aspect("equal"); ax.set_xlim(-8,D+R+8); ax.set_ylim(-R-8,R+8); ax.axis("off"); ax.legend(frameon=False,fontsize=9,loc="upper center",bbox_to_anchor=(.5,0),ncol=2)
ax.set_title("Option 4 — side-view geometry (arrows ×%d Gpc per c; distances to scale)\n"%S+lab,fontsize=10); f.savefig("option4_schematic.png",dpi=120,bbox_inches="tight"); plt.close(f)
# 5 Hubble diagram
Rs=np.linspace(0.2,25,120); f,ax=plt.subplots(figsize=(8,5)); cols=plt.cm.coolwarm(np.linspace(0,1,5))
for an,cl in zip((0,45,90,135,180),cols[::-1]):
    nn=np.array([np.cos(np.radians(an)),np.sin(np.radians(an))]); vv=[]
    for Rr in Rs:
        xp=xo+Rr*nn; rp=np.linalg.norm(xp); vv.append(((vf(rp)*xp/rp-vo)@nn))
    ax.plot(Rs,vv,color=cl,lw=2.2,label="θ = %d°"%an)
ax.axhline(0,color="k",lw=.8); ax.axvline(R,color="0.6",ls=":",lw=1); ax.text(R,ax.get_ylim()[1]*.92," R used above",fontsize=8,color="0.4")
ax.set_xlabel("galaxy distance R (Gpc)"); ax.set_ylabel("v_LOS / c  (+ = receding)"); ax.set_xlim(0,25)
ax.set_title("Option 5 — Hubble diagram by direction (0° = away from CoE)\ncollapsing, stall 1.05·D, D = 58.2 Gpc, M = 7.5×10⁵² kg",fontsize=10); ax.legend(frameon=False,fontsize=9)
f.savefig("option5_hubble.png",dpi=120,bbox_inches="tight"); plt.close(f)
# 6 region map
X=np.linspace(0.02,0.9,50); Y=np.linspace(0.02,0.99,50); F=np.zeros((50,50)); Sm=F.copy()
for a_,y_ in enumerate(Y):
    for b_,x_ in enumerate(X):
        hh=curve(M,D,x_*D,D/(1-y_**2))[:i180]; F[a_,b_],Sm[a_,b_]=stats(hh)
f,ax=plt.subplots(figsize=(7.5,5.5)); ok=((F>50)&(Sm>0)).astype(float)
ax.contourf(X,Y,ok,[.5,1.5],colors=["#f6b8b8"]); ax.contour(X,Y,ok,[.5],colors="#c0392b",linewidths=1.5)
ax.plot(RD,yv,"k*",ms=16); ax.annotate("this case",(RD,yv),(RD+.1,yv+.3),arrowprops=dict(arrowstyle="->"),fontsize=10)
ax.text(.45,.08,"looks expanding\n(>50% redshifted AND\nsky-mean redshift)",fontsize=10,color="#7b1f1f"); ax.text(.45,.8,"blueshift dominated",fontsize=10,color="#2e5fa8")
ax.set_xlabel("R / D  (galaxy distance / our distance to CoE)"); ax.set_ylabel("|v| / v_esc  (our infall speed / escape speed)")
ax.set_title("Option 6 — where a collapsing universe looks expanding",fontsize=10); f.savefig("option6_region.png",dpi=120,bbox_inches="tight"); plt.close(f)
