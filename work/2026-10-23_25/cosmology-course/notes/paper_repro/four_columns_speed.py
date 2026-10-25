import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Wedge, Circle
from scipy.optimize import brentq
c=1.0; p=1.9684                        # 12-shell 3%-lock content law from density_r2 (m ∝ r^p)
cases=[("1) observer 1.5c, R_rec fixed at 3 → gap",1.5,0),("2) observer 0.5c, R_rec fixed at 3 → overlap",0.5,0),("3) observer exactly c, R_rec fixed at 3 = D−R",1.0,2),("4) observer 0.5c, R_rec shrinks at (c−v) = 0.5c per step → R_rec = D−R",0.5,1)]
ROY=["red","orange","gold","green","blue","indigo","violet"]
def geom(v,t):
    d=3+v*t; R=c*t; n=2*t; L=np.linspace(d+R,d-R,n+1); return d,R,n,L
P={}
for _,v,_m in cases:
    for t in (1,2,3):
        d,R,n,L=geom(v,t); m=(L[:-1]+L[1:])/2; V=(4/3)*np.pi*(L[:-1]**3-L[1:]**3)
        P[(v,t)]=brentq(lambda q:(m[-1]**q/V[-1])/(m[0]**q/V[0])-1.03,0.1,6)
norm=mcolors.Normalize(0.95,1); cmap=plt.cm.Greys
fig,axs=plt.subplots(3,4,figsize=(34,24))
log=[]
for col,(title,v,mode) in enumerate(cases):
    for row,t in enumerate((1,2,3)):
        ax=axs[row,col]; d,R,n,L=geom(v,t); inn=d-R; Rrec=(3.0-(c-v)*t) if mode==1 else 3.0
        assert abs(L[-1]-(d-R))<1e-9 and abs(L[0]-(d+R))<1e-9 and abs(R-c*t)<1e-12
        p=P[(v,t)]; mm=(L[:-1]+L[1:])/2; VV=(4/3)*np.pi*(L[:-1]**3-L[1:]**3); rr=mm**p/VV; rmax=rr[-1]
        assert len(L)-1==n==2*t and np.allclose(np.diff(L),-c) and abs(rr[-1]/rr[0]-1.03)<1e-10 and np.all(np.diff(rr)>0)
        assert (abs(Rrec-inn)<1e-12 and abs((3.0-Rrec)/t-(c-v))<1e-12) if mode==1 else ((abs(Rrec-3)<1e-12 and abs(inn-Rrec)<1e-12) if mode==2 else ((inn>Rrec) if v>c else (inn<Rrec)))
        log.append((mode,v,t,d,inn,d+R,p,rr[0]/rr[-1],(mm[0]/mm[-1])**p))
        for i in range(1,n+1):
            ro,ri=L[i-1],L[i]; m=(ro+ri)/2; rho=m**p/((4/3)*np.pi*(ro**3-ri**3))/rmax
            ax.add_patch(Wedge((0,0),ro,-90,90,width=ro-ri,fc=cmap(norm(rho)),alpha=0.55,ec="none",zorder=0))
            ax.text(-0.2,-m,str(i),ha="right",va="center",fontsize=11,weight="bold")
            if i%2==1:
                for ang in np.deg2rad([-50,0,50]):
                    u=np.array([np.cos(ang),np.sin(ang)]); ln=0.09*m
                    ax.annotate("",tuple((m+ln/2)*u),tuple((m-ln/2)*u),arrowprops=dict(arrowstyle="->",color="navy",lw=1.6),zorder=2)
        for s in L: ax.add_patch(Wedge((0,0),s,-90,90,width=1e-4,ec="k",lw=1))
        ax.add_patch(Wedge((0,0),min(Rrec,inn),-90,90,fc="white",ec="none",zorder=0.5))
        ax.add_patch(Wedge((0,0),Rrec,-90,90,fc="none",ec="crimson",hatch="////",lw=2.5,zorder=4))
        if mode in (1,2): pass
        elif inn>Rrec:
            ax.add_patch(Wedge((0,0),inn,-90,90,width=inn-Rrec,fc="white",ec="none",zorder=0.6))
            u=np.array([np.cos(-1.25),np.sin(-1.25)]); ax.annotate("",tuple(inn*u),tuple(Rrec*u),arrowprops=dict(arrowstyle="<->",color="crimson",lw=1.8),zorder=6)
            pass
        else:
            ax.add_patch(Wedge((0,0),Rrec,-90,90,width=Rrec-inn,fc="crimson",alpha=0.25,ec="none",zorder=3.5))
            pass
        pass
        cen=[]
        for k in range(t,0,-1):
            ax.add_patch(Circle((d,0),c*k,fc=ROY[k-1],alpha=0.3,ec="none",zorder=3)); cen.append(d)
            ax.add_patch(Circle((d,0),c*k,fc="none",ec=ROY[k-1],lw=2.5,zorder=3))
        assert all(abs(x-d)<1e-12 for x in cen)
        pass
        pass
        ax.plot(d,0,"ko",ms=8,zorder=7)
        ax.annotate("",(d+0.5*v,0),(d,0),arrowprops=dict(arrowstyle="->",lw=2.5),zorder=7)
        ax.plot([0,13],[0,0],color="gray",lw=0.7); ax.axvline(0,color="k",lw=2); ax.text(-0.2,0.2,"CoE",ha="right")
        y=-(d+R)-0.6; ax.annotate("",(d,y),(0,y),arrowprops=dict(arrowstyle="<->")); ax.plot([d,d],[y,0],"k:",lw=1)
        pass
        ax.annotate("",(d+R,-R-0.35),(d-R,-R-0.35),arrowprops=dict(arrowstyle="<->",color="tab:blue"))
        ax.set_xlim(-1.5,11.5); ax.set_ylim(-12.5,11); ax.set_aspect("equal"); ax.axis("off")
        rel=("R_rec = D−R" if mode in (1,2) else (f"gap {inn-Rrec:.1f}" if inn>Rrec else f"overlap {Rrec-inn:.1f}"))
        ax.set_title((title+"\n" if row==0 else "")+f"t = {t}: {n} shells · D = {d:.1f} · D−R = {inn:.1f} · R_rec = {Rrec:.1f} · {rel}",fontsize=12)
cax=fig.add_axes([0.94,0.3,0.012,0.4]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label("density ÷ innermost band (content ∝ r^p, p solved per panel for 3%)",fontsize=12)
fig.suptitle("Red hatched = recombination line R_rec (fixed at 3 in columns 1–3; shrinking at 0.5c per step in column 4). Same scale in all panels.\nR = c·t; visible circle touches D−R (youngest shell) and D+R (shell 1). Arrows radial from CoE, length ∝ distance.",fontsize=14)
plt.subplots_adjust(left=0.03,right=0.93,top=0.93,bottom=0.07,hspace=0.08,wspace=0.02)
fig.text(0.47,0.012,"Key: black dot + arrow = observer (OBS) and its speed; red hatched = recombination region; red tint = shells inside R_rec (column 2); white band between hatch and shells = gap (column 1);\n"
 "red/orange/gold circles around OBS = light seen at t=1/2/3 (inner = already seen, outermost = seeing now); blue double arrow = 2R across the visible circle; black double arrow = D; navy arrows = radial shell speed; numbers at left = shell number (1 outermost).",ha="center",fontsize=13)
plt.savefig("four_columns_speed.png",dpi=150)
for l in log: print("mode=%d v=%.1f t=%d D=%.1f D-R=%.1f D+R=%.1f p=%.4f rhoout/in=%.4f content out/in=%.3f"%l)
