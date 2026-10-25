import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Wedge, Circle
from scipy.optimize import brentq
c=1.0; p=1.9684                        # 12-shell 3%-lock content law from density_r2 (m ∝ r^p)
cases=[("observer faster than c: D = 3 + 1.5·t",1.5),("observer slower than c: D = 3 + 0.5·t",0.5)]
ROY=["red","orange","gold","green","blue","indigo","violet"]
def geom(v,t):
    d=3+v*t; R=c*t; n=2*t; L=np.linspace(d+R,d-R,n+1); return d,R,n,L
P={}
for _,v in cases:
    for t in (1,2,3):
        d,R,n,L=geom(v,t); m=(L[:-1]+L[1:])/2; V=(4/3)*np.pi*(L[:-1]**3-L[1:]**3)
        P[(v,t)]=brentq(lambda q:(m[-1]**q/V[-1])/(m[0]**q/V[0])-1.03,0.1,6)
norm=mcolors.Normalize(0.95,1); cmap=plt.cm.Greys
fig,axs=plt.subplots(3,2,figsize=(18,24))
log=[]
for col,(title,v) in enumerate(cases):
    for row,t in enumerate((1,2,3)):
        ax=axs[row,col]; d,R,n,L=geom(v,t); Rrec=3.0; inn=d-R
        assert abs(L[-1]-(d-R))<1e-9 and abs(L[0]-(d+R))<1e-9 and abs(R-c*t)<1e-12
        p=P[(v,t)]; mm=(L[:-1]+L[1:])/2; VV=(4/3)*np.pi*(L[:-1]**3-L[1:]**3); rr=mm**p/VV; rmax=rr[-1]
        assert len(L)-1==n==2*t and np.allclose(np.diff(L),-c) and abs(rr[-1]/rr[0]-1.03)<1e-10 and np.all(np.diff(rr)>0)
        assert (inn>Rrec) if v>c else (inn<Rrec)
        log.append((v,t,d,inn,d+R,p,rr[0]/rr[-1],(mm[0]/mm[-1])**p))
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
        if inn>Rrec:
            ax.add_patch(Wedge((0,0),inn,-90,90,width=inn-Rrec,fc="white",ec="none",zorder=0.6))
            u=np.array([np.cos(-1.25),np.sin(-1.25)]); ax.annotate("",tuple(inn*u),tuple(Rrec*u),arrowprops=dict(arrowstyle="<->",color="crimson",lw=1.8),zorder=6)
            ax.text(*(1.02*inn*u+np.array([0.2,-0.3])),f"gap {inn-Rrec:.1f}\n(D−R > R_rec)",color="crimson",fontsize=10,zorder=6,bbox=dict(fc="white",ec="none",alpha=0.85))
        else:
            ax.add_patch(Wedge((0,0),Rrec,-90,90,width=Rrec-inn,fc="crimson",alpha=0.25,ec="none",zorder=3.5))
            ax.text(Rrec*0.5,Rrec*0.75,f"overlap {Rrec-inn:.1f}\nshells inside R_rec\n(D−R < R_rec)",color="crimson",fontsize=10,ha="center",bbox=dict(fc="white",ec="none",alpha=0.8),zorder=6)
        ax.text(Rrec*0.45,-Rrec*0.55,f"R_rec = {Rrec:.0f}\n(fixed)",ha="center",fontsize=10,bbox=dict(fc="white",ec="none"),zorder=6)
        cen=[]
        for k in range(t,0,-1):
            ax.add_patch(Circle((d,0),c*k,fc=ROY[k-1],alpha=0.3,ec="none",zorder=3)); cen.append(d)
            ax.add_patch(Circle((d,0),c*k,fc="none",ec=ROY[k-1],lw=2.5,zorder=3))
        assert all(abs(x-d)<1e-12 for x in cen)
        ax.text(d,R+0.15,"seeing now",ha="center",fontsize=11,bbox=dict(fc="white",ec="none",alpha=0.85),zorder=6)
        if t>1: ax.text(d,-(t-1)+0.3,"already seen",ha="center",fontsize=9,weight="bold",zorder=6)
        ax.plot(d,0,"ko",ms=8,zorder=7); ax.text(d+0.1,0.25,"OBS",fontsize=12,weight="bold",zorder=7)
        ax.annotate("",(d+0.5*v,0),(d,0),arrowprops=dict(arrowstyle="->",lw=2.5),zorder=7)
        ax.plot([0,13],[0,0],color="gray",lw=0.7); ax.axvline(0,color="k",lw=2); ax.text(-0.2,0.2,"CoE",ha="right")
        y=-(d+R)-0.6; ax.annotate("",(d,y),(0,y),arrowprops=dict(arrowstyle="<->")); ax.plot([d,d],[y,0],"k:",lw=1)
        ax.text(d/2,y-0.45,f"D_{t} = {d:.1f}",ha="center",fontsize=12)
        ax.annotate("",(d+R,-R-0.35),(d-R,-R-0.35),arrowprops=dict(arrowstyle="<->",color="tab:blue")); ax.text(d,-R-0.8,"2R",ha="center",color="tab:blue",weight="bold")
        ax.set_xlim(-1.5,11.5); ax.set_ylim(-12.5,11); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title((title+"\n" if row==0 else "")+f"t = {t}  ({n} shells, width c·Δt)",fontsize=14)
cax=fig.add_axes([0.93,0.3,0.012,0.4]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label("density ÷ innermost band (content ∝ r^p, p solved per panel for 3%)",fontsize=12)
fig.suptitle("Red hatched = fixed recombination line R_rec = 3. Same scale in all panels.\nR = c·t; visible circle touches D−R (youngest shell) and D+R (shell 1). Arrows radial from CoE, length ∝ distance.",fontsize=14)
plt.subplots_adjust(left=0.03,right=0.91,top=0.93,bottom=0.02,hspace=0.08,wspace=0.02)
plt.savefig("two_columns_speed.png",dpi=150)
for l in log: print("v=%.1f t=%d D=%.1f D-R=%.1f D+R=%.1f p=%.4f rhoout/in=%.4f content out/in=%.3f"%l)
