import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("heat_combo.py").read().split("N=42")[0])          # calc, vesc, P, curve, ang
exec("th="+open("cone_smooth.py").read().split("\nth=",1)[1].split("names=[]")[0])
def calc2(M,D,R,v):
    ve=vesc(M,D)
    if not (0<v<ve) or R>=0.95*D: return (np.nan,)*3
    y=v/ve; rst=D/(1-y*y); h=curve(M,D,R,rst)[:i180]
    return stats(h)[0],stats(h)[1],stats(K60@h)[0]
use={p[0]:p for p in P}
N=60
for name,lg in (("RD_vvesc",False),("DM_Rv",True),("MD_RD",True)):
    _,xl,yl,xr,yr,mp,held=use[name]
    if name=="RD_vvesc": yr=(0.02,0.99)
    sp=np.geomspace if lg else np.linspace; X=sp(*xr,N); Y=sp(*yr,N)
    F=np.full((N,N),np.nan); S=F.copy(); C=F.copy()
    for a,yv in enumerate(Y):
        for b,xv in enumerate(X): F[a,b],S[a,b],C[a,b]=calc2(*mp(xv,yv))
    both=(F>50)&(S>0); one=((F>50)^(S>0))&~np.isnan(F)
    f,ax=plt.subplots(figsize=(8,6))
    reg=np.where(np.isnan(F),np.nan,np.where(both,2,np.where(one,1,0)))
    ax.pcolormesh(X,Y,np.where(np.isnan(F),1,np.nan),cmap=matplotlib.colors.ListedColormap(["0.85"]),shading="auto")
    ax.pcolormesh(X,Y,reg,cmap=matplotlib.colors.ListedColormap(["white","#fdf0d5","#f6b8b8"]),vmin=-.5,vmax=2.5,shading="auto")
    c1=ax.contour(X,Y,F,[50],colors="#c0392b",linewidths=2); c2=ax.contour(X,Y,S,[0],colors="#2e5fa8",linewidths=2)
    c3=ax.contour(X,Y,C,[50],colors="#c0392b",linewidths=1.5,linestyles="--")
    if lg: ax.set_xscale("log"); ax.set_yscale("log")
    from matplotlib.patches import Patch; from matplotlib.lines import Line2D
    ax.legend(handles=[Patch(fc="#f6b8b8",label="BOTH: >50% receding AND sky mean > 0"),Patch(fc="#fdf0d5",label="only one of the two"),Patch(fc="white",ec="0.7",label="neither"),Patch(fc="0.85",label="excluded (v ≥ v_esc or R ≥ 0.95 D)"),
        Line2D([],[],color="#c0392b",lw=2,label="50% receding (raw)"),Line2D([],[],color="#2e5fa8",lw=2,label="sky mean = 0"),Line2D([],[],color="#c0392b",ls="--",label="50% receding, 60° Migkas cone")],
        frameon=False,fontsize=8.5,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.12))
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    if name=="RD_vvesc":
        sec=ax.secondary_yaxis("right",functions=(lambda y:y**2,lambda e:np.sqrt(np.clip(e,0,None)))); sec.set_ylabel("v² D / (2GM)  = (|v|/v_esc)²")
    ax.set_title("Collapse (stall-radius free fall): where >50%% recedes AND sky mean > 0\nheld: %s"%held,fontsize=10)
    f.savefig("region_%s.png"%name,dpi=120,bbox_inches="tight"); plt.close(f)
    if name=="RD_vvesc":
        for xv in (0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9):
            b=np.argmin(abs(X-xv)); col=lambda Z,t:(Y[np.where(Z[:,b]>t)[0]].max() if (Z[:,b]>t).any() else np.nan)
            print("R/D=%.2f  y50=%.3f  ymean0=%.3f  yboth=%.3f  y50_cone60=%.3f"%(X[b],col(F,50),col(S,0),col(np.where(both,1,0),.5),col(C,50)))
    else:
        print(name,"both-region cells",both.sum(),"of",(~np.isnan(F)).sum())
