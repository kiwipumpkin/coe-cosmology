import numpy as np, pickle, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D; from scipy.interpolate import griddata
plt.rcParams.update({"font.size":17,"axes.titlesize":17,"axes.labelsize":18,"legend.fontsize":12,"xtick.labelsize":15,"ytick.labelsize":15})
A=np.arange(0,181.); w=np.sin(np.radians(A)); Ds=[10,15,20,30,40,50,60]; Ms=[0,0.1,0.3,1,3,10,30,100]
S={(D,m):pickle.load(open("s_D%g_rs%g.pkl"%(D,m),"rb")) for D in Ds for m in Ms}
RL=[1.,2.,3.,4.5,6.,8.,10.,12.,14.3]
PL,PLe,SH,SHe=67.4,0.5,73.04,1.04; MID=(PL+SH)/2; HW=(SH-PL)/2; TF=100*(SH-PL)/MID
def stats(y):
    ok=~np.isnan(y); mn=(w[ok]*y[ok]).sum()/w[ok].sum(); return mn,100*(np.nanmax(y)-np.nanmin(y))/mn,y/mn,ok
okcase=lambda D,m: S[(D,m)]['Mfrac']>=0.95
def tension(ax,X,Y,Z,kind):
    h=[]
    if kind=="H":
        for lv,ls,lw,col in ((PL,"-",3,"cyan"),(PL-PLe,"--",1.5,"cyan"),(PL+PLe,"--",1.5,"cyan"),(SH,"-",3,"orange"),(SH-SHe,"--",1.5,"orange"),(SH+SHe,"--",1.5,"orange"),(MID,":",3,"white")):
            if np.nanmin(Z)<lv<np.nanmax(Z): ax.contour(X,Y,Z,levels=[lv],colors=col,linestyles=ls,linewidths=lw)
        h=[Line2D([],[],color="cyan",lw=3,label="Planck 2018: 67.4±0.5"),Line2D([],[],color="orange",lw=3,label="SH0ES 2022: 73.04±1.04"),Line2D([],[],color="white",lw=3,ls=":",label="midpoint %.1f±%.1f"%(MID,HW))]
    elif kind=="AA":
        for lv,col in ((PL/MID,"cyan"),(SH/MID,"orange"),(1.0,"white")):
            if np.nanmin(Z)<lv<np.nanmax(Z): ax.contour(X,Y,Z,levels=[lv],colors=col,linewidths=3,linestyles="-" if lv!=1 else ":")
        h=[Line2D([],[],color="cyan",lw=3,label="Planck/mid = %.3f"%(PL/MID)),Line2D([],[],color="orange",lw=3,label="SH0ES/mid = %.3f"%(SH/MID))]
    elif kind=="anis":
        if np.nanmin(Z)<TF<np.nanmax(Z): ax.contour(X,Y,Z,levels=[TF],colors="white",linewidths=3.5)
        h=[Line2D([],[],color="white",lw=3.5,label="(73.04−67.4)/70.2 = %.1f%%"%TF)]
    h.append(Line2D([],[],color="red",lw=3,label="invalid / unreachable edge"))
    ax.legend(handles=h,loc="upper center",bbox_to_anchor=(.5,-.17),ncol=2,frameon=False)
def fillnan(Z,xs,ys):
    Zf=Z.copy()
    for i in range(Zf.shape[0]):
        ok=~np.isnan(Zf[i])
        if ok.any(): Zf[i,~ok]=np.interp(np.array(xs)[~ok],np.array(xs)[ok],Zf[i,ok])
    for j in range(Zf.shape[1]):
        ok=~np.isnan(Zf[:,j])
        if ok.any() and (~ok).any(): Zf[~ok,j]=np.interp(np.array(ys)[~ok],np.array(ys)[ok],Zf[ok,j])
    return Zf
def drawmap(X,Y,Z,V,title,fn,cbl,kind,xl,yl,yt=None):
    fig,ax=plt.subplots(figsize=(8,7.2)); Zf=fillnan(Z,X[0],Y[:,0])
    cf=ax.contourf(X,Y,Zf,levels=24,cmap="viridis"); cb=fig.colorbar(cf); cb.set_label(cbl)
    if (~V).any() and V.any(): ax.contour(X,Y,V.astype(float),levels=[0.5],colors="red",linewidths=3)
    tension(ax,X,Y,Zf,kind)
    if yt is not None: ax.set_yticks(yt[0]); ax.set_yticklabels(yt[1])
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title,fontsize=16); fig.tight_layout(); fig.savefig(fn,dpi=110); plt.close(fig)
Mx=np.arange(len(Ms)); X,Y=np.meshgrid(Ds,Mx); yt=(Mx,[str(m) for m in Ms])
tab=open("sweep_table.txt","w")
for R in (14.3,3.0):
    Z={k:np.full(X.shape,np.nan) for k in ("mean","anis","AA0","AA90","AA180","H0","H90","H180")}; V=np.ones(X.shape,bool); V180=V.copy()
    for i,m in enumerate(Ms):
        for j,D in enumerate(Ds):
            y=S[(D,m)]['H'][R]; mn,an,q,ok=stats(y)
            Z["mean"][i,j]=mn; Z["anis"][i,j]=an; 
            for a in (0,90,180): Z["AA%d"%a][i,j]=q[a]; Z["H%d"%a][i,j]=y[a]
            V[i,j]=ok.all() and okcase(D,m); V180[i,j]=ok[180] and okcase(D,m)
            yc=S[(D,m)]['cl'][R]; mc,ac,qc,_=stats(yc)
            tab.write("R=%-4g D=%2d m=%-5g exact: mean %.1f anis %.1f%% H 0/90/180 %.1f %.1f %s | A/Aall %.3f %.3f %s | classical mean %.1f anis %.1f%% | valid %d/181%s\n"%(R,D,m,mn,an,y[0],y[90],("%.1f"%y[180]) if ok[180] else "beyond",q[0],q[90],("%.3f"%q[180]) if ok[180] else "beyond",mc,ac,ok.sum(),"" if okcase(D,m) else " INVALID(observer inside core)"))
    t="R = %g Gpc"%R; sx="R%g"%R; xl="D, our areal distance from CoE (Gpc)"; yl="2GM/c² (Gpc)"
    drawmap(X,Y,Z["mean"],V,"sky-mean H (km/s/Mpc), "+t,"map_mean_%s.png"%sx,"km/s/Mpc","H",xl,yl,yt)
    drawmap(X,Y,Z["anis"],V,"anisotropy (max−min)/mean, "+t,"map_anis_%s.png"%sx,"%","anis",xl,yl,yt)
    for a in (0,90,180):
        VV=V180 if a==180 else V
        drawmap(X,Y,Z["AA%d"%a],VV,"A/A_all at %d°, %s"%(a,t),"map_AA%d_%s.png"%(a,sx),"A/A_all","AA",xl,yl,yt)
        drawmap(X,Y,Z["H%d"%a],VV,"H at %d° (km/s/Mpc), %s"%(a,t),"map_H%d_%s.png"%(a,sx),"km/s/Mpc","H",xl,yl,yt)
P=np.array([[S[(D,m)]['dprop'] for D in Ds] for m in Ms]); MR=np.array([[S[(D,m)]['maxd'][180] for D in Ds] for m in Ms]); VC=np.array([[okcase(D,m) for D in Ds] for m in Ms])
drawmap(X,Y,P,VC,"our proper distance to CoE today (Gpc)","map_dist_to_CoE.png","Gpc","none","D areal (Gpc)","2GM/c² (Gpc)",yt)
drawmap(X,Y,MR,VC,"max reachable R toward CoE (180°), Gpc","map_maxR180.png","Gpc","none","D areal (Gpc)","2GM/c² (Gpc)",yt)
with open("sweep_table.txt","a") as f:
    for i,m in enumerate(Ms): f.write("m=%-5g proper dist to CoE: %s | maxR180: %s\n"%(m," ".join("%.1f"%x for x in P[i])," ".join("%.1f"%x for x in MR[i])))
tab.close()
# ---- ratio maps: x=R/Dprop, y=(2GM/c^2)/Dprop  (log axes), scattered data from all D x R
pts=[];vals={k:[] for k in ("mean","anis","AA0","AA90","AA180","H0","H90","H180")};vv=[];v180=[]
for D in Ds:
  for m in Ms[1:]:
    if not okcase(D,m): continue
    Dp=S[(D,m)]['dprop']
    for R in RL:
        y=S[(D,m)]['H'][R]; mn,an,q,ok=stats(y); pts.append((np.log10(R/Dp),np.log10(m/Dp)))
        vals["mean"].append(mn); vals["anis"].append(an)
        for a in (0,90,180): vals["AA%d"%a].append(q[a]); vals["H%d"%a].append(y[a])
        vv.append(float(ok.all())); v180.append(float(ok[180]))
pts=np.array(pts); gx=np.linspace(pts[:,0].min(),pts[:,0].max(),120); gy=np.linspace(pts[:,1].min(),pts[:,1].max(),120); GX,GY=np.meshgrid(gx,gy)
def grid(v):
    v=np.array(v); ok=~np.isnan(v); Zl=griddata(pts[ok],v[ok],(GX,GY),method="linear"); Zn=griddata(pts[ok],v[ok],(GX,GY),method="nearest"); return np.where(np.isnan(Zl),Zn,Zl)
Vg=grid(vv)>0.5; V180g=grid(v180)>0.5
xt=[0.03,0.1,0.3,1]; ytk=[0.003,0.01,0.03,0.1,0.3,1,3,10]
def rmap(Z,V,title,fn,cbl,kind):
    fig,ax=plt.subplots(figsize=(8,7.2)); cf=ax.contourf(GX,GY,Z,levels=24,cmap="viridis"); cb=fig.colorbar(cf); cb.set_label(cbl)
    if (~V).any() and V.any(): ax.contour(GX,GY,V.astype(float),levels=[0.5],colors="red",linewidths=3)
    tension(ax,GX,GY,Z,kind); ax.plot(pts[:,0],pts[:,1],"k.",ms=3,alpha=.5)
    ax.set_xticks(np.log10([x for x in xt if gx[0]<=np.log10(x)<=gx[-1]])); ax.set_xticklabels([str(x) for x in xt if gx[0]<=np.log10(x)<=gx[-1]])
    yy=[y for y in ytk if gy[0]<=np.log10(y)<=gy[-1]]; ax.set_yticks(np.log10(yy)); ax.set_yticklabels([str(y) for y in yy])
    ax.set_xlabel("R / D  (both proper distances)"); ax.set_ylabel("(2GM/c²) / D"); ax.set_title(title,fontsize=16); fig.tight_layout(); fig.savefig(fn,dpi=110); plt.close(fig)
rmap(grid(vals["mean"]),Vg,"sky-mean H (km/s/Mpc)","ratio_mean.png","km/s/Mpc","H")
rmap(grid(vals["anis"]),Vg,"anisotropy (max−min)/mean","ratio_anis.png","%","anis")
for a in (0,90,180):
    rmap(grid(vals["AA%d"%a]),V180g if a==180 else Vg,"A/A_all at %d°"%a,"ratio_AA%d.png"%a,"A/A_all","AA")
    rmap(grid(vals["H%d"%a]),V180g if a==180 else Vg,"H at %d° (km/s/Mpc)"%a,"ratio_H%d.png"%a,"km/s/Mpc","H")
# angle line plots & HvsR
for R in (14.3,3.0):
  for D in (20,40,60):
    for m in (1,10,100):
        r=S[(D,m)]; fig,ax=plt.subplots(figsize=(8,6)); ax.axhline(1,color="k",ls="--",lw=2.5,label="no mass")
        for y,lab,st,col in ((r['H'][R],"exact GR","-","tab:blue"),(r['cl'][R],"classical",":","tab:orange")): ax.plot(A,stats(y)[2],ls=st,color=col,lw=3.5,label=lab)
        ax.set_xlim(0,180); ax.set_xticks(range(0,181,45)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("A / A_all")
        ax.set_title("D=%d, 2GM/c²=%g Gpc, R=%g Gpc"%(D,m,R),fontsize=17); ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=14)
        fig.tight_layout(); fig.savefig("AA_D%d_m%g_R%g.png"%(D,m,R),dpi=110); plt.close(fig)
for D,m in ((20,10),(40,1),(60,100)):
    r=S[(D,m)]['HvsR']; Rs=sorted(r); fig,ax=plt.subplots(figsize=(8,6)); ax.axhline(70,color="k",ls="--",lw=2.5,label="no mass")
    for k,(a,col) in enumerate(((0,"tab:blue"),(90,"tab:green"),(180,"tab:red"))): ax.plot(Rs,[r[x][k] for x in Rs],color=col,lw=3.5,label="%d°"%a)
    ax.set_xlabel("distance R (Gpc)"); ax.set_ylabel("H (km/s/Mpc)"); ax.set_title("exact GR: D=%d, 2GM/c²=%g Gpc"%(D,m),fontsize=17)
    ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=14); fig.tight_layout(); fig.savefig("HvsR_D%d_m%g.png"%(D,m),dpi=110); plt.close(fig)
print("ok")
