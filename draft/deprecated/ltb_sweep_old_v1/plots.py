import numpy as np, pickle, glob, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
A=np.arange(0,181.); w=np.sin(np.radians(A)); Ds=[10,15,20,30,40,50,60]; Ms=[0,0.1,0.3,1,3,10,30,100]
S={(D,m):pickle.load(open("s_D%g_rs%g.pkl"%(D,m),"rb")) for D in Ds for m in Ms}
def stats(y):
    ok=~np.isnan(y); mn=(w[ok]*y[ok]).sum()/w[ok].sum(); return mn,(np.nanmax(y)-np.nanmin(y))/mn,y/mn,ok
Mx=np.arange(len(Ms)); X,Y=np.meshgrid(Ds,Mx)
def cmap(Z,valid,title,fn,cbl):
    fig,ax=plt.subplots(figsize=(8,6.5))
    Zf=Z.copy(); 
    # fill invalid by nearest valid in same row for full-chart fill; boundary drawn
    for i in range(Zf.shape[0]):
        ok=~np.isnan(Zf[i]); 
        if ok.any(): Zf[i,~ok]=np.interp(np.array(Ds)[~ok],np.array(Ds)[ok],Zf[i,ok])
    for j in range(Zf.shape[1]):
        ok=~np.isnan(Zf[:,j])
        if ok.any() and (~ok).any(): Zf[~ok,j]=np.interp(Mx[~ok],Mx[ok],Zf[ok,j])
    cf=ax.contourf(X,Y,Zf,levels=20,cmap="viridis"); cb=fig.colorbar(cf); cb.set_label(cbl)
    if (~valid).any() and valid.any(): ax.contour(X,Y,valid.astype(float),levels=[0.5],colors="red",linewidths=3)
    ax.set_yticks(Mx); ax.set_yticklabels([str(m) for m in Ms]); ax.set_xlabel("D, our distance from CoE (Gpc)"); ax.set_ylabel("2GM/c² (Gpc)")
    ax.set_title(title,fontsize=16); fig.tight_layout(); fig.savefig(fn,dpi=110); plt.close(fig)
tab=open("sweep_table.txt","w")
for R in (14.3,3.0):
    Zm=np.full(X.shape,np.nan); Za=Zm.copy(); Z0=Zm.copy(); Z90=Zm.copy(); Z180=Zm.copy(); V=np.zeros(X.shape,bool)
    for i,m in enumerate(Ms):
        for j,D in enumerate(Ds):
            y=S[(D,m)]['H'][R]; mn,an,q,ok=stats(y); Zm[i,j]=mn; Za[i,j]=100*an; Z0[i,j]=q[0]; Z90[i,j]=q[90]; Z180[i,j]=q[180]; V[i,j]=ok.all()
            if R==14.3 or True:
                yc=S[(D,m)]['cl'][R]; mc,ac,qc,_=stats(yc)
                tab.write("R=%-4g D=%2d m=%-5g exact: mean %.1f anis %.1f%% A/Aall 0 %.3f 90 %.3f 180 %s | classical: mean %.1f anis %.1f%% 180 %.3f | valid angles %d/181\n"%(R,D,m,mn,100*an,q[0],q[90],("%.3f"%q[180]) if ok[180] else "beyond",mc,100*ac,qc[180],ok.sum()))
    t="R=%g Gpc"%R; sfx="R%g"%R
    cmap(Zm,V,"sky-mean H (km/s/Mpc), "+t,"map_mean_%s.png"%sfx,"km/s/Mpc")
    cmap(Za,V,"anisotropy (max−min)/mean %%, "+t,"map_anis_%s.png"%sfx,"%")
    for nm,Z in (("0",Z0),("90",Z90),("180",Z180)):
        cmap(Z,V if nm!="180" else ~np.isnan(Z180),"A/A_all at %s°, %s"%(nm,t),"map_AA%s_%s.png"%(nm,sfx),"A/A_all")
tab.close()
# proper distance & max reach 180
P=np.array([[S[(D,m)]['dprop'] for D in Ds] for m in Ms]); MR=np.array([[S[(D,m)]['maxd'][180] for D in Ds] for m in Ms]); MR0=np.array([[S[(D,m)]['maxd'][0] for D in Ds] for m in Ms])
cmap(P,np.ones(X.shape,bool),"our proper distance to CoE today (Gpc)","map_dist_to_CoE.png","Gpc")
cmap(MR,np.ones(X.shape,bool),"max reachable R toward CoE (180°), Gpc","map_maxR180.png","Gpc")
with open("sweep_table.txt","a") as f:
    for i,m in enumerate(Ms): f.write("m=%-5g proper dist to CoE: %s | maxR180: %s | maxR0: %s\n"%(m," ".join("%.1f"%x for x in P[i])," ".join("%.1f"%x for x in MR[i])," ".join("%.1f"%x for x in MR0[i])))
# line plots
for R in (14.3,3.0):
  for D in (20,40,60):
    for m in (1,10,100):
        r=S[(D,m)]; fig,ax=plt.subplots(figsize=(8,6)); ax.axhline(1,color="k",ls="--",lw=2.5,label="no mass")
        for y,lab,st,col in ((r['H'][R],"exact GR","-","tab:blue"),(r['cl'][R],"classical",":","tab:orange")):
            ax.plot(A,stats(y)[2],ls=st,color=col,lw=3.5,label=lab)
        ax.set_xlim(0,180); ax.set_xticks(range(0,181,45)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("A / A_all")
        ax.set_title("D=%d, 2GM/c²=%g Gpc, R=%g Gpc"%(D,m,R),fontsize=17); ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,-.18))
        fig.tight_layout(); fig.savefig("AA_D%d_m%g_R%g.png"%(D,m,R),dpi=110); plt.close(fig)
for D,m in ((20,10),(40,1),(60,100)):
    r=S[(D,m)]['HvsR']; Rs=sorted(r); fig,ax=plt.subplots(figsize=(8,6)); ax.axhline(70,color="k",ls="--",lw=2.5,label="no mass")
    for k,(a,col) in enumerate(((0,"tab:blue"),(90,"tab:green"),(180,"tab:red"))): ax.plot(Rs,[r[x][k] for x in Rs],color=col,lw=3.5,label="%d°"%a)
    ax.set_xlabel("distance R (Gpc)"); ax.set_ylabel("H (km/s/Mpc)"); ax.set_title("exact GR: D=%d, 2GM/c²=%g Gpc"%(D,m),fontsize=17)
    ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18)); fig.tight_layout(); fig.savefig("HvsR_D%d_m%g.png"%(D,m),dpi=110); plt.close(fig)
print("ok")
