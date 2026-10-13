import numpy as np, pickle, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":17,"axes.labelsize":18,"xtick.labelsize":15,"ytick.labelsize":15})
A=np.arange(0,181.); w=np.sin(np.radians(A)); P="../sweep/s_D%g_rs%g.pkl"
bad={(10,30),(10,100),(15,100)}
def aa(y):
    ok=~np.isnan(y); return y/((w[ok]*y[ok]).sum()/w[ok].sum())
def plot(cases,lab,title,fn,R,cmap):
    fig,ax=plt.subplots(figsize=(8,6.5)); ax.axhline(1,color="k",ls="--",lw=2.5,label="no mass")
    cs=plt.get_cmap(cmap)(np.linspace(.35,1,len(cases))); info=[]
    for (D,m),c in zip(cases,cs):
        if (D,m) in bad: continue
        r=pickle.load(open(P%(D,m),"rb")); y=aa(r['H'][R]); ok=~np.isnan(y)
        L=lab(D,m)+("" if ok.all() else " (ends %d°)"%A[ok].max())
        ax.plot(A,y,color=c,lw=3.5,label=L); info.append((D,m,y[0],A[np.nanargmin(y)],np.nanmin(y),y[ok][-1],A[ok].max()))
    ax.set_xlim(0,180); ax.set_xticks(range(0,181,45)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("A / A_all")
    ax.set_title(title+", R=%g Gpc"%R,fontsize=17); ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.5,-.17),fontsize=14)
    fig.tight_layout(); fig.savefig(fn,dpi=110); plt.close(fig)
    for i in info: print(fn,"D=%g m=%g AA0=%.3f min@%d=%.3f end=%.3f@%d"%i)
for R in (3.0,14.3):
    for m in (1,10,100): plot([(D,m) for D in (15,20,30,40,60)],lambda D,m:"D=%d"%D,"exact GR: 2GM/c²=%g Gpc"%m,"fam_m%g_R%g.png"%(m,R),R,"Blues")
    for D in (20,40,60): plot([(D,m) for m in (1,3,10,30,100)],lambda D,m:"2GM/c²=%g"%m,"exact GR: D=%d Gpc"%D,"fam_D%d_R%g.png"%(D,R),R,"Reds")
