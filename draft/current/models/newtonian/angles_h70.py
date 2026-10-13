"""v_LOS and local slope dv/dR vs R per angle, D=3.2 Gpc, H=70 runs (expanding & collapsing)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open('levelD_h70.py').read().split("plt.rcParams")[0].replace("DS=[0.8,1.6,2.4,3.2]","DS=[3.2]"))
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
D=3.2; RM=0.246*D; Rs=np.linspace(0.01,RM,150); A=[0,45,90,135,180,270]
cols=["#b2182b","#ef8a62","#4d4d4d","#67a9cf","#2166ac","#1b7837"]
for mode in MODES:
    s=MODES[mode][2]; T="Expanding" if s>0 else "Collapsing"
    for m in Ms:
        rf,v=field(mode,m); Fm=lambda x:np.interp(x,rf,v); xo=np.array([D,0.])
        def vl(R,a):
            nn=np.array([np.cos(np.radians(a)),np.sin(np.radians(a))]); xp=xo+np.outer(R,nn); rp=np.linalg.norm(xp,axis=1)
            return ((Fm(rp)[:,None]*xp/rp[:,None]-Fm(D)*xo/D)@nn)/1e3
        for kind in ("v","slope"):
            f,ax=plt.subplots(figsize=(8,6)); out=[]
            for a,cl in zip(A,cols):
                y=vl(Rs,a); 
                if kind=="slope": y=np.gradient(y,Rs*1000)
                ax.plot(Rs*1000,y,color=cl,lw=3.5,ls=":" if a==270 else "-",label="%d°"%a); out.append((a,y))
            ax.plot(Rs*1000,s*70*Rs*1000 if kind=="v" else s*70+0*Rs,"k--",lw=2,label="no mass")
            ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_xlim(0,RM*1000)
            ax.set_ylabel("v_LOS (km/s)" if kind=="v" else "dv/dR (km/s/Mpc)")
            ax.set_title("%s, %s"%(T,"no mass" if m==0 else "M = %g×"%m)); ax.legend(frameon=False,ncol=2,fontsize=13,loc="best")
            f.tight_layout(); f.savefig("h70_%s_angles_%s_M%g.png"%(mode,kind,m),dpi=110); plt.close(f)
            if kind=="slope": print(mode,m," ".join("%d°:%.1f→%.1f"%(a,y[0],y[-1]) for a,y in out), "| sign-flip:",any((np.sign(y)!=s).any() for a,y in out))
