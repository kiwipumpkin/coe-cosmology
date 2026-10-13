import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
M,R,RST=2.0,5.0,56.0; Ds=np.linspace(15,50.9,150); Lv=[20,30,40,45,48,50]
def st(D):
    y=curve(M,D,R,RST); return (w*y).sum()/w.sum(), 100*(w*(y>0)).sum()/w.sum()
S=np.array([st(D) for D in Ds]); L=np.array([st(D) for D in Lv])
f,(a1,a2)=plt.subplots(2,1,figsize=(8.5,7),sharex=True,gridspec_kw=dict(height_ratios=[1.6,1]))
for sg,col,nm in ((1,"#c0392b","expanding (mirror, outward)"),(-1,"#2e5fa8","collapsing (free fall from r_stall = 56 Gpc)")):
    a1.plot(Ds,(S[:,0] if sg==-1 else -S[:,0]),color=col,lw=2.2,label=nm)
    a1.plot(Lv,(L[:,0] if sg==-1 else -L[:,0]),"o",color=col,mec="k")
    a2.plot(Ds,(S[:,1] if sg==-1 else 100-S[:,1]),color=col,lw=2.2); a2.plot(Lv,(L[:,1] if sg==-1 else 100-L[:,1]),"o",color=col,mec="k")
a1.axhline(0,color="k",lw=1); a2.axhline(50,color="k",lw=1,ls=":")
a1.set_ylabel("sky-average v_LOS / c\n(solid-angle weighted)"); a2.set_ylabel("% of sky receding"); a2.set_xlabel("observer distance to CoE, D (Gpc)")
a1.legend(frameon=False,fontsize=9); a1.set_title("Sky-average LOS velocity and % receding vs D (dots = level-set D values)\nM = 2× 7.5×10⁵² kg, R = 5 Gpc, r_stall = 56 Gpc (fixed absolute)",fontsize=10)
f.tight_layout(); f.savefig("skymean_vs_D_c.png",dpi=130)
i=np.where(np.diff(np.sign(S[:,0])))[0]; j=np.where(np.diff(np.sign(S[:,1]-50)))[0]
print("collapse mean=0 at D≈%.1f; 50%% at D≈%.1f"%(Ds[i[0]],Ds[j[0]]) if len(i) and len(j) else "no crossing")
