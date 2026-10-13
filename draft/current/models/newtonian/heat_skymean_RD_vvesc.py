import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
D=58.2; X=np.linspace(0.02,0.9,45); Y=np.linspace(0,1,41)
def frac(x,y):
    rst=np.inf if y>=1 else D/(1-y**2); yy=curve(1,D,x*D,rst); return (w*yy).sum()/w.sum()
Z=np.array([[frac(x,y) for x in X] for y in Y]); print("range %.4f – %.4f"%(Z.min(),Z.max()))
f,ax=plt.subplots(figsize=(7,5.5)); pc=ax.pcolormesh(X,Y,Z,norm=matplotlib.colors.TwoSlopeNorm(0),cmap="RdBu_r",shading="auto")
cs=ax.contour(X,Y,Z,[0],colors="k",linewidths=1.6); ax.clabel(cs,fmt="0")
for p in cs.allsegs[0]:
    for x in (0.05,0.25,0.5,0.9):
        i=np.argmin(abs(p[:,0]-x)); print("x=%.2f y0=%.3f"%(p[i,0],p[i,1]))
ax.set_xlabel("R / D   (galaxy distance / our distance to CoE)"); ax.set_ylabel("|v| / v_esc(D)   (our infall speed / escape speed at D)")
plt.colorbar(pc,label="sky-average v_LOS / c  (solid-angle weighted, + = receding)")
ax.set_title("Collapse (free fall from stall radius, all shells same law): sky-average v_LOS / c\nstall/D = 1/(1−y²); M = 1× 7.5×10⁵² kg, D = 58.2 Gpc (speeds scale as √(M/D))",fontsize=10)
f.savefig("heat_skymean_RD_vvesc.png",dpi=120,bbox_inches="tight")
