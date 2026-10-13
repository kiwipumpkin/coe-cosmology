import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
D=58.2; X=np.linspace(0.02,0.9,45); Y=np.linspace(0,1,41)
def frac(x,y):
    rst=np.inf if y>=1 else D/(1-y**2); yy=curve(1,D,x*D,rst); return 100*(w*(yy>0)).sum()/w.sum()
Z=np.array([[frac(x,y) for x in X] for y in Y]); print("range %.1f–%.1f"%(Z.min(),Z.max()))
f,ax=plt.subplots(figsize=(7,5.5)); pc=ax.pcolormesh(X,Y,Z,vmin=0,vmax=100,cmap="RdBu_r",shading="auto")
cs=ax.contour(X,Y,Z,[50],colors="k",linewidths=1.6); ax.clabel(cs,fmt="50%%")
for p in cs.allsegs[0]:
    for x in (0.05,0.25,0.5,0.9):
        i=np.argmin(abs(p[:,0]-x)); print("x=%.2f y50=%.3f"%(p[i,0],p[i,1]))
ax.set_xlabel("R / D   (galaxy distance / our distance to CoE)"); ax.set_ylabel("|v| / v_esc(D)   (our infall speed / escape speed at D)")
plt.colorbar(pc,label="% of sky receding (solid-angle weighted)")
ax.set_title("Collapse (free fall from stall radius, all shells same law): % sky receding\nstall/D = 1/(1−y²); independent of M and of D at fixed R/D",fontsize=10)
f.savefig("heat_recede_RD_vvesc.png",dpi=120,bbox_inches="tight")
