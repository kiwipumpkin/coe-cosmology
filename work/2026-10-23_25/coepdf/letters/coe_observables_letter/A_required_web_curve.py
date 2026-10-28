# (b) A/A_all before/after adding strip_avg(delta), with Migkas data and errors. Data from A_required_web_solve.py.
import numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
R=np.load("A_required_web_results.npz",allow_pickle=True)
l,y=R["l"],R["y"];plt.figure(figsize=(11,5))
plt.errorbar(l,y,yerr=[R["elo"],R["ehi"]],fmt="k.",ms=4,lw=.7,label="Migkas Fig.6 A/A_all")
plt.plot(l,y,"k-",lw=.8,label="before (connect the dots)");plt.plot(l,y+R["D"]-R["D"].mean(),"r.-",lw=1,label="after: y + strip_avg(delta(l)) (mean removed)")
plt.plot(l,y+R["Dfit"]-R["Dfit"].mean(),"g--",lw=1,label="after: cluster/void version (mean removed)")
plt.xlabel("l [deg]");plt.ylabel("A/A_all");plt.xlim(0,360);plt.grid(alpha=.3);plt.legend(fontsize=8)
plt.title(f"Required change: rms {R['D'].std():.3f} in A/A_all");plt.tight_layout();plt.savefig("A_required_web_curve.png",dpi=140)
