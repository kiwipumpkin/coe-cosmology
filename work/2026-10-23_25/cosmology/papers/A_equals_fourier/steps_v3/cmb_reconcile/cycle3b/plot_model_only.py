import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, sys
sys.path.insert(0,'.'); import compare as C
d=json.load(open('af_tables.json'))
fig,axs=plt.subplots(1,2,figsize=(13,4.8))
for ax,name,title in zip(axs,('omb','free'),('Omega_s = 0.049 fit (D = 49 Gpc)','Omega_s free fit (D = 94 Gpc)')):
    for Rs,c in [(0.07,'C0'),(0.21,'C1'),(0.47,'C2'),(2.0,'C3'),(6.0,'C4'),(14.3,'k')]:
        A=np.asarray(C.A_curve(name,Rs)); A=A/np.nanmean(A)
        ax.plot(C.TH1,A,color=c,lw=1.8,label=f'sampled at {Rs} Gpc')
    ax.set_title(title); ax.set_xlabel('Galactic longitude (deg)'); ax.set_ylabel('model A / mean'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig('model_only_vs_distance.png',dpi=120)
