import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=json.load(open("cmb-cancel-story.json")); c=np.array(d["c_m"]); phi=d["params"][4]
w=json.load(open("migkas-geometry-fit.json"))["wls_fit"]; l0=w["l0"]
b=np.loadtxt("migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
f=np.loadtxt("migkas-best-fit-curve-0-360-full.csv",delimiter=",",skiprows=1)
bmax=b[np.argmax(b[:,1]),0]; bmin=b[np.argmin(b[:,1]),0]
fmax=f[np.argmax(f[:,1]),0]
plt.rcParams.update({"font.size":15})
fig=plt.figure(figsize=(15,9)); ax=fig.add_subplot(121,projection="polar")
L=np.radians(np.arange(361)); ax.plot(L,0.75+(f[np.arange(361)%360,1]-1)*3,"k",lw=1.5,alpha=.5)
ax.scatter(np.radians(b[:,0]),0.75+(b[:,1]-1)*3,s=14,color="gray",zorder=3)
def ray(a,col,lab,ls="-",lw=3): ax.plot([np.radians(a)]*2,[0,1.35],color=col,ls=ls,lw=lw,label=lab)
ray(l0,"k",f"model shell axis ℓ₀ = {l0:.2f}°")
ray((l0+180)%360,"k","model anti-axis "+f"{(l0+180)%360:.2f}°",ls="--")
ray(phi,"C3",f"SR fit + all Fourier modes φ = {phi:.2f}°",lw=2)
ray(fmax,"C1",f"model curve peak {fmax:.0f}°",ls=":",lw=2.5)
ray(bmax,"C0",f"Migkas bins highest {bmax:.0f}°",ls="-.")
ray(bmin,"C2",f"Migkas bins lowest {bmin:.0f}°",ls="-.")
ax.set_ylim(0,1.35); ax.set_yticklabels([]); ax.set_title("Directions (galactic longitude)",pad=20)
fig.legend(loc="upper center",ncol=2,frameon=False,fontsize=14,bbox_to_anchor=(0.5,1.0))
ax2=fig.add_subplot(122); M=12; m=np.arange(1,M+1)
pk=np.array([(phi+(0 if c[k-1]>0 else 180/k))%(360/k) for k in m])
for k in m:
    peaks=(pk[k-1]+np.arange(k)*360/k)%360
    ax2.scatter(peaks,[k]*k,s=30+4e4*abs(c[k-1])**0.7,color="C3" if c[k-1]>0 else "C0",zorder=3)
ax2.axvline(phi,color="k",lw=2); ax2.axvline((phi+180)%360,color="k",ls="--",lw=2)
ax2.set_xlim(0,360); ax2.set_yticks(m); ax2.set_ylabel("harmonic m"); ax2.set_xlabel("ℓ of each mode's peaks (deg)")
ax2.set_title("Peak angles of each mode\n(red: peak on axis; blue: flipped by 180°/m)",fontsize=15); ax2.grid(alpha=.3); ax2.invert_yaxis()
fig.tight_layout(rect=[0,0,1,0.86]); fig.savefig("cmb-angles.png",dpi=100)
print(l0,phi,fmax,bmax,bmin)
