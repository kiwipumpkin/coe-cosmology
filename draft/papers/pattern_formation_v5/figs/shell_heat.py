"""Heat map of mass density in the wrinkled shell (shell_fig.py geometry: 45 deg wedge, R=3.2 Gpc, lambda=3*Delta, ell~131 for Delta/R=0.016).
Shell thickness drawn x4 for visibility. Density in shell relative to unwrinkled mean: 1+d*cos(ell*theta); ripple dr=amp*Delta*cos(ell*theta). Schematic."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":16})
R=3.2; fr=0.016; Dl=fr*R*1.5; ell=2*np.pi*R/(3*fr*R); SPAN=45
th=np.radians(np.linspace(0,SPAN,2500)); rr=np.linspace(R-3*Dl,R+3*Dl,500); T,Rg=np.meshgrid(th,rr)
for n,(amp,d,title) in enumerate([(0,0,"1. Flat shell: uniform density"),(0.25,0.15,"2. Wrinkled shell: slight gathering"),(0.45,0.8,"3. Concentrated: crests hot, troughs cold")],1):
    c=np.cos(ell*T); ctr=R+amp*Dl*c; w=Dl*(1+0.5*d*c)      # crests thicken and fill
    rho=np.where(np.abs(Rg-ctr)<w/2,(1+d*c),np.nan)
    fig=plt.figure(figsize=(8,9)); fig.text(0.5,0.965,title,ha="center",fontsize=18,weight="bold"); fig.text(0.5,0.93,"Δ/R = 0.016, λ = 3Δ (ℓ ≈ %d); thickness drawn ×1.5; schematic"%ell,ha="center",fontsize=13)
    ax=fig.add_axes([0.03,0.2,0.94,0.7]); ax.set_aspect("equal"); ax.axis("off")
    m=ax.pcolormesh(Rg*np.cos(T),Rg*np.sin(T),rho,cmap="inferno",vmin=0,vmax=2,shading="auto")
    Rm=R*1.15; ax.plot([0,Rm],[0,0],color="0.6",lw=1); ax.plot([0,Rm*np.cos(th[-1])],[0,Rm*np.sin(th[-1])],color="0.6",lw=1)
    ax.plot(0,0,"k*",ms=16); ax.text(0.08,-0.15,"CoE",fontsize=14); ax.set_xlim(-0.05,Rm); ax.set_ylim(-0.25,Rm*0.72)
    cb=fig.add_axes([0.12,0.1,0.76,0.035]); b=fig.colorbar(m,cax=cb,orientation="horizontal"); b.set_label("density / unwrinkled mean",fontsize=16); b.ax.tick_params(labelsize=14)
    fig.savefig(f"shell_heat_{n}.png",dpi=150); plt.close(fig)
