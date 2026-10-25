import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=json.load(open("cmb-cancel-story.json")); c=np.array(d["c_m"]); phi=d["params"][4]
cG=np.array(d["lnG_m"]); cS=np.array(d["lnInvS_m"])
src=open("cmb-cancel.py").read().split('if __name__')[0]; g={}; exec(src,g)
S=g["S_model"](*d["params"])[0]; R=g["resid"](g["G_FULL"],S)
L=np.radians(np.arange(0,360.5,0.5)); th=L-np.radians(phi); M=12
A=np.abs(c); pk=[(phi+(0 if c[m-1]>0 else 180/m))%(360/m) for m in range(1,M+1)]
var=c**2/2; vf=var/np.sum(R**2)*len(R)/len(R)  # fraction of residual variance
vf=var/np.mean(R**2)
cols=list(plt.cm.tab10.colors)+["k","goldenrod"]
plt.rcParams.update({"font.size":14})
fig=plt.figure(figsize=(22,30))
fig.suptitle("Fourier modes we must solve for: residual of the SR + central-mass correction",fontsize=26,y=0.995)
fig.text(0.5,0.978,rf"$R(\ell)=\Sigma_m\,c_m\cos m(\ell-\varphi)$,  $\varphi$ = {phi:.2f}° (shell axis), run $s\leq3$: $\beta_1$=0.587, $q$=0.230, $\mu$=0.0933",ha="center",fontsize=18)
def polar(ax,f,col,title,k=0.6):
    r0=1.0; rad=r0+f/np.max(np.abs(f))*k
    ax.plot(L,rad,color=col,lw=2.6); ax.plot(L,np.full_like(L,r0),color="gray",lw=0.8,ls=":")
    ax.plot([np.radians(phi)]*2,[0,1.7],color="k",lw=1,ls="--")
    ax.set_ylim(0,1.7); ax.set_yticklabels([]); ax.set_thetagrids([0,90,180,270]); ax.set_title(title,fontsize=16,pad=14)
for m in range(1,M+1):
    ax=fig.add_subplot(6,4,m,projection="polar"); f=c[m-1]*np.cos(m*th)
    polar(ax,f,cols[m-1],f"m = {m}   period {360/m:.1f}°\n$|c_{{{m}}}|$ = {A[m-1]:.2e}, peak {pk[m-1]:.2f}°")
ax=fig.add_subplot(6,4,(13,18),projection="polar"); f=sum(c[m-1]*np.cos(m*th) for m in range(1,M+1))
polar(ax,f,"k",rf"Sum of m = 1–12 (radius = 1 + scaled $R$)"+f"\nRMS {np.sqrt(np.mean(R**2)):.3e}, max {np.abs(R).max():.2e}",k=0.3)
ax=fig.add_subplot(6,4,(15,20),projection="polar")
for m in range(1,7): ax.plot(L,1+c[m-1]*np.cos(m*th)/A.max()*0.6,color=cols[m-1],lw=2,label=f"m={m}")
ax.set_ylim(0,1.7); ax.set_yticklabels([]); ax.set_thetagrids([0,90,180,270]); ax.set_title("m = 1–6 overlaid, same scale",fontsize=16,pad=14)
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.08),ncol=6,frameon=False,fontsize=13)
ax=fig.add_subplot(6,1,6); ax.axis("off")
rows=[[m,f"{360/m:.1f}",f"{c[m-1]:+.2e}",f"{pk[m-1]:.2f}",f"{100*vf[m-1]:.2f}",f"{cG[m-1]:+.2e}",f"{-cS[m-1]:+.2e}" if False else f"{cS[m-1]:+.2e}",
       "fit matches (cancels)" if m<=2 else ("shell vs Doppler mismatch" if m<=6 else "Doppler tail only")] for m in range(1,M+1)]
t=ax.table(cellText=rows,colLabels=["m","period (°)",r"$c_m$","peak (°)","% of var.",r"$(\ln G)_m$ shell",r"$(\ln 1/S)_m$ SR+mass","origin"],loc="center",cellLoc="center")
t.auto_set_font_size(False); t.set_fontsize(14); t.scale(1,2.0)
for (i,j),cell in t.get_celld().items():
    cell.set_width([0.05,0.09,0.11,0.08,0.08,0.14,0.15,0.22][j])
fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig("cmb-modes-sheet.png",dpi=80)
print(np.round(100*vf,2), sum(vf[2:5]))
