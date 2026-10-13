"""All-sky acceleration-sign test for CoE shell (Paper I), Newtonian, simultaneous (no light-cone).
Geometry: CoE at origin, observer at D=58.2 Gpc on +z; sightline n at angle th from the outward (CoE->us) axis;
th=180 deg = CoE direction. Galaxy at r_g = D z + s n; |r_g| always within 43.9..72.5 Gpc for s<=14.3.
Velocity law (Paper I eq. v2, run 2): v^2 = b1^2 r0^2 + 2 mu (1/r - 1/r0), r = s_exp*r0, in units c and d1 = D/s_exp.
Mass profiles: P = central mass only (GM = mu c^2 d1; the mass that sets the paper's v-law)
               PS = central mass + exp-lock shell rho = 2.9e-27 exp(-(r-D)/922.65 Gpc) for 43.9<r<72.5
               U = uniform shell rho = 2.9e-27 in 43.9..72.5, empty interior (contrast)
Expanding: v outward; collapsing: same |v| inward. Gravity identical in both."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
G=6.674e-11; c=2.998e8; GPC=3.0857e25
D=58.2; RIN,ROUT=43.9,72.5; RV=14.3; LEXP=922.65; RHO=2.9e-27
b1,mu,sexp=0.5866,0.09332,3.0; d1=D/sexp
GMc=mu*c**2*d1*GPC
def vmag(r):                       # m/s, r in Gpc
    r0=r/sexp/d1; x=r/d1
    return c*np.sqrt(np.maximum(b1**2*r0**2+2*mu*(1/x-1/r0),0))
rr=np.linspace(RIN,ROUT,4001)
def shellmass(dens):
    w=4*np.pi*rr**2*dens(rr)*GPC**3; M=np.concatenate([[0],np.cumsum(0.5*(w[1:]+w[:-1])*np.diff(rr))]); return M
Mexp=shellmass(lambda r:RHO*np.exp(-(r-D)/LEXP)); Muni=shellmass(lambda r:RHO+0*r)
PROF={"P":lambda r:GMc+0*r,
      "PS":lambda r:GMc+G*np.interp(r,rr,Mexp),
      "U":lambda r:G*np.interp(r,rr,Muni)}
NAMES={"P":"central mass only (paper v-law mass)","PS":"central mass + exp-lock shell","U":"uniform shell, empty interior"}
# grid: th 0..180 (axisymmetric), phi for maps
th=np.radians(np.linspace(0,180,1801)); ct,st=np.cos(th),np.sin(th)
def compute(prof,sgn,s):
    GM=PROF[prof]
    xo=np.array([0,0,D]); n=np.stack([st,0*th,ct],1); xg=xo+s*n; rg=np.linalg.norm(xg,axis=1)
    ag=-(GM(rg)/(rg*GPC)**3)[:,None]*xg*GPC; ao=-(GM(np.array(D))/(D*GPC)**3)*xo*GPC
    arel=ag-ao; alos=(arel*n).sum(1)
    vg=sgn*vmag(rg)[:,None]*xg/rg[:,None]; vo=sgn*vmag(D)*xo/D
    vrel=vg-vo; vlos=(vrel*n).sum(1); vt2=(vrel**2).sum(1)-vlos**2
    sdd=alos+vt2/(s*GPC)
    q=-s*GPC*sdd/np.where(vlos==0,np.nan,vlos**2)
    return dict(alos=alos,vlos=vlos,sdd=sdd,q=q,trans=vt2/(s*GPC))
def skyfrac(mask): return float(np.trapezoid(mask*st,th)/2)
def flips(y):
    i=np.where(np.diff(np.sign(y))!=0)[0]; return [round(float(180-np.degrees(th[k]+0.5*(th[k+1]-th[k]))),2) for k in i]  # angle from CoE dir
# Mollweide with CoE at map centre: angle from CoE = 180-th
lon=np.radians(np.linspace(-180,180,721)); lat=np.radians(np.linspace(-90,90,361)); LO,LA=np.meshgrid(lon,lat)
ang_coe=np.degrees(np.arccos(np.clip(np.cos(LA)*np.cos(LO),-1,1)))
def skymap(y,title,unit,fn,signonly=False):
    Y=np.interp(ang_coe,180-np.degrees(th[::-1]),y[::-1])
    fig=plt.figure(figsize=(9,5.2)); ax=fig.add_subplot(111,projection="mollweide")
    if signonly: Y=np.sign(Y); norm=TwoSlopeNorm(0,-1,1)
    else:
        m=np.nanmax(np.abs(Y)); norm=TwoSlopeNorm(0,-m,m)
    pc=ax.pcolormesh(LO,LA,Y,cmap="RdBu_r",norm=norm,shading="auto")
    ax.contour(LO,LA,Y,levels=[0],colors="k",linewidths=1)
    ax.plot(0,0,"k*",ms=14); ax.text(0.06,0.06,"CoE",fontsize=11)
    ax.set_xticklabels([]); ax.set_yticklabels([]); ax.grid(alpha=0.3)
    ax.set_title(title,fontsize=12,pad=12)
    cb=fig.colorbar(pc,ax=ax,orientation="horizontal",fraction=0.05,pad=0.06); cb.set_label(unit,fontsize=10)
    fig.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(fig)
res={}; S=[1,3,7,14]
import os; os.makedirs("maps",exist_ok=True)
for p in PROF:
  for case,sg in (("expanding",1),("collapsing",-1)):
    for s in S:
      d=compute(p,sg,s); k=f"{p}|{case}|s={s}"
      acc=d["sdd"]>0; rec=d["vlos"]>0
      res[k]=dict(frac_alos_pos=skyfrac(d["alos"]>0),frac_sdd_pos=skyfrac(acc),frac_receding=skyfrac(rec),
        frac_receding_and_accel=skyfrac(acc&rec),frac_speed_increasing=skyfrac(np.sign(d["sdd"])==np.sign(d["vlos"])),
        flip_alos_deg_from_CoE=flips(d["alos"]),flip_sdd_deg_from_CoE=flips(d["sdd"]),flip_vlos_deg_from_CoE=flips(d["vlos"]),
        alos_range=[float(d["alos"].min()),float(d["alos"].max())],sdd_range=[float(d["sdd"].min()),float(d["sdd"].max())],
        vlos_kms_range=[float(d["vlos"].min()/1e3),float(d["vlos"].max()/1e3)],trans_max=float(d["trans"].max()))
      lab=f"{NAMES[p]}, {case}, s = {s} Gpc"
      skymap(d["sdd"],f"d²s/dt² (gravity + transverse term)\n{lab}","m s⁻²",f"maps/sdd_{p}_{case}_s{s}.png")
      skymap(d["vlos"]/1e3,f"Line-of-sight relative velocity\n{lab}","km s⁻¹",f"maps/vlos_{p}_{case}_s{s}.png") if p=="P" else None
      if case=="expanding": skymap(d["alos"],f"Gravitational a_los (same for expanding/collapsing)\n{NAMES[p]}, s = {s} Gpc","m s⁻²",f"maps/alos_{p}_s{s}.png")
      skymap(np.where(np.sign(d["sdd"])==np.sign(d["vlos"]),1.,-1.),f"Recession speed |ds/dt| growing (red) / shrinking (blue)\n{lab}","sign",f"maps/speedup_{p}_{case}_s{s}.png",signonly=True)
json.dump(res,open("results.json","w"),indent=1)
# formula check vs reviewer closed form
d=compute("U",1,7); xo=D; s=7; rg=np.sqrt(D**2+s**2+2*D*s*ct); f=lambda r:PROF["U"](r)/(r*GPC)**3
alt=-f(rg)*(s+D*ct)*GPC+f(np.array(D))*D*ct*GPC; print("formula max rel diff",np.max(np.abs(alt-d["alos"]))/np.max(np.abs(d["alos"])))
print("v_obs/c",vmag(D)/c,"vmin/vmax/c",vmag(RIN)/c,vmag(ROUT)/c,"Mc",GMc/G,"Mshell exp",Mexp[-1],"uni",Muni[-1])
for k,v in res.items(): print(k,{a:(round(b,3) if isinstance(b,float) else b) for a,b in v.items() if a.startswith(("frac","flip"))})
