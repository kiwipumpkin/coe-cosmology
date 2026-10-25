"""CMB sky for Chris's CoE shell model, exactly-c geometry, t=3. All assumptions in ASSUMPTIONS.md."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import quad
c=1.0; t_now=3.0; v_obs=1.0; t0=3.0            # D = v_obs (t+t0) = 3+t  -> t0 = 3
D=v_obs*(t_now+t0); R=c*t_now; Rrec=3.0; T_rec=3000.0
assert abs(D-6)<1e-12 and abs(D-R-Rrec)<1e-12
th=np.radians(np.linspace(0,180,361))           # angle from the CoE direction (0 = looking at CoE)
# --- ray trace: emission point on visible circle, distance R from OBS, emitted at t_e = t_now - R/c
t_e=t_now-R/c
P=np.stack([D-R*np.cos(th),R*np.sin(th)],1)     # emission point at t_now coordinates? no: positions AT t_e (static-space light path)
r_e=np.hypot(*P.T)                               # radius of emission point at t_e
assert np.allclose(np.hypot(P[:,0]-D,P[:,1]),R)  # (observer position used at reception; see A5)
# observer at t_e is at D_e = v(t_e+t0)=3; light path in static space from P(t_e) to O(t_now) at D=6 -- A5
O=np.array([D,0.0]); k=(O-P)/np.linalg.norm(O-P,axis=1)[:,None]; L=np.linalg.norm(O-P,axis=1)
# matter at r_e at t_e coasts: v = r_e/(t_e+t0). It crossed R_rec at t_rec: r=v(t+t0)=Rrec
v_src=r_e/(t_e+t0); t_rec=Rrec/v_src-t0
S_full=(t_now+t0)/(t_rec+t0)                    # cube root of shell volume growth t_rec -> now (homologous coasting)
G=6.674e-11; cms=2.998e8; Msun=1.989e30; ell=3.086e25   # length unit = 1 Gpc (A9)
Mell=cms**2*ell/G/Msun                            # mass whose GM/(c^2 ell)=1
def laws():
    k3=np.log(1.03)/np.log(3)
    yield "power_law",  "ρ ∝ r^-0.0269 on 3–9", lambda r:np.where((r>=3)&(r<=9),r**-k3,0)
    BA1=0.03/(3**-1-1.03*9**-1); yield "A+B_r","A + B/r on 3–9",lambda r:np.where((r>=3)&(r<=9),1+BA1/r,0)
    BA3=0.03/(3**-3-1.03*9**-3); yield "A+B_r3","A + B/r³ on 3–9",lambda r:np.where((r>=3)&(r<=9),1+BA3/r**3,0)
    bA=-0.03/(3*2.03); yield "quadratic","A[1+β(r−6)] (γ=0) on 3–9",lambda r:np.where((r>=3)&(r<=9),1+bA*(r-6),0)
    yield "step","flat on 3–9, 5× flat inside R_rec (0–3)",lambda r:np.where(r<3,5.0,np.where(r<=9,1.0,0))
def potential(rho,M):                            # Newtonian Φ/c² at t_now positions, total mass M (code units GM/c²ℓ)
    if M==0: return lambda r: 0*r
    Z=quad(lambda x:4*np.pi*x*x*rho(x),0,9,points=[3],limit=200)[0]
    grid=np.linspace(1e-3,40,8001); dens=rho(grid)/Z
    Menc=np.cumsum(4*np.pi*grid**2*dens)*(grid[1]-grid[0])
    outer=np.cumsum((4*np.pi*grid*dens)[::-1])[::-1]*(grid[1]-grid[0])
    phi=-M*(Menc/grid+outer)
    return lambda r: np.interp(r,grid,phi)
fs=[0,0.25,0.5,0.75,1.0]; Ms=[0,1e18,1e20,1e21,1e22]   # solar masses
def stats(T):
    ok=np.isfinite(T)&(T>0)
    note="" if ok.all() else f"undefined in {int((~ok).sum())}/{len(T)} directions (excluded from stats)"
    w=np.sin(th)[ok]; x=np.cos(th)[ok]; T=T[ok]; thv=th[ok]
    mean=np.trapezoid(T*w,thv)/np.trapezoid(w,thv); a1=np.trapezoid(T*x*w,thv)/np.trapezoid(x*x*w,thv)
    res=T-mean-a1*x; rms=np.sqrt(np.trapezoid(res**2*w,thv)/np.trapezoid(w,thv))
    return dict(mean=mean,dipole=a1/mean,resid=rms/mean,note=note)
rows=[]
for key,lab,rho in laws():
    fig,ax=plt.subplots(2,len(fs),figsize=(26,10),sharex=True)
    for j,f in enumerate(fs):
        u_s=(1-f)*v_src[:,None]*(P/r_e[:,None]); u_o=(1-f)*v_obs*np.array([1.0,0])
        num=1-(k@u_o); den=1-np.sum(u_s*k,1)
        dop=np.where(den>0,num/den,np.nan)       # classical Doppler ν_o/ν_e; undefined when source outruns its own light toward us
        dop=np.where(num>0,dop,np.nan)          # observer outruns the light: never received
        stretch=S_full**f
        for M,col in zip(Ms,plt.cm.viridis(np.linspace(0,0.9,len(Ms)))):
            phi=potential(rho,M/Mell)
            phi_e=phi(r_e*(t_now+t0)/(t_e+t0))*(t_now+t0)/(t_e+t0)   # homologous: Φ_te(r)=λΦ_now(λr), λ=a_now/a_te
            phi_o=phi(np.array([D]))[0]
            grav=1+(phi_e-phi_o)                 # weak-field ν_o/ν_e
            weak_ok=bool(max(np.max(np.abs(phi_e)),abs(phi_o))<0.1)
            T=T_rec*dop*grav/stretch
            st=stats(T); rows.append(dict(law=key,f=f,M=M,weak_field_ok=weak_ok,max_abs_phi=float(max(np.max(np.abs(phi_e)),abs(phi_o))),
                 stretch_min=float(stretch.min()),stretch_max=float(stretch.max()),
                 dop_min=float(np.nanmin(dop)) if np.isfinite(dop).any() else None,dop_max=float(np.nanmax(dop)) if np.isfinite(dop).any() else None,**st))
            ls="-" if weak_ok else ":"
            ax[0,j].plot(np.degrees(th),T,color=col,ls=ls,label=f"M={M:.0e} M☉"+("" if weak_ok else " (weak field FAILS)"))
            if st["mean"]: ax[1,j].plot(np.degrees(th),T/st["mean"]-1,color=col,ls=ls)
            else: ax[1,j].plot(np.degrees(th),T/np.nanmean(T)-1,color=col,ls=ls)
        ax[0,j].set_title(f"f = {f} (expansion share)\nclassical Doppler"+("\nSR undefined: speeds ≥ c" if (1-f)*v_src.max()>=1 else ""),fontsize=10)
        ax[0,j].set_yscale("log"); ax[1,j].set_xlabel("angle from CoE direction (deg)")
    ax[0,0].set_ylabel("T observed (K), log"); ax[1,0].set_ylabel("ΔT/T (full curve, dipole NOT removed)")
    ax[0,-1].legend(fontsize=8,loc="center left",bbox_to_anchor=(1.02,0.5))
    fig.suptitle(f"CoE shell model, exactly-c geometry (t=3, D=6, R=3, R_rec=3). Density: {lab}. Gaps = undefined (Doppler singular).",fontsize=13)
    plt.tight_layout(rect=[0,0,0.92,0.94]); plt.savefig(f"cmb_{key}.png",dpi=110); plt.close()
# pieces (law-independent)
fig,ax=plt.subplots(1,4,figsize=(24,5.5)); d=np.degrees(th)
ax[0].plot(d,r_e); ax[0].set_title("emission radius r_e at t_e=0 (c·Δt)")
ax[1].plot(d,t_rec,label="t_rec (shell crossed R_rec)"); ax[1].axhline(t_e,color="k",label="t_e (light leaves)"); ax[1].axhline(t_now,color="r",label="t_now (received)"); ax[1].legend(); ax[1].set_title("times: recombination, emission, reception")
ax[2].plot(d,S_full); ax[2].set_title("full stretch S = a(now)/a(t_rec) (raised to f in T)")
for f in fs:
    u_s=(1-f)*v_src[:,None]*(P/r_e[:,None]); num=1-(k@((1-f)*np.array([1.0,0]))); den=1-np.sum(u_s*k,1)
    ax[3].plot(d,np.where((den>0)&(num>0),num/den,np.nan),label=f"f={f}")
ax[3].set_yscale("log"); ax[3].legend(); ax[3].set_title("classical Doppler factor ν_o/ν_e")
for a in ax: a.set_xlabel("angle from CoE direction (deg)"); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig("cmb_pieces.png",dpi=110)
json.dump(rows,open("cmb_table.json","w"),indent=1)
print("Mell(Msun)=%.3e"%Mell, "r_e range",r_e.min(),r_e.max(),"t_rec range",t_rec.min(),t_rec.max(),"S range",S_full.min(),S_full.max(),"v_src",v_src.min(),v_src.max())
for r in rows:
    print(f"{r['law']:10s} f={r['f']:<4} M={r['M']:.0e} weak={str(r['weak_field_ok']):5s} |Φ|max={r['max_abs_phi']:.2e} "+
          f"mean={r['mean']:.4g}K dip={r['dipole']:+.4f} resid={r['resid']:.4f} {r['note']}")
