"""Isotropic-emission CMB for Chris's CoE model at paper scale. See ASSUMPTIONS.md."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import quad
c=0.306601            # Gpc/Gyr
D=58.2; Rv=14.3; T=13.97; Rrec=0.0245; Trec=3000.0   # paper lines 285-287, 164, 1121; R_rec from chat (not in paper)
Gc2=4.785e-20         # G/c^2 in Gpc per Msun
lo,hi=D-Rv,D+Rv
k=np.log(1.03)/np.log(hi/lo)                  # 3% lock between annulus edges, in r
BA3=0.03/(lo**-3-1.03*hi**-3); bq=-0.03/(lo*1.03-hi+ (hi-lo)*0) if False else None
# quadratic (gamma=0 -> linear) in r about D with rho(lo)/rho(hi)=1.03:  (1+b(lo-D))/(1+b(hi-D))=1.03
b_lin=0.03/((lo-D)-1.03*(hi-D))
LAWS={"power_law":lambda r:np.maximum(r,1e-9)**-k,
      "A+B_r3":lambda r:1+BA3/np.maximum(r,1e-9)**3,
      "quadratic_clipped":lambda r:np.maximum(1+b_lin*(r-D),0.0),
      "step":lambda r:np.where(r<lo,5.0,1.0),
      "uniform":lambda r:np.ones_like(np.asarray(r,float))}
for nm,fn in LAWS.items():
    assert abs(fn(np.array([lo]))[0]/fn(np.array([hi]))[0]-(1.0 if nm in("step","uniform") else 1.03))<1e-9, nm
def potential(rho,M):
    if M==0: return lambda r:0*r
    Zan=quad(lambda x:4*np.pi*x*x*rho(np.array([x]))[0],lo,hi)[0]
    def phi(r):
        r=np.atleast_1d(r).astype(float); out=np.empty_like(r)
        for i,x in enumerate(r):
            pts=np.logspace(-6,np.log10(max(x,1e-5)),400); Me=np.trapezoid(4*np.pi*pts**2*rho(pts),pts)
            out[i]=-Gc2*M*Me/Zan/x
        return out
    return phi
fs=[0,0.25,0.5,0.75,1.0]; Ms=[0,1e12,1e15,1e18,1e20]
th=np.radians(np.linspace(0,180,181)); xo=np.array([D,0.0]); u_o_tot=D/T
tc=T*np.logspace(-60,-1e-9,60000)
rows=[]; curves={}
for f in fs:
    chi=(c*T*np.log(T/tc)) if f==1 else c*T**f*(T**(1-f)-tc**(1-f))/(1-f)   # comoving distance light covers tc->T (a(T)=1)
    rs=Rrec*(T/tc)**f                       # comoving radius of R_rec at time tc
    stretch=(T/tc)**f; v=Rrec/tc; rnow=v*T  # matter crossing at tc: speed, radius now
    ev=[]                                   # per direction list of events
    for t_ in th:
        s=np.array([-np.cos(t_),np.sin(t_)]); P=xo[None,:]+chi[:,None]*s[None,:]
        g=np.hypot(*P.T)-rs; idx=np.where(np.sign(g[:-1])!=np.sign(g[1:]))[0]
        E=[]
        for i in idx:
            w_=g[i]/(g[i]-g[i+1]); j=lambda a: a[i]+w_*(a[i+1]-a[i])
            Pe=P[i]+w_*(P[i+1]-P[i]); n_s=Pe/np.hypot(*Pe); dgdt=abs((g[i+1]-g[i])/(tc[i+1]-tc[i]))
            E.append(dict(tc=j(tc),stretch=j(stretch),v=j(v),rnow=j(rnow),n_s=n_s,k=-s,dgdt=dgdt,delay=T-j(tc)))
        ev.append(E)
    for law,rho in LAWS.items():
        for M in Ms:
            phi=potential(rho,M); phio=float(phi(np.array([D]))[0]); Tdir=np.full(len(th),np.nan)
            nund=0; ndop=0; maxphi=abs(phio)/c**2*0+abs(phio); npieces=[]
            for a,E in enumerate(ev):
                num=0; den=0; pieces=[]
                for e in E:
                    lam=e["rnow"]/Rrec                     # homologous factor now/tc
                    wt=rho(np.array([e["rnow"]]))[0]*lam**3*e["v"]/e["dgdt"]
                    if wt<=0: continue
                    us=(1-f)*e["v"]/c; uo=(1-f)*u_o_tot/c
                    dn=1-uo*np.dot(np.array([1.0,0]),e["k"]); dd=1-us*np.dot(e["n_s"],e["k"])
                    if dn<=0 or dd<=0: ndop+=1; continue
                    phie=lam*float(phi(np.array([e["rnow"]]))[0]); maxphi=max(maxphi,abs(phie))
                    grav=1+(phie-phio)
                    if grav<=0: ndop+=1; continue
                    Tobs=Trec*(dn/dd)*grav/e["stretch"]
                    num+=wt*Tobs**4; den+=wt; pieces.append((e["delay"],e["stretch"],dn/dd,grav))
                if den>0: Tdir[a]=(num/den)**0.25
                if pieces: npieces.append(pieces[0])
            ok=np.isfinite(Tdir); st=dict(mean=None,dipole=None,resid=None)
            if ok.all():
                w=np.sin(th[ok]); x=np.cos(th[ok]); Tv=Tdir[ok]; tv=th[ok]
                mean=np.trapezoid(Tv*w,tv)/np.trapezoid(w,tv); a1=np.trapezoid(Tv*x*w,tv)/np.trapezoid(x*x*w,tv)
                res=Tv-mean-a1*x; st=dict(mean=float(mean),dipole=float(a1/mean),resid=float(np.sqrt(np.trapezoid(res**2*w,tv)/np.trapezoid(w,tv))/mean))
            P=np.array(npieces) if npieces else None
            rows.append(dict(note=('' if ok.all() else ('no directions defined' if not ok.any() else 'partial sky: stats not computed')),f=f,law=law,M=M,directions_defined=int(ok.sum()),of=len(th),events_rejected_doppler_or_grav=ndop,
                weak_field_ok=bool(maxphi<0.1),max_abs_phi=float(maxphi),
                stretch_range=[float(P[:,1].min()),float(P[:,1].max())] if P is not None else None,
                doppler_range=[float(P[:,2].min()),float(P[:,2].max())] if P is not None else None,
                grav_range=[float(P[:,3].min()),float(P[:,3].max())] if P is not None else None,
                delay_Gyr=[float(P[:,0].min()),float(P[:,0].max())] if P is not None else None,**st))
            curves[(f,law,M)]=Tdir
    nhit=sum(1 for E in ev if E)
    if f>0: assert nhit==len(th), f"photon field missing in {len(th)-nhit} directions at f={f}"
    else:
        reach=c*T+Rrec; assert reach<D-Rv*0 and nhit==0
        print(f"f=0 FRAMEWORK FAILURE (not a bug): photons travel at most c*T+R_rec={reach:.3f} Gpc from the CoE; observer at D={D} Gpc -> none can arrive")
    print(f"f={f}: directions with any crossing on ray: {nhit}/{len(th)}; stretch available range {stretch.min():.3g}-{stretch.max():.3g}")
json.dump(rows,open("cmb_iso2_table.json","w"),indent=1)
d=np.degrees(th); cols=plt.cm.viridis(np.linspace(0,0.9,len(Ms)))
for law in LAWS:
    fig,ax=plt.subplots(2,len(fs),figsize=(26,10))
    for j,f in enumerate(fs):
        for M,cl in zip(Ms,cols):
            Tdir=curves[(f,law,M)]; r=[x for x in rows if x["f"]==f and x["law"]==law and x["M"]==M][0]
            ls="-" if r["weak_field_ok"] else ":"
            ax[0,j].plot(d,Tdir,color=cl,ls=ls,label=f"M={M:.0e} M☉"+("" if r["weak_field_ok"] else " (weak field FAILS)"))
            if r["mean"]: ax[1,j].plot(d,Tdir/r["mean"]-1,color=cl,ls=ls)
            elif np.isfinite(Tdir).any(): ax[1,j].plot(d,Tdir/np.nanmean(Tdir)-1,color=cl,ls="--")
        for a_ in ax[:,j]: a_.set_xlim(0,180)
        ax[0,j].set_title(f"f = {f}\nclassical Doppler"+("; SR undefined (>c)" if f<1 else "; no Doppler"),fontsize=10)
        ax[1,j].set_xlabel("angle from CoE direction (deg)")
        if any(np.isfinite(curves[(f,law,M)]).any() for M in Ms): ax[0,j].set_yscale("log")
        if not np.isfinite(curves[(f,law,0)]).any(): ax[0,j].text(0.5,0.5,("f=0: FRAMEWORK FAILURE\nno photon can reach D in time\n(c·T = 4.28 Gpc ≪ D = 58.2 Gpc)" if f==0 else "photon field exists in every direction,\nbut every event undefined\n(classical Doppler or gravity ≤ 0)"),transform=ax[0,j].transAxes,ha="center")
    ax[0,0].set_ylabel("T (K), log"); ax[1,0].set_ylabel("ΔT/T, dipole kept\n(dashed: partial sky, rel. to mean of defined part)")
    ax[0,-1].legend(fontsize=8,loc="center left",bbox_to_anchor=(1.02,0.5))
    fig.suptitle(f"Isotropic emission at R_rec = 24.5 Mpc crossing, D = 58.2 Gpc, age 13.97 Gyr. Density law: {law}. Gaps = undefined.",fontsize=12)
    plt.tight_layout(rect=[0,0,0.92,0.94]); plt.savefig(f"cmb_iso2_{law}.png",dpi=110); plt.close()
for r in rows:
    print(f"f={r['f']:<4} {r['law']:9s} M={r['M']:.0e} dirs={r['directions_defined']}/{r['of']} rej={r['events_rejected_doppler_or_grav']} weak={r['weak_field_ok']} |Φ|={r['max_abs_phi']:.2e} "
          +(f"mean={r['mean']:.4g}K dip={r['dipole']:+.4f} res={r['resid']:.4f}" if r['mean'] else "no stats: "+r["note"])+f" S={r['stretch_range']} Dop={r['doppler_range']} delay={r['delay_Gyr']}")
