import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c=0.306601; D=58.2; Rv=14.3; T=13.97; Rrec=0.0245; Trec=3000.0; Gc2=4.785e-20
lo,hi=D-Rv,D+Rv; k=np.log(1.03)/np.log(hi/lo)
rho=lambda r: np.maximum(r,1e-9)**-k                 # power law, all radii (cmb_iso3 law)
Zan=np.trapezoid(4*np.pi*np.linspace(lo,hi,4001)**2*rho(np.linspace(lo,hi,4001)),np.linspace(lo,hi,4001))
def Menc1(r):                                          # enclosed mass per unit M (M = mass in annulus)
    x=np.logspace(-6,np.log10(r),2000); return np.trapezoid(4*np.pi*x*x*rho(x),x)/Zan
th=np.radians(np.linspace(0,180,181)); mu=np.cos(th); W=np.sin(th); d=np.degrees(th)
tc=T*np.logspace(-60,-1e-9,60000); xo=np.array([D,0.0]); vobs=D/T
def events(beta):
    F=1-beta*c/vobs; chi=(c*T*np.log(T/tc)) if F==1 else c*T**F*(T**(1-F)-tc**(1-F))/(1-F); rs=Rrec*(T/tc)**F
    out=[]
    for t_ in th:
        s=np.array([-np.cos(t_),np.sin(t_)]); P=xo+chi[:,None]*s; g=np.hypot(*P.T)-rs
        i=np.where(np.sign(g[:-1])!=np.sign(g[1:]))[0]; assert len(i)==1, "expected one crossing per direction"
        i=i[0]; w=g[i]/(g[i]-g[i+1]); t=tc[i]+w*(tc[i+1]-tc[i]); Pe=P[i]+w*(P[i+1]-P[i]); n=Pe/np.hypot(*Pe); kk=-s
        v=Rrec/t; uo=beta; us=beta*min(1,v/vobs); go=1/np.sqrt(1-uo**2); gs=1/np.sqrt(1-us**2)
        dop=go*(1-uo*kk[0])/(gs*(1-us*np.dot(n,kk)))
        out.append(dict(stretch=(T/t)**F,dop=dop,rnow=v*T,delay=T-t))
    return out
def leg(y):
    return [float((2*l+1)/2*np.trapezoid(y*np.polynomial.legendre.legval(mu,[0]*l+[1])*W,th)) for l in range(5)]
def norm_aniso(Tv):
    m=np.trapezoid(Tv*W,th)/2; return Tv/m-1, m
Ms=np.r_[0,np.logspace(10,24,57)]
res={}
for beta in (0,0.1,0.5):
    E=events(beta); base=np.array([Trec*e["dop"]/e["stretch"] for e in E])
    me=np.array([Menc1(e["rnow"]) for e in E]); mo=Menc1(D)
    phi_e1=-Gc2*me/Rrec; phi_o1=-Gc2*mo/D               # per M☉ of M; Φ_e = λΦ_now(r_now) = −G M_enc(r_now)/R_rec
    a_base=leg(norm_aniso(base)[0])
    rows=[]
    for M in Ms:
        gw=1+M*(phi_e1-phi_o1); weak_ok=bool(M*max(np.abs(phi_e1).max(),abs(phi_o1))<0.1)
        ce=2*M*Gc2*me/Rrec; co=2*M*Gc2*mo/D; ex_def=bool((ce<1).all() and co<1)
        gx=np.sqrt((1-ce)/(1-co)) if ex_def else np.full_like(base,np.nan)
        r=dict(M=float(M),weak_ok=weak_ok,max_phi=float(M*np.abs(phi_e1).max()),exact_defined=ex_def,max_compactness_emit=float(ce.max()))
        for nm,g in (("weak",gw),("exact",gx)):
            if np.all(np.isfinite(g)) and np.all(g>0):
                comb=base*g; ac,mc=norm_aniso(comb); A=leg(ac); Ag=leg(norm_aniso(g)[0])
                r[nm]=dict(defined=True,mean=float(mc),a=A[1:],grav_a=Ag[1:],cancel=1-A[1]/a_base[1],rms=float(np.sqrt(np.trapezoid(ac**2*W,th)/2)),
                           rms_after_l1=float(np.sqrt(np.trapezoid((ac-A[1]*mu)**2*W,th)/2)))
            else: r[nm]=dict(defined=False)
        rows.append(r)
    res[beta]=dict(base_a=a_base[1:],base_mean=float(norm_aniso(base)[1]),rows=rows,base=base,phi_e1=phi_e1,phi_o1=phi_o1,me=me,mo=mo)
    cw=[(r["M"],r["weak"]["cancel"]) for r in rows if r["weak"].get("defined")]
    print(f"beta={beta}: base mean={res[beta]['base_mean']:.4g}K a1..a4={np.round(a_base[1:],4)}")
    for r in rows[::4]:
        w=r["weak"]; x=r["exact"]
        print(f"  M={r['M']:.2e} weak_ok={r['weak_ok']} |Φ|max={r['max_phi']:.2e} "+
              (f"WEAK a1..a4={np.round(w['a'],4)} grav-only a1={w['grav_a'][0]:+.2e} cancel={w['cancel']:+.4f}" if w.get("defined") else "WEAK undefined")+" | "+
              (f"EXACT a1={x['a'][0]:+.4f} cancel={x['cancel']:+.4f}" if x.get("defined") else f"EXACT undefined (max 2GM/rc² at emitter={r['max_compactness_emit']:.2e})"))
    cancels=[r["M"] for r in rows for nm in ("weak","exact") if r[nm].get("defined") and r[nm]["cancel"]>=1]
    print("  exact dipole cancellation mass:", cancels[0] if cancels else "NONE in range (cancel fraction never reaches 1)")
    print("  sign check: grav-only a1 at smallest nonzero M:", rows[1]["weak"]["grav_a"][0] if rows[1]["weak"].get("defined") else None)
# plots
fig,ax=plt.subplots(1,1,figsize=(9,5.5))
for beta,col in zip(res,("k","tab:blue","tab:red")):
    rr=res[beta]["rows"]
    for nm,ls in (("weak","-"),("exact","--")):
        x=[r["M"] for r in rr if r[nm].get("defined") and r["M"]>0]; y=[r[nm]["cancel"] for r in rr if r[nm].get("defined") and r["M"]>0]
        ax.plot(x,y,ls=ls,color=col,marker=".",label=f"β={beta}, {nm}")
ax.axhline(1,color="g",lw=0.8); ax.axhline(0,color="gray",lw=0.8); ax.set_xscale("log"); ax.set_xlabel("M in annulus (M☉)"); ax.set_ylabel("dipole cancelled = 1 − a1(comb)/a1(no grav)")
ax.legend(fontsize=8); ax.grid(alpha=.3); ax.set_title("Dipole cancellation fraction vs mass (power-law density)")
plt.tight_layout(); plt.savefig("cancel_vs_mass.png",dpi=130); plt.close()
key=[1e12,1e14,1e15,1e16]
fig,ax=plt.subplots(3,len(key),figsize=(22,12))
for j,M in enumerate(key):
    for beta,col in zip(res,("k","tab:blue","tab:red")):
        R=res[beta]; gw=1+M*(R["phi_e1"]-R["phi_o1"]); ok=np.all(gw>0)
        ax[0,j].plot(d,R["base"],color=col,label=f"β={beta}")
        ax[1,j].plot(d,gw if ok else np.where(gw>0,gw,np.nan),color=col)
        ax[2,j].plot(d,R["base"]*np.where(gw>0,gw,np.nan),color=col)
    ax[0,j].set_title(f"M={M:.0e} M☉\n(a) stretch+Doppler only"); ax[1,j].set_title("(b) gravity factor alone (weak field)"); ax[2,j].set_title("(c) combined")
    for a_ in ax[:,j]: a_.grid(alpha=.3)
    ax[2,j].set_xlabel("angle from CoE direction (deg)")
ax[0,0].legend(); ax[0,0].set_ylabel("T (K)"); ax[1,0].set_ylabel("ν_o/ν_e gravity"); ax[2,0].set_ylabel("T (K)")
plt.tight_layout(); plt.savefig("overlay_key_masses.png",dpi=110); plt.close()
json.dump({str(b):dict(base_a=v["base_a"],base_mean=v["base_mean"],rows=v["rows"]) for b,v in res.items()},open("cmb_grav_table.json","w"),indent=1)
