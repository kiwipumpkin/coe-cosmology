import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("cmb-cancel.py").read().split('if __name__')[0]; g={}; exec(src,g)
out=json.load(open("cmb-cancel.json")); r=out["run2"]; P=[r[k] for k in ("beta1","q","mu","s","phi_deg")]
ELL=g["ELL"]; G=g["G_FULL"]; full=g["full"]; S,Dm,Dg=g["S_model"](*P); phi=P[4]; th=np.radians(ELL-phi)
Sreq=(1/G)/np.mean(1/G); R=g["resid"](G,S)
M=12; c=np.array([2*np.mean(R*np.cos(m*th)) for m in range(1,M+1)])
modes=np.array([c[m-1]*np.cos(m*th) for m in range(1,M+1)]); series=modes.sum(0); left=R-R.mean()-series
def cm(f): f=f-f.mean(); return np.array([2*np.mean(f*np.cos(m*th)) for m in range(1,M+1)])
cG, cS, cDm, cDg = cm(np.log(G)), cm(-np.log(S)), cm(-np.log(Dm)), cm(-np.log(Dg))
b=np.loadtxt("migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":15})
def leg(a,n=2): a.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=n,frameon=False,fontsize=14)
fig,ax=plt.subplots(6,1,figsize=(12,31))
a=ax[0]; a.errorbar(b[:,0],b[:,1],yerr=[b[:,2],b[:,3]],fmt="o",color="gray",ms=4,capsize=2,label="Migkas Fig. 6 bins")
a.plot(ELL,full,"k",lw=2.4,label=r"shell geometry: $D$=58.2, $R$=14.3 Gpc, $\ell_0$=15.2°")
a.set_ylabel(r"$A/A_{\rm all}$"); leg(a,1); a.set_title("1. Our Migkas fit (no sine), nothing removed",pad=70)
a=ax[1]; a.plot(ELL,Sreq,"k",lw=3,label=r"required $S_{\rm req}=(1/G)/\langle 1/G\rangle$")
a.plot(ELL,S/np.mean(S),"C1--",lw=2.2,label=r"model $S=D_{\rm motion}D_{\rm grav}$")
a.plot(ELL,Dm/np.mean(Dm),"C0",lw=1.5,label=r"$D_{\rm motion}$ (SR Doppler)")
a.plot(ELL,Dg/np.mean(Dg),"C2",lw=1.5,label=r"$D_{\rm grav}$ (central mass)")
a.set_ylabel("factor (mean 1)"); leg(a,2)
a.set_title(rf"2. Before removal: SR velocity + mass correction: $\beta_1$={P[0]:.3f}, $q$={P[1]:.3f}, $\mu$={P[2]:.4f}, $s$={P[3]:.0f}, $\varphi$={phi:.2f}°",pad=70,fontsize=14)
a=ax[2]; a.plot(ELL,R*1e3,"C3",lw=2.4,label=rf"residual $R=GS/\langle GS\rangle-1$: RMS {np.sqrt(np.mean(R**2)):.2e}")
a.plot(ELL,series*1e3,"k--",lw=1.6,label=r"solved series $\Sigma_{m=1}^{12}\,c_m\cos m(\ell-\varphi)$")
a.set_ylabel(r"$R\times10^{3}$"); leg(a,1); a.set_title("3. After SR + mass correction, no modes removed",pad=70)
a=ax[3]; cols=list(plt.cm.tab10.colors)+["k","gold"]
for m in range(1,M+1):
    if abs(c[m-1])>2e-5: a.plot(ELL,modes[m-1]*1e3,color=cols[m-1],lw=2.0,label=f"m={m}")
a.set_ylabel(r"$c_m\cos m(\ell-\varphi)\times10^{3}$"); leg(a,5); a.set_title(r"4. The modes $m=1$–$12$ to be removed (cosines about $\varphi$=15.16°)",pad=70)
a=ax[4]; a.plot(ELL,left*1e5,"C0",lw=2,label=rf"$R-\Sigma_{{m=1}}^{{12}}$: RMS {np.sqrt(np.mean(left**2)):.2e}")
a.axhline(1,ls="--",color="gray"); a.axhline(-1,ls="--",color="gray",label=r"$\pm10^{-5}$")
a.set_ylabel(r"leftover $\times10^{5}$"); leg(a,2); a.set_title("5. Full-precision curve, $m=1$–$12$ removed",pad=48)
rr=out["run2_rounded"]; Pr=[rr[k] for k in ("beta1","q","mu","s","phi_deg")]
Rr=g["resid"](g["G_RND"],g["S_model"](*Pr)[0]); thr=np.radians(ELL-Pr[4]); Lr=Rr-Rr.mean()
for mm in range(1,M+1):
    Lr=Lr-2*np.mean(Lr*np.cos(mm*thr))*np.cos(mm*thr)-2*np.mean(Lr*np.sin(mm*thr))*np.sin(mm*thr)
a=ax[5]; a.plot(ELL,Lr*1e5,"C3",lw=1.3,label=rf"same, on 4-decimal CSV: RMS {np.sqrt(np.mean(Lr**2)):.2e}, max {np.abs(Lr).max():.2e}")
a.axhline(1,ls="--",color="gray"); a.axhline(-1,ls="--",color="gray",label=r"$\pm10^{-5}$"); a.axhline(0,color="k",lw=.6)
a.set_ylabel(r"leftover $\times10^{5}$"); leg(a,1); a.set_title("6. Rounded 4-decimal curve, $m=1$–$12$ removed (rounding noise)",pad=70)
for a in ax: a.set_xlim(0,360); a.grid(alpha=.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)")
fig.tight_layout(h_pad=2.5); fig.savefig("cmb-cancel-story.png",dpi=110)
# mode-origin bars
fig,a=plt.subplots(figsize=(12,7)); m=np.arange(1,M+1); w=0.27
a.bar(m-w,np.abs(cG),w,color="k",label=r"$|\ln G|_m$ shell geometry")
a.bar(m,np.abs(cS),w,color="C1",label=r"$|\ln(1/S)|_m$ SR Doppler + mass")
a.bar(m+w,np.abs(c),w,color="C3",label=r"$|c_m|$ residual $\approx \ln G_m-\ln(1/S)_m$")
a.set_yscale("log"); a.set_ylim(1e-8,0.2); a.set_xticks(m); a.set_xlabel("harmonic m"); a.set_ylabel("cosine amplitude")
a.grid(alpha=.3,axis="y"); leg(a,1); a.set_title("Where the residual modes come from",pad=95)
fig.tight_layout(); fig.savefig("cmb-cancel-mode-origin.png",dpi=120)
json.dump(dict(params=P,c_m=c.tolist(),lnG_m=cG.tolist(),lnInvS_m=cS.tolist(),lnInvDm_m=cDm.tolist(),lnInvDg_m=cDg.tolist(),
  leftover_rms=float(np.sqrt(np.mean(left**2)))),open("cmb-cancel-story.json","w"),indent=1)
print(np.round(cS[1:]/cS[:-1],3))
