"""Brief (Oct 11): what the leftover looks like once the large scales are cancelled. Saved results only, no reruns.
Residual = the saved final ring residual 'resid' (360 points, 1 deg) of each run. Modes: m = 1..179 about the CoE axis (l = 15.16 deg).
Randomness checks: symmetric share over all modes; 2000 random-phase rings with the same mode amplitudes (Gaussian-noise rings with the
same spectrum and rms) give the expected range; KS test of the values against a Gaussian of the same rms; autocorrelation first zero."""
import json,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({"font.size":12})
SO="/workspace/clump_scan/sequential_open"
EX=[("Compensated random web, seed 3 (entropy rule)",f"{SO}/minedit_comp/comp_s3_entropy.json","#1a3a8f"),
    ("Foam web, seed 3 (no deletion, 3 Gpc backstop)",f"{SO}/minedit2/b_foam_rim_s3_tidalcore_nodel_reloc.json","#1a3a8f"),
    ("Foam web, seed 1 (no deletion, 3 Gpc backstop)",f"{SO}/minedit2/b_foam_rim_s1_tidalcore_nodel_reloc.json","#1a3a8f"),
    ("Contrast: compensated web, seed 2 (stalled just above)",f"{SO}/minedit_comp/comp_s2_entropy.json","#7a3a8f")]
AX=np.radians(15.16); th=np.radians(np.arange(360)); L=np.arange(360); M=np.arange(1,180); FL=1e-5
def modes(v):
    a=2*np.mean(v[None]*np.cos(M[:,None]*(th-AX)),1); b=2*np.mean(v[None]*np.sin(M[:,None]*(th-AX)),1); return a,b
rng=np.random.default_rng(0); OUT={}
fig,axs=plt.subplots(4,2,figsize=(12,10.5),gridspec_kw=dict(width_ratios=[1.75,1],hspace=0.75,wspace=0.18))
fh,ah=plt.subplots(1,4,figsize=(12,3.3),sharey=True); fh.subplots_adjust(wspace=0.08,bottom=0.3,top=0.82)
for i,(name,f,col) in enumerate(EX):
    d=json.load(open(f)); v=np.array(d["resid"]); v=v-v.mean(); x=v/FL
    a,b=modes(v); A=np.hypot(a,b)/FL; P=A**2; sym=float((a**2).sum()/(a**2+b**2).sum())
    S=[];Kk=[]
    for _ in range(2000):
        ph=rng.uniform(0,2*np.pi,len(M)); w=(A[:,None]*np.cos(M[:,None]*(th-AX)+ph[:,None])).sum(0)
        aa,bb=modes(w); S.append((aa**2).sum()/(aa**2+bb**2).sum()); Kk.append(stats.kurtosis(w))
    c=np.array([np.mean(v*np.roll(v,k)) for k in range(180)]); c/=c[0]
    cum=np.cumsum(P)/P.sum()
    R=dict(rms=float(x.std()),maxabs=float(abs(x).max()),frac_outside_band=float(np.mean(abs(x)>1)),
           share_m1_12=float(P[:12].sum()/P.sum()),rms_m1_12=float(np.sqrt(P[:12].sum()/2)),m10=int(M[np.searchsorted(cum,0.1)]),m90=int(M[np.searchsorted(cum,0.9)]),
           largest_mode=int(M[P.argmax()]),largest_mode_share=float(P.max()/P.sum()),sym_all=sym,sym_m1_12_saved=d["after"]["sym"],
           sym_noise_5_95=[float(np.percentile(S,5)),float(np.percentile(S,95))],sym_percentile=float(100*np.mean(np.array(S)<sym)),
           kurtosis=float(stats.kurtosis(x)),kurt_noise_5_95=[float(np.percentile(Kk,5)),float(np.percentile(Kk,95))],skew=float(stats.skew(x)),
           ks_p=float(stats.kstest(x/x.std(),"norm").pvalue),acl_first_zero_deg=int(np.argmax(c<0)),leftover_saved=float(d["after"]["total"]))
    OUT[name]=R
    ax=axs[i,0]; ax.axhspan(-1,1,color="#bff0d4",alpha=0.8,lw=0,label="±1× CMB band"); ax.plot(L,x,color=col,lw=1.1)
    ax.axhline(0,color="0.5",lw=0.5); ax.set_xlim(0,360); ax.set_xticks([0,90,180,270,360]); ax.set_ylim(-4,4)
    ax.set_title(f"{name}: leftover {R['rms']:.3f}×",fontsize=12,loc="left"); ax.set_ylabel("residual (×10⁻⁵)")
    if i==3: ax.set_xlabel("position around the ring, ℓ (°)")
    bx=axs[i,1]; mm=M[:120]
    bx.bar(mm[:12],A[:12],width=0.9,color="#d1495b",label="m = 1–12 (large scales)"); bx.bar(mm[12:],A[12:120],width=0.9,color="#3a6fb0",label="m ≥ 13 (fine texture)")
    bx.set_xlim(0,121); bx.set_ylim(0,max(A[:120])*1.25); bx.set_ylabel("amplitude (×)")
    bx.set_title(f"m = 1–12 carry {100*R['share_m1_12']:.2g}% of the power",fontsize=11,loc="left")
    if i==3: bx.set_xlabel("mode m (waves around the ring)")
    hx=ah[i]; hx.hist(x/x.std(),bins=np.linspace(-3.5,3.5,29),density=True,color="#9db7d9",edgecolor="white")
    g=np.linspace(-3.5,3.5,200); hx.plot(g,stats.norm.pdf(g),color="k",lw=1.4,label="Gaussian noise, same rms")
    hx.set_title(name.split(" (")[0].replace("Contrast: c","C").replace("Compensated random web","Compensated web"),fontsize=10.5)
    hx.text(0.03,0.97,f"symmetric share {100*sym:.0f}%\n(noise: {100*R['sym_noise_5_95'][0]:.0f}–{100*R['sym_noise_5_95'][1]:.0f}%)\nKS p = {R['ks_p']:.2f}\ncorr. length {R['acl_first_zero_deg']}°",transform=hx.transAxes,va="top",fontsize=8.5,bbox=dict(facecolor="white",alpha=0.85,edgecolor="none",pad=1.5))
    hx.set_xlabel("residual ÷ its rms")
ah[0].set_ylabel("density of values")
h,l=axs[0,0].get_legend_handles_labels(); h2,l2=axs[0,1].get_legend_handles_labels()
fig.legend(h+h2,l+l2,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(0.5,0.0)); fig.subplots_adjust(bottom=0.09,top=0.96,left=0.07,right=0.98)
fig.savefig("figures/new/noise_examples_rows.png",dpi=150); plt.close(fig)
fh.legend(*ah[0].get_legend_handles_labels(),loc="lower center",frameon=False,ncol=1,bbox_to_anchor=(0.5,0.0)); [q.set_ylim(0,0.62) for q in ah]; fh.savefig("figures/new/noise_examples_hist.png",dpi=150,bbox_inches="tight")
json.dump(OUT,open("figures/new/noise_examples_stats.json","w"),indent=1)
for k,R in OUT.items(): print(k,{q:(round(w,3) if isinstance(w,float) else w) for q,w in R.items()})
