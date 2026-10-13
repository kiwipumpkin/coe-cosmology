import json, numpy as np
F=1e-5; o=json.load(open("ensemble.json")); base=o["base_rms_over_floor"]; A0=np.array(o["base_amps_over_floor"])*F
q35_0=np.linalg.norm(A0[2:5]); q12_0=np.linalg.norm(A0[0:2])
LBL={"a":"(a) 1e48–1e49","b":"(b) 3e49–3e50","c":"(c) 1e51–1e52"}
def cnt(d,ka,kr):   # identical to ../counts.py
    r=d[kr]/F; A=d[ka]; q35=np.linalg.norm(A[:,2:5],axis=1); q12=np.linalg.norm(A[:,0:2],axis=1)
    lo=int((r<base).sum()); lo35=int((q35<q35_0).sum()); hi12=int((q12>q12_0).sum()); both=int(((q35<q35_0)&(q12<=q12_0)).sum())
    strict=int((np.all(A[:,2:5]<A0[2:5],1)&np.all(A[:,0:2]<=A0[0:2],1)).sum())
    return dict(n=len(r),rms_lower=lo,rms_higher=len(r)-lo,m35_lower=lo35,m12_higher=hi12,m35_lower_and_m12_not_higher=both,
                each_m3_m4_m5_lower_and_each_m1_m2_not_higher=strict,
                median=float(np.median(r)),p5=float(np.percentile(r,5)),min=float(r.min()),
                median_q35=float(np.median(q35)/F),median_q12=float(np.median(q12)/F))
C={}
for gtag in ("grad","unif"):
    for tag in ("bare","buf"):
        for N in (200,500):
            for lev in "abc":
                k=f"{gtag}_{tag}_N{N}_{lev}"; d=np.load(f"raw_{k}.npz")
                C[k]=dict(placed=cnt(d,"amps","rms"),best=cnt(d,"amps_best","rms_best"),
                          frac_far_hemisphere=float(d["frac_back"].mean()),frac_within_28p6_gpc=float(d["frac_near"].mean()),frac_beyond=float(d["frac_beyond"].mean()))
O={}   # original numbers recomputed from the original raw files with the same function
for tag in ("bare","buf"):
    for N in (200,500):
        for lev in "abc":
            k=f"{tag}_N{N}_{lev}"; d=np.load(f"../raw_{k}.npz"); O[k]=dict(placed=cnt(d,"amps","rms"),best=cnt(d,"amps_best","rms_best"))
json.dump(dict(base_rms_over_floor=base,base_q35_over_floor=q35_0/F,base_q12_over_floor=q12_0/F,geometry=o["geometry"],
               counts=C,original_counts_recomputed=O),open("counts.json","w"),indent=1)
PL={"grad":"surround + density gradient","unif":"surround, no gradient","orig":"original (shell only)"}
with open("summary_table.md","w") as f:
    G=o["geometry"]
    f.write(f"Baseline (no web): RMS = {base:.1f}× floor; |m3–5| = {q35_0/F:.1f}×, |m1–2| = {q12_0/F:.1f}×. Four fitted clumps: 83×. 1000 realizations per setting.\n\n")
    fb=C["grad_bare_N500_a"]["frac_beyond"]
    f.write(f"Surround placement: {G['region']}. Shell {G['shell_gpc'][0]:.1f}–{G['shell_gpc'][1]:.1f} Gpc from CoE; beyond D+d = {G['beyond_gpc'][0]:.1f}–{G['beyond_gpc'][1]:.1f} Gpc from CoE "
            f"(≈{100*fb:.0f}% of objects, set by volume × density, not forced). Density: {G['density_law']}; {G['gradient_rule']}.\n")
    f.write("Masses, cluster/void split, zero net contrast, residual, best-scaling and scoring identical to the original.\n")
    for tag,t in (("bare","objects with centre >14.3 Gpc from observer (spec)"),("buf","whole object outside vision sphere (centre > 14.3 Gpc + own radius)")):
        f.write(f"\n### {t}\n\n|N|level|placement|median RMS/floor|5th pct|min|best-scaled median|best 5th pct|best min|RMS lower / higher (placed)|m3–5 lower|m1–2 higher|3–5 lower & 1–2 not higher|best: 3–5 lower|best: 1–2 higher|best: both|\n|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|\n")
        for N in (200,500):
            for lev in "abc":
                for pl in ("grad","unif","orig"):
                    X=O[f"{tag}_N{N}_{lev}"] if pl=="orig" else C[f"{pl}_{tag}_N{N}_{lev}"]; p=X["placed"]; b=X["best"]
                    f.write(f"|{N}|{LBL[lev]}|{PL[pl]}|{p['median']:.1f}|{p['p5']:.1f}|{p['min']:.1f}|{b['median']:.1f}|{b['p5']:.1f}|{b['min']:.1f}|{p['rms_lower']} / {p['rms_higher']}|{p['m35_lower']}|{p['m12_higher']}|{p['m35_lower_and_m12_not_higher']}|{b['m35_lower']}|{b['m12_higher']}|{b['m35_lower_and_m12_not_higher']}|\n")
    f.write("\n### Where the objects sit (mean over realizations)\n\n|setting|share beyond D+d|share in far hemisphere (away from CoE)|share within 28.6 Gpc of us|\n|-|-|-|-|\n")
    for k in ("grad_bare_N500_c","unif_bare_N500_c"): f.write(f"|{k}|{C[k]['frac_beyond']:.2f}|{C[k]['frac_far_hemisphere']:.2f}|{C[k]['frac_within_28p6_gpc']:.2f}|\n")
    f.write("|original shell placement (100,000 draws)|0.00|0.03|0.05|\n")
print(open("summary_table.md").read())
