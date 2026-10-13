import json, numpy as np
o=json.load(open("ensemble.json")); F=1e-5; base=o["base_rms_over_floor"]; A0=np.array(o["base_amps_over_floor"])*F
q35_0=np.linalg.norm(A0[2:5]); q12_0=np.linalg.norm(A0[0:2])
LBL={"a":"(a) 1e48–1e49","b":"(b) 3e49–3e50","c":"(c) 1e51–1e52"}
rows=[]; C={}
def cnt(d,ka,kr):
    r=d[kr]/F; A=d[ka]; q35=np.linalg.norm(A[:,2:5],axis=1); q12=np.linalg.norm(A[:,0:2],axis=1)
    lo=int((r<base).sum()); lo35=int((q35<q35_0).sum()); hi12=int((q12>q12_0).sum()); both=int(((q35<q35_0)&(q12<=q12_0)).sum())
    return dict(n=len(r),rms_lower=lo,rms_higher=len(r)-lo,m35_lower=lo35,m12_higher=hi12,m35_lower_and_m12_not_higher=both,
                median=float(np.median(r)),p5=float(np.percentile(r,5)),min=float(r.min()),
                median_q35=float(np.median(q35)/F),median_q12=float(np.median(q12)/F))
for tag in ("bare","buf"):
    for N in (200,500):
        for lev in "abc":
            k=f"{tag}_N{N}_{lev}"; d=np.load(f"raw_{k}.npz")
            C[k]=dict(placed=cnt(d,"amps","rms"),best=cnt(d,"amps_best","rms_best"))
json.dump(dict(base_rms_over_floor=base,base_q35_over_floor=q35_0/F,base_q12_over_floor=q12_0/F,counts=C),open("counts.json","w"),indent=1)
R=o["radii"]
with open("summary_table.md","w") as f:
    f.write(f"Baseline (no web): RMS = {base:.1f}× floor; |m3–5| = {q35_0/F:.1f}×, |m1–2| = {q12_0/F:.1f}× (quadrature sums of A_m). 1000 realizations per setting.\n\n")
    f.write("Mass scale → structure: radius of a δ=+1 cluster / δ=−0.5 void holding |M| (ρ̄ = 2.9e-27 kg/m³). Check: 7.6e51 kg → %.2f / %.2f Gpc.\n\n"%tuple(o["radius_check_7.6e51"]))
    f.write("|level|kg per object|cluster radius δ=+1 (Gpc)|void radius δ=−0.5 (Gpc)|\n|-|-|-|-|\n")
    for lev in "abc":
        r=R[lev]; f.write(f"|{LBL[lev]}|{r['M_lo']:.1e}–{r['M_hi']:.1e}|{r['R_cluster_d+1_gpc'][0]:.2f}–{r['R_cluster_d+1_gpc'][1]:.2f}|{r['R_void_d-0.5_gpc'][0]:.2f}–{r['R_void_d-0.5_gpc'][1]:.2f}|\n")
    for tag,t in (("bare","objects with centre >14.3 Gpc from observer (spec)"),("buf","whole object outside vision sphere (centre > 14.3 Gpc + own radius)")):
        f.write(f"\n### {t}\n\n|N|level|median RMS/floor|5th pct|min|best-scaled median|best 5th pct|best min|RMS lower / higher (placed)|m3–5 lower|m1–2 higher|3–5 lower & 1–2 not higher|best: 3–5 lower|best: 1–2 higher|best: both|\n|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|\n")
        for N in (200,500):
            for lev in "abc":
                k=f"{tag}_N{N}_{lev}"; p=C[k]["placed"]; b=C[k]["best"]
                f.write(f"|{N}|{LBL[lev]}|{p['median']:.1f}|{p['p5']:.1f}|{p['min']:.1f}|{b['median']:.1f}|{b['p5']:.1f}|{b['min']:.1f}|{p['rms_lower']} / {p['rms_higher']}|{p['m35_lower']}|{p['m12_higher']}|{p['m35_lower_and_m12_not_higher']}|{b['m35_lower']}|{b['m12_higher']}|{b['m35_lower_and_m12_not_higher']}|\n")
print(open("summary_table.md").read())
