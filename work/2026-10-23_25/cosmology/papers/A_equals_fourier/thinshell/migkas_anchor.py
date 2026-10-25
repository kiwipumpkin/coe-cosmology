"""Anchor <H> on Migkas et al. absolute H0 values: H(l) = <H> * A(l)/<A>, A = model af with fitted sine (planar, galactic l)."""
from tm_core import *
OUT = "/workspace/cosmology/papers/A_equals_fourier/thinshell/"
S_ = {}
for key in ("free", "omb"):
    a, d = newt_params(key); D = a["DR"] * RV; b = dict(a, DR=D / RV); rr, _ = newt_run(b)
    As = A_from(shifted(rr["af"], d["l0"]), TH1 - d["l0"], d["amp"], d["phi"]); S_[key] = As / As.mean()
RL = json.load(open(BUB + "/universe_rel_mcfree.json")); Ar = np.array(RL["A"]); S_["rel"] = Ar / Ar.mean()
for k, s in S_.items():
    print(k, "R2 vs Migkas pts", round(r2(MA, s[IMIG] * (MA.mean() / s[IMIG].mean())), 3) if True else "", "min/max l", int(np.argmin(s)), int(np.argmax(s)), round(s.min(), 3), round(s.max(), 3))
def pt(s, l): return s[int(round(l)) % 360]
def outside(s, l, th):
    dl = np.abs((TH1 - l + 180) % 360 - 180); return s[dl > th].mean()
# (label, paper, [(H, sig, kind, l, theta)]) kind: 'pt' = cone centre, 'out' = rest of sky outside cone / opposite hemisphere
CH = [
 ("M20 our sample LX-T, max/min", "M20 Sect.7", [(75.17, 1.81, "pt", 17, 0), (66.20, 1.72, "pt", 274, 0)]),
 ("M20 ACC LX-T, max/min", "M20 Sect.8.1", [(78.76, 4.15, "pt", 77, 0), (58.12, 2.68, "pt", 317, 0)]),
 ("M20 XCS-DR1 LX-T, max/min", "M20 Sect.8.1", [(77.91, 2.20, "pt", 31, 0), (63.56, 2.32, "pt", 281, 0)]),
 ("M20 joint (3 samples), max/min", "M20 Sect.8.2", [(76.64, 1.41, "pt", 34, 0), (65.20, 1.48, "pt", 303, 0)]),
 ("M21 our LX-T, region/rest", "M21 Sect.6.1.1", [(66.0, 1.7, "pt", 274, 0), (72.1, 1.4, "out", 274, 75)]),
 ("M21 ours+ACC LX-T, region/rest", "M21 Sect.6.1.2", [(66.2, 1.6, "pt", 284, 0), (72.7, 1.5, "out", 284, 75)]),
 ("M21 our YSZ-T, region/rest", "M21 Sect.6.2", [(62.3, 2.2, "pt", 268, 0), (72.1, 0.9, "out", 268, 60)]),
 ("M21 LBCG-T, region/opposite hemisphere", "M21 Sect.6.3", [(66.2, 2.5, "pt", 257, 0), (72.6, 2.3, "out", 257, 90)]),
 ("M21 ACC LX-T, region/rest", "M21 App.B.2.1", [(61.3, 4.2, "pt", 318, 0), (72.2, 2.4, "out", 318, 75)]),
 ("M21 ACC YSZ-T, region/opposite hemisphere", "M21 App.B.3.1", [(60.6, 3.6, "pt", 311, 0), (73.4, 1.9, "out", 311, 90)]),
 ("M21 ALL COMBINED (joint), region/rest", "M21 Sect.6.4.1", [(66.5, 1.0, "pt", 273, 0), (72.8, 0.6, "out", 273, 75)]),
]
res = []
for lab, src, P in CH:
    row = dict(label=lab, src=src, pts=P)
    for key, s in S_.items():
        m = np.array([pt(s, l) if kd == "pt" else outside(s, l, th) for H, sg, kd, l, th in P]); H = np.array([p[0] for p in P]); w = 1 / np.array([p[1] for p in P])**2
        Hm = (w * H * m).sum() / (w * m * m).sum(); eH = 1 / np.sqrt((w * m * m).sum())
        row[key] = dict(H=Hm, err=eH, model_ratio=m[0] / m[1], obs_ratio=H[0] / H[1], pred=list(Hm * m))
    res.append(row)
for r in res:
    print("%-44s" % r["label"], " ".join("%s %.2f±%.2f (ratio mod %.3f obs %.3f)" % (k, r[k]["H"], r[k]["err"], r[k]["model_ratio"], r[k]["obs_ratio"]) for k in S_))
summ = {}
for k in S_:
    v = np.array([r[k]["H"] for r in res]); summ[k] = dict(joint=res[-1][k]["H"], joint_err=res[-1][k]["err"], mean=v.mean(), sd=v.std(ddof=1), min=v.min(), max=v.max())
print(json.dumps(summ, indent=1))
json.dump(dict(choices=res, summary=summ, shape_minmax={k: [int(np.argmin(s)), float(s.min()), int(np.argmax(s)), float(s.max())] for k, s in S_.items()}), open(OUT + "migkas_anchor.json", "w"), indent=1, default=float)
