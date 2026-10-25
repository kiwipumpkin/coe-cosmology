"""Recompute Chris's per-angle columns literally from his cell formulas
(workbook chris_model.xlsx = attachment b3737b04...; tabs 'density w rotation'
and its copies share identical row formulas, only parameter cells differ)."""
import json, math, numpy as np, openpyxl, csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
FN = os.path.join(HERE, "chris_model.xlsx")
wv = openpyxl.load_workbook(FN, data_only=True)
SHEETS = ["density w rotation", "Copy of density w rotation",
          "Copy of Copy of density w rotat", "Copy of Copy of Copy of density"]
PRIMARY = "Copy of Copy of Copy of density"
COLS = "B C D E F G H I J K L N O P Q R S T U X Y AB AF AN AO BA BE BH BL AG AM".split()

def migkas_table():
    ws = wv["migkas"]
    n = sum(1 for r in range(1, ws.max_row + 1) if isinstance(ws.cell(r, 6).value, (int, float)))
    E = []; e = 360.0 / n
    for i in range(n):                      # E3 = 1/n*360 ; E(k)=E(k-1)+360/n
        E.append(e); e = e + 360.0 / n
    Ecached = [ws.cell(3 + i, 5).value for i in range(n)]   # values as stored by the sheet
    print("migkas!E recomputed vs cached max diff:", max(abs(a - b) for a, b in zip(E, Ecached)))

    F = [ws.cell(3 + i, 6).value for i in range(n)]
    G = [ws.cell(3 + i, 7).value for i in range(n)]
    H = [ws.cell(3 + i, 8).value for i in range(n)]
    return np.array(E), np.array(F, float), np.array(G, float), np.array(H, float)

ME, MF, MG, MH = migkas_table()
def vlookup_idx(key):                        # VLOOKUP(...,TRUE): largest E <= key
    i = np.searchsorted(ME, key + 1e-9, side="right") - 1   # sheet grid hits nodes exactly
    return i

def params(ws):
    g = lambda a: ws[a].value
    k = g("B25")

    if ws.title == "Copy of Copy of density w rotat":   # cached B25 exported as 0; formula =65000/(3.1E+22)
        k = 65000 / 3.1e22
    return dict(k=k, G=g("B26"), M=g("B27"), D=g("B28"), R=g("B29"), L=g("B30"),
                Vo=g("B32"), dt=g("B33"), Axo=g("B34"), vstar=g("B36"), frac=g("B55"),
                Mtot=g("B56"), E39=g("E39"), ampA=g("BA57"), w=g("BA58"), phi=g("BA59"),
                ampB=g("BE57"))

def compute(p):
    th = np.arange(0, 361, 1.0)
    B = np.radians(th)
    C = p["D"] + p["R"] * np.cos(B)
    Dy = p["R"] * np.sin(B)
    E = (C * C + Dy * Dy) ** 0.5
    F = p["k"] * E
    G = np.arctan(Dy / C)                    # Excel ATAN (row 66 adds RADIANS(B58)=0)
    H = F * np.cos(G); I = F * np.sin(G)
    J = H - p["Vo"]; K = I
    L = np.sqrt(J * J + K * K)
    x = E / p["L"]
    AO = (1 - np.exp(-x) * (1 + x + x * x / 2)) * p["Mtot"]
    N = -p["G"] * AO / (E * E)
    O = N * np.cos(G); P = N * np.sin(G)
    Q = O - p["Axo"]; R = P
    S = J + Q * p["dt"]; T = R + K * p["dt"]
    U = np.sqrt(S * S + T * T)
    X = np.arctan2(T, S)                     # Excel ATAN2(x=S, y=T)
    Y = np.arctan(P / O)
    AB = L; AF = U / AB
    AN = AO / p["M"]
    BA = p["ampA"] * np.sin(p["w"] * th + p["phi"])
    BE = p["ampB"] * np.sin(p["w"] * th + p["phi"])
    BH = (1 - BA) * AF; BL = (1 - BE) * AN
    # Migkas lookups: row n uses key AD(n+10); last 10 rows wrap to AD(n-349)
    key = np.where(th <= 350, th + 10, th - 349 + 0.0)
    key = np.where(th <= 350, key, th - 350 + 1)   # rows 417..426 -> AD68..AD77 = 2..11
    idx = vlookup_idx(key)
    AG = MG[idx]; AM = MH[idx] - p["E39"]; lmap = MF[idx]
    return dict(th=th, B=B, C=C, D=Dy, E=E, F=F, G=G, H=H, I=I, J=J, K=K, L=L, N=N, O=O, P=P,
                Q=Q, R=R, S=S, T=T, U=U, X=X, Y=Y, AB=AB, AF=AF, AN=AN, AO=AO, BA=BA, BE=BE,
                BH=BH, BL=BL, AG=AG, AM=AM, lmap=lmap)

def cached(ws, col):
    ci = openpyxl.utils.column_index_from_string(col)
    return np.array([ws.cell(r, ci).value for r in range(66, 427)], float)

def AF_fourier(l):
    l = np.radians(l)
    return (0.9965 + 0.0296*np.cos(l) + 0.0434*np.sin(l) + 0.0636*np.cos(2*l) + 0.0366*np.sin(2*l)
            - 0.0299*np.cos(3*l) - 0.0114*np.sin(3*l))

def r2(y, m):
    return 1 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2)
def pearson2(y, m):
    return np.corrcoef(y, m)[0, 1] ** 2

mig = list(csv.DictReader(open(os.path.join(HERE, "../../migkas/fig6_digitized.csv"))))
ml = np.array([float(r["l_deg"]) for r in mig]); mA = np.array([float(r["A_over_Aall"]) for r in mig])

out = {}
for s in SHEETS:
    ws = wv[s]; p = params(ws); c = compute(p)
    chk = {}
    for col in COLS:
        ca = cached(ws, col); d = np.abs(c[col] - ca) / np.nanmax(np.abs(ca))

        chk[col] = float(np.nanmax(d))
    lk_mis = {col: int(np.sum(np.abs(c[col] - cached(ws, col)) > 1e-9)) for col in ("AG", "AM")}
    th = c["th"]; BH = c["BH"]; AFm = c["AF"]
    # (i) identity: theta == l
    Af = AF_fourier(th)
    res = BH - Af
    BH_at_m = np.interp(ml, th, BH)
    # (ii) his lookup mapping l(theta): average BH over theta mapped to each Migkas l
    hm = np.array([BH[c["lmap"] == l].mean() if np.any(c["lmap"] == l) else np.nan for l in ml])
    ok = ~np.isnan(hm)
    # his own sheet statistic: LINEST(BI,BH) R^2
    out[s] = dict(params=p, maxabsdiff=chk,
        BH_min=float(BH.min()), BH_max=float(BH.max()), BH_mean=float(BH[:-1].mean()),
        th_min=float(th[np.argmin(BH)]), th_max=float(th[np.argmax(BH)]),
        vsAF=dict(rms=float(np.sqrt(np.mean(res[:-1]**2))), maxabs=float(np.abs(res).max()),
                  R2=float(r2(Af[:-1], BH[:-1])), r2_pearson=float(pearson2(Af[:-1], BH[:-1]))),
        vsAF_noRot=dict(rms=float(np.sqrt(np.mean((AFm-Af)[:-1]**2))), R2=float(r2(Af[:-1], AFm[:-1]))),
        vsMig_identity=dict(n=int(len(ml)), rms=float(np.sqrt(np.mean((mA-BH_at_m)**2))),
                            R2=float(r2(mA, BH_at_m)), r2_pearson=float(pearson2(mA, BH_at_m))),
        vsMig_hismap=dict(n=int(ok.sum()), rms=float(np.sqrt(np.mean((mA[ok]-hm[ok])**2))),
                          R2=float(r2(mA[ok], hm[ok])), r2_pearson=float(pearson2(mA[ok], hm[ok]))),
        lookup_mismatch_rows=lk_mis,
        model_cols_maxreldiff=float(max(v for k, v in chk.items() if k not in ("AG", "AM"))),
        sheet_LINEST_R2_A=float(pearson2(cached(ws, "AG"), BH)),
        sheet_LINEST_R2_A_cachedcell=float(ws["BH58"].value),
        sheet_LINEST_R2_noRot=float(pearson2(cached(ws, "AG")[2:], AFm[2:])),
        sheet_LINEST_R2_noRot_cachedcell=float(ws["AF51"].value),
        AF_vs_Mig=dict(rms=float(np.sqrt(np.mean((mA-AF_fourier(ml))**2))), R2=float(r2(mA, AF_fourier(ml)))),
        resid_by_l={int(t): float(res[int(t)]) for t in range(0, 360, 30)})
    if s == PRIMARY:
        import csv as _c
        with open(os.path.join(HERE, "primary_columns.csv"), "w", newline="") as fh:
            wr = _c.writer(fh); keys = ["th"] + [k for k in COLS if k in c] + ["lmap"]
            wr.writerow(keys); [wr.writerow([c[k][i] for k in keys]) for i in range(len(th))]
        np.savez(os.path.join(HERE, "primary_series.npz"), **c, Af=Af, ml=ml, mA=mA, hm=hm)
json.dump(out, open(os.path.join(HERE, "recompute_results.json"), "w"), indent=1)
for s in SHEETS:
    o = out[s]; print("==", s)
    print(" model cols max rel diff:", o["model_cols_maxreldiff"], "lookup mismatches", o["lookup_mismatch_rows"])
    print(" BH min/max/mean", o["BH_min"], o["th_min"], o["BH_max"], o["th_max"], o["BH_mean"])
    print(" vsAF", o["vsAF"], "noRot", o["vsAF_noRot"])
    print(" vsMig id", o["vsMig_identity"]); print(" vsMig map", o["vsMig_hismap"])
    print(" sheet LINEST R2", o["sheet_LINEST_R2_A"], o["sheet_LINEST_R2_A_cachedcell"], o["sheet_LINEST_R2_noRot"], o["sheet_LINEST_R2_noRot_cachedcell"], " AF vs Mig", o["AF_vs_Mig"])
    print(" resid", o["resid_by_l"])
