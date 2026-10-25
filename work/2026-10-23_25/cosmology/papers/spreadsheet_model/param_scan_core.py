"""Chris's BH column as a function of his parameter cells (formulas identical to the sheet)."""
import numpy as np, csv, os, openpyxl
HERE = os.path.dirname(os.path.abspath(__file__))
# fixed cells of working tab 'Copy of Copy of Copy of density'
R = float(os.environ.get("PSCAN_R", "4.4127e23")); k0 = 2.18e-18; G = 6.67e-17; dt = 1.0; E39 = 0.09; W = 0.01745329252
CUR = dict(DR=5.0, LD=0.18, M=4.44e68, amp=-0.06, phi=0.35, l0=0.0)

def BH(theta_deg, DR, LD, M, amp, phi, k=k0):
    D = R * DR; L = D * LD
    B = np.radians(theta_deg)
    C = D + R * np.cos(B); Dy = R * np.sin(B)
    E = (C * C + Dy * Dy) ** 0.5
    F = k * E
    Gc = np.arctan(Dy / C)
    H = F * np.cos(Gc); I = F * np.sin(Gc)
    Vo = k * D; J = H - Vo; K = I
    Lc = np.sqrt(J * J + K * K)
    xD = D / L; frac = 1 - np.exp(-xD) * (1 + xD + xD * xD / 2); Mtot = M / frac
    x = E / L; AO = (1 - np.exp(-x) * (1 + x + x * x / 2)) * Mtot
    N = -G * AO / (E * E)
    O = N * np.cos(Gc); P = N * np.sin(Gc)
    Axo = -G * M / (D * D)
    Q = O - Axo; Rr = P
    S = J + Q * dt; T = Rr + K * dt          # T exactly as in the sheet
    U = np.sqrt(S * S + T * T)
    AF = U / Lc
    BA = amp * np.sin(W * theta_deg + phi)
    return (1 - BA) * AF

def AFfour(l):
    l = np.radians(l)
    return (0.9965 + 0.0296*np.cos(l) + 0.0434*np.sin(l) + 0.0636*np.cos(2*l) + 0.0366*np.sin(2*l)
            - 0.0299*np.cos(3*l) - 0.0114*np.sin(3*l))

mig = list(csv.DictReader(open(os.path.join(HERE, "../../migkas/fig6_digitized.csv"))))
ML = np.array([float(r["l_deg"]) for r in mig]); MA = np.array([float(r["A_over_Aall"]) for r in mig])
MLO = np.array([float(r["A_err_low"]) for r in mig]); MHI = np.array([float(r["A_err_high"]) for r in mig])
LG = np.arange(0, 360, 1.0); AFG = AFfour(LG)
# sheet's own lookup-mapped Migkas values (cached column AG, rows theta=0..360)
_wv = openpyxl.load_workbook(os.path.join(HERE, "chris_model.xlsx"), data_only=True)["Copy of Copy of Copy of density"]
AGS = np.array([_wv.cell(r, 33).value for r in range(66, 427)], float); THS = np.arange(0, 361, 1.0)

def r2(y, m): return 1 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2)
TARGETS = {
  "mig":   lambda p: r2(MA,  BH(ML - p.get("l0", 0.0), p["DR"], p["LD"], p["M"], p["amp"], p["phi"])),
  "four":  lambda p: r2(AFG, BH(LG - p.get("l0", 0.0), p["DR"], p["LD"], p["M"], p["amp"], p["phi"])),
  "sheet": lambda p: r2(AGS, BH(THS, p["DR"], p["LD"], p["M"], p["amp"], p["phi"])),
}
