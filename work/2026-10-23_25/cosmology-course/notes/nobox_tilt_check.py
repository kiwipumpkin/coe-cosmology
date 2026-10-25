"""Accuracy check of the no-boxcar + sine-tilt extremum formulas against the exact
no-boxcar curve with tilt (model functions unchanged). No fitting."""
import numpy as np
from nobox_check import exact_ext, h_exact, K0
FU = 0.94; AMP = -0.06; PH = 0.35

def tilt_formulas(D, L, d, M, al=1.0, G=1.0, H=1.0, a=AMP, ph=PH, half=True):
    S = (1-FU)*(6*L**3*np.exp(D/L) - 6*L**3 - 6*D*L**2 - 3*D**2*L) + FU*D**3   # written out on the card
    g = al*G*M/H**2; ig = H**2/(al*G*M)
    Aall = 1 + g/2*(1/D**3 - 3/S)
    hpk = lambda sg: 1 + g*(2/D**3 - 3/S + sg*3*d/D*(1/D**3 - (2*L + (1-FU)*D)/(2*L*S)))
    hmin = 1 - g/D**3
    def delta(sg):  # sg=-1 -> delta_1 (d), +1 -> delta_2 (d -> -d)
        num = ig + 2/D**3 - 3/S + sg*3*d/D*(1/D**3 - (2*L + (1-FU)*D)/(2*L*S))
        den = 6/D**3 - 6/S + sg*9*d/(2*D)*(4/D**3 - (4*L + (1-FU)*D)/(L*S))
        return -a*np.cos(ph)*num/den
    d1, d2 = delta(-1), delta(+1)
    b = 3*d/(4*D)
    c = (ig - 1/D**3)/(6/D**3 - 6/S)
    dp1 = -a*(np.sin(ph) + b*np.cos(ph))*c
    dp2 = -a*(np.sin(ph) - b*np.cos(ph))*c
    l = (np.degrees(d1), 180 - np.degrees(d2), 90 + np.degrees(b + dp1), 270 - np.degrees(b + dp2))
    hf = 0.5 if half else 0.0
    A = (hpk(-1)/Aall*(1 - a*np.sin(ph) - hf*a*np.cos(ph)*d1),
         hpk(+1)/Aall*(1 + a*np.sin(ph) - hf*a*np.cos(ph)*d2),
         hmin/Aall*(1 - a*np.cos(ph) + a*b*np.sin(ph) + hf*a*np.sin(ph)*dp1),
         hmin/Aall*(1 + a*np.cos(ph) + a*b*np.sin(ph) + hf*a*np.sin(ph)*dp2))
    meanT = Aall - 9/16*a*np.sin(ph)*g*d/D*((L + (1-FU)*D)/(L*S) - 1/D**3)  # mean of T*H_obs/H
    return l, A, Aall, meanT, (d1, d2, dp1, dp2)

def cd(u, v): return abs((u - v + 180) % 360 - 180)

def exact_sorted(x, e, al, ka):
    ex, hbar, _ = exact_ext(x, e, al, ka, amp=AMP)
    return ex, hbar

def compare(x, e, al, ka, half=True):
    ex, hbar = exact_sorted(x, e, al, ka)
    l, A, *_ = tilt_formulas(1.0, 1/x, e, ka, al=al, half=half)
    mx = [(ll, AA) for t, ll, AA in ex if t == "max"]; mn = [(ll, AA) for t, ll, AA in ex if t == "min"]
    if len(mx) != 2 or len(mn) != 2: return None
    dl, dA = [], []
    for i in range(4):
        pool = mx if i < 2 else mn
        j = min(pool, key=lambda z: cd(z[0], l[i]))
        dl.append(cd(j[0], l[i])); dA.append(abs(j[1] - A[i]))
    return dl, dA

if __name__ == "__main__":
    x, e, al, ka = 5.9, 0.1, 1.0, K0
    ex, hbar = exact_sorted(x, e, al, ka)
    th = np.radians(np.arange(0, 360, 0.01)); h, _ = h_exact(th, x, e, al, ka)
    l, A, Aall, meanT, ds = tilt_formulas(1.0, 1/x, e, ka)
    print("exact:", [(t, round(ll, 2), round(AA, 4)) for t, ll, AA in ex], "hbar", hbar,
          "mean(T h)", np.mean((1 - AMP*np.sin(th + PH))*h))
    print("formula l", np.round(l, 2), "A", np.round(A, 4), "Aall", Aall, "meanT", meanT, "deltas deg", np.degrees(ds))
    print("no-half A", np.round(tilt_formulas(1.0, 1/x, e, ka, half=False)[1], 4))
    SW = [("x", 1.01, 15, {}), ("al", 0.01, 8, {}), ("ka", 0.001, 2, {}), ("e", 0.005, 0.999, {"al": 1.3})]
    for key, lo, hi, over in SW:
        print("==", key)
        for v in np.linspace(lo, hi, 40):
            kw = dict(x=5.9, e=0.1, al=1.0, ka=K0); kw.update(over); kw[key] = v
            r = compare(**kw); r0 = compare(**kw, half=False)
            if r is None: print("  %.4g  extrema count != 2 max / 2 min" % v); continue
            print("  %.4g  dl=%s  dA=%s  dA(noHalf)=%s" % (v, "/".join("%.2f" % t for t in r[0]),
                  "/".join("%.4f" % t for t in r[1]), "/".join("%.4f" % t for t in r0[1])))
