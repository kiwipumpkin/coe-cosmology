"""Can the full numerical model (dq_model_lib) give the measured a3/a2 with d/D <= 0.022?"""
import sys, json, numpy as np
sys.path.insert(0, "/workspace/cosmology/figures")
import dq_model_lib as M
meas = {"m1": 0.0398, "m2": 0.0733, "m3": -0.0292}
def proj(h, l0, m):  # coefficient of cos m(l-l0)
    return h[f"c{m}"]*np.cos(m*np.radians(l0)) + h[f"s{m}"]*np.sin(m*np.radians(l0))
out = []
for eps in (0.005, 0.0225, 0.1, 0.3, 0.76):
    for x in (3, 5, 7.62, 10, 14, 20):
        for b0 in (0.0, 15.0):
            l, A = M.model_fine(0.458, eps, x, 15.0, b0)
            h = M.harm_fine(l, A)
            p = {m: proj(h, 15.0, m) for m in (1, 2, 3)}
            out.append(dict(eps=eps, x=x, b0=b0, m1=p[1], m2=p[2], m3=p[3], r31=p[3]/p[2], r12=p[1]/p[2]))
            print(f"eps={eps:<6} D/L={x:<5} b0={b0:<4} m1={p[1]:+.5f} m2={p[2]:+.5f} m3={p[3]:+.5f}  m3/m2={p[3]/p[2]:+.4f} m1/m2={p[1]/p[2]:+.4f}")
print("measured m3/m2 =", meas["m3"]/meas["m2"], " m1/m2 =", meas["m1"]/meas["m2"])
json.dump(out, open("/workspace/cosmology/papers/dq_equality/scan_small_dD.json", "w"), indent=1)
