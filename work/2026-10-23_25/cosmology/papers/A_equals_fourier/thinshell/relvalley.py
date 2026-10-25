from tm_core import *
import relfree as RF, rm_core as C, time
J = json.load(open("gr_ltb/refit_mass/relfree2.json")); F = [o for o in J["refined"] if "R2_s" in o]
S = sorted([o for o in J["screen"] if "R2_s" in o], key=lambda o: -o["R2_s"])
pts = [F[1], F[3], S[2], F[4], F[2]]
out = []
for o in pts:
    t = time.time(); e = RF.ev(o["Mc_kg"], o["D"], o["kin"], o["kout"], full=True)
    if "Ad" not in e: print("fail", o["D"]); continue
    g = C.curve(np.array(e["Ad"]), 0.0); h = [harm(g, m)[0] for m in range(4)]
    out.append(dict(D=o["D"], R_in=o["R_in"], R_out=o["R_out"], mcc=e["mcc"], R2s=e["R2_s"], R2ns=e["R2_ns"], beta=e["H_r"] / e["H_perp"] - 1, a=h, eps=1.19 / o["D"]))
    print("D %.4g mcc %.4f beta %.4f a1 %.4f a2 %.4f a3 %.4f R2s %.4f  (%.0fs)" % (o["D"], e["mcc"], out[-1]["beta"], h[1], h[2], h[3], e["R2_s"], time.time() - t), flush=True)
json.dump(out, open("/workspace/cosmology/papers/A_equals_fourier/thinshell/relvalley.json", "w"), indent=1, default=float)
