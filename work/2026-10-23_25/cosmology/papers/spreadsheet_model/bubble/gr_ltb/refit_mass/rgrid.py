"""R = Migkas sample distance. Chris's layout (core + empty interior + uniform shell D-R..D+R + vacuum), model A linear
(A propto H(psi) at fixed z, i.e. D_L^-1), exact LTB GR-A, H = 70, t_B = 0. Sampling z matched to R:
 headline R = 0.316 Gpc (median z = 0.075, Migkas 2020), secondary R = 1.19 Gpc (z = 0.3). Coarse grid; l0, Amp, phi fitted."""
import json, sys, numpy as np, multiprocessing as mp
import rm_core as C
CASES = {"median": (0.316, 0.075), "zmax": (1.19, 0.3)}
DRS = [1.5, 3.0, 6.55, 12.0]; OSS = [0.02, 0.05, 0.2, 1.0]; MU = np.geomspace(0.003, 1.0, 9)
def one(a):
    case, DR, Os, mu = a; R, z = CASES[case]; C.ZL = (z,); D = DR * R
    Mc = mu * 4 * np.pi / 3 * C.RHOC * D**3 * C.KG
    try: o = C.evaluate(Mc, D, R, Os=Os)
    except Exception as e: o = dict(fail=str(e))
    o.update(case=case, DR=DR, Os=Os, mu=mu, z=z, D=D, R=R, Mc_kg=Mc); return o
if __name__ == "__main__":
    T = [(c, d, o, m) for c in CASES for d in DRS for o in OSS for m in MU]
    with mp.Pool(7) as P: res = P.map(one, T, chunksize=2)
    json.dump(res, open("rgrid.json", "w"), default=float)
    for c in CASES:
        V = [o for o in res if o["case"] == c and "R2_s" in o and o["twoE_min"] > -1]
        print(c, "valid", len(V), "of", sum(o["case"] == c for o in res))
        for o in sorted(V, key=lambda o: -o["R2_s"])[:6]:
            f = o["fits"][str(o["z"])]
            print("  DR %.2f D %.3f Os %.2f mu %.4f Mc %.2e | R2s %.3f R2ns %.3f l0s %.1f amp %.3f phi %.2f | along %s 2E %.3f age %.2f Hp %.1f Hr %.1f m2 %.3f" % (
                o["DR"], o["D"], o["Os"], o["mu"], o["Mc_kg"], o["R2_s"], o["R2_ns"], f["l0_s"], f["amp"], f["phi"], o["along"], o["twoE"], o["age_Gyr"], o["H_perp"], o["H_r"], f["m2"]))
