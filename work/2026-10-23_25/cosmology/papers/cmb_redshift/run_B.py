import numpy as np, json, time, models as m, analyze as an, ltb
t1 = time.time()
p, t0, ro = m.model_B()
info = m.describe(p, t0, ro, "B")
c = an.cmb(p, t0, ro, N=96)
nb_ = an.nearby(p, t0, ro)
print("time", time.time() - t1)
np.savez("res_B.npz", p=p, t0=t0, ro=ro, **{k: v for k, v in c.items() if isinstance(v, np.ndarray)},
         **{"nb_" + k: np.array(v) for k, v in nb_.items()})
S = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in c.items() if k not in ("mu", "psi", "z_a", "z_b", "dT_a", "dT_b")}
S["info"] = info; S["nearby_dH_leg"] = [x.tolist() for x in nb_["dH_leg"]]
S["nearby_Hperp_Hr"] = [nb_["Hperp"] / m.H0, nb_["Hr"] / m.H0]
S["nearby_A_fourier_z"] = {str(z): an.fourier(nb_["ell"], nb_["A"][:, k]).tolist() for k, z in enumerate(nb_["zlist"])}
json.dump(S, open("res_B.json", "w"), indent=1, default=float)
print(json.dumps(S, indent=1, default=float))
