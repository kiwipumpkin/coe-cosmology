import numpy as np, json, models as m, analyze as an, warnings
warnings.filterwarnings("ignore")
p, t0, ro = m.model_C(0.1 * m.M)
c = an.cmb(p, t0, ro, N=64); nb_ = an.nearby(p, t0, ro)
np.savez("res_C.npz", t0=t0, beta_a=c["beta_a"], **{k: v for k, v in c.items() if isinstance(v, np.ndarray)}, **{"nb_" + k: np.array(v) for k, v in nb_.items()})
print(c["c_a"][:4], c["cb_a"][:4], c["beta_a"], an.fourier(nb_["ell"], nb_["A"][:, 1]))
