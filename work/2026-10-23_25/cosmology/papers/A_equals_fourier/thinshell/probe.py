from tm_core import *
for key in ("free", "omb"):
    a, d = newt_params(key); rr, m1 = newt_run(a); eps = 1 / a["DR"]
    af = rr["af"] / rr["af"].mean()
    qs = np.linspace(1 - 1.2 * eps, 1 + 1.2 * eps, 41); h, v1 = radial_h(a, rr, qs)
    cf = np.polyfit(qs - 1, h, 6); h1 = np.polyval(cf, 0); beta = np.polyval(np.polyder(cf), 0) / h1; gam = np.polyval(np.polyder(cf, 2), 0) / h1
    Ar = A_radial(lambda q: np.polyval(cf, q - 1), eps)
    print(key, "eps", eps, "Ht", rr["Ht"], "mu", rr["mu"], "beta", beta, "gamma", gam, "v1/h1", v1, h1)
    for m in range(4): print("  m", m, "code", np.round(harm(af, m), 5), "radial-exact", np.round(harm(Ar, m), 5))
