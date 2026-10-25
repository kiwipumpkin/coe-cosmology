import numpy as np, models as m, rays, ltb, time
for name, mk in (("A", m.model_A), ("B", m.model_B)):
    p, t0, ro = mk()
    for psi in (0.02, np.pi / 2, np.pi - 0.02):
        t1 = time.time(); s = rays.trace(p, t0, ro, psi, t0 * 1e-7)
        print(name, "psi=%.2f nfev=%d  %.1fs null=%.1e" % (psi, s.nfev, time.time() - t1, rays.null_err(p, t0, s)))
        for f in (0.5, 0.1, 0.03, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
            y = s.sol(np.log(1 / f)); R = ltb.geom(y[0], t0 * f, p)[0]
            print("   t/t0=%g z=%.6g r=%.4g R=%.4g th=%.3f rho*6pi t^2=%.4g" % (f, -y[2] - 1, y[0], R, y[1], ltb.density(y[0], t0 * f, p) * 6 * np.pi * (t0 * f)**2))
