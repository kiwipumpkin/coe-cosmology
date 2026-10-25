import sys, time, numpy as np
from ltb_common import *
F = fitparams(); k = sys.argv[1]; cav = sys.argv[2] == "cav"; f = F[k]
p = make_p(f["Mc"], f["D"], f["Os"], f["Rin"], 5000.0, cavity=cav); t0 = setup(p, f["D"]); print(k, "t0", t0, flush=True)
mu, w = leggauss(24)
for m in mu:
    t = time.time()
    try:
        s = rays.trace(p, t0, f["D"], np.arccos(m), t0 * 5e-5)
        rmin = s.y[0].min(); z = rays.z_at_t(s, t0, t0 * 1e-4) if s.status == 0 else np.nan
        print("psi %.1f deg status %d rmin %.3f z(1e-4 t0) %.1f sec %.1f" % (np.degrees(np.arccos(m)), s.status, rmin, z, time.time() - t), flush=True)
    except Exception as e: print("psi %.1f fail %s" % (np.degrees(np.arccos(m)), str(e)[:80]), flush=True)
