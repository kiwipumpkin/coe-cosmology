import setup as S, ltb_thin as LT, numpy as np
from scipy.optimize import brentq
p = S.base()
def tpar(X): M = LT.massfun(X, p)[0]; return (2/3)*np.sqrt(X**3/(2*M))
for X in (S.Rin+1, S.D, S.Rout-1, S.Rout, 120, 150):
    print("X %.2f  M %.1f  2M/X %.1f  t_parabolic %.3f Gpc/c = %.2f Gyr" % (X, LT.massfun(X,p)[0], 2*LT.massfun(X,p)[0]/X, tpar(X), tpar(X)*S.GYR))
# GR-A
for t0 in (3.4, 3.6, 3.65, 3.7, 3.9, 4.24):
    p = S.base(); p[3] = 1; p[5] = t0
    try:
        Hp, Hr = LT.local_rates(S.D, t0, p); E = LT.geom(S.D, t0, p)[5]
        es = []
        for X in np.linspace(S.Rin, S.Rout, 60):
            es.append(LT.geom(X, t0, p)[5])
        es = np.array(es)
        print("t0 %.3f (%.2f Gyr): H_perp %.1f H_r %.1f Hloc %.1f km/s/Mpc, 2E(D) %.1f, min 2E in shell %.2f" % (t0, t0*S.GYR, Hp/S.H0*70, Hr/S.H0*70, (Hr+2*Hp)/3/S.H0*70, 2*E, 2*es.min()))
    except Exception as e: print(t0, "fail", e)
