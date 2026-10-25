import setup as S, ltb_thin as LT, numpy as np
k = 70 / S.H0
for t0 in (2.0, 2.4, 2.8, 3.0, 3.2):
    p = S.base(); p[3] = 1; p[5] = t0
    Hp, Hr = LT.local_rates(S.D, t0, p); E = LT.geom(S.D, t0, p)[5]
    print("A t0 %.3f (%.2f Gyr): H_perp %.1f H_r %.1f Hloc %.1f, 2E(D) %.1f" % (t0, t0*S.GYR, Hp*k, Hr*k, (Hr+2*Hp)/3*k, 2*E))
p, t0, ro = S.model_areal(); print("GR-A solved t0", t0, t0*S.GYR, [x*k for x in LT.local_rates(ro, t0, p)])
p, t0, ro = S.model_lagr()
R, Rr, Rrr, Rt, Rtr, E, E1, M, M1 = LT.geom(ro, t0, p)
print("GR-B: ro %.3f t0 %.4f (%.3f Gyr)  H_perp %.2f H_r %.2f  2E %.2f  2M/R %.2f Rdot %.3f  rho/crit %.4f" % (ro, t0, t0*S.GYR, Rt/R*k, Rtr/Rr*k, 2*E, 2*M/R, Rt, LT.density(ro,t0,p)/S.RHOC))
for lab, r in (("inner edge", p[9]), ("outer edge", p[10]), ("core Lc", p[2])):
    g = LT.geom(r, t0, p); print(lab, "label %.2f -> areal today %.2f" % (r, g[0]))
