from bubble_core import *
import time
# isotropy check: background only
for Ob in (0.1, 0.3, 1.0):
    rr = run(2.9, 0.07, 1e-9, Ob); a = rr["af"]; print("bg only Ob", Ob, "aniso (max-min)/mean", (a.max()-a.min())/a.mean(), "mean", a.mean())
p0 = dict(DR=2.904374576719836, LD=0.07027609488588729, mc=0.02395, Ob=0.0, l0=13.166, )
for Ob in (0.0, 0.3):
    p = dict(p0, Ob=Ob)
    for norm in ("kR", "mean"):
        v, a, ph, rr = score(p, "mig", norm)
        print("prev best Ob", Ob, norm, "R2", round(v, 4), "contrast", rr["contrast"], "c34", rr["c34"], "vesc2", rr["vesc2"], "age", rr["age"], "Ht", rr["Ht"])
t = time.time(); run(2.9, 0.07, 0.02, 0.3); print("t/run", time.time() - t)
