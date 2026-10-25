"""GR-A family (today's exact mass distribution, uniform bang): largest age with 1 + 2E > 0 on every matter shell,
and Earth-shell numbers at the Newtonian age 13.83 Gyr and for a bound Earth shell."""
import setup as S, ltb_thin as LT, numpy as np, json
from scipy.optimize import brentq
K = 70 / S.H0
def min2E(t0):
    p = S.base(); p[3] = 1; p[5] = t0
    return min(2 * LT.geom(r, t0, p)[5] for r in np.linspace(S.Rin - S.W, S.Rout + S.W, 400))
tmax = brentq(lambda t: min2E(t) + 1, 3.5, 3.8, xtol=1e-10)
p = S.base(); p[3] = 1; p[5] = tmax; Hp, Hr = LT.local_rates(S.D, tmax, p)
tN = S.N["age"] / S.GYR; p2 = S.base(); p2[3] = 1; p2[5] = tN; Hp2, Hr2 = LT.local_rates(S.D, tN, p2)
M = S.M_D; D = S.D
out = dict(t0_max_Gyr=tmax * S.GYR, at_tmax=dict(H_perp=Hp * K, H_r=Hr * K, H_loc=(Hr + 2 * Hp) / 3 * K, twoE_D=2 * LT.geom(D, tmax, p)[5]),
           at_newtonian_age=dict(age_Gyr=S.N["age"], min2E_matter=min2E(tN), H_perp=Hp2 * K, H_r=Hr2 * K, H_loc=(Hr2 + 2 * Hp2) / 3 * K, twoE_D=2 * LT.geom(D, tN, p2)[5]),
           earth_bound_window=dict(H_perp_min=np.sqrt(2 * M / D - 1) / D * K, H_perp_max=np.sqrt(2 * M / D) / D * K,
                                   age_parabolic_Gyr=(2 / 3) * np.sqrt(D**3 / (2 * M)) * S.GYR, Rmax_min_Gpc=2 * M, t_turn_min_Gyr=np.pi * M * S.GYR),
           newtonian_comoving_check=dict(twoE_if_V_eq_H0D=(S.H0 * D)**2 - 2 * M / D, V2_over_vesc2_enclosed=(S.H0 * D)**2 / (2 * M / D),
                                         V2_over_vesc2_enclosed_VO=S.N["vE"]**2 / (2 * M / D)))
print(json.dumps(out, indent=1, default=float)); json.dump(out, open("agelimit.json", "w"), indent=1, default=float)
