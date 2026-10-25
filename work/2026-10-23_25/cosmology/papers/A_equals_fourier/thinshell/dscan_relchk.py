from tm_core import *
import relfree as RF, json, multiprocessing as mp
RHO = 9.2042e-27
def one(a):
    DR, mcc, ki, ko = a; D = DR * RV; Mc = mcc * 4 * np.pi / 3 * RHO * (D * 3.0857e25)**3
    o = RF.ev(Mc, D, ki, ko); Msh = 4 * np.pi / 3 * 0.049 * RHO * (((D + ko * RV) * 3.0857e25)**3 - ((D - ki * RV) * 3.0857e25)**3)
    return dict(DR=DR, mcc=o.get("mcc"), kin=ki, kout=ko, R2=o.get("R2_s"), share=Mc / (Mc + Msh))
if __name__ == "__main__":
    J = [(10, m, 1.03, 1.12) for m in (0.065, 0.075, 0.09)] + [(20, m, 1.03, 1.12) for m in (0.06, 0.065, 0.08)] + [(20, 0.098, 1.0, 25.7), (5, 0.098, 1.0, 25.7)]
    with mp.Pool(8) as P: res = P.map(one, J)
    for r in res: print(r)
    json.dump(res, open("/workspace/cosmology/papers/A_equals_fourier/steps_v3/dscan/dscan_rel_checks.json", "w"), indent=1, default=float)
