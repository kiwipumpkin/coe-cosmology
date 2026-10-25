"""Accuracy check of the no-boxcar extremum formulas against the exact
no-boxcar curve built from the unchanged model functions. No fitting."""
import sys, importlib.util, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "/workspace/cosmology/density_local_match"); sys.path.insert(0, "/workspace/cosmology")
from run_density_local_match import params_twoscale_local_match, peculiar_u, m_enc_at_r, D, H, G
spec = importlib.util.spec_from_file_location("hm", "/workspace/figs_letter_exact_python/makers/make_amp_phi_and_geometry_heatmaps.py")
hm = importlib.util.module_from_spec(spec); spec.loader.exec_module(hm)
K0 = 0.23058828928336284; FU = 0.94

def h_exact(th, x, eps, alpha, kappa):
    p = hm.scale_kappa(params_twoscale_local_match(x, FU), kappa)
    n = np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], -1)
    d = eps*D; Dv = np.array([D, 0, 0.]); rv = Dv + d*n
    ur = peculiar_u(rv, p, alpha=alpha); uD = peculiar_u(Dv.reshape(1, 3), p, alpha=alpha)[0]
    return np.maximum(H + np.sum((ur-uD)*n, -1)/d, 1e-6*H), p  # same clamp as af_an_alpha

def exact_ext(x, eps, alpha, kappa, amp=0.0, ph=0.35):
    th = np.radians(np.arange(0, 360, 0.01))
    h, p = h_exact(th, x, eps, alpha, kappa); hbar = h.mean()
    A = (1 - amp*np.sin(th+ph))*h/hbar
    out = []
    n = len(A)
    for i in range(n):
        a, b, c = A[i-1], A[i], A[(i+1) % n]
        if b > a and b >= c: out.append(("max", np.degrees(th[i]), b))
        if b < a and b <= c: out.append(("min", np.degrees(th[i]), b))
    return out, hbar, p

def formula(x, eps, alpha, kappa, L=None):
    # physical form, D=H=G=1 units; M = kappa*D^3*H^2/G
    Dd = 1.0; L = Dd/x; d = eps*Dd; M = kappa; aGM = alpha*G*M; H2 = H**2
    S = (1-FU)*(6*L**3*np.exp(Dd/L) - 6*L**3 - 6*Dd*L**2 - 3*Dd**2*L) + FU*Dd**3
    den = (2*Dd**3*H2 + aGM)*S - 3*aGM*Dd**3
    Aall = den/(2*Dd**3*H2*S)
    def amax(sg):  # sg=-1 at 0 deg, +1 at 180 deg
        return ((2*Dd**3*H2 + 4*aGM + sg*6*aGM*d/Dd)*S - 6*aGM*Dd**3
                - sg*3*aGM*d*(2*Dd**2 + (1-FU)*Dd**3/L))/den
    amin = 2*(Dd**3*H2 - aGM)*S/den
    lmin = 90 + np.degrees(3*d/(4*Dd))
    return dict(l=(0.0, 180.0, lmin, 360-lmin), A=(amax(-1), amax(+1), amin, amin), Aall=Aall)

if __name__ == "__main__":
    # sanity: code A_raw vs closed h form
    x, e, a, k = 5.9, 0.1, 1.0, K0
    ex, hbar, p = exact_ext(x, e, a, k)
    print("kappa from p", p["kappa"], "rho_local/mean", p["rho_local_over_mean"])
    f = formula(x, e, a, k)
    print("exact no-box:", [(t, round(l, 2), round(v, 4)) for t, l, v in ex], "hbar", hbar)
    print("formula:", [round(v, 2) for v in f["l"]], [round(v, 4) for v in f["A"]], "Aall", f["Aall"])
    ext_t, _, _ = exact_ext(x, e, a, k, amp=-0.06)
    print("exact no-box WITH tilt:", [(t, round(l, 2), round(v, 4)) for t, l, v in ext_t])
