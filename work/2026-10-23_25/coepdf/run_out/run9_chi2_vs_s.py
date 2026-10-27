"""run9_chi2_vs_s.png -- REMADE generating script (original plotting code not found on box, git, or GitHub).
Data: run9/log9.txt, written by run9/run9.py (line 1 header dict; lines 'red_surface {...}' and 'R0.6 {...}'
with chi^2 vs Planck TT (l=2..30) for each shear scale s). Plot = chi^2 vs s, log-log, as in the original.
Run: python3 run9_chi2_vs_s.py -> run9_chi2_vs_s_regen.png (original untouched)."""
import os, ast
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
L = open(os.path.join(H, 'run9/log9.txt')).read().splitlines()
hdr = ast.literal_eval(L[0]); curves = {}
for ln in L[1:]:
    tag, d = ln.split(' ', 1)
    if tag in ('red_surface', 'R0.6'):
        d = ast.literal_eval(d); pts = sorted((float(k.rsplit('_', 1)[1]), v) for k, v in d.items() if v is not None)
        curves[tag] = pts
fig = plt.figure(figsize=(6.5, 13/3))
for tag in ('red_surface', 'R0.6'):
    s, c = zip(*curves[tag]); plt.loglog(s, c, 'o-', label=tag)
plt.xlabel('shear scale s'); plt.ylabel('chi^2 vs Planck'); plt.legend()
plt.title(f"run9: n_bodies={hdr['n_bodies']}, rms/floor={hdr['pattern_rms_over_floor']:.2f}")
plt.tight_layout(); plt.savefig(os.path.join(H, 'run9_chi2_vs_s_regen.png'), dpi=120)
