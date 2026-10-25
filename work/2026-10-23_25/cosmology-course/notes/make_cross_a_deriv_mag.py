#!/usr/bin/env python3
"""Make the log-magnitude derivative plot from the recorded analytic formulas."""
import math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

OUT = Path('/workspace/cosmology-course/notes/cross-a-deriv-mag.png')
REPORT = Path('/workspace/cosmology-course/notes/cross-a-deriv-mag-report.txt')
T0 = 0.550
T_TODAY = 61.2
T_TODAY_S = 61.75
D0 = 2.5
T_KILL = 1.10
R_TODAY = 14.3
C = 0.3066
P = math.log(1000.0) / math.log(T_TODAY_S / T0)
T_CUT = 1.0979391162917367  # exact last plotted point recorded in cross-a-deriv.json


def quantities(t):
    t = np.asarray(t, dtype=float)
    ts = t + T0
    D = D0 * (ts / T0) ** (2.0 / 3.0)
    v = (2.0 / 3.0) * D / ts
    a = (ts / T_TODAY_S) ** P
    da = P * a / ts
    R = R_TODAY * t / T_TODAY
    dR = R_TODAY / T_TODAY
    out = {
        'd(Da)/dt': v * a + D * da,
        'd((D-R)a)/dt': (v - dR) * a + (D - R) * da,
        'd((D+R)a)/dt': (v + dR) * a + (D + R) * da,
        'd(Ra)/dt': dR * a + R * da,
        'd(v/c)/dt': (-v / (3.0 * ts)) / C,
        'da/dt': da,
    }
    u = 1.0 - t / T_KILL
    uu = np.clip(u, 1e-30, None)
    drec = D0 * (1.0 / 3.0) * uu ** (-2.0 / 3.0) * (-1.0 / T_KILL)
    Drec = D0 * uu ** (1.0 / 3.0)
    drec_a = drec * a + Drec * da
    drec_a = np.where(u > 0.0, drec_a, np.nan)
    out['d(D_rec a)/dt'] = drec_a
    return out


# Use the same chart window as cross-a-deriv.png: 10^-3 through 74 Gyr.
t_log = np.logspace(-3, math.log10(74.0), 3000)
t_rec = np.logspace(-3, math.log10(T_CUT), 3000)
s = quantities(t_log)
sr = quantities(t_rec)

plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 15,
    'axes.titlesize': 17,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
})
fig, ax = plt.subplots(figsize=(14.5, 7.2), layout='constrained')
# Plot positive magnitudes on a logarithmic y axis (equivalent to log10|derivative|).
styles = [
    ('d(Da)/dt', 'tab:red', '-', 2.2, r'$|d(Da)/dt|$ (Gpc/Gyr)'),
    ('d(D_rec a)/dt', 'tab:blue', '-', 2.2, r'$|d(D_{\mathrm{rec}}a)/dt|$ (Gpc/Gyr)'),
    ('d((D+R)a)/dt', 'tab:green', '--', 2.0, r'$|d((D+R)a)/dt|$ (Gpc/Gyr)'),
    ('d((D-R)a)/dt', 'tab:green', ':', 2.0, r'$|d((D-R)a)/dt|$ (Gpc/Gyr)'),
    ('d(Ra)/dt', 'tab:green', '-', 2.2, r'$|d(Ra)/dt|$ (Gpc/Gyr)'),
    ('d(v/c)/dt', 'tab:purple', '-', 2.2, r'$|d(v/c)/dt|$ (1/Gyr)'),
    ('da/dt', 'black', '-', 2.2, r'$|da/dt|$ (1/Gyr)'),
]
for key, color, ls, lw, label in styles:
    if key == 'd(D_rec a)/dt':
        ax.plot(t_rec, np.abs(sr[key]), color=color, lw=lw, ls=ls, label=label)
    else:
        ax.plot(t_log, np.abs(s[key]), color=color, lw=lw, ls=ls, label=label)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-3, 74.0)
ax.set_xlabel('time after the crossing (Gyr)')
ax.set_ylabel(r'magnitude $|d(\mathrm{quantity})/dt|$ (log$_{10}$ scale)')
ax.set_title('Magnitude of time derivatives')
ax.axvline(T_KILL, color='0.45', lw=1.0, zorder=0)
ax.axvline(T_TODAY, color='0.2', lw=1.0, ls='--', zorder=0)
ax.legend(loc='upper left', framealpha=0.95)
ax.annotate(
    r'$|d(D_{\mathrm{rec}}a)/dt|$ stops at the same finite point before $t=1.10$ Gyr; no curve is drawn after',
    xy=(0.99, 0.02), xycoords='axes fraction', ha='right', va='bottom', fontsize=14,
    color='tab:blue',
)
fig.savefig(OUT, dpi=140)
plt.close(fig)

# Requested velocity-magnitude samples, evaluated directly from the analytic formula.
vals = {}
for t in (0.1, 1.0, 5.0, 13.8, 61.2):
    vals[t] = float(abs(quantities(np.array([t]))['d(v/c)/dt'][0]))
REPORT.write_text(
    'cross-a derivative magnitude report\n'
    'PNG: /workspace/cosmology-course/notes/cross-a-deriv-mag.png\n'
    'The plot uses the exact analytic formulas recorded in cross-a-deriv.json, with absolute values on a logarithmic y axis (equivalent to log10 magnitudes).\n'
    'The recombination derivative is plotted only through t = 1.0979391162917367 Gyr, the recorded last finite point (~0.25 Gpc/Gyr magnitude), and not at/after 1.10 Gyr.\n\n'
    '|d(v/c)/dt| samples (1/Gyr):\n' +
    ''.join(f'  t = {t:g} Gyr: {v:.12g}\n' for t, v in vals.items())
)
print(f'wrote {OUT} ({OUT.stat().st_size} bytes)')
print(REPORT.read_text())
