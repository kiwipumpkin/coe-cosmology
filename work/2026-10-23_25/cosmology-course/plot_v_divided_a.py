import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Exact model specified
t = np.logspace(-3, np.log10(80.0), 2000)  # Gyr after crossing
t_sing = t + 0.550
D = 2.5 * (t_sing / 0.550)**(2/3)  # Gpc
c = 0.3066  # Gpc/Gyr
p = np.log(1000.0) / np.log(61.75 / 0.550)
a = (t_sing / 61.75)**p
old_shell_speed = ((2/3) * D / t_sing) / c
old_speed_divided_by_a = old_shell_speed / a

plt.rcParams.update({"font.size": 14})
fig, ax = plt.subplots(figsize=(12, 7.5), constrained_layout=True)
ax.plot(t, old_speed_divided_by_a, color="#c23b22", lw=2.8, label="old speed divided by a")
ax.plot(t, old_shell_speed, color="#2166ac", lw=2.8, label="old shell speed")
ax.axhline(1, color="black", lw=1.8, ls="--", label="c")
ax.axvline(13.8, color="#555555", lw=1.5, ls=":")
ax.axvline(61.2, color="#555555", lw=1.5, ls=":")
ax.text(13.8, 0.97, "13.8 Gyr", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=14)
ax.text(61.2, 0.97, "61.2 Gyr", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=14)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-3, 80)
ax.set_xlabel("Time after crossing, t (Gyr)")
ax.set_ylabel("Speed / c")
ax.set_title("Shell speed divided by a(t)", fontsize=18)
ax.grid(True, which="both", alpha=0.25)
ax.legend(fontsize=14, loc="best")
path = Path("/workspace/cosmology-course/notes/v-divided-by-a.png")
path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(path, dpi=180)
plt.close(fig)

# Requested values, including t=0 by direct evaluation of the exact formula.
def values(t0):
    ts = t0 + 0.550
    d = 2.5 * (ts / 0.550)**(2/3)
    aa = (ts / 61.75)**p
    old = ((2/3) * d / ts) / c
    return old / aa, old, aa
print(f"p={p:.15g}")
for t0 in (0.0, 13.8, 61.2):
    div, old, aa = values(t0)
    print(f"t={t0:g} Gyr: divided={div:.12g}, old={old:.12g}, a={aa:.12g}")
print(path)
