import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

out_dir = '/workspace/cosmology-course/notes'

# Exact formulas requested.
t = np.logspace(-3, np.log10(80.0), 1600)
t_sing = t + 0.550
D = 2.5 * (t_sing / 0.550) ** (2/3)
c = 0.3066
p = np.log(1000) / np.log(61.75 / 0.550)
a = (t_sing / 61.75) ** p
v_without_a = ((2/3) * D / t_sing) / c
v_scale_factor = (a * p * D / t_sing) / c
v_with_a = v_without_a + v_scale_factor

plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 22,
    'axes.labelsize': 17,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
})

def style_log_axes(ax):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-3, 80)
    ax.set_ylim(bottom=1e-2)
    ax.grid(True, which='major', color='#b0b0b0', alpha=0.35)
    ax.grid(True, which='minor', color='#d0d0d0', alpha=0.25)
    ax.tick_params(which='both', direction='out')
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())

# Keep a generous top margin: title first, then the centered one-row legend.
def finish(fig, ax, title, legend_handles, legend_labels, filename, ncol):
    ax.set_title(title, pad=62)
    ax.legend(
        legend_handles, legend_labels,
        loc='lower center', bbox_to_anchor=(0.5, 1.13),
        ncol=ncol, frameon=True, borderaxespad=0.0,
        handlelength=2.4, columnspacing=1.2,
    )
    # tight_layout reserves the top region; bbox_inches='tight' retains the
    # title and legend in the output rather than clipping them.
    fig.tight_layout(rect=(0.03, 0.06, 0.99, 0.70))
    fig.savefig(f'{out_dir}/{filename}', dpi=150, bbox_inches='tight')
    plt.close(fig)

# 1) Velocity with and without scale factor.
fig, ax = plt.subplots(figsize=(10.24, 6.60), dpi=150)
style_log_axes(ax)
line_without, = ax.plot(t, v_without_a, lw=2.5, color='#1f77b4', label='without scale factor (through space only)')
line_scale, = ax.plot(t, v_scale_factor, lw=2.5, color='#ff7f0e', label='scale-factor piece only')
line_with, = ax.plot(t, v_with_a, lw=2.5, color='#2ca02c', label='with scale factor (sum)')
line_c = ax.axhline(1.0, color='black', lw=1.8, ls='--', label='c')
for xv in (13.8, 61.2):
    ax.axvline(xv, color='#666666', lw=1.6, ls=':')
ax.set_ylim(1e-2, float(np.max(v_with_a) * 1.15))
ax.set_xlabel('Time after crossing, t (Gyr)')
ax.set_ylabel('Speed (c)')
finish(
    fig, ax, 'Velocity with and without the scale factor',
    [line_without, line_scale, line_with, line_c],
    ['without scale factor', 'scale-factor piece', 'with scale factor', 'c'],
    'v-with-without-a-logy-1e-2.png', 4,
)

# 2) Scale factor only.
fig, ax = plt.subplots(figsize=(10.24, 6.60), dpi=150)
style_log_axes(ax)
line_a, = ax.plot(t, a, lw=2.5, color='#9467bd', label='a(t)')
for xv in (13.8, 61.2):
    ax.axvline(xv, color='#666666', lw=1.6, ls=':')
ax.set_ylim(1e-2, float(np.max(a) * 1.15))
ax.set_xlabel('Time after crossing, t (Gyr)')
ax.set_ylabel('Scale factor, a(t)')
finish(fig, ax, 'Scale factor', [line_a], ['a(t)'], 'scale-factor-logy-1e-2.png', 1)
