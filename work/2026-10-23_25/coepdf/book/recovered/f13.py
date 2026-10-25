# Distributed spherical clumps and voids: wide/narrow ell=220 spectra
# Original model code used to generate distributed_spherical_clumps_ell220.png
# Requires numpy, scipy, matplotlib
# Note: inverse-designed gravitational-potential spectrum, NOT CMB temperature.

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import eval_legendre
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

rng = np.random.default_rng(20261021)
ell = np.arange(1, 351)
shell_radius = 1.012
n_axes = 5
n_theta, n_phi = 130, 260
theta = np.linspace(.003, np.pi-.003, n_theta)
phi = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
tt, pp = np.meshgrid(theta, phi, indexing='ij')
directions = np.stack([np.sin(tt)*np.cos(pp), np.sin(tt)*np.sin(pp), np.cos(tt)], axis=-1)
xyz = directions.reshape(-1, 3)

def unit(a):
    return a / np.linalg.norm(a, axis=-1, keepdims=True)

def build_case(width, seed):
    rg = np.random.default_rng(seed)
    target = np.exp(-.5*((ell-220)/width)**2)
    # Draw distinct orientations for each degree; multiple m components.
    density = np.zeros(len(xyz))
    recovered = np.zeros(len(ell))
    for i, l in enumerate(ell):
        if target[i] < 1e-9:
            continue
        axes = unit(rg.normal(size=(n_axes, 3)))
        weights = rg.choice([-1., 1.], size=n_axes)
        gram = eval_legendre(int(l), np.clip(axes @ axes.T, -1, 1))
        S = max(float(weights @ gram @ weights), 1e-12)
        # Choose each harmonic amplitude to enforce the prescribed target.
        A = shell_radius**(l+1) * np.sqrt(target[i]*(2*l+1)/(l*(l+1)*S))
        basis = sum(weights[j]*eval_legendre(int(l), xyz @ axes[j])
                    for j in range(n_axes))
        density += A*basis
        recovered[i] = l*(l+1)/(2*l+1)*(A/shell_radius**(l+1))**2*S
    density = density.reshape(n_theta, n_phi)
    return target, recovered/recovered.max(), density

cases = [(32, build_case(32, 32)), (8, build_case(8, 8))]
fig = plt.figure(figsize=(14, 11), facecolor='white', layout='constrained')
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15])
for row, (width, (target, power, density)) in enumerate(cases):
    ax = fig.add_subplot(gs[row, 0], projection='3d')
    # Display clipping only: actual density values are not changed.
    clip = np.percentile(np.abs(density), 97)
    norm = TwoSlopeNorm(vmin=-clip, vcenter=0, vmax=clip)
    colors = plt.get_cmap('RdBu_r')(norm(density))
    X = np.sin(tt)*np.cos(pp)
    Y = np.sin(tt)*np.sin(pp)
    Z = np.cos(tt)
    ax.plot_surface(X, Y, Z, facecolors=colors, shade=False,
                    linewidth=0, rstride=1, cstride=1, antialiased=False)
    ax.set_box_aspect((1, 1, 1))
    ax.set(xlim=(-1, 1), ylim=(-1, 1), zlim=(-1, 1))
    ax.view_init(elev=24, azim=35)
    ax.set_axis_off()
    ax.set_title(f'Distributed 3D clumps and voids\n'
                 f'sigma_l={width}; 97th-percentile color clip', fontsize=12)

    ax2 = fig.add_subplot(gs[row, 1])
    ax2.plot(ell, target, 'k', lw=2, label='Prescribed target')
    ax2.plot(ell, power, color='#d44b40', ls='--', lw=1.7,
             label='Calculated spherical power')
    ax2.axvline(220, color='.5', ls=':')
    ax2.set(xlim=(110, 330), ylim=(0, 1.08), xlabel='Spherical multipole l',
            ylabel='Normalized potential D_l', title=f'Peak l=220, sigma_l={width}')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=.2)

fig.suptitle('Non-axisymmetric source shell: density varies across latitude AND longitude\n'
             'Exterior shell r=1.012R | gravitational potential, not CMB temperature',
             fontsize=14)
fig.savefig('distributed_spherical_clumps_ell220.png', dpi=170, bbox_inches='tight')
plt.show()
