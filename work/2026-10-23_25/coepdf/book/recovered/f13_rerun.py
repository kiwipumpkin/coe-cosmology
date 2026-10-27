# AUTO-WRAP generator for f13_rerun.png
"""f13_rerun.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: byte-identical to fig7/R14_inverse_ell220.png (script R14_inverse_ell220.py = Chris's f13 code)
Embedded source: /workspace/coepdf/fig7/R14_inverse_ell220.py  (verbatim below, unmodified)
Run:  python3 f13_rerun.py   -> writes f13_rerun_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/fig7 so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'f13_rerun' + '_regen.png')
SRC_CWD = '/workspace/coepdf/fig7'
MATCH = ['R14_inverse_ell220.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_f13_rerun_')
import matplotlib; matplotlib.use('Agg')
import matplotlib.figure as _mf
_orig_open = builtins.open
def _divert(p):
    p = os.fspath(p)
    if os.path.abspath(p) == OUT: return p
    return os.path.join(SCRATCH, os.path.basename(p))
def _open(file, mode='r', *a, **k):
    if isinstance(file, (str, bytes, os.PathLike)) and any(c in mode for c in 'wax+'):
        file = _divert(file)
    return _orig_open(file, mode, *a, **k)
builtins.open = _open
_orig_save = _mf.Figure.savefig
def _savefig(self, fname, *a, **k):
    if isinstance(fname, (str, os.PathLike)):
        b = os.path.basename(os.fspath(fname))
        fname = OUT if b in MATCH else os.path.join(SCRATCH, b)
    return _orig_save(self, fname, *a, **k)
_mf.Figure.savefig = _savefig
try:
    from PIL import Image as _I
    _ois = _I.Image.save
    def _isave(self, fp, *a, **k):
        if isinstance(fp, (str, os.PathLike)):
            b = os.path.basename(os.fspath(fp)); fp = OUT if (b in MATCH or os.path.abspath(os.fspath(fp)) == OUT) else os.path.join(SCRATCH, b)
        return _ois(self, fp, *a, **k)
    _I.Image.save = _isave
except Exception: pass
os.chdir(SRC_CWD); sys.path.insert(0, SRC_CWD)
sys.argv = ['/workspace/coepdf/fig7/R14_inverse_ell220.py'] + []

SOURCE = "# rev7: original inverse-design code, run unchanged. Dimensionless in R (shell_radius=1.012 R_ls), so it holds at R_ls=0.18957 Gpc without change.\n# Distributed spherical clumps and voids: wide/narrow ell=220 spectra\n# Original model code used to generate distributed_spherical_clumps_ell220.png\n# Requires numpy, scipy, matplotlib\n# Note: inverse-designed gravitational-potential spectrum, NOT CMB temperature.\n\nimport numpy as np\nimport matplotlib.pyplot as plt\nfrom scipy.special import eval_legendre\nfrom matplotlib.colors import TwoSlopeNorm\nfrom pathlib import Path\n\nrng = np.random.default_rng(20261021)\nell = np.arange(1, 351)\nshell_radius = 1.012\nn_axes = 5\nn_theta, n_phi = 130, 260\ntheta = np.linspace(.003, np.pi-.003, n_theta)\nphi = np.linspace(0, 2*np.pi, n_phi, endpoint=False)\ntt, pp = np.meshgrid(theta, phi, indexing='ij')\ndirections = np.stack([np.sin(tt)*np.cos(pp), np.sin(tt)*np.sin(pp), np.cos(tt)], axis=-1)\nxyz = directions.reshape(-1, 3)\n\ndef unit(a):\n    return a / np.linalg.norm(a, axis=-1, keepdims=True)\n\ndef build_case(width, seed):\n    rg = np.random.default_rng(seed)\n    target = np.exp(-.5*((ell-220)/width)**2)\n    # Draw distinct orientations for each degree; multiple m components.\n    density = np.zeros(len(xyz))\n    recovered = np.zeros(len(ell))\n    for i, l in enumerate(ell):\n        if target[i] < 1e-9:\n            continue\n        axes = unit(rg.normal(size=(n_axes, 3)))\n        weights = rg.choice([-1., 1.], size=n_axes)\n        gram = eval_legendre(int(l), np.clip(axes @ axes.T, -1, 1))\n        S = max(float(weights @ gram @ weights), 1e-12)\n        # Choose each harmonic amplitude to enforce the prescribed target.\n        A = shell_radius**(l+1) * np.sqrt(target[i]*(2*l+1)/(l*(l+1)*S))\n        basis = sum(weights[j]*eval_legendre(int(l), xyz @ axes[j])\n                    for j in range(n_axes))\n        density += A*basis\n        recovered[i] = l*(l+1)/(2*l+1)*(A/shell_radius**(l+1))**2*S\n    density = density.reshape(n_theta, n_phi)\n    return target, recovered/recovered.max(), density\n\ncases = [(32, build_case(32, 32)), (8, build_case(8, 8))]\nfig = plt.figure(figsize=(14, 11), facecolor='white', layout='constrained')\ngs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15])\nfor row, (width, (target, power, density)) in enumerate(cases):\n    ax = fig.add_subplot(gs[row, 0], projection='3d')\n    # Display clipping only: actual density values are not changed.\n    clip = np.percentile(np.abs(density), 97)\n    norm = TwoSlopeNorm(vmin=-clip, vcenter=0, vmax=clip)\n    colors = plt.get_cmap('RdBu_r')(norm(density))\n    X = np.sin(tt)*np.cos(pp)\n    Y = np.sin(tt)*np.sin(pp)\n    Z = np.cos(tt)\n    ax.plot_surface(X, Y, Z, facecolors=colors, shade=False,\n                    linewidth=0, rstride=1, cstride=1, antialiased=False)\n    ax.set_box_aspect((1, 1, 1))\n    ax.set(xlim=(-1, 1), ylim=(-1, 1), zlim=(-1, 1))\n    ax.view_init(elev=24, azim=35)\n    ax.set_axis_off()\n    ax.set_title(f'Distributed 3D clumps and voids\\n'\n                 f'sigma_l={width}; 97th-percentile color clip', fontsize=12)\n\n    ax2 = fig.add_subplot(gs[row, 1])\n    ax2.plot(ell, target, 'k', lw=2, label='Prescribed target')\n    ax2.plot(ell, power, color='#d44b40', ls='--', lw=1.7,\n             label='Calculated spherical power')\n    ax2.axvline(220, color='.5', ls=':')\n    ax2.set(xlim=(110, 330), ylim=(0, 1.08), xlabel='Spherical multipole l',\n            ylabel='Normalized potential D_l', title=f'Peak l=220, sigma_l={width}')\n    ax2.legend(fontsize=9)\n    ax2.grid(alpha=.2)\n\nfig.suptitle('Non-axisymmetric source shell: density varies across latitude AND longitude\\n'\n             'Exterior shell r=1.012R | gravitational potential, not CMB temperature',\n             fontsize=14)\nfig.savefig('R14_inverse_ell220.png', dpi=170, bbox_inches='tight')\nplt.show()\n"
exec(compile(SOURCE, '/workspace/coepdf/fig7/R14_inverse_ell220.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/fig7/R14_inverse_ell220.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
