> ## ⚠️ DRAFT — work in progress: results and figures are preliminary and may change

# CoE Cosmology

Center-of-Expansion (CoE) shell cosmology: an off-centre observer in an anisotropic Hubble flow, with exact LTB ray tracing.

## The model
- **Big bang from a point:** every shell starts at the CoE at t = 0 and moves radially outward.
- **Central mass:** a mass M_c at the CoE slows the shells; its Schwarzschild radius may extend beyond the observer.
- **Observer offset:** we sit on a shell a distance D (10–60 Gpc) from the CoE, so the Hubble flow depends on direction (0° = away from the CoE, 180° = toward it).
- **Tuning:** launch speeds and age are set so the local Hubble rate today is H0 = 70 km/s/Mpc.
- **Structure:** thin decelerating shells wrinkle, seeding clumps, walls and voids.

## Layout (everything lives under `draft/`)
| Path | Contents |
|---|---|
| `draft/models/ltb_exact/` | **Current** exact GR (Lemaître–Tolman–Bondi) ray tracer `ltbx.py`, sweep driver `prod.py`, sweep results and figures |
| `draft/models/newtonian/` | Newtonian shell model, classical H(θ), A/A_all plots, older cosmological-redshift scripts |
| `draft/wrinkle/` | Thin-shell wrinkle instability scripts and figures |
| `draft/analyses/` | Clump/void/network structure studies, shell ensembles, figures, reports |
| `draft/papers/` | Pattern-formation paper v5 (PDF + LaTeX) and wrinkle PDFs |

## Running
Python 3 with `numpy`, `scipy`, `matplotlib`:
```
cd draft/models/ltb_exact && python prod.py       # exact LTB sweep
cd draft/models/newtonian && python classical14.py # classical H(θ)
```
Large intermediate data (`*.pkl`, `*.npz`, `*.npy`) is not included; the scripts regenerate it.

## Note on superseded code
The Newtonian/Doppler cosmological-redshift scripts (`cosmoz*`, `ourframe14*` in `draft/models/newtonian/`) are **superseded** by the exact LTB code (`draft/models/ltb_exact/ltbx.py`), which is validated against the classical limit (it matches classical H(θ) at short range and low speed, and gives H = 70 at every angle with no mass). Earlier results showing a low reading toward the CoE at 14.3 Gpc came from the older scripts and are not correct.

## License
MIT
