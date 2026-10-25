# Cycle 6: Celestial sphere between two emitters

Chris's ask: observer between **two parallel plates** OR **two concentric shells**, both emitting uniform light at constant rest-frame temperature \(T\) — what does the sky look like?

## Setup

- Specific intensity along each ray = blackbody \(B_\nu(T_{\rm emitter})\) of the surface the ray hits (vacuum between emitters; no absorption). Brightness temperature \(T_{\rm obs}=T_{\rm emitter}\) in the rest frame.
- Lambertian / blackbody: \(I_\nu\) independent of viewing angle on the surface → foreshortening cancels; no extra \(\cos\theta\) on \(T_{\rm obs}\).
- **Plates:** infinite at \(z=\pm d\); observer at \(z\in(-d,d)\). Every upward ray hits \(+\)plate, every downward ray hits \(-\)plate (hemispheres fixed — offset never changes solid angles).
- **Shells:** \(R_{\rm in}=1\), \(R_{\rm out}=3\), observer at \(D=1.8\) (off-centre). Inner shell subtends angular radius \(\alpha=\arcsin(R_{\rm in}/D)\approx 33.7^\circ\). Rays with \(|\theta|<\alpha\) (toward CoE) hit the inner shell; all others hit the outer shell.
- Optional second step: Doppler from observer motion (plates) or Hubble-like \(V_{\rm rec}\propto r\) from CoE (shells).

## Key bullets

### Rest-frame equal \(T\)
- **Parallel plates (mid or offset):** sky is **perfectly isotropic** at \(T\). Offset does nothing — both plates still fill exact hemispheres at the same \(T\).
- **Concentric shells (any off-centre \(D\)):** sky is **perfectly isotropic** at \(T\). Geometry only decides *which* shell a ray hits; same \(T\) ⇒ same \(T_{\rm obs}\) everywhere.
- Scanning the sky: flat \(T_{\rm obs}(\theta)=T\) on the full \(0^\circ\)–\(360^\circ\) great circle; Mollweide is a blank uniform map.

### When anisotropy appears
1. **Different \(T\) on the two emitters**
   - Plates: hard **hemisphere step** at the horizon (\(\theta=90^\circ/270^\circ\)).
   - Shells: hard **disk / annulus step** at the inner-shell limb \(\alpha\).
2. **Relative motion**
   - Plates + observer \(\beta\): classic smooth **dipole** on an otherwise flat sky.
   - Shells + recession Doppler \(V_{\rm rec}\): smooth angle dependence from line-of-sight velocity; thin-shell model also jumps at \(\alpha\) where hit radius (hence \(\beta\)) discontinuous. Comoving vs CoE-rest observers differ in dipole orientation/amplitude.
3. **Offset alone does *not* create anisotropy** when both emitters share the same rest-frame \(T\).

### Pedagogical takeaway
Uniform rest-frame \(T\) filling every line of sight ⇒ **isotropic celestial sphere**, independent of plate offset or shell decentring. Geometry without a \(T\) contrast or velocity field is invisible on the sky. Anisotropy requires unequal emitter temperatures or Doppler (\(V_{\rm obs}\) / \(V_{\rm rec}\)).

## Files
- `compute_two_emitters.py` → `two_emitters.json`
- `plot_two_emitters.py` →
  - `plates_T_vs_theta.png` — plates, full great circle
  - `shells_T_vs_theta.png` — shells, full great circle
  - `plates_mollweide.png` — plates sky maps
  - `shells_mollweide.png` — shells sky maps

## Caveats
- Thin shells / infinite plates idealization; continuum last-scattering surface in FLRW is isotropic for comoving observers even with expansion (toy \(V_{\rm rec}\) here is discontinuous in \(r\)).
- Non-relativistic Doppler for shells (\(\beta_{\rm rec}=0.02\)); plates use exact special-relativistic boost.
- No gravitational redshift / SW / ISW in this cycle — pure emitter geometry + optional kinematic Doppler.
