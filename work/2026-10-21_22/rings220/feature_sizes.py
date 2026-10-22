"""Physical size of ℓ≈220 stochastic-field features at observation-surface distance R.
Half-wavelength scale θ = π/ℓ; physical size s = R·θ. Compared with observed clump/void sizes."""
import numpy as np, json
ell=np.array([200,220,240]); R=np.array([4,8,14])  # Gpc
th=np.pi/ell
out={f"R={r}Gpc":{f"l={l}":{"theta_deg":round(np.degrees(t),3),"size_Mpc":round(r*1e3*t,1)} for l,t in zip(ell,th)} for r in R}
out["observed_reference_Mpc"]={"voids_typical":"10-100 (up to ~300 supervoids)","superclusters":"~100-200","BAO_scale":147}
json.dump(out,open("feature_sizes.json","w"),indent=1)
for k,v in out.items(): print(k,v)
