# Progress notes, 2026-10-19: from "silly" to positive

## 1. Low point: tilting the sky did not work
We tried to line up the model's CMB peaks with Planck by pinching or tilting the sky, spreading the curve from -10 to +10 degrees (the `rotfit/` and `rotfit/free/spread/` runs). It failed. The model's cold peak sits about 90 degrees away from Planck's cold spot, so no small tilt can line them up.
Chris: "That didn't work very well huh."

## 2. Turnaround: direction is fine, amplitude is the gap
Chris reframed the problem. The peaks are in the right place on the sky and have the right sign, so the direction works and only the amplitude is off.
Chris: "The peaks are just fine in the sky... If they were on the opposite side it would be a different story."

## 3. Positive turn: one coherent geometry
A single geometry does three things at once. It hits the Migkas cluster anisotropy, fits the CatWISE density gradient, and reproduces the hot and cold sides of the CMB.
Chris: "Right now we've got a geometry that can hit Migkas. It can fit CatWISE density gradient. And it can get the hot and cold side of the CMB. We just need to keep developing the model."

## 4. Second positive: the dipole amplitude
The ~89x dipole comes from the spread in emission times (sequential recombination). If recombination is close to simultaneous, the dipole moves toward the Planck value. A recombination radius that shrinks quickly with a critical-density threshold does this without any new knobs (~1.72x; see `dipole_amp/` and `shrink_rec/`). Chris was happy with this. He said a fixed recombination radius never made sense and that it should shrink. He first had the sign backwards, said so, and corrected it: simultaneous recombination lowers the CMB dipole, and sequential recombination gives the large dipole.

## 5. Third positive: accepting the geometric floor
The 1.73x geometric floor is the minimum set by mass, because mass makes the expansion uneven (see `shrink_rec/floor/`). We accept it as the honest minimum.
Chris: "That's perfectly acceptable."

## Still open
The non-dipole structure (m = 3,4,5) is still too large under our density-at-recombination rule. See `clusters_voids/`.
