# Spherical small circles and Mollweide projection — geometry audit

## 1. Exact spherical locus

For Galactic longitude/latitude `(l,b)`,

\[
\mathbf n(l,b)=(\cos b\cos l,\;\cos b\sin l,\;\sin b).
\]

For center `(l0,b0)` with unit vector `n0` and angular radius
`r in [0,180 deg]`, the small circle is exactly

\[
\boxed{\mathbf n(l,b)\cdot\mathbf n_0=\cos r}.
\]

Choose a temporary vector `q` not parallel to `n0` and form an orthonormal
basis

\[
\mathbf e_1={\mathbf q\times\mathbf n_0\over\lVert\mathbf q\times\mathbf n_0\rVert},
\qquad
\mathbf e_2=\mathbf n_0\times\mathbf e_1.
\]

A stable parameterization is

\[
\boxed{\mathbf n(t)=\cos r\,\mathbf n_0+\sin r
(\cos t\,\mathbf e_1+\sin t\,\mathbf e_2)},\quad 0\le t<2\pi.
\]

It has unit norm and `n(t).n0=cos(r)` to roundoff. Reversing either basis
orientation only reverses traversal and does not change the locus.

## 2. Standard Mollweide

For display longitude `lambda in [-pi,pi]` and latitude `phi`, solve the
monotone implicit equation

\[
2\psi+\sin(2\psi)=\pi\sin\phi,
\qquad \psi\in[-\pi/2,\pi/2].
\]

Then

\[
\boxed{x={2\sqrt2\over\pi}\lambda\cos\psi,\qquad
 y=\sqrt2\sin\psi}.
\]

The Migkas Galactic convention is

\[
\lambda=-\operatorname{wrap}(l),
\]

so `l` increases left. Pixel normalization in the overlay is correct:
`Xn=x/(2 sqrt(2))`, `Yn=y/sqrt(2)`, then `px=cx+Rx*Xn`,
`py=cy-Ry*Yn`.

The implementation now uses a bracketed bisection for `psi`; maximum residual
of the implicit equation on a dense grid was `8.9e-16`. It was compared to
Matplotlib's trusted `MollweideAxes.MollweideTransform`: maximum coordinate
difference was `1.30e-3` map units (about `0.023%` of full map width), caused
by Matplotlib's intentionally loose Newton stopping threshold / polar Taylor
approximation; both use the same standard formula. Cardinal points and the
boundary agree.

## 3. Antimeridian behavior

A small circle is one closed curve on the sphere. A single Mollweide chart is
cut at `lambda=+/-pi`. If the curve crosses that cut, its chart image is
multiple polyline segments whose endpoints occur on both the left and right
oval boundaries. The two boundaries represent the same spherical meridian.
The pieces are closed only after those identified edges are glued; they are
not a closed Euclidean ellipse in the chart interior.

Therefore:

- split at longitude jumps before projection;
- merge the artificial parameter start/end if that join is not itself a seam;
- put cut endpoints on both `+pi` and `-pi` boundaries;
- never draw an interior chord, ellipse fit, or oval-rim bridge.

In the supplied infographic, case 5 `(l0,b0)=(180,0), r=20 deg` is wrong if
shown as a closed interior oval near `l=90 deg`: its center is on the seam, so
it must appear as matching arcs at both map edges. Case 8 `(200,-60), r=20`
also crosses the chosen seam and must split left/right. Case 1 `(0,+60),r=20`
does not reach the pole (`b_max=80 deg`) and is a closed, compressed interior
curve, not a genuinely open spherical arc.

## 4. Code audit

Correct before and after:

- `D_HAT=(1,0,0)` is anti-CoE `(0,0)`; CoE is `(-1,0,0)` at `(180,0)`.
- Cartesian-to-Galactic conversion is correct.
- The small-circle ONB in the old code was orthonormal; its orientation only
  changed point order.
- `lon_display=-wrap(l)` correctly makes Galactic longitude increase left.
- Mollweide x/y normalization to the detected oval was correct.

Actual bugs / misleading behavior found and corrected in
`make_migkas_fig8_coe_rings_overlay.py`:

1. **Fabricated 185-deg rim loop.** For `mu < -0.95`, the old code discarded
   the exact projection, snapped points to the Mollweide rim, and added north
   and south rim bridges. Those bridges are not the projection of the cone.
2. **Inward coordinate scaling.** Curves/markers were multiplied by an inset,
   changing their positions. They are now projected at exact coordinates and
   clipped by the registered oval path only for stroke width.
3. **Potential seam gaps / false third segment.** Linear splitting did not
   merge the cyclic parameterization start/end. The implementation now merges
   that artificial break and leaves only the real chart segments.
4. **Near-pole Newton instability.** The 12-step unbracketed auxiliary-angle
   solve could be inaccurate very near a pole. It is now a stable monotone
   bisection.

## 5. Meaning of the four scan labels

Spherical angular separation has range `[0,180 deg]`. If a 0–360 deg strip
scan angle is converted only through `mu=cos(theta)`, its directional / branch
information is discarded. The unique radius is
`r=arccos(cos(theta))`:

| label | mu | radius from +x | equivalent radius from CoE (-x) |
|---|---:|---:|---:|
| 30 deg | +0.866025 | 30 deg | 150 deg |
| 185 deg | -0.996195 | 175 deg | 5 deg |
| 115 deg | -0.422618 | 115 deg | 65 deg |
| 265 deg | -0.087156 | 95 deg | 85 deg |

Thus `mu=cos(185 deg)` is mathematically a valid cone value but **not** a
185-deg spherical distance: it is exactly the same cone as 175 deg from +x,
or a 5-deg circle about CoE. Likewise 265 deg reduces to 95 deg from +x (85
deg about CoE). If the original 0–360 strip scan distinguishes theta from
`360-theta`, an axisymmetric `A(mu)` ring model cannot preserve that
distinction; additional azimuth/branch information is required.
