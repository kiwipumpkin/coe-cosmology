# Distance to recombination vs time and offset from the CoE

Guided by Chris Alameddin, *Centre-of-Expansion shell cosmology: from Migkas anisotropy to a flat CMB* (draft preprint, 1 October 2026, not for submission). The draft is a **guide, not a rigid template**. Where it is unfinished, this note says so. No cheerleading of Lambda cold dark matter (ΛCDM) and no cheerleading of the center-of-expansion (CoE) sketch.

**Working picture as of 3 October 2026, from Chris, overrides the draft's fixed recombination radius.** Mass leaving a surface drops the pressure, the temperature falls through 3000 K, and the recombination distance shrinks. The draft's "fixed physical radius" sentence stays below as what the draft wrote, not as the equation to use.

## Equations: our position and \(D_{\mathrm{recomb}}(t)\)

Clock starts at \(t = 0\) (the earlier "\(T = 0\)" is this start, not zero temperature). Temperature is \(T\). We start a small distance \(D_0\) from the center of expansion (CoE) with outward speed \(v_0\). Gravity pulls toward the CoE, so the motion is always decelerating. Deceleration is not the same as turning around.

**Our distance from the CoE**

\[
\ddot D = -\frac{G M(<D, t)}{D^2}, \qquad D(0) = D_0, \qquad \dot D(0) = v_0.
\]

\(M(<D, t)\) is the mass still interior to us. While that mass is constant,

\[
\dot D^2 = v_0^2 - 2GM\left(\frac{1}{D_0} - \frac{1}{D}\right).
\]

Write \(A = v_0^2 - 2GM/D_0\) and \(B = 2GM\). Then \(\dot D = \sqrt{A + B/D}\), and the time to reach distance \(D\) is the integral of \(ds/\sqrt{A + B/s}\).

If \(A > 0\) we keep going outward forever:

\[
t = \left[\frac{\sqrt{s(As+B)}}{A} - \frac{B}{A^{3/2}}\ln\left(\sqrt{As} + \sqrt{As+B}\right)\right]_{D_0}^{D}.
\]

If \(A = 0\) (exactly the speed that only dies at infinite distance):

\[
D(t) = \left[D_0^{3/2} + \frac{3}{2}\sqrt{2GM}\, t\right]^{2/3}.
\]

If \(A < 0\) we turn around once the speed hits zero:

\[
t = \left[-\frac{\sqrt{s(B-|A|s)}}{|A|} + \frac{B}{|A|^{3/2}}\arcsin\sqrt{\frac{|A|s}{B}}\right]_{D_0}^{D}.
\]

These three antiderivatives were checked by differentiating them back to \(1/\dot D\). They are Newtonian, constant interior mass. They are not the LTB shell equation in the draft, and they are not a fit to Table 1.

**Recombination distance**

\(D_{\mathrm{recomb}}(t)\) is the radius at which the temperature crosses 3000 K. Inside that radius the gas is still hotter. Mass leaving the surface drops the pressure, the temperature falls through 3000 K, and the radius moves in:

\[
T\big(D_{\mathrm{recomb}}(t), t\big) = 3000\,\mathrm{K}, \qquad \dot D_{\mathrm{recomb}} < 0
\]

while mass is still leaving. Let \(M_{\mathrm{rec}}(t)\) be the mass still inside that surface, so mass leaving means \(\dot M_{\mathrm{rec}} < 0\). Holding the surface on 3000 K,

\[
\dot D_{\mathrm{recomb}} = -\frac{\partial T/\partial M}{\partial T/\partial r}\, \dot M_{\mathrm{rec}}.
\]

Temperature rises with the enclosed mass and falls with radius, so the fraction is positive and a loss of mass shrinks \(D_{\mathrm{recomb}}\). Do not set this radius equal to the 14.3 Gpc vision radius.

**Ratio, constant-thickness shell, Migkas mass (later the same day).** \(D_{\mathrm{rec}}\) is a named parameter until a fit pins it. Position is the decelerating \(D(t)\) above. Crossing is \(D(t_\star)=D_{\mathrm{rec}}(t_\star)\). After that we see only \(c(t-t_\star)\).

\[
\frac{D_{\mathrm{rec}}(t)}{D_{\mathrm{rec}}(t_0)}=\frac{M_{\mathrm{left}}(t)}{M_{\mathrm{left}}(t_0)}.
\]

\(M(<D)\) comes from the Migkas \(A/A_{\mathrm{all}}\) versus angle fit. That is the mass inside our distance, not automatically the mass still inside \(D_{\mathrm{rec}}\).

The mass still inside the recombination surface is a shell of constant thickness \(\Delta\), and \(\Delta\) is seen. Sample density \(\rho\):

The thin-shell volume \(4\pi D_{\mathrm{rec}}^2\Delta\) is the wrong truncation. The shell runs from \(D_{\mathrm{rec}}\) to \(D_{\mathrm{rec}}+\Delta\) at constant thickness \(\Delta\):

\[
V=\frac{4}{3}\pi\left[(D_{\mathrm{rec}}+\Delta)^3-D_{\mathrm{rec}}^3\right]
=\frac{4}{3}\pi\left(3D_{\mathrm{rec}}^2\Delta+3D_{\mathrm{rec}}\Delta^2+\Delta^3\right).
\]

\(M_{\mathrm{left}}=\rho V\) and \(D_{\mathrm{rec}}=C M_{\mathrm{left}}\) with \(C=\mu m_u G/(k_B\cdot 3000\,\mathrm{K})\) give the quadratic

\[
3\Delta\, D_{\mathrm{rec}}^2+\left(3\Delta^2-\frac{3}{4\pi C\rho}\right)D_{\mathrm{rec}}+\Delta^3=0.
\]

Both roots are positive when the linear coefficient is negative. Keep the root that sits a thickness \(\Delta\) inside the mass that is seen. Density enters only through that linear coefficient: higher \(\rho\) gives a smaller \(D_{\mathrm{rec}}\). The Migkas \(A/A_{\mathrm{all}}\) mass inside \(D\) is what later has to match \(\rho V\). Not fit here.

**Ideal gas closure (3 October 2026).** Chris asked to close this with the ideal gas law. Mean particle mass \(\mu m_u\), with \(m_u\) the hydrogen mass:

\[
P = \frac{\rho}{\mu m_u} k_B T.
\]

On the recombination surface \(T = 3000\,\mathrm{K}\), so \(P/\rho = k_B (3000\,\mathrm{K})/(\mu m_u)\). The pressure that supports the gas against gravity scales as the weight of the mass still inside, \(P \sim G M_{\mathrm{rec}}\rho / D_{\mathrm{recomb}}\). Then

\[
D_{\mathrm{recomb}}(t) = \frac{\mu m_u G\, M_{\mathrm{rec}}(t)}{k_B \cdot 3000\,\mathrm{K}},
\qquad
\frac{\dot D_{\mathrm{recomb}}}{D_{\mathrm{recomb}}} = \frac{\dot M_{\mathrm{rec}}}{M_{\mathrm{rec}}}.
\]

The factor \(1/2\) from the central pressure of a uniform ball is the wrong object. Chris rejected it. Also rejected: a vision sphere reaching a distance \(R\) on either side of the recombination line. We pass the line. We do not straddle it.

After the crossing \(D(t_\star)=D_{\mathrm{recomb}}(t_\star)\):

- If \(\dot D > c\), we never see the recombination line again.
- If \(\dot D < c\), we see the place in space where the line was, and we keep seeing it only while that radius does not fall faster than light, \(|\dot D_{\mathrm{recomb}}| < c\). If the 3000 K radius shrinks faster than \(c\), we lose it.

The ideal-gas ratio \(P/\rho = k_B(3000\,\mathrm{K})/(\mu m_u)\) is unchanged. \(\mu\) stays a symbol. The shell mass of thickness \(2R\) is not the geometry. \(R=ct\) is not a radius measured on both sides of the line.

**Vision radius**

Let \(t_\star\) be the crossing, \(D(t_\star) = D_{\mathrm{recomb}}(t_\star)\). The vision clock starts there. Before that, nothing. After that, the only distance we can see is \(c\,dt\) with \(dt = t - t_\star\):

\[
R_v = 0 \quad (t \le t_\star), \qquad R_v = c\,dt = c(t - t_\star) \quad (t > t_\star).
\]

Not \(ct\) measured from the start of the motion.

Once we are outside, the gap back to the 3000 K surface is \(\ell(t) = D(t) - D_{\mathrm{recomb}}(t)\). We see that surface again when \(c(t-t_\star) = \ell(t)\). The gap opens at

\[
\dot\ell = \dot D - \dot D_{\mathrm{recomb}} = \dot D + |\dot D_{\mathrm{recomb}}|.
\]

Two scenarios, from Chris, not from a gap that straddles the line:

- \(\dot D > c\): we never see the recombination line again.
- \(\dot D < c\): we see the place where the line was, and we keep seeing it only while \(|\dot D_{\mathrm{recomb}}| < c\). If the radius falls faster than light, we lose it.

Coordinate speed above \(c\) is allowed for us, as in the draft.

---


### Acronyms

| Acronym | Meaning |
|---------|---------|
| **CoE** | Center of expansion: the preferred center in the draft's spherical shell picture. |
| **CMB** | Cosmic microwave background. |
| **FLRW** | Friedmann–Lemaître–Robertson–Walker: homogeneous and isotropic metric family. |
| **LTB** | Lemaître–Tolman–Bondi: spherically symmetric dust cosmology, written in the draft about the CoE, not about us. |

**Do not mix these radii.** The draft uses several lengths that an earlier course note collapsed into one \(R\).

| Symbol in the draft | What it is |
|---------------------|------------|
| \(D\) | Distance from the CoE to the observer. Fitted values are tens of Gpc, not a small offset inside a 14 Gpc last-scattering shell. |
| \(R\) | **Vision** radius used for today's-field diagnostics, \(R \approx 14.3\,\mathrm{Gpc}\). Not the recombination radius. |
| \(R_{\mathrm{in}}, R_{\mathrm{out}}\) | Inner and outer edges of the thick ordinary-matter shell. Observer is meant to sit mid-shell. Shell rule in the draft: \(D \gtrsim 14.7\,\mathrm{Gpc}\) so the observer is not outside the shell and not in an empty cavity that would spoil the fitted anisotropy. One Newtonian row sets \(R_{\mathrm{in}}=D-R\), \(R_{\mathrm{out}}=D+R\). |
| Recombination surface | A spherical surface of **fixed physical radius** around the CoE. Call that radius \(R_{\mathrm{rec}}\). The draft does not quote a number for it. |
| \(D_{\mathrm{recomb}}\) | The CoE-distance at which "vision turns on." Hubble-bubble volume about the observer is zero until the observer passes this distance. |

---

## Plain-language answer

In the draft, recombination is not a comoving shell that grows with the scale factor. It is a sphere of fixed physical radius around the CoE. The thick shells expand outward. The sphere of events we can see grows at the speed of light, and only after the observer's distance from the CoE passes \(D_{\mathrm{recomb}}\).

So the distance from us to that recombination surface changes with time mainly because **we move** (the shell we ride expands, so \(D(t)\) changes) while **\(R_{\mathrm{rec}}\) stays put**. It is not the homologous case in which \(D\) and \(R\) stretch together and \(D/R\) is frozen.

While we are still inside that sphere (\(D < R_{\mathrm{rec}}\)), the near side is about \(R_{\mathrm{rec}}-D\) and the far side about \(R_{\mathrm{rec}}+D\). As \(D\) grows, that split gets larger, not smaller. Crossing \(D_{\mathrm{recomb}}\) is the draft's switch: before it, the Hubble-bubble volume is zero; after it, the visible sphere grows at \(c\).

The draft's own CMB calculation is not finished. A kinematic–Doppler "kill function" is proposed to cancel the gravitational last-scattering multipole. The shapes are said to match; the amplitude is **not** yet shown to reach \(\sim 10^{-5}\). Treat that cancellation as an open problem, not as a result.

An earlier course note used a different picture: last-scattering radius \(\sim 14\,\mathrm{Gpc}\) and \(D/R \lesssim 0.02\). In this draft, \(\sim 14.3\,\mathrm{Gpc}\) is the **vision** radius, and the fitted \(D/R\) values are of order 2 to 7. Do not mix those two uses of "\(R\)."

---

## What the draft actually says about time

From the methods (step 4) and section 2.5, as a guide:

1. The CoE stays fixed.
2. Thick shells expand outward. The draft allows that expansion to be faster than light in a coordinate sense.
3. The sphere of events visible to us grows at \(c\).
4. Recombination is a spherical surface of fixed physical radius around the CoE.
5. The observer stays mid-shell, between the inner edge \(S_1\) and the outer edge \(S_2\).
6. The Hubble-bubble volume about the observer, \(V_{\mathrm{HB}}\), is zero until the observer passes \(D_{\mathrm{recomb}}\). After that crossing, the visible sphere grows at \(c\).
7. The draft suggests tracking the time-dependent volumes \(V_{\mathrm{HB}}(t)\), \(V_{S_1}(t)\), \(V_{S_2}(t)\) and the fractions of the bubble relative to the shell annulus and the outer sphere. Figure 9 is that sequence (frames A–D): shells cross the recombination surface while the visible sphere grows more slowly.

That is the dynamics. It is not "everything homologous, \(D/R\) constant."

---

## Distance from the observer to the fixed recombination sphere

Euclidean chord, same geometry as a point inside or outside a sphere. This is the draft's ruler, not a finished null geodesic in LTB. The draft itself says the full light-cone integral is still open.

**Inside** (\(D < R_{\mathrm{rec}}\)), with \(\theta = 0\) looking away from the CoE:

\[
\ell(\theta) = -D\cos\theta + \sqrt{R_{\mathrm{rec}}^2 - D^2\sin^2\theta}.
\]

Axis limits: \(\ell \approx R_{\mathrm{rec}}-D\) away from the CoE, \(\ell \approx R_{\mathrm{rec}}+D\) toward the CoE and through to the far side. Fractional split \(( \ell_{\mathrm{far}}-\ell_{\mathrm{near}} ) / (2 R_{\mathrm{rec}}) = D/R_{\mathrm{rec}}\).

**Time, if the observer rides an expanding shell and \(R_{\mathrm{rec}}\) is fixed:**

\[
\frac{d}{dt}\left(\frac{D}{R_{\mathrm{rec}}}\right) = \frac{\dot D}{R_{\mathrm{rec}}}.
\]

\(\dot D > 0\) while the shell expands, so \(D/R_{\mathrm{rec}}\) **grows**. The short path \(R_{\mathrm{rec}}-D\) shrinks. The long path grows. The angular contrast increases until the crossing. That is the opposite of a homologous freeze, and also the opposite of "the light cone runs outward so \(D/\chi\) falls."

**At the crossing** \(D = D_{\mathrm{recomb}}\). The draft identifies this with vision turning on, not with a calculated temperature. If \(D_{\mathrm{recomb}}\) means the same sphere as \(R_{\mathrm{rec}}\), the nearest distance to the surface hits zero as the observer meets it. The draft does not write that identification as an equation; do not force it.

**Outside** (\(D > R_{\mathrm{rec}}\)): the nearest point on the sphere, looking back toward the CoE, is \(D - R_{\mathrm{rec}}\). Rays that still hit the sphere are a cone, not a full sky of last-scattering in the inside-the-shell sense. Half-angle of a related eclipse in the appendix (inner shell, not recombination) is \(\alpha = \arcsin(R_{\mathrm{in}}/D)\). For the Newtonian \(\Omega_s = 0.049\) row, \(R_{\mathrm{in}}/D \approx 34.5/49.2 \approx 0.70\), and the appendix quotes \(\alpha \approx 33^\circ\) for its cartoon, so that cartoon is not the same ratio as the headline fit. Guide, not rigid.

The visible sphere is a **second** surface. After the crossing its radius grows at \(c\), independent of the fixed \(R_{\mathrm{rec}}\). Distance to recombination and radius of the vision bubble are different functions of time.

---

## Numbers in the draft (fits, not measurements of \(R_{\mathrm{rec}}\))

Table 1, best-fit sheet cells. Vision radius \(R \approx 14.30\,\mathrm{Gpc}\) in all three columns.

| Quantity | Newtonian, \(\Omega_s\) free | Newtonian, \(\Omega_s = 0.049\) | Relativistic LTB, \(\Omega_s = 0.049\) |
|----------|------------------------------|--------------------------------|----------------------------------------|
| \(D/R\) | 6.554 | 3.439 | 2.303 |
| \(D\) | 93.72 Gpc | 49.17 Gpc | 32.94 Gpc |
| \(R_{\mathrm{in}}\) | 79.4 Gpc (\(D-R\)) | 34.5 Gpc | 18.1 Gpc |
| \(R_{\mathrm{out}}\) | 108.0 Gpc (\(D+R\)) | 70.3 Gpc | 49.7 Gpc |
| \(R^2\) vs Migkas | 0.909 | 0.893 (0.54 geometry only) | 0.885 (0.57 geometry only) |

Headline the draft prefers for the relativistic case: \(D \approx 32.9\,\mathrm{Gpc}\), \(R_{\mathrm{in}} \approx 18.1\,\mathrm{Gpc}\), \(R_{\mathrm{out}} \approx 49.7\,\mathrm{Gpc}\). Newtonian \(\Omega_s = 0.049\): \(D \approx 49.17\,\mathrm{Gpc}\), \(R_{\mathrm{in}} \approx 34.48\,\mathrm{Gpc}\), \(R_{\mathrm{out}} \approx 70.27\,\mathrm{Gpc}\), \(\ell_0 = 13^\circ\).

No row gives \(R_{\mathrm{rec}}\) or \(D_{\mathrm{recomb}}\) as a number. Do not borrow 14.3 Gpc for that radius.

Parentheses on \(R^2\): geometry only, no sine. The draft treats the sine piece as part of the fit, not as the physics claim by itself. Same caution as before: \(R^2\) is a diagnostic against the Migkas \(A/A_{\mathrm{all}}\) curve, not a proof.

---

## What this does to the old dipole-budget note

The earlier note (`distance-to-recombination-vs-time-and-D.md` before this rewrite, and Module 6) assumed:

- \(R_{\mathrm{CMB}} \sim 14\,\mathrm{Gpc}\) as a last-scattering shell,
- homologous stretch, \(D/R\) constant in time,
- \((\delta T/T)_{\mathrm{geom}} \approx \alpha (D/R)(\Delta H/H)\) forcing \(D/R \lesssim 0.02\) at \(\Delta H/H \sim 0.06\).

The draft's guide replaces the first two:

- 14.3 Gpc is the vision radius around the **observer**.
- Fitted \(D\) is several times that radius, so the observer is not a small offset inside a Gpc-scale last-scattering sphere of that size.
- The draft says recombination's physical radius is fixed. Chris's working picture (above) replaces that: \(D_{\mathrm{recomb}}\) shrinks as mass leaves and \(T\) falls through 3000 K. \(D/R\) is not frozen either way.

The draft does **not** replace the CMB bound. It still treats Planck isotropy at \(\sim 10^{-5}\) as a fact to be survived. Its proposed mechanism is a kill function: an electromagnetic temperature pattern, kinematic–Doppler in character, opposite in sign to the gravitational multipole and matching in angular shape (even multipoles for concentric cells that share a bulk velocity). Figures 11–12: morphology at large \(D\), most demeaned power through \(\ell \le 4\) removed in exploratory partners, **not** yet \(10^{-5}\) on the full light cone.

Keep the \(D/R \lesssim 0.02\) arithmetic as a constraint on the **old** small-offset sketch. Do not apply it unchanged to Table 1. Those fits have \(D/R \sim 2\)–\(7\) with a different \(R\). Whether they survive a real last-scattering integral is the draft's open problem, not a solved cancellation.

---

## Sketch (guide geometry, two times)

```
CoE fixed. Recombination sphere: fixed physical radius R_rec.
Shell edges S1, S2 move out. Observer stays mid-shell, so D grows.

t early, D < D_recomb          t later, D > D_recomb
V_HB = 0                       visible sphere grows at c

        R_rec (fixed)                  R_rec (fixed)
     ················               ················
   ·                  ·           ·                  ·
  ·        CoE         ·         ·        CoE         ·
  ·         *          ·         ·         *          ·
  ·          \         ·         ·                    ·
  ·           * us     ·         ·          S1   us * ·  S2
  ·          D(t)      ·         ·         D(t) larger·
   ·                  ·           ·                  ·
     ················               ················
                                      \ vision bubble
                                        radius ~ c * time
                                        after the crossing
```

Not to scale. Table 1 has \(D\) several times the 14.3 Gpc vision radius. \(R_{\mathrm{rec}}\) is deliberately not drawn equal to that vision radius, because the draft does not set them equal.

---

## Caveats (where the guide is unfinished)

- Draft, 1 October 2026, marked not for submission. Acknowledgments say the author may be wrong. Use the geometry, not the conclusions, as settled.
- "Apparent late-time acceleration is an off-centre sampling artefact" and "no dark energy is required" are claims in the draft, not results of this note.
- The kill function is a proposal. Closing it to \(\sim 10^{-5}\) along null geodesics is stated as open.
- Chord lengths are not LTB null geodesics. The draft says the same.
- \(R_{\mathrm{rec}}\) and \(D_{\mathrm{recomb}}\) are not given as numbers. Any plot of \(\ell(t)\) needs an assumption the draft does not lock in.
- Two Newtonian columns differ by orders of magnitude in mass (\(M(<D) \sim 10^{56}\,\mathrm{kg}\) vs \(\sim 10^{54}\,\mathrm{kg}\)). The draft shows a \(D\)–core-mass degeneracy. Do not pick one column as "the" geometry.
- Alnes-style "observer within \(\sim 15\,\mathrm{Mpc}\)" and the old \(D/R \lesssim 0.02\) gate describe different models than Table 1. State them when they apply; do not paste them onto \(D \sim 30\)–\(90\,\mathrm{Gpc}\) without a new calculation.

---

## What would be a next calculation, still following the guide loosely

1. Choose whether \(D_{\mathrm{recomb}}\) equals the fixed recombination radius. The draft's prose suggests the crossing is when vision turns on; it does not write \(D_{\mathrm{recomb}} = R_{\mathrm{rec}}\).
2. Give \(D(t)\) from one shell worldline in the LTB equation the draft already writes, \(\dot R^2 = 2E(r) + 2GM(r)/R\), with the observer's shell label held fixed.
3. Then \(\ell(\theta, t)\) is the chord at that \(D(t)\) until the crossing, and a different exterior distance after it.
4. Only then compare the resulting temperature pattern with \(\sim 10^{-5}\). The draft has not done step 4.
