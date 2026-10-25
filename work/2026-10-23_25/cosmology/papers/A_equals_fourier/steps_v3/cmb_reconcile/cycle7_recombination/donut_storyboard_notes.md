# Donut-universe storyboard notes

**Status:** working notes for Chris’s Centre-of-Expansion (CoE) / donut-universe geometry, tied to the stacked storyboard in this directory.  
**Figure:** `donut_universe_storyboard.png` (built by `make_donut_universe_storyboard.py`).  
**Drawing units:** radii below are pedagogical drawing units; relative volumes follow \(R^3\).
**Layout:** the PNG is **geometry only** (vertical A→E stack). Stage prose and volume discussion live in the **right column of `letter_brief/letter.tex`**, not baked into the figure.

These notes record the storyboard rules that are easy to lose when regenerating panels: when the Hubble bubble appears, how volumes are defined, and what the geometry implies for apparent acceleration.

---

## 1. Illusory accelerated expansion

In this geometry the **apparent** accelerated expansion arises **without dark energy** and **without any actual outward acceleration of the shell**. The shells coast (or otherwise expand without a late-time \(\ddot a > 0\) driven by a dark component). What looks like acceleration is an **illusion** produced by geometry and by **off-centre sampling of the expansion field**: how an observer at distance \(D\) from the CoE measures \(H\) across the sky and along the past light cone. The storyboard therefore must not be read as a Friedmann \(\ddot a > 0\) epoch; it is a picture of concentric shells about a fixed CoE, with US mid-shell, in which directional and light-cone averages can mimic acceleration even though the material shell is not itself accelerating outward under dark energy.

---

## 2. Contracting universe as a distinguisher

A **contracting** universe would also show, on average, the **appearance** of accelerated expansion under analogous off-centre / light-cone sampling. The angular signature on the sky would differ from the expanding-shell case. That difference is a **distinguisher** and a natural **future check**: compare the predicted angular pattern of apparent acceleration (and of related Hubble anisotropy) in expanding versus contracting CoE geometries against directional \(H\) and supernova data, rather than relying on a single sky-averaged \(q\).

---

## 3. Hubble bubble volume until \(D_{\mathrm{recomb}}\)

Define the Hubble bubble as the observer-centred sphere of current vision radius \(R_{\mathrm{vis}}\) around US. In the storyboard,

\[
V_{\mathrm{Hubble\,bubble}}(t) = 0
\]

until the observer passes \(D_{\mathrm{recomb}}\) (the CoE distance at which vision turns on / recombination is crossed by US). Equivalently, \(R_{\mathrm{vis}} = 0\) on all panels with \(D \le D_{\mathrm{recomb}}\): no Hubble bubble is drawn (or only a point at US), and the annotation reads \(V_{\mathrm{HB}} = 0\).

Only after \(D > D_{\mathrm{recomb}}\) does the visible sphere around US grow at \(c\). Shells continue to expand faster than light in the coordinate sense used here, so \(R_{\mathrm{vis}}\) must not appear to outrun \(R_{\mathrm{in}}\) and \(R_{\mathrm{out}}\).

In the figure, \(D_{\mathrm{recomb}}\) coincides with the fixed CoE-centred recombination front radius \(R_{\mathrm{rec}}\) (purple dashed): US “passes recombination” when \(D = R_{\mathrm{rec}}\).

---

## 4. Scale, volumes, and volume fractions

The storyboard is scaled and annotated using volumes as functions of time (drawing units):

| Symbol | Meaning | Definition |
|--------|---------|------------|
| \(V_{\mathrm{HB}}(t)\) | Hubble bubble | \(\frac{4}{3}\pi R_{\mathrm{vis}}^{3}\) (zero while \(R_{\mathrm{vis}}=0\)) |
| \(V_{\mathrm{S1}}(t)\) | Outer sphere (through S1) | \(\frac{4}{3}\pi R_{\mathrm{out}}^{3}\) |
| \(V_{\mathrm{S2}}(t)\) | Inner sphere (through S2) | \(\frac{4}{3}\pi R_{\mathrm{in}}^{3}\) |
| \(V_{\mathrm{ann}}(t)\) | Shell material (annulus) | \(V_{\mathrm{S1}} - V_{\mathrm{S2}}\) |

**Volume fractions** (conceptual; discussed in the letter right column, not drawn on the PNG):

- \(V_{\mathrm{HB}} / V_{\mathrm{ann}}\) — bubble as a fraction of shell material volume;
- \(V_{\mathrm{HB}} / V_{\mathrm{S1}}\) — bubble as a fraction of the outer sphere;
- \(V_{\mathrm{S2}} / V_{\mathrm{S1}}\) — inner-to-outer sphere ratio (geometry of the cavity vs full ball).

Numerical stage table below is for author notes; the published figure carries only geometry labels (CoE, US, S1, S2).

---

## 5. Geometry constraints (kept on every panel)

- CoE fixed at the origin; US mid-shell: \(D = (R_{\mathrm{in}} + R_{\mathrm{out}})/2\).
- Recombination front: CoE-centred, purple dashed, fixed \(R_{\mathrm{rec}}\) in drawing units.
- After the shells have crossed recombination: \(R_{\mathrm{rec}} < R_{\mathrm{in}} < R_{\mathrm{out}}\).
- Shells expand faster than light; the visible sphere grows at \(c\) and must not appear to outrun the shells.
- Layout: vertical stack (A→E), geometry-only PNG for the letter left column; stage prose in `letter.tex` right column.
- Colours: red S2 (filled), blue-ish S1, purple dashed recomb., green/cyan visible sphere.
- Captions describe **what the panel shows**, not how the PNG was built.

---

## 6. Panel sequence (conceptual)

| Stage | Role |
|-------|------|
| **A** | Before \(D_{\mathrm{recomb}}\): shells still near/inside the recomb. front; \(R_{\mathrm{vis}}=0\), \(V_{\mathrm{HB}}=0\). |
| **B** | US at \(D = D_{\mathrm{recomb}}\): vision turns on but bubble volume still zero; \(V_{\mathrm{HB}}=0\). |
| **C** | Early post-recomb.: shells have crossed \(R_{\mathrm{rec}}\); visible sphere grows at \(c\). |
| **D** | Later expansion: mid-shell US; bubble a growing fraction of the annulus. |
| **E** | Today-like donut: thick shell, large \(D\), Hubble bubble still slower than shell growth. |

## 7. Numerical stages (drawing units)

Coast with \(c=1\), \(v_{\mathrm{in}}=2\), \(v_{\mathrm{out}}=4\), so \(D=(R_{\mathrm{in}}+R_{\mathrm{out}})/2\) and half-thickness grows at \(c\). Fixed \(R_{\mathrm{rec}}=D_{\mathrm{recomb}}=1.20\). Vision: \(R_{\mathrm{vis}}=\max(0,\,t-t_{\mathrm{rec}})\) with \(t_{\mathrm{rec}}=0.40\).

| Stage | \(t\) | \(D\) | \(R_{\mathrm{in}}\) | \(R_{\mathrm{out}}\) | \(R_{\mathrm{vis}}\) | \(R_{\mathrm{rec}}\) | \(V_{\mathrm{HB}}\) | \(V_{\mathrm{S1}}\) | \(V_{\mathrm{S2}}\) | \(V_{\mathrm{ann}}\) | \(V_{\mathrm{HB}}/V_{\mathrm{ann}}\) | \(V_{\mathrm{HB}}/V_{\mathrm{S1}}\) | \(V_{\mathrm{S2}}/V_{\mathrm{S1}}\) |
|-------|------:|------:|-------------------:|--------------------:|--------------------:|--------------------:|-------------------:|-------------------:|-------------------:|--------------------:|----------------------------------:|--------------------------------:|--------------------------------:|
| A | 0.25 | 0.75 | 0.50 | 1.00 | **0** | 1.20 | **0** | 4.19 | 0.52 | 3.67 | 0 | 0 | 0.125 |
| B | 0.40 | 1.20 | 0.80 | 1.60 | **0** | 1.20 | **0** | 17.16 | 2.14 | 15.01 | 0 | 0 | 0.125 |
| C | 0.90 | 2.70 | 1.80 | 3.60 | 0.50 | 1.20 | 0.52 | 195.4 | 24.4 | 171.0 | 0.0031 | 0.0027 | 0.125 |
| D | 1.60 | 4.80 | 3.20 | 6.40 | 1.20 | 1.20 | 7.24 | 1098 | 137.3 | 960.8 | 0.0075 | 0.0066 | 0.125 |
| E | 2.40 | 7.20 | 4.80 | 9.60 | 2.00 | 1.20 | 33.51 | 3706 | 463.2 | 3243 | 0.0103 | 0.0090 | 0.125 |

Confirmed: \(V_{\mathrm{HB}}=0\) on stages A and B (\(D\le D_{\mathrm{recomb}}\)). For C–E, \(R_{\mathrm{rec}}<R_{\mathrm{in}}<R_{\mathrm{out}}\) and \(R_{\mathrm{vis}}\) stays inside the shell half-thickness.
