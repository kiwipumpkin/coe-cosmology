# Original exponential density profile (pre-ptshell)

**Date:** 2026-10-06 (America/Chicago)  
**Goal:** Recover Chris’s **original** CoE model — total mass + characteristic length with \(\rho\propto e^{-r/L}\)** — not the later ptshell constant-\(\Omega_s\) thick shell. No Migkas refit.

---

## 1. Exact original \(\rho(r)\)

### Physical form (what the spreadsheet encodes)

\[
\boxed{\rho(r)=\rho_0\,e^{-r/L}}
\qquad\text{with enclosed mass}\qquad
\boxed{M(<r)=M_\infty\,f(r/L),\quad
f(u)=1-e^{-u}\bigl(1+u+\tfrac12 u^2\bigr)}
\]

Normalization: \(M_\infty=\int_0^\infty 4\pi r^2\rho_0 e^{-r/L}\,dr=8\pi\rho_0 L^3\), so

\[
\rho_0=\frac{M_\infty}{8\pi L^3}.
\]

Chris’s sheet parameter \(M\) is the mass **inside Earth–CoE distance \(D\)**, not \(M_\infty\):

\[
M\equiv M(<D)=M_\infty\,f(D/L)
\quad\Rightarrow\quad
M_\infty=\frac{M}{f(D/L)},\qquad
\rho_0=\frac{M}{8\pi L^3\,f(D/L)}.
\]

Hence the working formula used everywhere in code:

\[
\boxed{
\rho(r)=\frac{M}{8\pi L^3\,f(D/L)}\,e^{-r/L}
}
\qquad
f(u)=1-e^{-u}(1+u+\tfrac12 u^2).
\]

*(Optional later add-on — same \(L\):)* a uniform floor \(\rho_b\) (bubble `Ob`) or analytic floor share \(f_u\) (no-boxcar notes):

\[
\rho_{\rm tot}(r)=\rho_0 e^{-r/L}+\rho_b.
\]

### Symbols

| Symbol | Meaning |
|--------|---------|
| \(r\) | Distance from CoE |
| \(L\) | Characteristic / scale length of \(e^{-r/L}\) (sheet `L`; code `LD` \(= L/D\)) |
| \(D\) | Earth–CoE distance today (sheet `D`; code unit of length) |
| \(d\) or \(R\) | Observer / vision-sphere radius (sheet `R` = Earth→P; code \(d=R\), often 14.3 Gpc) |
| \(M\) | Mass inside radius \(D\) (sheet `M` / `total mass` after undoing \(f(D/L)\)) |
| \(M_\infty\) | Asymptotic total mass of the pure exponential |
| \(f(u)\) | Incomplete-gamma / Fx enclosed-mass fraction |
| \(\mu\) | Code strength \(GM t_0^2/D^3\) (same \(M(<D)\)) |

Core ball (Chris’s wording): the sphere from the CoE out to the **near** side of the vision sphere, radius \(D-d\). The exponential fills that ball smoothly; there is **no** flat \(\Omega_s\) shell.

---

## 2. Where it lives (files & cells)

### Spreadsheet (source of truth for the original model)

**File:** `/workspace/cosmology/papers/spreadsheet_model/chris_model.xlsx`

| Sheet | Cells | Content |
|-------|-------|---------|
| `density` | A21–A22 | Label: `m(R)/M = 1 - exp(-X)*(1+X+{X²/2})`, `x=R/L` |
| `density` | B23 / **C61** | Characteristic length \(L\) (C61 = 2.0 in the current toy copy) |
| `density` | B27 / A27 | \(M\) (= 50 toy) — mass used as \(M(<D)\) |
| `density` | B28 / A28 | \(D\) (= 10 toy) — Earth to CoE |
| `density` | B29 / A29 | \(R\) (= 2 toy) — Earth to P (= \(d\)) |
| `density` | B55 | \(f(D/L)=(1-\exp(-B28/C61)*(1+B28/C61+B28^2/(2 C61^2)))\) |
| `density` | B56 | `total mass` \(=B27/B55\) \(= M_\infty\) |
| `density` | AN66… | \(M(<r_P)=f(E_{\rm row}/C\$61)*B\$56\) along sky rows |
| `density w rotation` | **A30/B30/C30** | `L` / 1.5 / **“density distribution”** — same \(L\); B55 uses `B30` |

Same \(f(u)\) appears on `density error fx` and the Copy-of-density sheets.

### Bubble / integrated code (analytic \(\rho\) + contrast)

| Path | Role |
|------|------|
| `/workspace/cosmology/papers/spreadsheet_model/bubble/bubble_core.py` L51–52 | `rho(r, mu, LD, mub) = mu*exp(-r/LD)/(8π LD³ Fx(1/LD)) + 3 mub/(4π)` |
| same file L70–71 | `contrast = (ρ(1−rs)−ρ(1+rs))/ρ(1)`, `c34` at ±3.4 Gpc |
| `/workspace/cosmology/papers/spreadsheet_model/integrated/integ_core.py` L1–14 | Docstring: “his … \(e^{-r/L}\) profile”; `Fx(X)=1-exp(-X)(1+X+X²/2)` |
| `/workspace/cosmology/papers/spreadsheet_model/bubble/prof_core.py` L14 | `kind=="exp"` → `np.exp(-r/L)` |
| `/workspace/cosmology/papers/spreadsheet_model/integrated.tex` | Paper writeup of Eq. (enclosed mass of \(e^{-r/L}\)) |
| `/workspace/cosmology/papers/spreadsheet_model/universe_bestfit.tex` | Drawing caption: grey shading \(e^{-r/L}\) |
| `/workspace/cosmology/papers/spreadsheet_model/uniform_bubble.tex` | Contrast limit \(\delta=[\rho(D-R)-\rho(D+R)]/\rho(D)\); exp as main profile |

### Later ptshell (NOT the original — for contrast)

`/workspace/cosmology/papers/spreadsheet_model/bubble/ptshell.py`: compact `exp(-RG/Lc)` **core only** + **constant** \(\Omega_s\) between `rin`–`rout`. That flat shell is what the matter-dipole note found has \(\delta\approx0\) at \(D\approx58.2\) Gpc.

### Same \(L\) in no-boxcar analytics

`/workspace/cosmology-course/notes/extrema-equations.md`: \(\rho(r)=\rho_{c0}e^{-(D/L)r}+\rho_u\), and the bracket

\[
S=(1-f_u)\bigl(6L^3 e^{D/L}-6L^3-6DL^2-3D^2 L\bigr)+f_u D^3
\]

is exactly \(L^3 N(D/L)\) from integrating that exponential+floor. Chris’s memory of \((1-f_u)(6L^3 e^{D/L}-\cdots)\) is this same characteristic length.

---

## 3. Signed contrast \((\rho(D-d)-\rho(D+d))/\bar\rho\)

Code definition (`bubble_core.run`), with \(\bar\rho=\rho(D)\) and \(d=R\):

\[
\delta(R)=\frac{\rho(D-R)-\rho(D+R)}{\rho(D)}.
\]

**Sign:** \(\delta>0\) means higher density on the **CoE-near** side (\(D-R\)) than the anti-CoE side (\(D+R\)) — the CatWISE-friendly sign if the count dipole is CoE-aligned.

### Pure exponential (no floor / `Ob=0`)

\[
\boxed{\delta(R)=2\sinh(R/L)}
\]

Independent of \(M\) and of \(D\) except through \(L\). (Huge unless \(L\) is a large fraction of \(D\), or a floor dilutes it.)

### With uniform floor / background

Numerator only sees the exponential piece; denominator is total density at \(D\):

\[
\delta(R)=f_{\rm CoE}\cdot 2\sinh(R/L),\qquad
f_{\rm CoE}=\frac{\rho_0 e^{-D/L}}{\rho_0 e^{-D/L}+\rho_b}
\]

(`fcoe` in `bubble_core`). Verified on the \(\delta\le2\%\) Migkas+\(\ell_0\) fit: `fcoe * 2 sinh(R/L) = contrast` exactly.

### Numbers

**A. Spreadsheet “your geometry” baseline** (`integrated.tex` / `his_geom`):  
\(D/R=5\), \(L/D=0.18\) ⇒ \(D\approx71.5\) Gpc, \(L\approx12.9\) Gpc (with \(R=14.3\) Gpc).

| \(R\) (Gpc) | pure-exp \(\delta\) |
|-------------|---------------------|
| 3.0 | \(2\sinh(3/12.9)\approx0.470\) (**47%**) |
| 14.3 | \(2\sinh(14.3/12.9)\approx2.71\) (**271%**) |

**B. Contrast as a function of \(L/D\) at \(D=58.2\) Gpc** (ptshell’s \(D\), for apples-to-apples with the flat-shell note; pure exp, no floor):

| \(L/D\) | \(L\) (Gpc) | \(\delta(R=3)\) | \(\delta(R=14.3)\) |
|--------:|------------:|----------------:|-------------------:|
| 0.02 | 1.16 | 13.09 | \(2.16\times10^5\) |
| 0.05 | 2.91 | 2.45 | 136 |
| 0.07 | 4.07 | 1.61 | 33.4 |
| 0.10 | 5.82 | 1.08 | 11.6 |
| 0.18 | 10.5 | 0.581 | 3.66 |
| 0.30 | 17.5 | 0.345 | 1.83 |
| 0.50 | 29.1 | 0.207 | 1.02 |
| 1.00 | 58.2 | 0.103 | 0.496 |

Few-percent CatWISE-scale \(\delta\) at \(R\sim3\)–14 needs either **large \(L/D\)** or (as in `uniform_bubble`) a **floor that drives \(f_{\rm CoE}\ll1\)** so \(\delta\approx f_{\rm CoE}\cdot 2\sinh(R/L)\) sits at 0.5–2%. Example: bubble fit `mean|c=0.02|a|Ob=0.3|mig_l0` has \(D\approx40.5\) Gpc, \(L\approx3.55\) Gpc, \(f_{\rm CoE}\approx3.6\times10^{-4}\), \(\delta(R=14.3)=2\%\) by construction, \(\delta(3.4\,\mathrm{Gpc})\approx0.08\%\).

**C. ptshell at \(D\approx58.2\) Gpc** (for comparison only): \(\delta(3)=\delta(14.3)\approx0\) — flat \(\Omega_s\) slab; see `matter-dipole-notes.md`.

---

## 4. How the exponential replaces the flat shell (CatWISE) while keeping the core ball

- **ptshell today:** empty-ish interior + constant-\(\Omega_s\) thick shell with Earth inside \(R_{\rm in}\ll D<R_{\rm out}\). Vision poles at \(D\pm R\) both sit in the same flat slab ⇒ **no** native matter dipole; Migkas \(A\) is kinematic.
- **Original exp:** one smooth \(\rho\propto e^{-r/L}\) (plus optional floor). Density **falls** with \(r\), so  
  \(\rho(D-d)>\rho(D+d)\) automatically ⇒ **signed** contrast for CatWISE without a hard shell edge.
- **Core ball kept:** still identify the ball of radius \(D-d\) (CoE → near vision surface) as the “core” region. For the exponential that ball is simply the inner part of the same profile (most of \(M(<D)\) already sits well inside \(D\) when \(L\ll D\)). No need for a separate constant-\(\Omega_s\) annulus; the “shell” for contrast is the vision sphere sampling two radii on the **same** falloff.
- **Practical CatWISE lever:** choose \(L\) (and optional \(\rho_b\) / \(f_u\)) so \(\delta(R)\) is a few percent at the relevant \(R\), while Migkas geometry still uses the same \(D\), \(d\), and enclosed-mass Gauss pull. Do **not** equate \(\delta\rho/\rho\) blindly to \(2\mathcal{D}_{\rm CatWISE}\) (tracer bias).

---

## 5. Pointer index

| Path | What |
|------|------|
| `chris_model.xlsx` → sheets `density`, `density w rotation` | Original sheet: \(M\), \(L\), \(D\), \(f(u)\) |
| `bubble/bubble_core.py` | Analytic \(\rho\), contrast |
| `integrated/integ_core.py` | `Fx`, integrated \(e^{-r/L}\) EOM |
| `bubble/prof_core.py` | `shape("exp")` |
| `bubble/ptshell.py` | **Later** flat-shell model (not original) |
| `integrated.tex`, `uniform_bubble.tex` | Written equations + \(\delta\) limit |
| `extrema-equations.md` / nobox notes | Same \(L\) in \(f_u\) bracket |
| `matter-dipole-notes.md` | ptshell \(\delta\approx0\) at \(D=58.2\); why CatWISE wants radial structure |

**Success check:** Chris can open the spreadsheet `density` sheet (A21, C61/B30, B55–B56) or `bubble_core.rho` and see \(\rho\propto e^{-r/L}\) with total mass + characteristic length — not the ptshell constant-\(\Omega_s\) shell.
