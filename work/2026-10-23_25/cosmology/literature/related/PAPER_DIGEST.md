# Literature digest — Chris’s Sedov–Taylor + Gauss centered-universe program

**Work dir:** `/workspace/cosmology/literature/related/`  
**Digest date:** 2026-09-19 (America/Chicago)  
**Chris thesis context:** real Center of Expansion (CoE) beyond horizon; off-center observer in Sedov–Taylor (ST) blast-wave expansion; enclosed-mass Gauss \(g\propto 1/r^2\) with decreasing \(\rho\); match Migkas 2020 \(A\) vs Galactic longitude; prefer ideal time-integrated gravity over bolt-on kick; challenge cosmological principle; honest tensions with \(\Lambda\)CDM/CMB.

**File integrity notes (important):**
- `alnes_amarzguioui_2006_sn.pdf` is **misfiled**: content is Capuzzo-Dolcetta et al. (astro-ph/0611205) on triaxial galaxies + DM haloes — **not** Alnes & Amarzguioui SN. Digested as-found; intended SN companion is missing.
- `smoller_temple_2003_shock_bh.pdf` is **misfiled**: content is Nipoti et al. (astro-ph/0306082) on galactic cannibalism in C0337-2522 — **not** Smoller–Temple. Digested as-found; intended Smoller–Temple 2003 is missing.
- `fennelly_nasa_cosmological_bh.pdf` is Debney (NASA TM X-73378, Feb 1977), not Fennelly — still on-topic for “cosmological black hole.”

Numbers below are taken from the extracted PDF text only; where OCR/layout garbled a figure, we say so rather than invent.

---

## 1. Fennelly 1977 — Hubble anisotropy + large-scale inhomogeneity

### Citation / year
A. J. Fennelly, *Mon. Not. R. astr. Soc.* **181**, 121–130 (1977). File: `fennelly_1977_hubble_anisotropy.pdf`.

- **One-sentence gist:** Fits reported large-scale Hubble-parameter anisotropy to a tilted self-similar cosmology with a Gaussian density inhomogeneity, placing us ~120 Mpc off-center, and argues this removes the conflict with microwave isotropy that a pure peculiar-velocity interpretation would create.
- **Setup:** Homogeneous background + large-scale Gaussian density distribution in a tilted self-similar spacetime; Galaxy 120 Mpc from centre toward \(\alpha=12^{\mathrm{h}}42^{\mathrm{m}}\), \(\delta=+20^\circ\); inhomogeneous region diameter ~540 Mpc.
- **What they compute:** Direction-dependent Hubble expansion; gravitational velocity shifts vs expansion contours; rough estimates of number-count anisotropy and luminosity-function anisotropy between hemispheres; qualitative CMB/microwave compatibility.
- **Key quantitative claims:** Diameter ~540 Mpc; offset ~120 Mpc toward \(\alpha=12^{\mathrm{h}}42^{\mathrm{m}}\), \(\delta=+20^\circ\); notes ~500±125 km/s peculiar-velocity interpretation conflicts with microwave limits ~250 km/s (as cited then); dimensional number-count estimate \(|N_1-N_2|/(N_1+N_2)\) over ~\(10^2\) Mpc giving density contrast order ~0.46 and \(N_1/N_2\sim 2.7\) (compared to de Vaucouleurs ratios ~3.2–4.1).
- **Overlap with Chris’s ST+Gauss+Migkas program:** Strong conceptual cousin — **off-center observer**, large-scale inhomogeneity, Hubble anisotropy vs longitude/sky direction, and explicit desire to keep CMB isotropy without discarding the anisotropy signal. Does **not** use Sedov–Taylor blast, Gauss enclosed-mass \(g\propto1/r^2\) cosmology, or X-ray \(L_X{-}T\) \(A(\ell)\). Metric is tilted self-similar + Gaussian density, not ST shell.
- **Tension / adversary note:** Copernican cost already acknowledged (special place relative to centre). 1970s Hubble-anisotropy data are historically contested; modern CMB dipole is kinematic-dominated in standard framing.
- **Quotable equations or phrases:** “This eliminates the conflict between the anomaly and the isotropy of the microwave background.” Fit to “tilted self-similar cosmology with a Gaussian density distribution superposed on a homogeneous Universe.”

---

## 2. Debney 1977 (NASA TM) — Techniques of global analysis: a cosmological black hole?

### Citation / year
George Debney, NASA TM X-73378 (February 1977). File: `fennelly_nasa_cosmological_bh.pdf` (filename misleading).

- **One-sentence gist:** Elementary model of freely falling observers/emitters inside a Schwarzschild black hole, computing redshift and luminosity distance to argue that an observer’s sky is bounded, with an exotic “curiosity” application to a cosmological-scale black hole.
- **Setup:** Eddington–Finkelstein / scaled Schwarzschild interior (\(R=r/2m\), horizon at \(R=1\)); test particles; no backreaction; spherical symmetry.
- **What they compute:** Redshift formulas; luminosity distance \(d\) for radial and non-radial emitters; maps of redshift/blueshift anisotropy vs luminosity distance inside the horizon.
- **Key quantitative claims:** Abstract states no luminosity distances greater than approximately [OCR-garbled bound] within the black hole’s radius (exact printed number unclear in extract — treat as UNKNOWN DETAILS if citing precisely). Horizon at \(R=1\) in mass-scaled coordinates.
- **Overlap with Chris’s program:** Transfers the **“universe inside a BH / beyond-horizon CoE”** language and bounded-sky intuition; does **not** provide ST dynamics, Gauss gravity for decreasing \(\rho\), or Migkas matching. Useful as historical scaffolding for Pathria/Debney–Farnsworth lineage.
- **Tension / adversary note:** Explicitly a “curiosity”; not an observational cosmology model. Smooth FRW–Schwarzschild identification later contested (see Pathria note §3).
- **Quotable equations or phrases:** Scaled metric with horizon at \(R=1\); “An exotic application can be made, as a curiosity, to a black hole…” of cosmological mass.

---

## 3. Khakshournia 2014 — A note on Pathria’s model of the universe as a black hole

### Citation / year
S. Khakshournia, arXiv:1412.0105 (2014). File: `pathria_note_1412.0105.pdf`.

- **One-sentence gist:** Pathria’s closed dust FRW interior identified with a Schwarzschild(+Λ) horizon at maximum expansion is **not** a smooth match: Barrabès–Israel null-shell formalism shows a null shell with **surface pressure**.
- **Setup:** Pressureless closed FRW (retarded coordinates) glued to Schwarzschild–de Sitter exterior along null hypersurface \(\Sigma\) at \(R_{\max}\equiv R_s\); Barrabès–Israel null shell.
- **What they compute:** Junction conditions across \(\Sigma\); jump of extrinsic curvature ⇒ surface stress; concludes surface pressure is required.
- **Key quantitative claims:** Qualitative: matching “is not smooth”; null hypersurface “is the history of a null shell admitting a surface pressure.” (No Migkas-style observational numbers.)
- **Overlap with Chris’s program:** Directly relevant to **beyond-horizon / universe-as-BH** matching honesty — if Chris uses an FRW- or ST-like interior with an exterior vacuum, this paper is a warning that null matching can force a shell (surface stress), not a passive horizon identification. Does not treat ST blast or off-center observers.
- **Tension / adversary note:** Undercuts naive “universe = BH interior” slogans; notes Knutsen (2009) critique of Pathria via timelike matching, and clarifies null-horizon limit still needs a shell.
- **Quotable equations or phrases:** “the matching is not smooth, and in fact, the null hypersurface is the history of a null shell admitting a surface pressure.”

---

## 4. Alnes & Amarzguioui 2006 — CMB anisotropies for off-center LTB observer

### Citation / year
H. Alnes & M. Amarzguioui, arXiv:astro-ph/0607334 (2006). File: `alnes_amarzguioui_2006.pdf`.

- **One-sentence gist:** In dust LTB void models that fit SNIa + first CMB peak with an on-center observer, an off-center observer induces a CMB dipole that forces the observer within ~15 Mpc of the centre — too close for the induced quadrupole/octopole to explain low-\(\ell\) alignment.
- **Setup:** LTB metric \(ds^2=-dt^2+[R']^2/(1+\beta(r))\,dr^2+R^2 d\Omega^2\); underdense spherical bubble in Einstein–de Sitter exterior; observer displaced from centre; Models I & II (I = prior centered fit).
- **What they compute:** Null geodesics for off-center observers; temperature map \(\delta T/T\); multipoles \(a_{\ell 0}\) vs displacement \(d_{\mathrm{obs}}\); compatibility with COBE-scale dipole.
- **Key quantitative claims:** Observer must lie within \(\sim 15\,\mathrm{Mpc}\) for induced dipole ≲ observed (\(|a_{10}|\sim 10^{-3}\)); at \(d_{\mathrm{obs}}=15\,\mathrm{Mpc}\) in Model I, exact \(a_{10}=1.4\times10^{-3}\), \(a_{20}=2.3\times10^{-7}\), \(a_{30}=1.3\times10^{-9}\) (quad/oct too small for alignment); illustrative maps also at 20 and 200 Mpc.
- **Overlap with Chris’s program:** Core **off-center observer + spherical CoE** lesson transfers: dipole is the killer constraint. Does **not** use ST self-similar blast or Gauss \(g\propto1/r^2\); LTB void is dust inhomogeneity, not an explosion shell. Alignment motivation fails once dipole is respected — Chris should not expect off-center ST alone to “fix” CMB low-\(\ell\) alignment without extra freedom.
- **Tension / adversary note:** Severe Copernican fine-tuning (~15 Mpc in Gpc-scale void). Intrinsic/ISW dipoles ~\(10^{-5}\) unlikely to cancel a large geometric dipole. Honest \(\Lambda\)CDM adversary: CMB still prefers near-isotropy about us.
- **Quotable equations or phrases:** \(d_{\mathrm{obs}}\lesssim 15\,\mathrm{Mpc}\); “for such small displacements… the induced quadru- and octopoles turn out to be insufficiently large to explain the alignment.”

---

## 5. `alnes_amarzguioui_2006_sn.pdf` — MISFILED (not Alnes SN)

### Citation / year (as found on disk)
R. Capuzzo-Dolcetta, L. Leccese, D. Merritt, A. Vicari, arXiv:astro-ph/0611205 (2006): *Self-consistent models of cuspy triaxial galaxies with dark matter haloes*. Filename claims Alnes SN paper — **content does not match**.

- **One-sentence gist:** Schwarzschild orbit-superposition models show triaxial ellipticals can coexist with triaxial DM haloes of several axis-ratio choices; chaos implies possible slow halo shape evolution.
- **Setup / compute / claims:** Galactic dynamics; not cosmology. Axis ratios e.g. luminous 0.7:0.86:1; halo variants more prolate/oblate.
- **Overlap with Chris’s program:** Essentially **none** for ST+Gauss+Migkas. Treat as accidental download.
- **Tension / adversary note:** N/A for cosmology.
- **Quotable:** “Our results demonstrate for the first time that triaxial dark matter haloes can co-exist with triaxial galaxies.”
- **Action for parent:** Re-fetch intended Alnes/Amarzguioui SN / off-center luminosity-distance paper (likely companion to astro-ph/0512006 / 0607334 lineage). Mark **UNKNOWN DETAILS** for the intended SN analysis until PDF is replaced.

---

## 6. Saito, Ishibashi & Kodama 2010 — Analytic CMB formulae in LTB

### Citation / year
K. Saito, A. Ishibashi, H. Kodama, arXiv:1002.3855 (2010). File: `kodama_etal_2010_cmb_ltb.pdf`. (4 pp. letter-style.)

- **One-sentence gist:** Derives analytic dipole and quadrupole CMB temperature formulae for LTB voids and recovers the ~15 Mpc off-center bound, consistent with Alnes & Amarzguioui numerics.
- **Setup:** LTB dust metric with curvature function \(k(r)\); observer at displacement \(\delta L\) from symmetry centre; last-scattering treated via photon distribution at centre + derivatives.
- **What they compute:** Analytic \((\delta T/T)^{(1)}\) and \((\delta T/T)^{(2)}\); numerical evaluation of \(a_{10}\) for an Alnes-like void profile.
- **Key quantitative claims:** For their fiducial void parameters (\(H_\perp(t_0,r_{\mathrm{out}})=51\,\mathrm{km\,s^{-1}Mpc^{-1}}\), etc.), \(\delta L\lesssim 15\,\mathrm{Mpc}\); consistency with Alnes & Amarzguioui [10].
- **Overlap with Chris’s program:** Transfers **analytic off-center dipole machinery** that could be adapted to any spherically symmetric ST-like metric (later generalized in Yoo/Kodama 2010). Not ST; no Migkas \(A(\ell)\).
- **Tension / adversary note:** Same Copernican dipole bound; BAO/kSZ cited as further pressure on LTB (footnote).
- **Quotable equations or phrases:** Dipole formula schematically involving \(\delta L_{nj}\), \(H_\parallel=Ṙ'/R'\), and radial derivatives of the distribution function; \(\delta L\lesssim 15\,\mathrm{Mpc}\).

---

## 7. Kodama, Saito & Ishibashi 2010 — Off-center CMB in general spherical symmetry

### Citation / year
H. Kodama, K. Saito, A. Ishibashi, arXiv:1004.3089 (2010). File: `yoo_etal_2010_offcenter_cmb.pdf` (filename “yoo” is misleading; authors are Kodama–Saito–Ishibashi).

- **One-sentence gist:** Generalizes analytic CMB dipole/quadrupole formulae to **any** spherically symmetric universe and updates the void off-center bound with WMAP to \(\delta L\lesssim 16\,\mathrm{Mpc}\) (example models \(\lesssim 14\,\mathrm{Mpc}\)).
- **Setup:** General spherical symmetry (not only LTB dust); photon distribution \(F\) near centre; observer displacement \(\delta L_n\); WMAP dipole normalization.
- **What they compute:** Model-independent dipole/quadrupole formulae; application to void models; consistency check vs prior numerics.
- **Key quantitative claims:** Induced \(a_{10}\lesssim 1.23\times10^{-3}\) (WMAP) ⇒ \(\delta L\lesssim 16\,\mathrm{Mpc}\); other LTB examples \(\delta L\lesssim 14\,\mathrm{Mpc}\); small (few×100 Mpc) voids already ruled out by other data (stated in intro).
- **Overlap with Chris’s program:** **Highest-transfer CMB toolkit** in this set for a spherical ST CoE: formulae do not require LTB dust specifically. Still: Gauss \(g\) and Migkas \(A(\ell)\) are outside scope. Off-center ST observers will face the same dipole ceiling unless compensated by peculiar velocity (fine-tuned).
- **Tension / adversary note:** Reinforces that spherical CoE + off-center us ⇒ CMB dipole dominates long before SN anisotropy does.
- **Quotable equations or phrases:** Dipole \(\propto \delta L_n\cdot\Omega\) times radial derivatives of \(F_0\); “update the constraints… with the latest WMAP data.”

---

## 8. Blomqvist & Mörtsell 2010 — SNe as seen by off-center void observers

### Citation / year
M. Blomqvist & E. Mörtsell, arXiv:0909.4723 (2010). File: `blomqvist_mortsell_2010_sn_void.pdf`.

- **One-sentence gist:** SNIa alone allow off-center observers out to ≳15% of the void scale radius, but combining with the CMB dipole forces the observer within ~1% of the void scale — strong Copernican fine-tuning against LTB voids.
- **Setup:** LTB local void (dust); SDSS-II first-year SN and Constitution SN sets; on- and off-center luminosity distances; CMB dipole prior.
- **What they compute:** Fits of void profile parameters; anisotropic SN Hubble diagram for displaced observers; joint SN+CMB location constraints.
- **Key quantitative claims:** On-center: SDSS-II best LTB \(\chi^2=229.3\) vs flat \(\Lambda\)CDM \(\chi^2=231.3\); Constitution best LTB \(\Omega_{\mathrm{in}}=0.13\), \(\Omega_{\mathrm{out}}=1\), \(r_s=3.5\,\mathrm{Gpc}\) (\(z\approx1.02\)), \(\chi^2=461.0\) vs \(\Lambda\)CDM \(465.5\); SN-only: displacement ≳15% of void scale still acceptable; with CMB dipole: within ~1% of void scale.
- **Overlap with Chris’s program:** Critical **methodology cousin** for off-center observables (anisotropic \(d_L(\hat n)\)) and for ranking SN vs CMB constraints. Transfers the warning that SN are weak on centering; CMB is ruthless. Not ST/Gauss/Migkas.
- **Tension / adversary note:** Explicit Copernican indictment of voids; Gpc voids improbable in standard structure formation (cited). Honest: SN fits can look competitive with \(\Lambda\)CDM in \(\chi^2\) while still being cosmologically baroque.
- **Quotable equations or phrases:** “the data compels the observer to be located within about one percent of the void scale radius”; “strongly disfavouring the model from a Copernican principle point of view.”

---

## 9. Bertschinger 1983 — Cosmological self-similar shock waves

### Citation / year
E. Bertschinger, *ApJ* **268**, 17–29 (1983). File: `bertschinger_1983_sedov_eds.pdf`.

- **One-sentence gist:** Finds cosmological self-similar blast-wave solutions in Einstein–de Sitter: early Sedov-like phase morphs into a thin dense shell that cools by expansion and can fragment into galaxies (Ostriker–Cowie explosive amplification).
- **Setup:** Point explosion in expanding, self-gravitating EdS (\(\Omega=1\)) or empty (\(\Omega=0\)) dust; adiabatic \(\gamma\); strong-shock jump conditions; similarity variable \(\lambda=r/r_s(t)\).
- **What they compute:** Similarity ODEs; shock trajectory; postshock profiles; application to explosive galaxy formation (cooling, Jeans mass, generational shells).
- **Key quantitative claims:** Initial Sedov behavior for \(t\sim t_i\); for \(\Omega=1\), shock decelerates until velocity ~\(6/5\) Hubble flow; dense shell thickness ~\(0.03\,r_s\); before \(1+z\sim9\), Compton/free-free cooling prevent galactic-scale shocks; one amplification generation ~\(3\times10^{11}\,M_\odot\) shell (group-like); possible second ~\(2\times10^{13}\,M_\odot\) (small-cluster-like). Exact OCR of some formula prefactors is messy — treat detailed coefficients with care when citing.
- **Overlap with Chris’s program:** **Closest classical ST-in-cosmology cousin.** Transfers: Sedov → cosmologically modified self-similar blast; thin shell; EdS background. Does **not**: place observer off-center for CMB/Hubble dipole; Gauss \(g\propto1/r^2\) from enclosed mass as a directional gravity program; Migkas \(A(\ell)\); beyond-horizon CoE/BH framing. Chris’s “ideal time-integrated gravity” is a different physical question than Bertschinger’s galaxy-formation cooling shells.
- **Tension / adversary note:** Applied to galaxy formation, not alternative cosmology; radiative cooling and structure-formation context assume standard FLRW background outside the blast.
- **Quotable equations or phrases:** “Initially the shock wave formed satisfies the Sedov solution… but it is modified by the cosmological expansion to form a self-similar thin dense shell”; shell thickness \(\sim0.03\,r_s\).

---

## 10. Smoller & Temple 1998 — Shock-wave cosmology (FRW–TOV matching)

### Citation / year
J. Smoller & B. Temple (1998 shock-wave cosmology paper in extract; file `smoller_temple_1998_shock_cosmo.pdf`). Long GR/PDE paper (~66 pp. PDF).

- **One-sentence gist:** Constructs an exact Einstein-equation shock matching a critically expanding (\(k=0\)) FRW interior to a TOV exterior, motivated by a leading-edge blast explaining Hubble expansion — with Copernican and information-loss consequences — and derives present-day shock-position bounds from \(H_0\) and \(T_{\mathrm{CMB}}\).
- **Setup:** \(k=0\) FRW inside shock; TOV static fluid outside; Rankine–Hugoniot / GR shock conditions; equation of state general enough to include matter+radiation discussion.
- **What they compute:** Shock ODEs; phase-plane analysis; constraints from \(\bar\rho>\bar p>0\); estimates of shock travel relative to freefall; bounds on present shock radius vs \(T_0\), \(H_0\), \(R_*\).
- **Key quantitative claims (from text):** Shock position at “present” (when \(H_0\) and \(T_0\) match observations) is **comparable to the Hubble length** in the model’s narrative summary; example with \(T_0\approx2.7\,\mathrm{K}\), \(R_*=2.7/4000\), \(h_0=0.55\): shock travel excess \(r^2-r_*^2\approx(0.019)^2/H_0^2\) (i.e. ~0.019 Hubble lengths of “extra” shock motion over freefall); pure-radiation comparison at \(T_0=2.7\,\mathrm{K}\) gives shock distance **significantly beyond** the Hubble length, with bounds of the schematic form \(36h_0/H_0\le r\le\cdots\). Use quoted formulae carefully — they are model-internal.
- **Overlap with Chris’s program:** Major cousin for **explosion / shock at the edge of the expanding region**, entropy selecting explosion over implosion, and willingness to break Copernican symmetry + early-universe reverse-engineering. Differs: Smoller–Temple use FRW–TOV shock matching, not Sedov–Taylor blast into EdS with Gauss \(g\); not aimed at Migkas \(A(\ell)\). Beyond-Hubble shock foreshadows 2002 “inside BH” paper.
- **Tension / adversary note:** Explicit Copernican violation; shock destroys information about pre-shock early universe (anti-inflation-reconstruction). CMB isotropy still needs explanation if we are not centred.
- **Quotable equations or phrases:** “the expansion of the universe… might be the result of a large scale localized explosion that generated a shock–wave at the leading edge”; shock “comparable to the Hubble distance”; entropy ⇒ loss of early-universe information.

---

## 11. Smoller & Temple 2002 — Shock-wave cosmology inside a black hole

### Citation / year
J. Smoller & B. Temple, arXiv:astro-ph/0210105 (2002). File: `smoller_temple_2002_shock_bh.pdf`.

- **One-sentence gist:** Extends FRW–TOV shock matching **beyond one Hubble length**, placing the exterior inside a black hole; Big Bang becomes a localized explosion inside an asymptotically flat Schwarzschild BH, with shock weakening to an Oppenheimer–Snyder interface.
- **Setup:** Critically expanding FRW matched to TOV-inside-BH metric across a shock with \(\bar r>\bar r_{\mathrm{crit}}=\mathrm{Hubble\ length}\) (for \(k=0\)); \(p=\rho/3\) distinguished for finite nonzero (light) emergence speed from \(t=0\).
- **What they compute:** Global exact solutions; mass behind shock decreasing; entropy condition selecting explosion; large-time OS asymptotics; qualitative thermalization beyond the central observer’s light cone without inflation.
- **Key quantitative claims:** For \(k=0\), \(\bar r_{\mathrm{crit}}\) equals the Hubble length — cannot match classical TOV beyond that without continuing into a BH. \(p=\rho/3\) ⇒ shock emerges at light speed. (Short paper; many details deferred to forthcoming work.)
- **Overlap with Chris’s program:** **Closest “beyond-horizon CoE / cosmological BH” GR cousin.** Transfers: CoE + blast origin of expansion; shock beyond Hubble scale; BH exterior. Does not deliver ST self-similarity, Gauss \(g(\rho)\), or Migkas longitude pattern. Ideal for contrasting Chris’s preferred **time-integrated gravity** vs Smoller–Temple’s continuum shock matching.
- **Tension / adversary note:** Maximal Copernican cost; standard early-universe + inflation narrative replaced. Must still face CMB/near-isotropy about us.
- **Quotable equations or phrases:** “the Big Bang arises from a localized explosion occurring inside the Black Hole of an asymptotically flat Schwarzschild spacetime”; “thermalize the radiation in a region well beyond the light cone… even though the model does not invoke inflation.”

---

## 12. `smoller_temple_2003_shock_bh.pdf` — MISFILED (not Smoller–Temple)

### Citation / year (as found on disk)
C. Nipoti, M. Stiavelli, L. Ciotti, T. Treu, P. Rosati, arXiv:astro-ph/0306082 (2003): *Galactic cannibalism in the galaxy cluster C0337-2522 at \(z=0.59\)*. Filename claims Smoller–Temple 2003 — **content does not match**.

- **One-sentence gist:** N-body study of five ellipticals in a \(z=0.59\) cluster core as a candidate early cD-formation (cannibalism) system.
- **Overlap with Chris’s program:** None for ST+Gauss cosmology.
- **Action for parent:** Re-fetch intended Smoller–Temple 2003 shock/BH follow-up if needed; mark UNKNOWN DETAILS until replaced.

---

## 13. Migkas et al. 2020 — Probing isotropy with \(L_X{-}T\) (focus: abstract, method, Fig. 6)

### Citation / year
K. Migkas et al., arXiv:2004.03305 / *A&A* (2020). File: `migkas2020.pdf` (44 pp.; digest uses pp. 1–20 focus; appendices skipped per brief).

- **One-sentence gist:** Directional \(L_X{-}T\) normalization \(A\) for 313 MCXC clusters (core-excised \(T\)) varies across the sky at \(\gtrsim4\sigma\) toward \((l,b)\sim(280^\circ,-20^\circ)\); joint with prior samples \(\sim5\sigma\) toward \((303^\circ,-27^\circ)\), with max \(D_L\) variation \(\sim16\pm3\%\).
- **Setup / method:** \(L_X\) depends on cosmology via luminosity distance \(D_L\); spectroscopic \(T\) does not. Fit \(L_X{-}T\) (normalization \(A\), slope \(B\), intrinsic scatter). Map anisotropy via (i) \(\Delta l=90^\circ\) longitude bins (Fig. 6) and (ii) scanning cones \(\theta=90^\circ,75^\circ,60^\circ,45^\circ\). Full-sky fit: \(A_{\mathrm{all}}=1.114^{+0.044}_{-0.040}\), \(B_{\mathrm{all}}=2.102\pm0.064\), \(\sigma_{\mathrm{int}}=0.242\,\mathrm{dex}\). Often fix \(B=B_{\mathrm{all}}\) and map \(A/A_{\mathrm{all}}\). Systematics battery: absorption, groups/low-\(z\), metallicity, instrument Chandra vs XMM, bootstraps (\(10^5\)), etc. — none explained away the anisotropy in the authors’ tests.
- **What they compute (for Chris):** Sky maps and longitude dependence of \(A\) (and \(B\)); significance vs rest-of-sky; joint analysis with ACC + XCS-DR1 (842 clusters).
- **Key quantitative claims:** \(\gtrsim4\sigma\) toward \((280^\circ,-20^\circ)\); joint \(\sim5\sigma\) toward \((303^\circ,-27^\circ)\); \(\Delta D_L\sim16\pm3\%\) max; Fig. 6: \(A\) varies ~31% with longitude while \(B\) varies ~18% with milder tension (largest independent-bin \(A\) tension ~3.2\(\sigma\); \(B\) ~1.8\(\sigma\)); lowest-\(A\) longitude sector \(\ell\in[210^\circ,300^\circ]\) (~2.9\(\sigma\) vs rest).
- **Overlap with Chris’s program:** **Primary observational target.** Chris wants ST+Gauss off-center gravity to **reproduce \(A(\ell)\)** (and related \(D_L\) anisotropy). Migkas supplies the empirical pattern and the clean observable (\(A\) mirrors \(D_L\) under isotropy assumptions). Does not provide a physical CoE model.
- **Tension / adversary note:** Could be cosmological **or** X-ray systematics; authors stress identification is crucial. Direction is **not** the CMB kinematic dipole apex \((264^\circ,+48^\circ)\); closer to some other reported late-universe anisotropy directions (authors compare to SNIa etc.). \(\Lambda\)CDM assumes isotropy of scaling relations — this is a direct challenge if cosmological.
- **Quotable equations or phrases:** “perfectly mirrored in the behavior of the normalization \(A\) of the \(L_X{-}T\) relation”; “maximum variation of \(D_L\) seems to be \(\sim16\pm3\%\)”; Fig. 6 = \(A/A_{\mathrm{all}}\) and \(B/B_{\mathrm{all}}\) vs Galactic longitude (\(\Delta l=90^\circ\) strips).

---

## 14. Heinesen 2020/2021 — Multipole decomposition of the general \(d_L\) Hubble law

### Citation / year
A. Heinesen, arXiv:2010.06534 (v2 2021). File: `heinesen_2020_multipole_hubble.pdf`.

- **One-sentence gist:** Model-independent series for luminosity distance vs redshift in a **general** spacetime, with coefficients as a finite set of multipoles — 9 / 25 / 61 degrees of freedom at \(O(z)\) / \(O(z^2)\) / \(O(z^3)\) — enabling anisotropic effective \(H_0\), \(q_0\), etc., and direct FLRW tests.
- **Setup:** No metric ansatz; geometric optics + observer 4-velocity congruence; Etherington reciprocity; series in redshift about the observer.
- **What they compute:** Multipole expansion of the generalized Hubble law; effective direction-dependent cosmographic parameters; discussion of peculiar motion.
- **Key quantitative claims:** Exact anisotropic Hubble law has **9, 25, 61** dof at \(O(z),O(z^2),O(z^3)\) respectively (finite because of tensor symmetries / traceless reductions).
- **Overlap with Chris’s program:** Ideal **language layer** for forecasting what an off-center ST+Gauss observer predicts for \(d_L(z,\hat n)\) and for comparing to Migkas-like directional \(D_L\) without committing to FLRW. Does not itself supply ST dynamics or CMB multipoles.
- **Tension / adversary note:** Framework is neutral; makes FLRW falsifiable in the neighbourhood of the observer. Peculiar-velocity accounting remains subtle (paper discusses assumptions).
- **Quotable equations or phrases:** “exact anisotropic Hubble law is given by 9, 25, 61 degrees of freedom”; “direct testing of the FLRW ansatz.”

---

## 15. Saadeh et al. 2016 — How isotropic is the Universe? (Bianchi / shear limits)

### Citation / year
D. Saadeh, S. M. Feeney, A. Pontzen, H. V. Peiris, J. D. McEwen, arXiv:1605.07178 (2016). File: `saadeh_etal_2016_bianchi.pdf`. (Focus: abstract + main limits.)

- **One-sentence gist:** First general Planck temperature+polarization test of **all** Bianchi VIIh shear modes finds anisotropic expansion strongly disfavoured (odds ~121,000:1 against), with \((\sigma_V/H)_0<4.7\times10^{-11}\) (95% CI) and weakest mode \((\sigma_{T,\mathrm{reg}}/H)_0<1.0\times10^{-6}\) (95% CI).
- **Setup:** Homogeneous but anisotropic Bianchi models (VIIh and limits VII0, V, I); shear decomposed into scalar/vector/regular-tensor/irregular-tensor modes; Planck CMB T+P.
- **What they compute:** Limits on \((\sigma/H)_0\) per mode; Bayesian evidence vs isotropic \(\Lambda\)CDM.
- **Key quantitative claims:** \((\sigma_V/H)_0<4.7\times10^{-11}\) (95% CI); \((\sigma_{T,\mathrm{reg}}/H)_0<1.0\times10^{-6}\) (95% CI); anisotropic expansion odds **121,000:1 against**; also vorticity \((\omega/H)_0<7.6\times10^{-10}\) (95% CI) in related discussion.
- **Overlap with Chris’s program:** **Adversarial constraint** on globally anisotropic (homogeneous) expansion. Chris’s model is **inhomogeneous** (spherical CoE + off-center observer), not Bianchi — so Saadeh does **not** directly kill ST+Gauss. But it shows CMB polarization+temperature leave almost no room for universe-wide shear; any Chris anisotropy must be **observer-geometric / inhomogeneous**, not Bianchi-like global anisotropy, and must still respect CMB.
- **Tension / adversary note:** Strongest modern statement that FLRW isotropy of expansion is an excellent CMB fit among homogeneous models. Parent should not claim “anisotropic expansion is open season” without distinguishing Bianchi vs LTB/ST-off-center.
- **Quotable equations or phrases:** “anisotropic expansion of the Universe is strongly disfavoured, with odds of 121,000:1 against”; \((\sigma_V/H)_0<4.7\times10^{-11}\) (95% CI).

---

## Paywalled stubs (no PDF here) — UNKNOWN DETAILS

### Pathria 1972, *Nature* **240**, 298
- **Known from title / secondary citations only:** Classic “universe as a black hole” note identifying closed FRW maximum radius with Schwarzschild horizon for suitable \(\Lambda\).
- **UNKNOWN DETAILS:** Full matching assumptions, equations, and numerical claims not verified from primary PDF here.
- **Chris relevance:** Ancestral slogan for beyond-horizon CoE; see Khakshournia 2014 for modern junction critique.

### Debney & Farnsworth 1981, *Nuovo Cimento B* **64**, 367
- **Known from title only (in this workspace):** Debney–Farnsworth collaboration in the cosmological / global-analysis orbit (Fennelly 1977 thanks Debney & Farnsworth).
- **UNKNOWN DETAILS:** Contents, equations, observational claims — not available.

### Debney & Farnsworth 1983, *Acta Cosmologica* **12**, 41
- **Known from title only:** Same lineage; likely further cosmological BH / global structure discussion.
- **UNKNOWN DETAILS:** Contents not available here.

---

## Synthesis for Chris

### Closest cousins ranked
1. **Bertschinger 1983** — actual **Sedov → cosmological self-similar blast** in EdS (dynamics DNA of ST program).
2. **Smoller & Temple 2002** (+1998) — **beyond-Hubble shock / Big Bang inside a BH**; Copernican+entropy honesty.
3. **Migkas et al. 2020** — **empirical \(A(\ell)\) / \(D_L\) anisotropy** to match.
4. **Kodama–Saito–Ishibashi 2010** (general spherical CMB formulae) + **Alnes & Amarzguioui 2006** — off-center **CMB dipole ceiling** (~15–16 Mpc in void models).
5. **Heinesen 2020** — model-independent **anisotropic \(d_L\) multipoles** language for forecasts.
6. **Blomqvist & Mörtsell 2010** — SN vs CMB ranking of off-center constraints (SN weak, CMB strong).
7. **Fennelly 1977** — historical Hubble-anisotropy + off-center Gaussian inhomogeneity.
8. **Pathria note / Debney 1977** — BH-universe matching caveats and curiosities.
9. **Saadeh et al. 2016** — Bianchi shear nearly excluded (adversary for *global* anisotropy, not identical to Chris’s inhomogeneous CoE).

### What nobody has done (gap Chris claims)
No paper in this set combines:
- **Sedov–Taylor / Bertschinger-like** blast-wave expansion with
- **Gauss enclosed-mass** gravity \(g\propto M(<r)/r^2\) and decreasing \(\rho\),
- a **real CoE beyond the horizon** (Smoller–Temple-like),
- an **off-center observer** predicting
- **Migkas \(A\) vs Galactic longitude** (and related late-universe anisotropy),
preferring **time-integrated gravity** over a bolt-on peculiar kick,
while staying honest about **CMB dipole / Copernican / Saadeh** tensions.

LTB void industry (Alnes, Kodama, Blomqvist, …) got off-center CMB+SN right but used dust voids, not ST+Gauss. Smoller–Temple got explosion+BH right but not Migkas phenomenology. Migkas got the sky pattern without a CoE mechanism. Bertschinger got ST-in-EdS without cosmology-as-blast-for-observers.

### Must-read shortlist of 5 for next session
1. **Bertschinger 1983** — ST cosmological self-similarity (mechanics).
2. **Smoller & Temple 2002** — beyond-horizon / BH explosion cosmology (ontology).
3. **Migkas et al. 2020** — abstract + method + Fig. 6 \(A(\ell)\) (target data).
4. **Kodama–Saito–Ishibashi 2010** (`yoo_etal_2010_offcenter_cmb.pdf`) — general spherical off-center CMB dipole formulae (adversary toolkit).
5. **Heinesen 2020** — anisotropic \(d_L\) multipole Hubble law (forecast language).

**Honorable immediate add-ons:** Alnes & Amarzguioui 2006 (numerical dipole lesson); Blomqvist & Mörtsell 2010 (SN vs CMB); Saadeh et al. 2016 (don’t confuse with Bianchi); Khakshournia 2014 (Pathria matching honesty).

### Data hygiene for parent
- Replace misfiled `alnes_amarzguioui_2006_sn.pdf` and `smoller_temple_2003_shock_bh.pdf`.
- Rename `fennelly_nasa_cosmological_bh.pdf` → Debney 1977 NASA TM in catalogs.
- Rename `yoo_etal_2010_offcenter_cmb.pdf` → Kodama–Saito–Ishibashi 2010 in catalogs.



## Hygiene fix (2026-09-19 evening)
- Replaced `alnes_amarzguioui_2006_sn.pdf` with Alnes & Amarzguioui astro-ph/0610331 (SN Hubble diagram for off-center observers).
- Renamed Debney NASA TM and Kodama–Saito–Ishibashi files.
- Removed misfiled Nipoti et al.; Smoller–Temple 2003 follow-up not uniquely pinned on arXiv (1998+2002 remain core).
- Alnes SN gist: SNIa alone allow offsets out to ~225 Mpc with only marginal χ² improvement; do not tighten the ~15 Mpc CMB dipole bound.
