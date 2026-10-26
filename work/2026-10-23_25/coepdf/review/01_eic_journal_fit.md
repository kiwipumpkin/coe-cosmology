# Peer Review Report — Journal-Fit Reviewer (EIC)

## Manuscript Information
Off-centre CoE cosmology: the mathematical workflow from geometry to the Planck TT spectrum (rev. 5). Reviewed as a submission to a PRD/JCAP-type journal.

## Reviewer Information
- **Role:** Journal-Fit Reviewer (card #1). **Identity:** editor handling non-standard cosmologies. **Focus:** originality, significance, framing, fit.

## Overall Assessment
- **Recommendation:** Major Revision (not publishable as a research article in current form; acceptable as an arXiv working note after reframing).
- **Confidence:** 4/5. **Calibration:** NOT_CALIBRATED.
- **Summary:** An honest, unusually self-critical working document that walks from an off-centre expanding-matter geometry to fits of the Migkas cluster anisotropy and the Planck TT spectrum. Its candour (self-flags, "the question is posed, not answered", CMB fit worse than ΛCDM) is a strength. But it is structured as lab notes, not a paper: it has no stated falsifiable claim, the central motivation (apparent acceleration from SNe) is never computed, and it does not engage the large literature on off-centre LTB observers that already constrains exactly this set-up.

## Strengths
- **S1 Transparency.** Every number is tied to a file; unverified items are flagged in-line (`\flag`).
- **S2 Honest CMB result.** States χ²=3,787 vs 2,558 for ΛCDM and that the half-scale mode has no physical origin (§Discussion).
- **S3 Clear derivations** of standard tools (Eqs. 6–10, 15–20) — useful pedagogically.

## Weaknesses
- **W1 [MAJOR] No engagement with prior off-centre/LTB literature.** Location: §1, References. Problem: Alnes & Amarzguioui (2006), Foreman et al. (2010), Clifton–Ferreira–Land (2008), Zibin–Moss–Scott kSZ tests (2008/2011), Bull–Clifton–Ferreira (2012), Caldwell & Stebbins (2008) all test off-centre/void models against the CMB. Fix: add a related-work section positioning the CoE model against these constraints.
- **W2 [MAJOR] Motivating claim untested.** Location: Abstract, §2.1, §5 next steps (5). Problem: the paper is motivated by apparent acceleration from SNe distances, yet D_L(z) is never computed. Fix: compute D_L(z, n̂) along off-centre lines of sight and compare with Pantheon+/Union3, or drop the motivation from the abstract.
- **W3 [MAJOR] Framing of the cluster result.** Location: Abstract, §4.1, Conclusions bullet 1. Problem: "χ² falls from 45.14 to about 10–12" is presented as the strongest point, but the bins are overlapping sliding cones; no significance follows. Fix: present with a proper likelihood (see methodology).
- **W4 [MINOR] Reference list incomplete** (no years/journal/volume for Ellis & Stoeger, Migkas, Planck, Eisenstein & Hu, Press–Schechter). Fix: full citations (Migkas et al. 2020 A&A 636 A15; Secrest et al. 2021 ApJL 908 L51; Planck 2018 VI; etc.).
- **W5 [MINOR] Title/Author.** "Working notes for C. A." and figures "restored from the original preprint" signal a draft; a journal needs a self-contained figure set with reproducible code.

## Detailed Comments
- **Title & Abstract:** the abstract reports 7 data numbers but no conclusion about viability; state plainly that the CMB stage disfavours the model at Δχ²≈1,200.
- **Introduction:** good statement of question; missing related work.
- **Methods/Results:** too long on textbook derivations (D_ℓ, cosmic variance, Limber) relative to model-specific content; move to an appendix.
- **Discussion:** the "Tensions" paragraph is the paper's real conclusion and should lead.

## Questions for Authors
1. What observation would falsify the CoE model in your view?
2. How does the model evade the Alnes/Foreman off-centre bounds?
3. Why are two different sky orientations (+80.8° roll; 45° rotation) both used?

## Minor Issues
- ℓ used both for Galactic longitude and multipole (§4.7). Two Planck peak lists disagree (§4.5 vs Table 2).

## Criterion-Bound Judgements
| Criterion | Judgement |
|---|---|
| Originality | Moderate (off-centre matter geometry + cluster anisotropy combination is uncommon; CoE/LTB idea is not new) |
| Significance | Low as presented (no viability claim survives the CMB stage) |
| Journal fit | Poor as a research article; fine as arXiv notes |
| Clarity | Good locally, poor globally (no throughline) |
