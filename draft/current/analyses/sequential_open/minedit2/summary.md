# Minimal edit with move / rescale / delete (9-10 Oct 2026)
Webs: B random s1-3 and B foam s1-3 (our sphere + shell to 30 Gpc incl. core cap). Target: SR + central mass + the web. No additions.
Edits: move, rescale (mass x0-10; void R >= 0.5, cluster >= 0.1 M0, below that = deleted), delete. Moved/kept bodies end >= 3 Gpc from the ring,
inside the region, no overlaps (checked independently: check_layout.py; all runs clean except one 0.12 Gpc overlap in edits-only random s2).
Fitter (minedit2.py): leverage ranking (best of linearised move/resize gain and exact deletion gain), joint bounded LSQ of each batch from 4 starts
(as found / deleted / two random relocations), joint re-fit of the last 60 edited bodies every 5 rounds, keep only improving rounds, stop at <=1x or plateau.
Two families: (a) edits only; (b) first delete every body whose surface is within D Gpc of the ring, then edit (D = 1-10).
Result: (a) never reaches 1x (43-137x, all at m>=13; m1-12 ~0.5-1.4x). (b) reaches <=1x for every web from D = 3 (the ring band itself):
  ~460-880 edits (2.7-5.2% of the web), of which 225-290 are deletions of the near-ring bodies. Larger D needs more deletions (D=10: ~4,300).
Earlier move/resize-only pass (minedit/) had an overlap bug (edited bodies checked against old positions); superseded. Invalid v2 outputs in minedit2/invalid_overlap/.
Figures: minedit2_curve.png, minedit2_changes.png, minedit2_changes_cleared.png, minedit2_side_views.png. Numbers: minedit2/summary.json.
D = 1 and D = 2 (finished): no web reaches 1x. D=1: 13-37x (1,480-2,350 edits); D=2: 2.3-4.7x (1,050-1,750 edits). D = 3 is the smallest clearing that works.
Figure style (default for edit figures from 10 Oct): plain before/after pairs, no arrows/markers (edit_before_after.py); minedit2_side_views.png.

## No deletion (10 Oct, 01:45-02:15): every body kept; move + resize only (minedit2.py NODEL=1)
Near-ring bodies (surface < 3 Gpc from the ring, 221-276 per web) must end >= 3 Gpc away; void R >= 0.5 Gpc, clusters 0.1-10 x M0; no overlaps.
reloc = move them to the nearest legal spot first, joint-refit them, then minimal move/resize of others. forced = fit them out of the band from their own spots.
Reached <= 1x: foam s1 reloc (614 edits, 3.6%, 0.954x) and foam s3 reloc (489, 2.8%, 0.995x). Others plateau at 1.18-2.86x; m1-12 is ~0.01x everywhere,
the leftover is all m>=13. Relocate-first beats forced-out in every web. All layouts legal (check_layout). Figures: nodel_curve.png, nodel_side_views.png; numbers: summary_nodel.json.
