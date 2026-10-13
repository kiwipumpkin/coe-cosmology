# Minimal-edit test (first pass, 9 Oct 2026)
Webs: B random and B foam, seeds 1-3, fill = our sphere + shell to 30 Gpc incl. core cap (tidalcore). Target = SR + central mass + web pattern.
No new bodies. Edited bodies end inside 30 Gpc, >=0.5 Gpc from us, surface >=3 Gpc from the ring, no overlaps; voids R 0.5-4 Gpc, clusters within x10 of grown mass.
Method: leverage ranking (linearised one-body gain), greedy batches + joint refit (full joint refit up to 150 bodies), stop at <=1x or plateau (<1% gain per 300 edits).
Variants: "frozen" = bodies grown within 3 Gpc of the ring cannot be edited; "movable" = they can, but must end >=3 Gpc from the ring.
No run reached 1x. Large scales m=1-12 fall to ~0.1-2x; leftover is almost all m>=13 (fine-scale texture from bodies near the ring).
See minedit/summary.json, minedit_curve.png, minedit_changes.png, minedit_side_views.png.
