# explanation column for the stacked loop flow (exec'd from loop_flow2.py)
import textwrap
FS = 11.5; LH = FS * 1.32 / 72
def para(top_in, text, width=66):
    lines = []
    for p in text.split("\n\n"): lines += textwrap.wrap(p, width) + [""]
    lines = lines[:-1]
    fig.text(XT / FW, 1 - top_in / FH, "\n".join(lines), fontsize=FS, va="top", ha="left", linespacing=1.32, family="serif")
    return top_in + LH * len(lines) + 0.18
def table(top_in, cap, head, rows, colw, rowcols=None, h=0.27):
    n = len(rows) + 1; H_ = h * n
    fig.text(XT / FW, 1 - top_in / FH, cap, fontsize=10, va="top", ha="left", family="serif", style="italic")
    ax_ = fig.add_axes(R(XT, top_in + 0.24, sum(colw), H_)); ax_.axis("off")
    t = ax_.table(cellText=rows, colLabels=head, colWidths=[c / sum(colw) for c in colw], loc="upper left", cellLoc="center", bbox=[0, 0, 1, 1])
    t.auto_set_font_size(False); t.set_fontsize(10)
    for (r_, c_), cell in t.get_celld().items():
        cell.set_edgecolor("0.6"); cell.set_linewidth(0.6); cell.get_text().set_family("serif")
        if r_ == 0: cell.set_facecolor("#eeeeee"); cell.get_text().set_fontweight("bold")
        elif rowcols and c_ == 0 and rowcols[r_ - 1]: cell.get_text().set_color(rowcols[r_ - 1]); cell.get_text().set_fontweight("bold")
    return top_in + 0.24 + H_ + 0.25
def step_head(top_in, s):
    fig.text(XT / FW, 1 - top_in / FH, s, fontsize=12.5, va="top", ha="left", fontweight="bold", family="serif")
    return top_in + 0.34
byk = {L["key"]: L for L in LOOPS}
ORD = ["MAX1", "MIN1", "MAX2", "MIN2"]           # step-1 order (increasing l)
nm = lambda k: "%s-%s" % (byk[k]["kind"], k[-1])
rc = [COLS[k] for k in ORD]
ring_u = dict(zip(keys, UNT["ring"])); ring_b = dict(zip(keys, FIT["ring"])); ray_b = dict(zip(keys, FIT["ray"]))
# step 1
y = step_head(TOP[1] - 0.75, "Step 1. Gravity breaks homology: your model A")
y = para(y, "A pure homologous blast (V = H R) would look isotropic from anywhere, so the CoE would drop out. Gravity decelerates matter at different positions relative to the CoE by different amounts, so the expansion rate seen from Earth depends on "
         "direction. Your model A (thin-shell best fit) turns this into A($\ell$) ∝ H0² along the Galactic plane, fitted to Migkas Fig. 6 with the parameters in "
         "Table 1a. Its two maxima and two minima are read off the curve (Table 1b); each coloured band is ±9° about an extremum, and the same colours are "
         "used for the loops on every map below.")
y = table(y, "Table 1a. Model A fit to Migkas Fig. 6", ["Parameter", "Value"],
          [["$\ell_0$ (anti-CoE longitude)", "%.2f°" % l0], ["Amp", "%.3f" % F["amp"]], ["φ", "%.3f rad" % F["phi"]], ["R²", "%.3f" % F["R2"]]], [3.4, 1.6])
y = table(y, "Table 1b. Extrema of model A", ["Loop", "$\ell_{\\rm ext}$ (°)", "A/A_all"],
          [[nm(k), "%.1f" % byk[k]["l"], "%.3f" % [p for p in pairs if abs(p["l_model"] - byk[k]["l"]) < 1e-6][0]["A"]] for k in ORD], [1.6, 1.7, 1.7], rc)
# step 2
y = step_head(TOP[2] - 0.6, "Step 2. Each extremum becomes a loop (untilted)")
y = para(y, "Model A is axisymmetric about the CoE axis, so a value read at longitude $\ell_{\\rm ext}$ in the plane holds on the whole cone of half-angle "
         "θ = $\ell_{\\rm ext}$ − $\ell_0$ about the anti-CoE direction. Drawing each cone on the sky turns the four extrema into four full loops, here with the axis kept in "
         "the Galactic plane (b0 = 0). The ring miss is the signed angular distance from each Migkas extremum to its loop (Table 2). MAX-2 is a small ring "
         "5.4° round the CoE (see inset), and the two MIN loops have almost equal radii, so they are drawn 2° apart for display.")
y = table(y, "Table 2. Loop radii and ring misses, untilted", ["Loop", "θ (°)", "cone radius (°)", "Migkas ($\ell$, b)", "miss (°)"],
          [[nm(k), "%.1f" % byk[k]["theta"], ("%.1f" % (180 - byk[k]["theta"]) + " (CoE)") if k == "MAX2" else "%.1f" % min(byk[k]["theta"], 360 - byk[k]["theta"]),
            "(%.0f°, %+.0f°)" % tuple(byk[k]["mig"]), "%+.1f" % ring_u[k]] for k in ORD] + [["RMS", "", "", "", "%.1f" % UNT["rms_ring"]]],
          [0.95, 0.8, 1.55, 1.45, 0.95], rc + [None])
# step 3
y = step_head(TOP[3] - 0.6, "Step 3. Tilt the CoE axis out of the plane")
y = para(y, "The anti-CoE longitude is fixed by $\ell_0$, but the axis may point out of the Galactic plane. Tilting it to latitude b0 moves all four loops "
         "rigidly on the sphere; b0 is fitted by least squares on the ring misses, giving b0 = +%.1f° and anti-CoE (%.1f°, %+.1f°). All four Migkas "
         "extrema then fall inside the ±9° bands (Table 3). Freeing $\ell_0$ as well only reaches %.1f° RMS, so the remaining ~6° comes from the offset between "
         "the model and Migkas extremum longitudes, not from the orientation." % (FIT["b"], l0, FIT["b"], F["fits"]["free"]["rms_ring"]))
y = table(y, "Table 3. Ring misses before and after the tilt (°)", ["Loop", "b0 = 0", "b0 = %+.1f°" % FIT["b"]],
          [[nm(k), "%+.1f" % ring_u[k], "%+.1f" % ring_b[k]] for k in ORD] + [["RMS", "%.1f" % UNT["rms_ring"], "%.1f" % FIT["rms_ring"]]], [1.6, 1.7, 1.7], rc + [None])
# step 4
y = step_head(TOP[4] - 0.6, "Step 4. Roll about the CoE axis: the gradient")
y = para(y, "The tilt fixes where the loops lie but not their orientation about the axis. A roll about the CoE axis sets each loop's θ-ray point, the "
         "point where the loop meets the great circle through the axis that plays the role of step 1's slice; there the loop value is the one actually "
         "read from A($\ell$). Each loop is drawn darkest at its θ-ray point and fades to 12%% opacity opposite, so the dark end marks the measured direction and "
         "the faded part is the axisymmetric extrapolation. Because sin θ has opposite signs for MIN-1 and MIN-2, their dark ends lie on opposite sides. "
         "The roll, %.1f°, minimises the distance from the θ-ray points to the Migkas extrema (Table 4)." % FIT["roll"])
y = table(y, "Table 4. θ-ray miss after roll %.1f° (°)" % FIT["roll"], ["Loop", "θ-ray miss"],
          [[nm(k), "%.1f" % ray_b[k]] for k in ORD] + [["RMS", "%.1f" % FIT["rms_ray"]]], [1.6, 1.7], rc + [None])
# step 5
y = step_head(TOP[5] - 0.6, "Step 5. Compare with Migkas Fig. 8")
y = para(y, "Finally the tilted and rolled loops are laid on Migkas et al.'s Fig. 8 map at θ = 75°, with the same registration as your earlier "
         "overlays, so they can be compared directly with the observed A/A_all pattern (higher = faster, A ∝ H0²). The white line is the CoE axis great "
         "circle and the black line is Migkas's celestial equator. No parameter is fitted in this step; it is the check of steps 1 to 4 against the data.")
