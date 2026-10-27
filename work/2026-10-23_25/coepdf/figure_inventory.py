"""figure_inventory.png -- REMADE generating script (the original renderer was not found; the table
content came from /tmp/build.py, which wrote figure_inventory.md). This script renders the table in
figure_inventory.md as a matplotlib table (wrapped cells, blue header row), in the same style as the original.
Run: python3 figure_inventory.py -> figure_inventory_regen.png (original untouched)."""
import os, textwrap
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
lines = [l for l in open(os.path.join(H, 'figure_inventory.md')).read().splitlines() if l.startswith('|')]
rows = [[c.strip() for c in l.strip('|').split('|')] for l in lines if not l.startswith('|---')]
hdr, body = rows[0], rows[1:]
widths = [3, 6, 22, 34, 30, 16, 38]   # wrap widths (chars)
cw = [w / sum(widths) for w in widths]
wrap = lambda r: ['\n'.join(textwrap.wrap(c, w)) or c for c, w in zip(r, widths)]
cells = [wrap(r) for r in body]; hdrw = ['\n'.join(textwrap.wrap(c, w)) for c, w in zip(hdr, widths)]
nl = [max(c.count('\n') + 1 for c in r) for r in cells]
fig = plt.figure(figsize=(24.2, 20.17)); ax = fig.add_axes([0.01, 0.01, 0.98, 0.98]); ax.axis('off')
t = ax.table(cellText=cells, colLabels=hdrw, colWidths=cw, loc='upper center', cellLoc='left')
t.auto_set_font_size(False); t.set_fontsize(9.5)
unit = 1 / (sum(nl) + 2 + len(nl) * 1.3)
for (i, j), c in t.get_celld().items():
    c.set_height(unit * (2 if i == 0 else nl[i - 1] + 1)); c.set_linewidth(0.5)
    if i == 0: c.set_facecolor('#dde6f3'); c.get_text().set_weight('bold'); c.get_text().set_ha('center')
fig.savefig(os.path.join(H, 'figure_inventory_regen.png'), dpi=100, facecolor='white')
