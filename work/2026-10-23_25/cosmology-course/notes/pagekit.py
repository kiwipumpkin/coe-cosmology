"""Shared page layout (wrapping, rendering) for the tilt-card note pages."""
import re, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})
EQ, NT, HD, TT = 18, 15, 17, 20
B = r"$\bullet$  "
PREFIX = "page"
WMAX = 9.25; IND = 0.55
fig0 = plt.figure(figsize=(30, 3)); R0 = fig0.canvas.get_renderer()
def width(t, sz):
    tx = fig0.text(0, .5, t, fontsize=sz); b = tx.get_window_extent(R0); tx.remove(); return b.width/fig0.dpi
def balance(body, carry):
    """Convert @L( / @R) markers; close/open across lines with \\left. / \\right."""
    toks = re.findall(r"@L(\\\{|\(|\[)|@R(\\\}|\)|\])", body)
    opens = carry
    for o, c in toks:
        opens += 1 if o else -1
    pre = r"\left." * carry
    body = re.sub(r"@L(\\\{|\(|\[)", r"\\left\1", body); body = re.sub(r"@R(\\\}|\)|\])", r"\\right\1", body)
    post = r"\right." * opens
    return pre + body + post, opens
def wrap(lhs, chunks, sz=EQ):
    lines, cur, carry, lvl = [], lhs, 0, 0
    def txt(body, c): return "$" + balance(body, c)[0] + "$"
    for ch, cont in chunks:
        trial = cur + " " + ch
        if cont != "!" and width(txt(trial, carry), sz) + lvl*IND <= WMAX:
            cur = trial
        else:
            t, carry_new = balance(cur, carry)
            lines.append(("\x01"*lvl + "$" + t + "$", sz)); carry = carry_new
            lvl = 1
            cur = (cont + " " if cont and cont != "!" else "") + ch
    t, _ = balance(cur, carry); lines.append(("\x01"*lvl + "$" + t + "$", sz))
    return lines
def render(n, items, dpi=250, colsfs=None):
    fig = plt.figure(figsize=(10, 40)); r = fig.canvas.get_renderer()
    H0 = 40; y = H0 - 0.35; x0 = 0.35; widest = 0; extra = []
    flat = []
    for it in items:
        if it[0] == "eq": flat += wrap(it[1], it[2], it[3] if len(it) > 3 else EQ)
        else: flat.append(it)
    for it in flat:
        if it[0] == "gap": y -= it[1]; continue
        if it[0] == "table":
            rows, cols, fs = it[1], it[2], it[3]
            for j, row in enumerate(rows):
                for c, t in zip(cols, row):
                    fig.text((x0 + c)/10, y/H0, t, fontsize=fs, va="top", fontweight="bold" if j == 0 else "normal")
                y -= 0.36 if j else 0.42
                if j == 0: extra.append(y + 0.05)
            continue
        t, sz = it
        ind = IND*(len(t) - len(t.lstrip("\x01"))); t = t.lstrip("\x01")
        tx = fig.text((x0 + ind)/10, y/H0, t, fontsize=sz, va="top")
        bb = tx.get_window_extent(r); wd = bb.width/fig.dpi + ind; h = bb.height/fig.dpi
        widest = max(widest, wd)
        if wd > 9.4: print(f"  p{n} TOO WIDE {wd:.2f}: {t[:80]}")
        y -= h + (0.26 if sz >= EQ else 0.13)
    used = H0 - y + 0.25
    for yy in extra:
        fig.add_artist(plt.Line2D([x0/10, 9.65/10], [yy/H0]*2, color="k", lw=0.8))
    fig.set_size_inches(10, used, forward=False)
    for t in fig.texts: t.set_y(1 - (1 - t.get_position()[1])*H0/used)
    for l in fig.artists: l.set_ydata([1 - (1 - v)*H0/used for v in l.get_ydata()])
    fig.savefig(f"{PREFIX}-p{n}.png", dpi=dpi); plt.close(fig)
    print(f"p{n}: 10 x {used:.2f} in, widest {widest:.2f} in" + ("   <-- TALL" if used > 12.2 else ""))
    return used


def note(text, sz=NT, bullet=True, maxw=9.2):
    """Greedy word-wrap of a note; $...$ groups are never split. Continuation lines indented."""
    toks, cur_t, inm = [], "", False
    for ch in text:                      # split on spaces outside $...$
        if ch == "$": inm = not inm
        if ch == " " and not inm:
            if cur_t: toks.append(cur_t)
            cur_t = ""
        else: cur_t += ch
    if cur_t: toks.append(cur_t)
    # merge tokens that split a $..$ group with spaces inside
    out, cur = [], (B if bullet else "")
    first = True
    for t in toks:
        trial = (cur + " " + t) if cur.strip() and cur != B else cur + t
        if width(trial, sz) > maxw and cur.strip() not in ("", B.strip()):
            out.append((cur, sz)); cur = r"$\qquad$ " + t
        else:
            cur = trial
    out.append((cur, sz))
    return out
