from PIL import Image, ImageDraw, ImageFont
import numpy as np
A = "/home/box/agent-data/agents/196ab413-0708-4f81-99fc-58a32c5366c5/attachments/"
P = "/workspace/cosmology/papers/spreadsheet_model/"
graph = Image.open(A + "f480c67c165d9485e7fd2b77af76a6260c8b8b7ccbb2b681e3a8d791bd706f5d.jpg").convert("RGB")
rows_src = [Image.open(f).convert("RGB") for f in (P + "loop_threepanel_notilt_web.jpg",
            A + "e2fab4ef4115016f7bc638a6ac1075749b6acbab6e4c2ca6b2f20e29f8735139.jpg",
            A + "fe905a57d5aaf4dede004a4cbdc50883bbf92c984372f9f260251bbfc057b757.jpg",
            A + "d6982109d332a918353df473def180356c133fb81196c90b6e0dea313238ef05.jpg")]
ovs = [Image.open(f).convert("RGB") for f in (A + "3b64499f8c49ab6537b8a73b1fd80965e8b5f8ae82d93d6c99b61dc4f5b8344d.jpg",
       A + "83d1ae71406713257f0623d2b732446ef1ee220de6167592d4a45c5cac55637c.jpg", P + "fig8_overlay_matched_grad_rot81_web.jpg")]
LAB = ["Untilted, no roll", "+15° tilt", "+15° tilt, gradient", "+15° tilt, gradient, +80.8° roll"]
OLAB = ["+15° tilt", "gradient", "+80.8° roll"]
FB = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 22)
FH = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 18)
FT = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 17)
W = 1700; LH = 40
def row_content(i):
    im = rows_src[i]
    if i == 0:
        return im.crop((0, 110, W, 628))            # keep A./B./C. headings + footnotes (first row only)
    c = im.crop((0, 222, W, 572)).copy()          # drops suptitle, A/C headings, footnotes, legend
    ImageDraw.Draw(c).rectangle((385, 234 - 222, 1212, 279 - 222), fill="white")   # B heading
    return c
legend = rows_src[3].crop((0, 712, W, 768))        # one legend strip (last row's)
def ov_crop(k, cbar):
    return ovs[k].crop((140, 355, 1445 if cbar else 1276, 1060))
def label(dr, x, y, t, f=FB, anchor="mm"):
    dr.text((x, y), t, font=f, fill="black", anchor=anchor)
def save(im, name):
    im.save(P + name + ".jpg", quality=92)
    w = 1400; im.resize((w, round(im.height * w / im.width)), Image.LANCZOS).save(P + name + "_web.jpg", quality=88)
    print(name, im.size)
# ---------------- 3-column version ----------------
parts = [("g", graph)]
blocks = []
H = graph.height + 20
for i in range(4):
    rc = row_content(i); blocks.append((i, rc)); H += LH + rc.height + 10
s3 = (W - 20 - 2 * 20) / (ov_crop(0, 0).width * 2 + ov_crop(2, 1).width)
oc = [ov_crop(0, 0), ov_crop(1, 0), ov_crop(2, 1)]
oc = [o.resize((round(o.width * s3), round(o.height * s3)), Image.LANCZOS) for o in oc]
OT = 26
H += LH + OT + max(o.height for o in oc) + 10 + legend.height + 10
im = Image.new("RGB", (W, H), "white"); dr = ImageDraw.Draw(im)
y = 0; im.paste(graph, ((W - graph.width) // 2, y)); y += graph.height + 20
for i, rc in blocks:
    label(dr, W // 2, y + LH // 2, LAB[i]); y += LH; im.paste(rc, (0, y)); y += rc.height + 10
label(dr, W // 2, y + LH // 2, "Migkas Fig. 8, θ = 75°, A/Aₐₗₗ"); y += LH
x = 10
for k, o in enumerate(oc):
    label(dr, x + round((1276 - 140) * s3) // 2, y + OT // 2, OLAB[k], FT)
    im.paste(o, (x, y + OT)); x += o.width + 20
y += OT + max(o.height for o in oc) + 10
im.paste(legend, (0, y)); y += legend.height + 10
save(im.crop((0, 0, W, y)), "stack_3col_clean")
# ---------------- 4-column version ----------------
s4 = 772 / 1304   # same overlay scale as stack_4col.jpg (no shrinking)
o4 = [ov_crop(0, 1), ov_crop(1, 0), ov_crop(2, 0)]
o4 = [o.resize((round(o.width * s4), round(o.height * s4)), Image.LANCZOS) for o in o4]
CW = max(o.width for o in o4) + 20; W4 = W + CW; HT = 30
Hs = [LH + blocks[0][1].height + 10]
for j in range(1, 4):
    rc = blocks[j][1]; Hs.append(max(LH + rc.height, (HT if j == 1 else 0) + o4[j - 1].height) + 12)
H4 = graph.height + 20 + sum(Hs) + legend.height + 10
im = Image.new("RGB", (W4, H4), "white"); dr = ImageDraw.Draw(im)
y = 0; im.paste(graph, ((W4 - graph.width) // 2, y)); y += graph.height + 20
for j in range(4):
    rc = blocks[j][1]
    label(dr, W // 2, y + LH // 2, LAB[j]); im.paste(rc, (0, y + LH))
    if j >= 1:
        oy = y
        if j == 1:
            label(dr, W + 10 + round((1276 - 140) * s4) // 2, y + HT // 2, "Migkas Fig. 8, θ = 75°, A/Aₐₗₗ", FH); oy += HT
        im.paste(o4[j - 1], (W + 10, oy))
    y += Hs[j]
im.paste(legend, (0, y)); y += legend.height + 10
save(im.crop((0, 0, W4, y)), "stack_4col_clean")
