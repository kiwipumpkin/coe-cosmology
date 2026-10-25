import subprocess, numpy as np
from PIL import Image
for f, out in (("termmatch_A_v2", "termmatch_A_equations_v2"), ("termmatch_B_v2", "termmatch_B_table_v2")):
    subprocess.run(["pdftoppm", "-r", "220", "-png", "-singlefile", f + ".pdf", "/tmp/" + f], check=True)
    im = Image.open("/tmp/" + f + ".png").convert("RGB"); a = np.asarray(im.convert("L")) < 245
    rows = np.where(a.any(1))[0]; y0, y1 = max(rows[0] - 40, 0), min(rows[-1] + 40, im.size[1])
    im = im.crop((0, y0, im.size[0], y1)); im.save(out + ".jpg", quality=92); print(out, im.size)
    if out.endswith("table_v2"):
        w = 1400; im.resize((w, round(im.size[1] * w / im.size[0])), Image.LANCZOS).save(out + "_web.jpg", quality=90)
