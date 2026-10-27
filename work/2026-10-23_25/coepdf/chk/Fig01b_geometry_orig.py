"""Fig01b_geometry_orig.png -- generating script (extraction).
This image is the original raster embedded in the preprint Word file
/workspace/cosmo_preprint/coe_shell_preprint-v5.docx as word/media/rId26.png (byte-identical, md5
f4139d461db6d61c6dbafc809fed1d51; also in chain/coe/coe_shell_preprint-v5.docx). No plotting code for it
exists on the box, in git, or on GitHub; the matplotlib redraws (coepdf/scripts/fig01b_geometry.py etc.)
trace this 481x371 image. This script extracts it.
Run: python3 Fig01b_geometry_orig.py -> Fig01b_geometry_orig_regen.png"""
import os, zipfile
H = os.path.dirname(os.path.abspath(__file__))
data = zipfile.ZipFile('/workspace/cosmo_preprint/coe_shell_preprint-v5.docx').read('word/media/rId26.png')
open(os.path.join(H, 'Fig01b_geometry_orig_regen.png'), 'wb').write(data)
