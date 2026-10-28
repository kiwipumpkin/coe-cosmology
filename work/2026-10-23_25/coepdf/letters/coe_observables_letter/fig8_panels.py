# Copies the book's original figures unchanged (no resampling) into this folder.
# Sky-orientation sequence: made by make_fig05_sky_orientation.py (book Fig 10 / original Fig 5).
# Migkas Fig 8 MAX/MIN overlays: made by make_migkas_fig8_overlay_matched.py (tilt; tilt+grad; tilt+grad+80.8 roll).
import os,shutil;H=os.path.dirname(os.path.abspath(__file__))
from PIL import Image
Image.open('/workspace/coepdf/fig/Fig05_sky_orientation_four_row.jpg').save(H+'/sky_orientation_sequence.png')  # lossless re-encode, same pixels
S='/workspace/cosmology/papers/spreadsheet_model/'
for src,dst in [('fig8_overlay_matched.png','migkas_fig8_tilt.png'),('fig8_overlay_matched_grad.png','migkas_fig8_tilt_grad.png'),('fig8_overlay_matched_grad_rot81.png','migkas_fig8_tilt_grad_roll81.png')]:
  shutil.copy(S+src,H+'/'+dst)
