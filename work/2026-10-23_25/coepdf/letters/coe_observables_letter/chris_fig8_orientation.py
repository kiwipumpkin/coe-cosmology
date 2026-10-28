# Book sky-orientation sequence (left) beside the stacked Migkas Fig 8 MAX/MIN overlays (right), as in the first draft.
import os;from PIL import Image;H=os.path.dirname(os.path.abspath(__file__))
L=Image.open('/workspace/coepdf/fig/Fig05_sky_orientation_four_row.jpg').convert('RGB')
R=Image.open('/workspace/coepdf/fig/Fig06_overlay_one.png').convert('RGB')
L=L.resize((round(L.width*R.height/L.height),R.height),Image.LANCZOS)
C=Image.new('RGB',(L.width+R.width,R.height),'white');C.paste(L,(0,0));C.paste(R,(L.width,0))
C.save(H+'/chris_fig8_orientation.png')
