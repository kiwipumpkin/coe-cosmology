# REMADE (not original): restack the 3x2 grid f19mc3b.png (from mc_br_plot.py, MC web #3) into one column.
from PIL import Image
g=Image.open('/workspace/coepdf/fig/f19mc3b.png').convert('RGB');W,H=g.size;w,h=W//3,H//2
tiles=[g.crop((c*w,r*h,(c+1)*w,(r+1)*h)) for r in range(2) for c in range(3)]
out=Image.new('RGB',(w,h*6),'white')
for i,t in enumerate(tiles):out.paste(t,(0,i*h))
out.save('f19mc3b_stack_remade.png')
