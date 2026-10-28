# Clean version of fig8_migkas_maps.png (letter figure 5). Removes crop/screenshot blemishes only; data, colours, layout unchanged.
# Source: /workspace/coepdf/fig/Fig06_overlay_one.png (the stack of the three make_migkas_fig8_overlay_matched.py panels).
# Fixes: (1) take the full stack instead of split_fig8.py's x=3043 cut, which pulled in 123 px of the left figure's legend text
#        ("(±9°", "of 360°)"); (2) blank the strip under each "Galactic longitude" label, which holds the next sub-panel's
#        clipped "-90°" tick, "1.2" colourbar label and red dotted map edge from the cropped Migkas raster; (3) replace the
#        half-clipped "Galactic latitude" label with a full one; (4) smooth the saw-tooth HEALPix pixel edge of the Mollweide map
#        by whitening only a thin annulus just outside a fitted ellipse.
import numpy as np; from PIL import Image, ImageDraw, ImageFont
SRC='/workspace/coepdf/fig/Fig06_overlay_one.png'; OUT='fig8_migkas_maps_clean.png'
im=Image.open(SRC).convert('RGB'); a=np.asarray(im).copy(); H,W,_=a.shape
OFFS=[0,1061,2121]          # vertical offsets of the 3 sub-panels (measured row structure)
INSET=float(__import__('os').environ.get('INSET','4'))
for o in OFFS:
    a[o+930:o+962,:]=255                       # stray next-panel fragments under the longitude label
    p=a[o+100:o+930].astype(int); s=p.max(2)-p.min(2); m=s>60; m[:,1450:]=False
    # least-squares ellipse through left/right edge points of each row (tooth noise averages out)
    P=[]
    for y in range(m.shape[0]):
        xs=np.where(m[y])[0]
        if len(xs)>50: P+= [(xs.min(),y+o+100),(xs.max(),y+o+100)]
    for x in range(200,1450,4):
        ys=np.where(m[:,x])[0]
        if len(ys)>50: P+= [(x,ys.min()+o+100),(x,ys.max()+o+100)]
    P=np.array(P,float); X,Y=P[:,0],P[:,1]
    # axis-aligned fit: x^2/Rx^2... solve a x^2 + b y^2 + d x + e y = 1
    c=np.linalg.lstsq(np.c_[X**2,Y**2,X,Y],np.ones_like(X),rcond=None)[0]
    cx=-c[2]/(2*c[0]); cy=-c[3]/(2*c[1]); k=1+c[0]*cx**2+c[1]*cy**2
    Rx=np.sqrt(k/c[0])-INSET; Ry=np.sqrt(k/c[1])-INSET; print('ellipse',o,round(cx,1),round(cy,1),round(Rx,1),round(Ry,1))
    Yg,Xg=np.mgrid[o+100:o+930,0:1450]; r=np.sqrt(((Xg-cx)/Rx)**2+((Yg-cy)/Ry)**2)
    sub=a[o+100:o+930,0:1450]; m2=(s>12)[:,:1450]; sub[(r>1.0)&(r<1.06)&m2]=255   # whiten only coloured map pixels just outside the ellipse (labels untouched)
    a[o+360:o+612,150:200]=255                 # half-clipped latitude label
im=Image.fromarray(a)
try: F=ImageFont.truetype('DejaVuSans-Bold.ttf',30)
except Exception: F=ImageFont.load_default()
for o in OFFS:
    t=Image.new('RGBA',(300,40),(255,255,255,0)); ImageDraw.Draw(t).text((150,20),'Galactic latitude',font=F,fill='black',anchor='mm')
    t=t.rotate(90,expand=True); im.paste(t,(158,o+486-150),t)
im.save(OUT); print(OUT,im.size)
