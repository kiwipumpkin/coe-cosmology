import sys,numpy as np
from PIL import Image
def load(p,size=None):
    im=Image.open(p).convert('L')
    if size: im=im.resize(size,Image.LANCZOS)
    return np.asarray(im,float)/255
def sim(a,b):
    A=load(a); B=load(b,(A.shape[1],A.shape[0]))
    try:
        from skimage.metrics import structural_similarity as ss
        s=ss(A,B,data_range=1)
    except Exception: s=float('nan')
    return s, 1-np.abs(A-B).mean(), Image.open(a).size, Image.open(b).size
if __name__=='__main__':
    o=sys.argv[1]
    r=[(sim(o,c),c) for c in sys.argv[2:]]
    for (s,m,so,sc),c in sorted(r,key=lambda x:-x[0][0])[:3]: print(f'{s:.3f} {m:.3f} {so} {sc} {c}')
