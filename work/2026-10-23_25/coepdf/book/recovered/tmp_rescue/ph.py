import sys,os,glob,numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
def h(p):
    try:
        im=Image.open(p).convert('L').resize((48,48))
        a=np.asarray(im,float); return (a-a.mean())/(a.std()+1e-6), Image.open(p).size
    except Exception: return None,None
T="f06 f02 f12 f11 f03 f17 f18 f13 f15 sky_angle sky_mapA_al sky_mapB_al sky_tab_only Fig06_overlay_one f19mc3b_stack Fig01a_geometry_top_test".split()
tg={}
for t in T:
    p=glob.glob('/workspace/coepdf/fig/'+t+'.*')[0]; tg[t]=(p,)+h(p)
cands=[l.strip() for l in open(sys.argv[1])]
res={t:[] for t in T}
for c in cands:
    v,sz=h(c)
    if v is None: continue
    for t,(p,tv,tsz) in tg.items():
        if c==p: continue
        r=(v*tv).mean()
        res[t].append((r,c,sz,tsz))
for t in T:
    print('##',t,tg[t][2])
    for r,c,sz,tsz in sorted(res[t],reverse=True)[:4]: print('  %.3f'%r,c,sz)
