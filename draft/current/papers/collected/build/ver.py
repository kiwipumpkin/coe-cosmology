import pymupdf as f, json, numpy as np
out=f.open('/workspace/combined/CoE_collected_papers.pdf'); cm=json.load(open('cmap.json')); srcs={}
bad=[];n=0
for k,v in cm.items():
    fi,i=k.rsplit('|',1); i=int(i)
    if fi in('supp.pdf','mich_c.pdf'): continue
    s=srcs.setdefault(fi,f.open(fi))[i]; o=out[v-1]
    a=np.frombuffer(s.get_pixmap(dpi=60).samples,np.uint8); b=np.frombuffer(o.get_pixmap(dpi=60).samples,np.uint8)
    n+=1
    if a.shape!=b.shape or (a!=b).any(): bad.append((fi.split('/')[-1],i+1,v, 'shape' if a.shape!=b.shape else int((a!=b).sum())))
    if s.rect!=o.rect or s.rotation!=o.rotation: bad.append(('size/rot',fi,i+1))
print(n,'compared; diffs:',bad)
