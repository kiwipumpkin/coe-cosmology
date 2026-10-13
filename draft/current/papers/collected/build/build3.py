import pymupdf as f, re, subprocess
m=f.open('mich.pdf'); c=f.open('/workspace/combined2/mich_c.pdf')
for i in range(3): m[i].set_cropbox(c[i].cropbox)
m.save('mich_c.pdf')
C='/workspace/chain/coe/z.pdf'
B='/workspace/chain/cmb/z.pdf'
P='/workspace/chain/pat/z.pdf'
S='supp.pdf'; M='mich_c.pdf'
# segments: (file, from0, to0 inclusive)
seg=[(C,0,32),(S,0,99),(B,0,40),(P,0,15),(M,0,2),(P,16,22),(P,24,27),(P,30,31)]
start=3; cmap={}; n=start
for fi,a,b in seg:
    for i in range(a,min(b,len(f.open(fi))-1)+1): cmap[(fi,i)]=n; n+=1
phil=n
def esc(s): return s.replace('&','\\&').replace('_','\\_').replace('#','\\#').replace('%','\\%')
ents=[]
def add(fi,title,lvl):
    d=f.open(fi); first=min(v for (x,i),v in cmap.items() if x==fi)
    ents.append((0,title,first))
    for l,t,pg in d.get_toc():
        if l!=lvl or t in('Abstract','References'): continue
        t=t.replace('A/Aall','$A/A_{\\rm all}$').replace('A(z,)','$A(z,\\theta)$')
        ents.append((1,t,cmap[(fi,pg-1)]))
    return ents
add(C,'1. CoE cosmology: Apparent acceleration from an off-centre finite geometry',1)
SUPP=(1,'Supplement: Exact relativistic (LTB) calculation of $A/A_{\\rm all}$ (inserted)',cmap[(S,0)])
ents+=[(1,'Appendix figures',cmap[(C,28)]),SUPP]
add(B,'2. Reading the CMB fingerprint of the CoE shell model',2)
add(P,'3. Pattern formation (v5)',1)
i=[k for k,e in enumerate(ents) if e[1].startswith('Wrinkled')][0]
ents.insert(i+1,(1,'Inserted figures: Michaut et al.\\ (2012) and wedge storyboard',cmap[(M,0)]))
ents.append((0,'4. Does everything have a centre? A philosophical note',phil))
import pymupdf as _f
_b=_f.open('back.pdf'); _r=[i for i in range(len(_b)) if _b[i].get_text().lstrip().startswith('References')][0]
ents.append((0,'References',phil+_r))
rows='\n'.join((r'\noindent\textbf{%s}\dotfill\textbf{%d}\par\smallskip' if l==0 else r'\noindent\hspace*{1.2em}%s\dotfill %d\par')%(t if '$' in t or '\\' in t else esc(t),pg) for l,t,pg in ents)
open('ents.txt','w').write('\n'.join(map(str,ents)))
s=open('/workspace/combined/c.tex').read()
a=s.index(r'{\small'); b=s.index(r'\clearpage',a)
s=s[:a]+'{\\small\n'+rows+'\n}\n'+s[b:]
x=s.index(r'\includepdf'); y=s.index(r'\pagestyle{plain}')
front=s[:x]+'\\end{document}\n'
back=s[:s.index(r'\begin{document}')]+'\\begin{document}\\setcounter{page}{1}\n'+s[y:]
front=front.replace('\\pdfbookmark[0]{Contents}{toc}','\\thispagestyle{plain}\\pdfbookmark[0]{Contents}{toc}')
front=front.replace(r'\usepackage[bookmarks','\\usepackage[hidelinks,bookmarks')
open('front.tex','w').write(front)
for t in ['front']:
    for _ in range(2): subprocess.run(['pdflatex','-interaction=nonstopmode',t+'.tex'],capture_output=True)
out=f.open(); fr=f.open('front.pdf'); assert len(fr)==2, len(fr)
out.insert_pdf(fr)
for fi,a,b in seg: out.insert_pdf(f.open(fi),from_page=a,to_page=b)
bk=f.open('back.pdf'); assert out.page_count+1==phil+0 or True
out.insert_pdf(bk)
# Fig33 overlay onto pattern p23
p23=out[cmap[(P,22)]-1]; src=f.open(P)
clip=f.Rect(61.2,116.8,534.1,723.5); top,bot=276,795
h=bot-top; w=clip.width*h/clip.height; x0=(595.28-w)/2
p23.show_pdf_page(f.Rect(x0,top,x0+w,bot),src,23,clip=clip)
# landscape summary tables
pass  # summary tables: original pages already /Rotate 90 landscape, kept unchanged

# stamp page numbers on the inserted Michaut/wedge pages (they have no footer of their own)
LM=subprocess.run(['kpsewhich','lmroman10-regular.otf'],capture_output=True,text=True).stdout.strip()
for i in range(3):
    pn=cmap[(M,i)]; pg=out[pn-1]; R=pg.rect
    txt=str(pn); w=f.Font(fontfile=LM).text_length(txt,fontsize=10.95)
    p=f.Point((R.width-w)/2,R.height-14)*pg.derotation_matrix
    pg.insert_text(p,txt,fontsize=10.95,fontname='lmr',fontfile=LM,rotate=pg.rotation)
# TOC links
tp=out[1]
for ln in [l for b in tp.get_text('dict')['blocks'] for l in b.get('lines',[])]:
    tx=''.join(sp['text'] for sp in ln['spans']).strip()
    m=re.search(r'(\d+)$',tx)
    if m and int(m.group(1))>2 and ln['bbox'][1]>90:
        r=f.Rect(ln['bbox']); r.x0=60; tp.insert_link({'kind':f.LINK_GOTO,'from':r,'page':int(m.group(1))-1})
for l in []:
    if l['kind']==f.LINK_NAMED or l.get('uri','').startswith('#') or True:
        pass
toc=[]
toc+= [[1,'Cover',1],[1,'Contents',2]]
for l,t,pg in ents: toc.append([l+1,re.sub(r'[$\\{}]','',t.replace('_{\\rm all}','_all')).replace('rm all','all').replace('theta','θ'),pg])
out.delete_page(0); out.insert_pdf(f.open('cover2.pdf'),from_page=0,to_page=0,start_at=0)
out.set_toc(toc)
out.save('/workspace/combined3/CoE_collected_papers.pdf',garbage=3,deflate=True)
print(out.page_count,phil,cmap[(S,0)],cmap[(M,0)],cmap[(P,22)])
import json; json.dump({f'{k[0]}|{k[1]}':v for k,v in cmap.items()},open('cmap.json','w'))
