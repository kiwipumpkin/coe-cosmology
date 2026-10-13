import pymupdf as f, re
P=[('CoE cosmology: Apparent acceleration from an off-centre finite geometry','/workspace/cosmo_preprint/coe_shell_preprint.pdf'),
   ('Reading the CMB fingerprint of the CoE shell model','/workspace/clump_scan/report/brief/main.pdf'),
   ('Pattern formation (v5)','/workspace/pattern_formation/v5/pattern_formation_v5.pdf')]
NTOC=1; start=2+NTOC; ents=[]
def esc(s): return s.replace('&','\\&').replace('_','\\_').replace('#','\\#').replace('%','\\%')
for i,(t,p) in enumerate(P):
    d=f.open(p); toc=d.get_toc(); lvl=2 if i==1 else 1
    secs=[(x[1],x[2]) for x in toc if x[0]==lvl]
    if i==0: secs+= [('References',41),('Appendix figures',42),('Acknowledgments',46)]
    if i==2: secs.insert(-1,('References',31))
    if i==1: secs=[s for s in secs if s[0]!='Abstract']
    ents.append((0,f'{i+1}. '+t,start)); n=0
    for s,pg in secs:
        s=s.replace('A/Aall','$A/A_{\\rm all}$').replace('A(z,)','$A(z,\\theta)$')
        ents.append((1,s,start+pg-1))
    start+=d.page_count
ents.append((0,'4. Does everything have a centre? A philosophical note',start))
rows='\n'.join((r'\noindent\textbf{%s}\dotfill\textbf{%d}\par\smallskip' if l==0 else r'\noindent\hspace*{1.2em}%s\dotfill %d\par') % (t if '$' in t else esc(t),pg) for l,t,pg in ents)
s=open('c.tex.bak_before_toc').read()
a=s.index(r'\noindent{\large\bfseries Contents}'); b=s.index(r'\clearpage',a)
toc=r'''\noindent\small Full table of contents on the next page. Page numbers refer to this combined PDF; each paper keeps its own original page numbering.
\clearpage
\pdfbookmark[0]{Contents}{toc}
{\noindent\Large\bfseries Contents\par}\medskip
{\small
'''+rows+'\n}\n'
s=s[:a]+toc+s[b:]
open('c.tex','w').write(s); open('ents.txt','w').write(repr(ents))
