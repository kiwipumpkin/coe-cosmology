import json,re
pat=json.load(open('pat_refs.json')); cmb=json.load(open('cmb_refs.json'))
cmbdoi={1:'10.1093/mnras/stv2215',2:'10.1086/591439',3:'10.1088/0004-637X/785/2/104',4:'10.1093/mnras/stw3250',5:'10.1051/0004-6361/201016288',6:'10.3847/1538-4357/835/2/161',7:'10.1103/PhysRevLett.112.251302',8:'10.1093/mnras/stv488',9:'10.1093/mnras/stw1752',10:'10.1088/0004-637X/775/1/62',11:'10.1093/mnras/stt2024',12:'10.3847/1538-4357/ab0ebf',13:'10.1093/mnras/stv547',14:'10.1111/j.1365-2966.2012.21402.x',15:'10.3847/2041-8213/abdd40',16:'10.1051/0004-6361/201936602',17:'10.1086/305424',18:'10.1038/385139a0',19:'10.1093/mnras/289.4.801',20:'10.1093/mnras/stt1169',21:'10.1111/j.1365-2966.2004.07661.x',22:'10.1086/148982',23:'10.1038/217511a0'}
def esc(s):
  s=s.replace('\\','\\textbackslash{}')
  for a,b in [('&','\\&'),('%','\\%'),('#','\\#'),('_','\\_'),('∼','$\\sim$'),('Liivam¨agi','Liivamägi'),('Kov´acs','Kovács'),('Garc´ıa','García'),('Gottl¨ober','Gottlöber'),('M¨uller','Müller'),('mass– richness','mass–richness'),('threedimensional','three-dimensional'),('LX–T','$L_X$–$T$'),('Lemaˆıtre','Lemaître')]: s=s.replace(a,b)
  return s
def lk(u,label=None): return '\\href{%s}{\\nolinkurl{%s}}'%(u,label or u)
D=lambda d:lk('https://doi.org/'+d,'doi:'+d)
A=lambda a:lk('https://arxiv.org/abs/'+a,'arXiv:'+a)
coe=[
('A. G. Riess et al., Astron. J. 116, 1009 (1998).',[D('10.1086/300499')]),
('S. Perlmutter et al., Astrophys. J. 517, 565 (1999).',[D('10.1086/307221')]),
('G. F. R. Ellis and W. R. Stoeger, Class. Quantum Grav. 4, 1697 (1987).',[D('10.1088/0264-9381/4/6/025')]),
('E. Hubble, Proc. Natl. Acad. Sci. USA 15, 168 (1929).',[D('10.1073/pnas.15.3.168')]),
('K. Migkas et al., Astron. Astrophys. 636, A15 (2020).',[D('10.1051/0004-6361/201936602'),A('2004.03305')]),
('L. I. Sedov, Prikl. Mat. Mekh. [J. Appl. Math. Mech.] 10, 241 (1946).',None),
('G. I. Taylor, Proc. R. Soc. Lond. A 201, 159 (1950).',[D('10.1098/rspa.1950.0049')]),
('G. Lemaître, Ann. Soc. Sci. Bruxelles A53, 51 (1933); English translation: Gen. Rel. Grav. 29, 641 (1997).',[D('10.1023/A:1018855621348')+' (translation)']),
('R. C. Tolman, Proc. Natl. Acad. Sci. USA 20, 169 (1934).',[D('10.1073/pnas.20.3.169')]),
('H. Bondi, Mon. Not. R. Astron. Soc. 107, 410 (1947).',[D('10.1093/mnras/107.5-6.410')]),
('Planck Collaboration, Astron. Astrophys. 641, A7 (2020).',[D('10.1051/0004-6361/201935201')]),
('N. J. Secrest, S. von Hausegger, M. Rameez, R. Mohayaee, S. Sarkar, and J. Colin, Astrophys. J. Lett. 908, L51 (2021).',[D('10.3847/2041-8213/abdd40'),A('2009.14826')])]
def items(lst):
  return '\n'.join('\\item[{[%d]}] %s %s'%(i,esc(t),', '.join(l) if l else '\\textit{[no online link found]}') for i,(t,l) in enumerate(lst,1))
cm=[]
for s in cmb:
  n=int(re.match(r'\[(\d+)\]',s).group(1)); body=re.sub(r'^\[\d+\]\s*','',s)
  ax=re.search(r',?\s*arXiv:(\S+?)\.?$',body); links=[D(cmbdoi[n])]
  if ax: body=body[:ax.start()]+'.'; links.append(A(ax.group(1)))
  cm.append((body,links))
pt=[]
for s in pat:
  n=int(re.match(r'\[(\d+)\]',s).group(1)); body=re.sub(r'^\[\d+\]\s*','',s)
  body=re.sub(r'\s*https://\S+?(,|$)','',body).strip().rstrip(',')
  body=re.sub(r';\s*arXiv:\S+?\.$','.',body)
  dois=re.findall(r'https://doi.org/(\S+?)(?:,|$)',s); axs=re.findall(r'arxiv.org/abs/(\S+?)(?:,|$)',s)
  if n==9: dois=['10.1080/00018739400101505']
  if n==29: dois=['10.1103/PhysRevE.76.066301']
  if n==55:
    body='C. Michaut, C. Cavet, S. E. Bouquet, F. Roy, H. C. Nguyen, Numerical study of the Vishniac instability in supernova remnants, Astrophys. J. 759, 78 (2012).'; dois=['10.1088/0004-637X/759/2/78']
  pt.append((body,[D(d) for d in dois]+[A(a) for a in axs]))
note=[('ellis79',"G.~F.~R. Ellis, ``The homogeneity of the universe,'' \\emph{Gen.\\ Rel.\\ Grav.} \\textbf{11}, 281 (1979).",[D('10.1007/BF00759562')]),
('ellis02',"G.~F.~R. Ellis, ``Cosmology and local physics,'' \\emph{New Astron.\\ Rev.} \\textbf{46}, 645 (2002).",[D('10.1016/S1387-6473(02)00234-8')]),
('ellisstoeger',"G.~F.~R. Ellis and W.~R. Stoeger, ``The `fitting problem' in cosmology,'' \\emph{Class.\\ Quantum Grav.} \\textbf{4}, 1697 (1987).",[D('10.1088/0264-9381/4/6/025')]),
('milne',"E.~A. Milne, \\emph{Relativity, Gravitation and World-Structure} (Oxford, 1935).",[lk('https://archive.org/details/dli.ernet.7948','archive.org/details/dli.ernet.7948')]),
('rindler',"W.~Rindler, \\emph{Relativity: Special, General, and Cosmological}, 2nd ed.\\ (Oxford, 2006).",[D('10.1093/oso/9780198567318.001.0001')]),
('migkas',"K.~Migkas et al., ``Probing cosmic isotropy with a new X-ray galaxy cluster sample through the $L_X$--$T$ scaling relation,'' \\emph{A\\&A} \\textbf{636}, A15 (2020).",[D('10.1051/0004-6361/201936602'),A('2004.03305')]),
('secrest',"N.~J. Secrest et al., ``A test of the cosmological principle with quasars,'' \\emph{ApJL} \\textbf{908}, L51 (2021).",[D('10.3847/2041-8213/abdd40'),A('2009.14826')])]
L=lambda lst:'\\begin{description}[leftmargin=2.6em,labelwidth=2.2em,labelsep=0.4em,itemsep=1pt,parsep=0pt,font=\\normalfont]\n'+items(lst)+'\n\\end{description}\n'
refs='\\clearpage\n\\phantomsection\\pdfbookmark[0]{References}{allrefs}\n\\begin{center}{\\Large\\bfseries References\\par}\\end{center}\n\\noindent\\small References are grouped by paper. Each group keeps that paper\'s own numbering, so the bracketed citations in the original pages refer to the matching group below.\n'
refs+='\\subsection*{1. CoE cosmology: Apparent acceleration from an off-centre finite geometry}\n'+L(coe)
refs+='\\noindent\\textit{The supplement on the exact relativistic (LTB) calculation cites no additional references.}\n'
refs+='\\subsection*{2. Reading the CMB fingerprint of the CoE shell model}\n'+L(cm)
refs+='\\subsection*{3. Pattern formation (v5)}\n'+L(pt)
refs+='\\noindent\\textit{Inserted Figs.~101--103} are reproduced from Michaut et al.\\ (2012), Astrophys. J. 759, 78, '+D('10.1088/0004-637X/759/2/78')+' (ref.~[55] above). Fig.~104 is original to this collection.\n'
refs+='\\subsection*{4. Does everything have a centre? A philosophical note}\n\\renewcommand{\\refname}{}\\vspace{-3.5em}\n\\begin{thebibliography}{9}\\setlength{\\itemsep}{1pt}\n'+'\n'.join('\\bibitem{%s} %s %s'%(k,t,', '.join(l)) for k,t,l in note)+'\n\\end{thebibliography}\n'
s=open('back.tex.bak_refs').read()
a=s.index('\\begin{thebibliography}'); b=s.index('\\end{thebibliography}')+len('\\end{thebibliography}')
s=s[:a]+refs+s[b:]
s=s.replace('\\usepackage[bookmarks','\\usepackage{lmodern}\\usepackage[T1]{fontenc}\\usepackage[utf8]{inputenc}\\usepackage{enumitem}\\usepackage{xurl}\\usepackage[colorlinks=true,urlcolor=blue!60!black,linkcolor=black,citecolor=black,bookmarks',1)
open('back.tex','w').write(s)
print('entries',len(coe)+len(cm)+len(pt)+len(note)+1)
