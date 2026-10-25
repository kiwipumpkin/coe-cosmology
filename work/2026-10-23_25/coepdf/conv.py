import re
t=open('main_rev4.tex').read()
t=t.replace(r'\documentclass[10pt]{article}',r'\documentclass[10pt,twocolumn]{article}')
t=t.replace(r'\usepackage{multicol}\setlength{\columnsep}{0.6cm}',r'\setlength{\columnsep}{0.6cm}'+'\n'+r'\setcounter{topnumber}{2}\setcounter{totalnumber}{3}\renewcommand{\topfraction}{0.85}\renewcommand{\textfraction}{0.12}\renewcommand{\floatpagefraction}{0.8}\renewcommand{\dbltopfraction}{0.75}\renewcommand{\dblfloatpagefraction}{0.8}\setcounter{dbltopnumber}{1}')
t=t.replace(t.splitlines()[6], r'\newcommand{\FIG}[4]{\begin{figure*}[tp]\centering\fig{\textwidth}{#2\textheight}{fig/#1}\caption{#3}\label{#4}\end{figure*}}'+'\n'+r'\newcommand{\FIGc}[4]{\begin{figure}[tbp]\centering\fig{\columnwidth}{#2\textheight}{fig/#1}\caption{#3}\label{#4}\end{figure}}')
t=t.replace(r'\begin{multicols}{2}'+'\n','').replace(r'\end{multicols}'+'\n','').replace(r'\begin{multicols}{2}','').replace(r'\end{multicols}','')
t=t.replace(r'(revision 4,',r'(revision 5,')
# column figures with height caps
col={'f02.png':0.36,'d4_radmat.png':0.2,'f07.jpg':0.3,'f11.png':0.4,'d5_mod.png':0.2,'f06.png':0.55,'f05.png':0.4,'f17.png':0.3,'f14.png':0.2,'f15.png':0.2,'f18.png':0.36,'f01.png':0.32}
def rep(m):
    f=m.group(1)
    if f in col: return r'\FIGc{%s}{%s}'%(f,col[f])
    return m.group(0)
t=re.sub(r'\\FIG\{([^}]+)\}\{[0-9.]+\}',rep,t)
# inverse problem block
t=re.sub(r'\\begin\{center\}\n(\\fig\{0\.6.*?)\\captionof\{figure\}(.*?)\n\\label\{fig:inv\}\\end\{center\}',
 lambda m:'\\begin{figure*}[tp]\\centering\n'+m.group(1)+'\\caption'+m.group(2)+'\n\\label{fig:inv}\\end{figure*}',t,flags=re.S)
t=t.replace('\\begin{center}\\small\n\\captionof{table}','\\begin{table*}[tp]\\centering\\small\n\\caption')
t=t.replace('Right: peak $\\ell$ ($D_\\ell$ in $\\mu$K$^2$).\n\\end{center}','Right: peak $\\ell$ ($D_\\ell$ in $\\mu$K$^2$).\n\\end{table*}')
# sky section
old=re.search(r'\\FIG\{sky_spots\.png\}.*?\{fig:skyspots\}',t).group(0)
new=(r'\FIG{sky_angle.png}{0.5}{A (top) and B (bottom) running-mean residuals vs model-frame angle ($\pm15/25/35^\circ$), with each spot drawn at its model-frame angle.}{fig:skyspots}'+'\n'
 r'\begin{figure*}[tp]\centering\fig{0.495\textwidth}{0.4\textheight}{fig/sky_mapA.png}\hfill\fig{0.495\textwidth}{0.4\textheight}{fig/sky_mapB.png}'
 r'\caption{Residual maps (Galactic Mollweide, $\ell=0$ at centre, increasing to the left) lifted $15^\circ$ and rotated $45^\circ$ about the CoE axis, with spots at published $(\ell,b)$. Left: A. Right: B.}\label{fig:skymaps}\end{figure*}')
t=t.replace(old,new)
t=t.replace(r'\FIG{sky_spots_tab.png}{0.6}',r'\FIG{sky_spots_tab.png}{0.55}')
open('main_rev5.tex','w').write(t)
