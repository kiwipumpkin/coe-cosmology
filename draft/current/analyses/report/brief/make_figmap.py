"""Writes figmap.tex (old -> new figure/table numbers) from the old-order aux numbers (/tmp/oldnums.json, taken from main.aux
before the reorder) and main_reordered.aux."""
import re,json
old=json.load(open('/tmp/oldnums.json')) if __import__('os').path.exists('/tmp/oldnums.json') else json.load(open('oldnums_before_reorder.json'))
json.dump(old,open('oldnums_before_reorder.json','w'))
new={m.group(1):m.group(2) for m in re.finditer(r'newlabel\{((?:fig|tab):[^}]*)\}\{\{([^}]*)\}',open('main_reordered.aux').read())}
D={'fig:geom':('Geometry of the model, side view','Introduction'),'fig:waterfall':('Starting logic: shell to 237\\fl\\ fingerprint','Introduction'),
'fig:source':('Where the fingerprint comes from (mode sizes)','Fingerprint'),'fig:fp':('The fingerprint mode by mode, polar panels','Fingerprint'),
'fig:overview':('Best leftover vs number of bodies, all scenarios','Fitting a few bodies'),'fig:families':('Same curves split by family','Fitting a few bodies'),
'fig:ex1':('Example universes 1--3','Fitting a few bodies'),'fig:ex4':('Core-open example universes 10--13','Fitting a few bodies'),
'fig:coreopen':('Best leftover vs $N$, CoE core fully open','Fitting a few bodies'),'fig:smooth':('Smooth-field solution and trade-off','Fitting a few bodies'),
'fig:counts':('Cluster/void counts and mixed ratio','Fitting a few bodies'),'fig:patch':('Over:under patch ratio domains','Fitting a few bodies'),
'fig:patchmap':('Map of the required patch ratio','Fitting a few bodies'),'fig:dist':('Extreme contrasts vs radius; 2M++ map','Observed structure'),
'fig:hist':('2M++ density-contrast histograms','Observed structure'),'fig:seq':('Sequential placement curves','Growing and editing webs'),
'fig:nodel':('Foam s3, move+resize under ring backstop','Growing and editing webs'),'fig:entrun':('Foam s2, entropy rule (2.40\\fl)','Growing and editing webs'),
'fig:ringdens':('Body density vs distance from the ring','Growing and editing webs'),'fig:fcc':('FCC cell; perfect-to-random ladder','Network rules'),
'fig:comprule':('Compensation rule: clusters per void','Network rules'),'fig:compmesh':('Compensated web s3 mesh (0.998\\fl)','Compensated web'),
'fig:comphist':('Compensated web mesh histograms','Compensated web'),'fig:noiserows':('Leftover alone + mode bars, four webs','What is left'),
'fig:noisehist':('Leftover histograms vs Gaussian noise','What is left'),'fig:trimesh':('Foam s2 Delaunay mesh before/after','Diagnostics'),
'fig:trihist':('Foam s2 mesh histograms','Diagnostics'),'fig:pairE':('Pair energies by link type','Diagnostics'),'fig:arrows':('Net-force arrows','Diagnostics'),
'fig:shellsA':('Equal-mass shells: constant, widening','Shell models'),'fig:shellsB':('Equal mass: narrowing; density compared','Shell models'),
'fig:shellsC':('Structure, equal mass: constant, widening','Shell models'),'fig:shellsD':('Structure, equal mass: narrowing; compared','Shell models'),
'fig:shellsE':('Mass rising outward: constant, widening','Shell models'),'fig:shellsF':('Mass out: narrowing; density compared','Shell models'),
'fig:shellsG':('Structure, mass out: constant, widening','Shell models'),'fig:shellsH':('Structure, mass out: narrowing; compared','Shell models'),
'fig:sequniv':('Sequential-placement universes 1--3','Appendix A'),'fig:sequniv2':('Sequential-placement universes 4--6','Appendix A'),
'fig:ex2':('Example universes 4--6','Appendix B'),'fig:ex3':('Example universes 7--9','Appendix B')}
rows=sorted(D,key=lambda l:int(new[l]))
assert set(D)=={k for k in old if k.startswith('fig:') and k!='fig:dvs'}
half=(len(rows)+1)//2
SEC={"Introduction":"I","Fingerprint":"II","Fitting a few bodies":"IV","Observed structure":"V","Growing and editing webs":"VI","Network rules":"VII","Compensated web":"VIII","What is left":"IX","Diagnostics":"X","Shell models":"XI","Appendix A":"App.\\,A","Appendix B":"App.\\,B"}
def line(l): o,n=old[l],new[l]; return f"{o} & \\textbf{{{n}}} & {D[l][0]} & {SEC[D[l][1]]}"
head="\\begin{center}\\scriptsize\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{rrll}\n\\toprule\nold & new & figure & sec.\\\\\n\\midrule\n"
foot="\\\\\n\\bottomrule\n\\end{tabular}\n\\end{center}\n"
with open('figmap.tex','w') as f:
    for k in range(0,len(rows),14):   # chunks of 14 rows so the table can break between columns
        f.write(head+"\\\\\n".join(line(l) for l in rows[k:k+14])+foot)
    f.write("Tables keep their numbers (I--VI). Figures~\\ref{fig:sequniv}--\\ref{fig:ex3} are in the appendices.\n")
print([(old[l],new[l]) for l in rows])
