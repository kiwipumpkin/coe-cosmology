import re
s = open("uniform_bubble_template.tex").read()
for n in ("ub_t1", "ub_t2", "ub_t3", "alt_rows", "thin_rows"):
    s = s.replace("\\input{%s.tex}" % n, open("../%s.tex" % n).read().rstrip("\n"))
open("../uniform_bubble.tex", "w").write(s)
