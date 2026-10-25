import sys,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
script,out=sys.argv[1],sys.argv[2]
_orig=plt.subplots
class FakeFig:
    def __init__(s,figs): s.figs=figs
    def tight_layout(s,**k): [f.tight_layout() for f in s.figs]
    def savefig(s,name,**k):
        _real_save(name,**k)
        with PdfPages('/workspace/cosmology-course/notes/paper_repro/'+out) as pdf:
            for f in s.figs: pdf.savefig(f)
def sub(n,m=1,figsize=None,gridspec_kw=None,**k):
    hr=(gridspec_kw or {}).get('height_ratios',[1]*n); figs=[];axs=[]
    for h in hr:
        f,a=_orig(1,1,figsize=(11,11) if h==1 else (14,8.5)); figs.append(f); axs.append(a)
    return FakeFig(figs),__import__('numpy').array(axs)
src=open(script).read()
# keep PNG from the original full-figure code path: run original once for PNG separately
plt.subplots=sub
_real_save=lambda name,**k: None
exec(src)
