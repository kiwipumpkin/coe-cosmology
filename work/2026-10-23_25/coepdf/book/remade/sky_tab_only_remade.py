# AUTO-WRAP generator for sky_tab_only_remade.png
"""sky_tab_only_remade.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: byte-identical to /tmp/scr/l50/tab120.png; table section of comb45c.py
Embedded source: /tmp/scr/l50/comb45c.py  (verbatim below, unmodified)
Run:  python3 sky_tab_only_remade.py   -> writes sky_tab_only_remade_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /tmp/scr/l50 so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'sky_tab_only_remade' + '_regen.png')
SRC_CWD = '/tmp/scr/l50'
MATCH = ['tab.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_sky_tab_only_remade_')
import matplotlib; matplotlib.use('Agg')
import matplotlib.figure as _mf
_orig_open = builtins.open
def _divert(p):
    p = os.fspath(p)
    if os.path.abspath(p) == OUT: return p
    return os.path.join(SCRATCH, os.path.basename(p))
def _open(file, mode='r', *a, **k):
    if isinstance(file, (str, bytes, os.PathLike)) and any(c in mode for c in 'wax+'):
        file = _divert(file)
    return _orig_open(file, mode, *a, **k)
builtins.open = _open
_orig_save = _mf.Figure.savefig
def _savefig(self, fname, *a, **k):
    if isinstance(fname, (str, os.PathLike)):
        b = os.path.basename(os.fspath(fname))
        fname = OUT if b in MATCH else os.path.join(SCRATCH, b)
        if b == 'tab.png': k['dpi'] = 120  # the shipped PNG (= /tmp/scr/l50/tab120.png) is this table at dpi 120 instead of the source's 110
    return _orig_save(self, fname, *a, **k)
_mf.Figure.savefig = _savefig
try:
    from PIL import Image as _I
    _ois = _I.Image.save
    def _isave(self, fp, *a, **k):
        if isinstance(fp, (str, os.PathLike)):
            b = os.path.basename(os.fspath(fp)); fp = OUT if (b in MATCH or os.path.abspath(os.fspath(fp)) == OUT) else os.path.join(SCRATCH, b)
        return _ois(self, fp, *a, **k)
    _I.Image.save = _isave
except Exception: pass
os.chdir(SRC_CWD); sys.path.insert(0, SRC_CWD)
sys.argv = ['/tmp/scr/l50/comb45c.py'] + []

SOURCE = "import matplotlib;matplotlib.use('Agg')\nimport matplotlib.pyplot as plt\nfrom PIL import Image\ntop=Image.open('combined_v5.png').crop((0,70,2006,1175))\ntop.save('top.png')\nrows=[('A','CMB dipole cold','84, −48','cold','< 0','+0.02','?','Our motion → distance/velocity'),\n('A','Bulk flow away','107, −8','cold','< 0','−2.27','✓','Recession velocity → distance'),\n('A','CMB Cold Spot','209, −57','cold','< 0','−2.33','✓','Eridanus supervoid → distance'),\n('A','CMB dipole hot','264, +48','hot','> 0','+1.54','✓','Our motion → distance/velocity'),\n('A','Bulk flow toward','287, +8','hot','> 0','+1.90','✓','Recession velocity → distance'),\n('B','NVSS cold','16, −46','cold','< 0','−4.87','✓','Source density → sample'),\n('B','Hemi. asym. cold','50, +20','cold','< 0','−2.11','✓','Power direction → sample'),\n('B','CatWISE cold','58, −31','cold','< 0','−4.87','✓','Matter dipole → sample'),\n('B','NVSS hot','196, +46','hot','> 0','+2.03','✓','Source density → sample'),\n('B','Hemi. asym. hot','230, −20','hot','> 0','+1.62','✓','Power direction → sample'),\n('B','CatWISE hot','238, +31','hot','> 0','+1.39','✓','Matter dipole → sample'),\n('B','Loop I / NPS','329, +17.5','hot','> 0','+1.60','✓','Soft X-ray foreground → slope')]\nfig,ax=plt.subplots(figsize=(18,5.2));ax.axis('off')\nt=ax.table(cellText=[r[:7]+(r[7],) for r in rows],colLabels=['Tied to','Spot','ℓ, b','Hot/cold','Needs','Residual at spot (%)','Match','Why A or B'],loc='center',cellLoc='center',colWidths=[.06,.15,.1,.07,.06,.13,.06,.3])\nt.auto_set_font_size(False);t.set_fontsize(12);t.scale(1,1.6)\nfor i,r in enumerate(rows,1):\n  t[i,3].get_text().set_color('red' if r[3]=='hot' else 'blue')\n  t[i,6].get_text().set_color({'✓':'green','?':'darkorange'}.get(r[6],'red'))\nax.set_title('Axis lifted 15°, then 45° counterclockwise about CoE axis  —  Score: A 4/5, B 7/7, total 11/12, 1 on the zero line (?)',fontsize=15)\nfig.savefig('tab.png',dpi=110,facecolor='white',bbox_inches='tight')\nims=[Image.open(p).convert('RGB') for p in ['top_new.png','overlay_lift15_r45_bluered.png','tab.png']]\nW=max(i.width for i in ims)\nims=[i.resize((W,int(i.height*W/i.width))) for i in ims]\nout=Image.new('RGB',(W,sum(i.height for i in ims)),'white');y=0\nfor i in ims: out.paste(i,(0,y));y+=i.height\nout.save('combined_r45_bluered.png')\n"
try:
    exec(compile(SOURCE, '/tmp/scr/l50/comb45c.py', 'exec'), {'__name__': '__main__', '__file__': '/tmp/scr/l50/comb45c.py'})
except FileNotFoundError as e:  # later combine step needs top_new.png, no longer in /tmp/scr/l50; tab.png is already written
    print('note:', e)
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
