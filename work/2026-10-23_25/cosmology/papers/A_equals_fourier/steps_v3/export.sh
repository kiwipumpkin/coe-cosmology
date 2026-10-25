set -e
cd /workspace/cosmology/papers/A_equals_fourier/steps_v3
J=/workspace/.nbvenv/bin/jupyter; PY=/workspace/.nbvenv/bin/python
[ "${EXEC:-1}" = 1 ] && $J nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=1800 coe_steps_v3.ipynb 2>&1 | grep -v "WARNING | Kernel"
$J nbconvert --to markdown coe_steps_v3.ipynb
$J nbconvert --to html coe_steps_v3.ipynb
# condensed print notebook: figure-script cells shortened to their first 30 lines (full verbatim code stays in .ipynb/.html/.md)
$PY - <<'P'
import nbformat
nb = nbformat.read("coe_steps_v3.ipynb", 4)
for c in nb.cells:
    if c.cell_type == "code" and c.source.startswith("%%runscript"):
        L = c.source.splitlines()
        if len(L) > 40:
            c.source = "\n".join(L[:30]) + "\n\n# ... %d more lines of the verbatim script: see coe_steps_v3.ipynb / .html / .md (or coe_steps_v3_fullcode.pdf)" % (len(L) - 30)
nbformat.write(nb, "_print.ipynb")
P
$J nbconvert --to html _print.ipynb
CSS="<style>@page{size:letter landscape;margin:0.35in} body{font-size:13px} .jp-Notebook{padding:0!important} pre,.highlight pre,.jp-OutputArea-output pre{white-space:pre-wrap!important;word-break:break-word;font-size:10.5px!important;line-height:1.25} .jp-RenderedMarkdown table{font-size:11px} .jp-InputPrompt,.jp-OutputPrompt{display:none!important} .jp-OutputArea-output img{max-width:100%!important;height:auto}</style>"
for pair in "_print.html:coe_steps_v3.pdf" "coe_steps_v3.html:coe_steps_v3_fullcode.pdf"; do
  in=${pair%%:*}; out=${pair##*:}
  $PY -c "import sys; s=open('$in').read(); open('_p.html','w').write(s.replace('</head>', sys.argv[1] + '</head>'))" "$CSS"
  timeout 900 google-chrome --headless=new --no-sandbox --disable-gpu --virtual-time-budget=30000 --no-pdf-header-footer --print-to-pdf=$out file://$PWD/_p.html 2>/dev/null
done
rm -f _p.html _print.html _print.ipynb
