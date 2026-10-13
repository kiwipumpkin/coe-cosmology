#!/usr/bin/env python3
"""DRAFT. Recreate the original working-directory layout of every figure script and input file
from this repository, so the figure scripts can be rerun.

    python3 rebuild_layout.py /tmp/coe_root
    cd /tmp/coe_root/workspace/<dir of script> && python3 <script>.py

PATHMAP.json maps each original absolute path (/workspace/...) to its location in the repo
(relative to draft/current/). Scripts were written with absolute /workspace/... paths; the copies
placed under <root> have that prefix rewritten to <root>/workspace/ so nothing outside <root> is touched."""
import json, os, shutil, sys
here = os.path.dirname(os.path.abspath(__file__)); cur = os.path.dirname(here)
root = os.path.abspath(sys.argv[1]); W = os.path.join(root, 'workspace') + '/'
for orig, rel in json.load(open(os.path.join(here, 'PATHMAP.json'))).items():
    src = os.path.normpath(os.path.join(cur, rel)); dst = W + orig[len('/workspace/'):]
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if orig.endswith('.py'):
        open(dst, 'w').write(open(src, errors='surrogateescape').read().replace('/workspace/', W))
    else:
        shutil.copy2(src, dst)
for d in json.load(open(os.path.join(here, 'OUTDIRS.json'))):
    os.makedirs(W + d[len('/workspace/'):], exist_ok=True)
print('layout rebuilt under', W)
