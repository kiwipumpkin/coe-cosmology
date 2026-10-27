"""Render script for book_pages/p-XX.png (89 pages, 850x1100 = 100 dpi).
Source: the 89-page build of the book, produced by /workspace/coepdf/book/build_prev.py (writes book.tex),
compiled with tectonic, then pdftoppm -r 100. (The current book/book.pdf is a later 133-page build.)
Everything runs in a temporary copy of /workspace/coepdf/book, so book/book.tex and book.pdf are not touched.
Output: book_pages_regen/ beside this folder (verified 89/89 pages pixel-identical on Oct 26 2026).
Run: python3 render_book_pages.py   (takes a few minutes)"""
import os, subprocess, shutil, tempfile, glob
H = os.path.dirname(os.path.abspath(__file__)); OUT = H + '_regen'
T = tempfile.mkdtemp(prefix='bookprev_'); B = os.path.join(T, 'book')
shutil.copytree('/workspace/coepdf/book', B)
subprocess.run(['python3', 'build_prev.py'], cwd=B, check=True, stdout=subprocess.DEVNULL)
subprocess.run(['tectonic', 'book.tex'], cwd=B, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
os.makedirs(OUT, exist_ok=True)
subprocess.run(['pdftoppm', '-r', '100', '-png', os.path.join(B, 'book.pdf'), os.path.join(OUT, 'p')], check=True)
print(len(glob.glob(OUT + '/p-*.png')), 'pages ->', OUT)
