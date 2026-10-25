import glob,os,json,numpy as np
from multiprocessing import Pool
import fitz
def w(f):
    r=fitz.fit(f);r['z0']=fitz.fit(f,np.full_like(fitz.zc,0.003));return os.path.basename(f)[:-4],r
R=json.load(open('fitz_results.json'))
todo=[f for f in sorted(glob.glob('runs/*.npz')) if os.path.basename(f)[:-4] not in R]
with Pool(6) as p:
    for k,r in p.imap_unordered(w,todo):
        R[k]=r;json.dump(R,open('fitz_results.json','w'));print(k,r['exp']['chi2'],r['full']['chi2'],flush=True)
