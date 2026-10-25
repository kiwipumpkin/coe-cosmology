import json,os
import interp_an as I
from multiprocessing import Pool
if __name__=='__main__':
    R=json.load(open('interp_results.json'))
    todo=[(k,v) for k,v in I.runs.items() if R.get(k,{}).get('missing')]
    print([t[0] for t in todo],[os.path.exists(v) for k,v in todo],flush=True)
    with Pool(3) as p:
        for k,A in p.map(I.work,todo):R[k]=A
    json.dump(R,open('interp_results.json','w'),indent=1);print('saved')
