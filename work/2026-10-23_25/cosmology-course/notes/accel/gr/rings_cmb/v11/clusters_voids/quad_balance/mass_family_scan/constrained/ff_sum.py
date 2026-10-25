import json,numpy as np
from ff_sim import pchi
I=json.load(open('interp_results.json'));S=json.load(open('fullfree_sims.json'));L0=json.load(open('l0scan_results.json'))
rows=[]
for k,v in I.items():
    if k=='_meta' or 'ratio' not in v or 'Dl' not in v:continue
    rows.append(dict(name=k,ratio=v['ratio'],Dl=v['Dl'][:31],src='interp',fn=None))
for d in S:
    if 'ratio' in d:rows.append(dict(name='sim%03d'%d['i'],ratio=d['ratio'],Dl=d['Dl'],C=d['C'],frac=d['frac'],DRV=d['DRV'],mig=d['mig'],mig360=d['mig360'],src='sim'))
json.dump(rows,open('/tmp/rows.json','w'))
print(len(rows))
