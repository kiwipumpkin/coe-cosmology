# Phase4: prof_core Ht (=H t0) vs mc at fixed geometry, to map age anchor -> mc.
import importlib.util,numpy as np,json
sp=importlib.util.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=importlib.util.module_from_spec(sp);sp.loader.exec_module(cw)
o=[]
for mc in [0.001,0.003,0.01,0.03,0.1,0.3,1,3]:
    for D,eta in [(34.75,0.0301),(25,1.0)]:
        r=cw.sky_locked(D,eta,mc,cw.RV)
        o.append(dict(mc=mc,D=D,eta=eta,Ht=None if r is None else r[1]['Ht']));print(o[-1],flush=True)
json.dump(o,open('/workspace/coepdf/joint_solve/phase4/pc_probe.json','w'),indent=1)
