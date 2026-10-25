import numpy as np
L=[l for l in open('/workspace/cosmology/migkas/migkas2020_tablec1_vizier.tsv') if not l.startswith('#') and l.strip()]
h=L[0].rstrip('\n').split('\t');rows=[r.rstrip('\n').split('\t') for r in L[3:]]
ix={k:h.index(k) for k in ['z','GLON','GLAT','T','LX']}
d=np.array([[float(r[ix[k]]) for k in ['z','GLON','GLAT','T','LX']] for r in rows]);print(d.shape,np.percentile(d[:,0],[0,33.3,50,66.7,100]))
np.save('clusters.npy',d)
