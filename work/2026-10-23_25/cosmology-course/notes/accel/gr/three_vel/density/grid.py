import numpy as np,pickle,multiprocessing as mp
from dens import ratio
DG=np.geomspace(15,400,40);PG={'exp':[None],'soft':list(np.geomspace(0.3,10,25)),'eta':list(np.geomspace(0.03,3,25))}
def job(a):k,i,j=a;return a,ratio(k,DG[i],PG[k][j])
jobs=[(k,i,j) for k in PG for i in range(len(DG)) for j in range(len(PG[k]))]
with mp.Pool(8) as P: res=dict(P.map(job,jobs))
out={k:np.array([[res[(k,i,j)] for j in range(len(PG[k]))] for i in range(len(DG))]) for k in PG}
pickle.dump(dict(DG=DG,PG=PG,R=out),open('grid.pkl','wb'))
for k in out: r=out[k][...,0];print(k,np.nanmin(r),np.nanmax(r),np.isnan(r).sum())
