exec(open('recon_n3.py').read().split("cols=['purple'")[0])
for nm,Y,sY,mb,k in (('A',A0,sA,AFl,0),('B',B0,sB,BLl,1)):
  m=mb.copy()
  for i in range(4):
    if i>0:n,p=steps2[i-1];m=m+T(n,p,th,k)
    r=Y-m;print(nm,i,f'{np.sqrt(np.mean(r**2)):.3f}',f'{np.sum(r**2/sY**2):.2f}')
