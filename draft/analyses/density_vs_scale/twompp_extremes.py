import numpy as np, json
d=np.load("data/twompp_density.npy").astype(np.float64)
print(d.shape, d.min(), d.max(), d.mean(), d[128,128,128])
n=257; dx=400/256
ax=(np.arange(n)-128)*dx
X,Y,Z=np.meshgrid(ax,ax,ax,indexing="ij"); r=np.sqrt(X**2+Y**2+Z**2)
F=np.fft.rfftn(d)
bstar=0.3**0.55/0.43   # b* = f/beta*, f=Om^0.55, Om=0.3
out=[]
for Rh in [5,10,20,30,50,70,100]:
    k=(r<=Rh).astype(float); k/=k.sum()
    ks=np.fft.ifftshift(k)   # kernel centred at index 0
    avg=np.fft.irfftn(F*np.fft.rfftn(ks),s=d.shape)
    ok=r+Rh<=150
    a=avg[ok]; imax=np.argmax(a); imin=np.argmin(a)
    pmax=(X[ok][imax],Y[ok][imax],Z[ok][imax]); pmin=(X[ok][imin],Y[ok][imin],Z[ok][imin])
    row=dict(R_hMpc=Rh,R_Mpc=Rh/0.7,max_dg=a.max(),min_dg=a.min(),rms_dg=a.std(),
             max_dm=a.max()/bstar,min_dm=a.min()/bstar,rms_dm=a.std()/bstar,
             pmax_gal_hMpc=[round(float(v),1) for v in pmax],pmin_gal_hMpc=[round(float(v),1) for v in pmin],ncentres=int(ok.sum()))
    out.append(row); print({k:(round(v,3) if isinstance(v,float) else v) for k,v in row.items()})
json.dump(dict(bstar=bstar,rows=out,notes="2M++ (Carrick+2015) luminosity-weighted delta_g*, 4 Mpc/h Gaussian pre-smoothed; top-hat sphere averages; centres with r+R<=150 Mpc/h; mass delta = delta_g*/b*, b*=Om^0.55/beta* with beta*=0.43, Om=0.3"),open("twompp_extremes.json","w"),indent=1,default=float)
