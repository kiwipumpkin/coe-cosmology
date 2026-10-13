# Counts-in-cells from the public 2M++ map (Carrick+2015). Top-hat spheres, mass delta = delta_g*/b*.
import numpy as np, json
d=np.load("../data/twompp_density.npy").astype(np.float64)
n=257; dx=400/256; h=0.7
ax=(np.arange(n)-128)*dx; X,Y,Z=np.meshgrid(ax,ax,ax,indexing="ij"); r=np.sqrt(X**2+Y**2+Z**2)
F=np.fft.rfftn(d); bstar=0.3**0.55/0.43; RMAX=150.0   # reliable radius, Mpc/h
out={}
for R_Mpc in [10,20,40,80,150]:
    Rh=R_Mpc*h; k=(r<=Rh).astype(float); k/=k.sum()
    avg=np.fft.irfftn(F*np.fft.rfftn(np.fft.ifftshift(k)),s=d.shape,axes=(0,1,2))/bstar
    ok=r+Rh<=RMAX; a=avg[ok]
    Vavail=4/3*np.pi*(RMAX-Rh)**3; Vs=4/3*np.pi*Rh**3
    nind=(4/3*np.pi*RMAX**3)/Vs if Rh<RMAX else 1   # independent spheres that fit in the map volume (rough)
    m=a.mean(); s=a.std()
    # lognormal parameters for 1+delta
    ln=np.log1p(np.clip(a,-0.999,None))
    out[R_Mpc]=dict(R_hMpc=Rh,ncentres=int(ok.sum()),n_independent_approx=float(nind),mean=float(m),std=float(s),
        p01=float(np.percentile(a,0.1)),p1=float(np.percentile(a,1)),p50=float(np.median(a)),p99=float(np.percentile(a,99)),p999=float(np.percentile(a,99.9)),
        min=float(a.min()),max=float(a.max()),skew=float(((a-m)**3).mean()/s**3),ln_mu=float(ln.mean()),ln_sig=float(ln.std()),
        frac_abs_gt_0055=float((np.abs(a-m)>0.055).mean()),frac_abs_gt_0068=float((np.abs(a-m)>0.068).mean()))
    np.save(f"cic_R{R_Mpc}.npy",a.astype(np.float32))
    print(R_Mpc,{k:(round(v,4) if isinstance(v,float) else v) for k,v in out[R_Mpc].items()})
json.dump(dict(bstar=bstar,RMAX_hMpc=RMAX,rows=out),open("cic_2mpp.json","w"),indent=1)
