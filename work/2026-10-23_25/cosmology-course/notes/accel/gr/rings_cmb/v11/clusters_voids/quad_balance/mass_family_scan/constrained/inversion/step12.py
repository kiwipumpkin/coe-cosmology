import json,numpy as np,healpy as hp
from core import *;from obs import *
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import importlib
def vec(L,b):L,b=np.radians(L),np.radians(b);return np.array([np.cos(b)*np.cos(L),np.cos(b)*np.sin(L),np.sin(b)])
def lb(v):return float(np.degrees(np.arctan2(v[1],v[0]))%360),float(np.degrees(np.arcsin(v[2])))
U={'red':('../fitDM_runs/red_D2.937_f20000.npz',shape_poly([-0.041,-0.035,0.1353,0],2.937)),'REF':('../fitDM_runs/ref_D2.937_f30.npz',shape_slc())}
G={k:obs_grid(f) for k,f in(('full','/workspace/inv/smica_alm_l30.fits'),('masked','/workspace/inv/smica_masked_alm_l30.fits'))}
iN=hp.vec2pix(NA,*vec(0,15));L2=slice(2,LM+1);w=1/ClP[L2]
R={'notes':'a_l0 in muK, pole=anti-CoE axis aG (model psi=0); sigma_l^2=C_l(Planck binned TT) cosmic variance; axis grid HEALPix nside 32'}
for k,(f,sh) in U.items():
    o=np.load(f)['o'];sm,gm,dist=terms(o,shape=sh);am,ag=al0(sm),al0(gm);np.savez(f'terms_{k}.npz',am=am,ag=ag,sm=sm,gm=gm,dist=dist)
    # verify vs full map in galactic frame (run.py style rotation)
    r=dict(a_motion=am.tolist(),a_grav=ag.tolist(),a_sum=(am+ag).tolist(),emission_dist_Gpc=[float(dist.min()),float(dist.max())])
    r['Dl_sum_from_al0']=Dl_al0(am+ag).tolist();r['Dl_map_poleframe']=Dl_map(sm+gm).tolist()
    aG=vec(0,15);Vg=V;psi=np.degrees(np.arccos(np.clip(Vg@aG,-1,1)))
    # galactic-frame map by re-evaluating the pole-frame map through interpolation in psi
    import scipy.interpolate as si;srt=np.argsort(PSI);
    ip=lambda m:np.interp(psi,PSI[srt],m[srt]);mg=ip(sm+gm)
    r['Dl_map_galactic_interp']=Dl_map(mg).tolist()
    r['max_rel_diff_al0_vs_poleframe_map_l2_30']=float(np.max(np.abs(Dl_al0(am+ag)[2:]/np.array(r['Dl_map_poleframe'])[2:]-1)))
    for gk,(Ga,Cl) in G.items():
        frac0=(Ga**2/((2*l+1)*Cl)[:,None]);ao=Ga[L2]
        def fit(A):return ((am[L2,None]+np.multiply.outer(ag[L2],np.atleast_1d(A))-ao)**2*w[:,None]).sum(0)
        c1=fit(1.0);num=((ao-am[L2,None])*ag[L2,None]*w[:,None]).sum(0);den=(ag[L2]**2*w).sum();Af=num/den;cf=fit(Af)
        c0=((ao)**2*w[:,None]).sum(0)
        i1=int(np.argmin(c1));i2=int(np.argmin(cf))
        r[gk]=dict(nominal_axis=dict(lb=[0,15],chi2_Ag1=float(c1[iN]),Ag_free=float(Af[iN]),chi2_Agfree=float(cf[iN]),a_obs=Ga[:,iN].tolist(),m0_power_frac=frac0[:,iN].tolist()),
          best_axis_Ag1=dict(lb=lb(AX[i1]),chi2=float(c1[i1]),a_obs=Ga[:,i1].tolist(),m0_power_frac=frac0[:,i1].tolist()),
          best_axis_Agfree=dict(lb=lb(AX[i2]),chi2=float(cf[i2]),Ag=float(Af[i2])),
          zero_model_chi2_min_over_axes=float(c0.min()),Cl_obs=Cl.tolist(),m0_frac_expected_isotropic=(1/(2*l+1)).tolist(),
          m0_frac_max_over_axes=frac0.max(1).tolist())
        np.save(f'chi2map_{k}_{gk}_Ag1.npy',c1)
    R[k]=r;json.dump(R,open('step12_results.json','w'),indent=1)
    print(k,'grav/motion rms',gm.std()/sm.std(),{g:(r[g]['nominal_axis']['chi2_Ag1'],r[g]['best_axis_Ag1']) for g in G} if False else '',flush=True)
    for g in G:print(g,r[g]['nominal_axis']['chi2_Ag1'],r[g]['nominal_axis']['Ag_free'],r[g]['nominal_axis']['chi2_Agfree'],r[g]['best_axis_Ag1']['lb'],r[g]['best_axis_Ag1']['chi2'],r[g]['best_axis_Agfree'])
    print('verify',r['max_rel_diff_al0_vs_poleframe_map_l2_30'],np.array(r['Dl_map_galactic_interp'])[2:6],np.array(r['Dl_sum_from_al0'])[2:6])
