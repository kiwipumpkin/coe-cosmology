import numpy as np,glob,json
from numpy.polynomial.legendre import legval
from scipy.interpolate import CubicSpline
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/'
def ml(f):
    o=np.load(f)['o'];o=o[np.argsort(o[:,0])];x=np.linspace(-1,1,4001);ps=np.degrees(np.arccos(x))
    T=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(ps));T=T/T.mean()-1   # axisymmetric: a_l0 -> rms amplitude per l
    out=[]
    for l in range(1,13):
        c=np.zeros(l+1);c[l]=1;al=(2*l+1)/2*np.trapezoid(T*legval(x,c),x)   # T=sum al P_l
        out.append(abs(al)/np.sqrt(2*l+1)/1e-5)                              # sqrt((2l+1)C_l/4pi) in 1e-5
    return out,float(np.mean(o[:,1]))
R={}
for f in sorted(glob.glob(B+'shrink_rec/*.npz')+glob.glob(B+'shrink_rec/floor/*.npz')+glob.glob(B+'dipole_amp/[sr]_*.npz')):
    try:L,z=ml(f);R[f.replace(B,'')]=dict(zmean=z,l1=L[0],l2=L[1],l3_5=float(np.sqrt(np.sum(np.square(L[2:5])))),l6_12=float(np.sqrt(np.sum(np.square(L[5:])))),per_l=L)
    except Exception as e:print(f,e)
for k,v in R.items():print(f"{k:38s} 1+z={v['zmean']:8.1f} l1={v['l1']:9.2f} l2={v['l2']:9.2f} l3-5={v['l3_5']:8.3f} l6-12={v['l6_12']:8.4f}")
json.dump(R,open('lever_scan.json','w'),indent=1)
