import numpy as np,json
from numpy.polynomial.legendre import legval
from scipy.interpolate import CubicSpline
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/shrink_rec/'
x=np.linspace(-1,1,8001);ps=np.degrees(np.arccos(x))   # x = cos(angle from anti-CoE)
def lnT(f):
    o=np.load((f if f.startswith('/') else B+f))['o'];o=o[np.argsort(o[:,0])];return -CubicSpline(o[:,0],np.log(o[:,1]))(ps),float(np.mean(o[:,1]))
def alm(y):  # signed a_l0-type amplitude: y=sum A_l P_l(x); report per-l rms = A_l/sqrt(2l+1) in 1e-5, sign kept
    out=[]
    for l in range(0,13):
        c=np.zeros(l+1);c[l]=1;A=(2*l+1)/2*np.trapezoid(y*legval(x,c),x);out.append(A/np.sqrt(2*l+1)/1e-5)
    return np.array(out)
R={}
for f in ['b_1100_k0.npz','none_1100_k0.npz','b_1100_k0.3.npz','none_1100_k0.3.npz','b_1100_k1.npz','none_1100_k1.npz','floor/sim_b_0.000909.npz','floor/sim_a_0.000909.npz','/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/b_1100_k0.45.npz','/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/b_1100_k0.5.npz']:
    y,z=lnT(f);R[f]=dict(zmean=z,a=alm(y-y.mean()).tolist())
    print(f"{f:28s} 1+z={z:7.1f} l1={R[f]['a'][1]:+9.2f} l2={R[f]['a'][2]:+9.2f} l3={R[f]['a'][3]:+8.2f} l4={R[f]['a'][4]:+7.2f} l5-12rms={np.sqrt(np.sum(np.square(R[f]['a'][5:]))):.3f}")
json.dump(R,open('signed_multipoles.json','w'),indent=1)
