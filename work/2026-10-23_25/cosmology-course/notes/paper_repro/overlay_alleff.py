import numpy as np,json,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
a1=np.load('alleff1.npy'); a0=np.load('alleff0.npy'); j1=json.load(open('alleff1.json')); j0=json.load(open('alleff0.json'))
assert a0.shape[0]==4 and a1.shape[0]==3   # no-mass: stretch,dop1,full ; mass: stretch,dop1 only (full failed)
assert 'error' in j1['full'] and j1['stretch']['gamma_defined_dirs']==0 and j1['stretch']['gravity_factor'] is None
b=a1[0]; l=(b+15.2)%360; o=np.argsort(l)
d=lambda T:(T-T.mean())/T.mean()/1e-5
m=d(a1[1]); n=d(a0[3]); assert abs(m.max()-19.675)<0.01 and abs(n).max()<1e-5
fig,ax=plt.subplots(figsize=(10,5.5))
ax.plot(l[o],m[o],'r',lw=2,label='With central mass: volume stretch only\n(speed and gravity effects CANNOT be computed: speeds are 40-43x light speed)')
ax.plot(l[o],n[o],'b',lw=2,label='No mass: stretch + speed effect + time-slowing + gravity (none) + travel delay')
ax.axhspan(-1,1,color='g',alpha=.15,label='Required flatness: within 1 part in 100,000')
ax.set_xlabel('Direction on the sky, galactic longitude (degrees); CoE is at 195'); ax.set_ylabel('Temperature change, parts per 100,000')
ax.set_xlim(0,360); ax.set_ylim(-25,25); ax.axvline(195.2,color='k',ls=':'); ax.legend(fontsize=8,loc='lower left'); ax.grid(alpha=.3)
ax.set_title('CMB temperature vs direction, Migkas-fit geometry (full curve, nothing removed)')
fig.tight_layout(); fig.savefig('overlay_cmb_alleffects.png',dpi=150)
