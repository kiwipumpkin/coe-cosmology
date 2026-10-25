import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
D=json.load(open('decomp.json'));P=json.load(open('paper_final.json'));M8=json.load(open('m2only_in30.json'))['8'];r=np.load('ring_patterns.npz');e=np.load('final_s3_edited.npz')
fig,ax=plt.subplots(1,2,figsize=(18,6.5));L=np.arange(360)
ax[0].plot(L,r['R']/1e-5,'k',label='smooth CoE ring R (237×)');ax[0].plot(L,(r['R']-r['w0'])/1e-5,'r',lw=.6,alpha=.6,label='+ seed-3 web, start (789×)');ax[0].plot(L,e['res']/1e-5,'b',label='+ edited web, paper (0.998×)')
ax[0].axhspan(-1,1,color='gold',alpha=.4,label='1× CMB');ax[0].set_xlabel('ring angle (deg, paper frame)');ax[0].set_ylabel('residual ΔT/T ×1e5');ax[0].legend()
m=np.arange(1,13);w=.2
for j,(y,lab) in enumerate(((D['smooth_R']['m'],'smooth R'),(D['web_start_residual']['m'],'web start'),(P['per_m'],'edited web (paper)'),(M8['per_m_after'],'R − 8 in-plane bodies fitted to m=2 only'))):ax[1].bar(m+(j-1.5)*w,np.maximum(y,1e-3),w,label=lab)
ax[1].axhline(1,c='gold',lw=3,label='1× CMB');ax[1].set_yscale('log');ax[1].set_xticks(m);ax[1].set_xlabel('ring harmonic m');ax[1].set_ylabel('per-m rms (× CMB 1e-5)');ax[1].legend(fontsize=9)
plt.savefig('ring_plane_cancel.png',dpi=110,bbox_inches='tight')
