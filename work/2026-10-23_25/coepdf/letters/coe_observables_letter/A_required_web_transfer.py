# (d) Transfer function of the 90 deg longitude strip average, H(m)=sin(m*45deg)/(m*45deg); zeros at m=4,8,12,...
import numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
m=np.arange(0,2001);H=np.sinc(m*45/180);fig,ax=plt.subplots(1,2,figsize=(13,4.5))
ax[0].plot(m[:41],H[:41],"k.-");ax[0].axhline(0,color="gray",lw=.6);[ax[0].axvline(k,color="r",lw=.5,ls=":") for k in range(4,41,4)]
ax[0].set(xlabel="m",ylabel="H(m)",title="Strip transfer function, m=0..40 (red: zeros m=4k)");ax[0].grid(alpha=.3)
ax[1].loglog(m[1:],np.abs(H[1:])+1e-20,"k-",lw=.6);ax[1].loglog(m[1:],1/(m[1:]*np.pi/4),"r--",lw=.8,label="envelope 4/(pi m)");ax[1].set(xlabel="m",ylabel="|H(m)|",ylim=(1e-4,1.2),title="|H(m)| to 2000");ax[1].legend();ax[1].grid(alpha=.3,which="both")
plt.tight_layout();plt.savefig("A_required_web_transfer.png",dpi=140)
