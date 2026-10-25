import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad
a,b=3.0,9.0; Ns=[6,12,24]; r=np.linspace(a,b,600)
qc=2-np.log(1.03)/np.log(b/a)                      # continuum: dm/dr ∝ r^qc -> rho ∝ r^(qc-2), rho(3)/rho(9)=1.03
def cont(kind):
    f={"eq":lambda x:1.0+0*x,"flat":lambda x:x**2,"3pct":lambda x:x**qc}[kind]
    Z=quad(f,a,b)[0]; return lambda x:f(x)/Z          # total mass = 1
def stairs(kind,N):
    e=np.linspace(a,b,N+1); ri,ro=e[:-1],e[1:]; m=(ri+ro)/2; w=ro-ri; V=(4/3)*np.pi*(ro**3-ri**3)
    if kind=="eq": M=np.ones(N); p=None
    elif kind=="flat": M=V.copy(); p=None
    else:
        p=brentq(lambda q:(m[0]**q/V[0])/(m[-1]**q/V[-1])-1.03,0.5,4); M=m**p
        assert abs((M[0]/V[0])/(M[-1]/V[-1])-1.03)<1e-12
    M=M/M.sum(); return e,m,M/w,M/V,p                 # total mass = 1
names={"eq":"equal mass per shell","flat":"flat density (mass ∝ volume)","3pct":"3% density drop D−R to D+R"}
fig,ax=plt.subplots(3,2,figsize=(16,14),sharex=True)
cols={6:"tab:red",12:"tab:green",24:"tab:blue"}; rep=[]
for i,k in enumerate(names):
    g=cont(k); rc=lambda x:g(x)/(4*np.pi*x**2)
    ax[i,0].plot(r,g(r),"k",lw=2.5,label="continuum"); ax[i,1].plot(r,rc(r),"k",lw=2.5,label="continuum")
    for N in Ns:
        e,m,mw,rho,p=stairs(k,N)
        ax[i,0].stairs(mw,e,color=cols[N],lw=1.6,baseline=None,label=f"{N} shells"+(f" (p={p:.4f})" if p else ""))
        ax[i,1].stairs(rho,e,color=cols[N],lw=1.6,baseline=None,label=f"{N} shells")
        # compare to continuum averaged over each shell
        cm=np.array([quad(g,e[j],e[j+1])[0]/(e[j+1]-e[j]) for j in range(N)])
        cr=np.array([quad(lambda x:g(x),e[j],e[j+1])[0] for j in range(N)])/((4/3)*np.pi*(e[1:]**3-e[:-1]**3))
        dmw=np.max(np.abs(mw/cm-1)); dr=np.max(np.abs(rho/cr-1)); rep.append((k,N,p,dmw,dr))
    ax[i,0].set_ylabel(f"{names[k]}\nmass per unit width"); ax[i,1].set_ylabel("density")
    for a_ in ax[i]: a_.legend(fontsize=9,loc="center left",bbox_to_anchor=(1.0,0.5)); a_.grid(alpha=.3)
for a_ in ax.flat:
    a_.set_xlim(a,b); a_.axvline(a,color="k",lw=2); a_.axvline(b,color="k",lw=2)
    a_.set_xticks(range(3,10)); a_.set_xticklabels(["D−R\n3","4","5","6","7","8","D+R\n9"])
for a_ in ax[2]: a_.set_xlabel("distance from CoE (c·Δt)")
ax[0,0].set_title("mass per unit width (mass per shell ÷ shell width)"); ax[0,1].set_title("density (mass ÷ shell volume)")
fig.suptitle(f"t = 3, exactly-c geometry. Total mass between D−R and D+R = 1 in every case. Continuum 3% law: dm/dr ∝ r^{qc:.4f}",fontsize=13)
plt.tight_layout(); plt.savefig("shell_refinement.png",dpi=170)
print("qc",qc)
for x in rep: print("%-5s N=%2d p=%s max dev mass/width=%.2e density=%.2e"%(x[0],x[1],"%.4f"%x[2] if x[2] else "-",x[3],x[4]))
d3=[x[3] for x in rep if x[0]=="3pct"]; assert d3[0]>d3[1]>d3[2] and d3[2]<2e-3
