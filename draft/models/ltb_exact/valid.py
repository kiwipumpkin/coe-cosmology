from ltbx import *; import sys
A=np.arange(0,181,15.)
print("(a) no mass, D=20, R=1,5,14.3")
M=Model(20,0.0); H,te,re,_=M.trace(A,[1.,5.,14.3])
for R in H: print("  R=%.1f H:"%R,np.round(H[R],2))
print("(b) small D=1 Gpc, small mass; cosmological vs classical vs local-limit (H_r cos^2+H_perp sin^2)")
for rs in (1e-3,1e-2):
    M=Model(1.0,rs,rc=0.1,rmax=6)
    RL=[0.025,0.05,0.1,0.2,0.5]; H,_,_,_=M.trace(A,RL,sub=4); loc=M.local(A)
    print(" rs=%g Gpc: v_obs=%.3fc  local-limit H 0/90/180: %.3f %.3f %.3f"%(rs,np.interp(1,M.R[-1],M.Rd[-1])/C,loc[0],loc[6],loc[12]))
    for R in RL:
        cl=M.classical(A,R); dc=H[R]-cl; an_c=H[R].max()-H[R].min(); an_l=cl.max()-cl.min()
        print("   R=%.3f cosmo 0/90/180 %.3f %.3f %.3f | classical %.3f %.3f %.3f | max|diff| %.4f  anis(cosmo,class) %.4f %.4f"%(R,H[R][0],H[R][6],H[R][12],cl[0],cl[6],cl[12],abs(dc).max(),an_c,an_l),flush=True)
