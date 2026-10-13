import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from matplotlib.collections import LineCollection,PatchCollection
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import sys; sys.argv=["x"]
G=6.674e-11; GPC=3.0857e25; MC=1.40e50; JU=G*MC**2/GPC; FS=12; FT=13; FL=11
Z=np.load("comp_layout.npz"); C=json.load(open("comp_rows.json")); L=float(Z["L"]); K=Z["K"]; M=Z["M"]; RAD=Z["RAD"]; OWN=Z["OWN"]; N=len(K)
rows=C["rows"]; a=np.sqrt(2*L*L/(np.sqrt(3)*N))
def geom(P):
    sh=np.array([[i,j] for i in (-1,0,1) for j in (-1,0,1)])*L; Q=np.vstack([P+s for s in sh]); idx=np.tile(np.arange(N),9)
    S=Delaunay(Q).simplices; cen=Q[S].mean(1); keep=np.all((cen>=0)&(cen<L),1); S=S[keep]; q=Q[S]
    A=0.5*np.abs((q[:,1,0]-q[:,0,0])*(q[:,2,1]-q[:,0,1])-(q[:,2,0]-q[:,0,0])*(q[:,1,1]-q[:,0,1]))
    E=np.unique(np.sort(np.vstack([idx[S][:,[0,1]],idx[S][:,[1,2]],idx[S][:,[0,2]]]),1),axis=0)
    Ev=np.vstack([S[:,[0,1]],S[:,[1,2]],S[:,[0,2]]]); Ev=np.unique(np.sort(Ev,1),axis=0); Lv=np.linalg.norm(Q[Ev[:,0]]-Q[Ev[:,1]],axis=1)
    # unique physical links with their length
    key={}
    for (u,v),l in zip(np.sort(idx[Ev],1),Lv): key[(u,v)]=l
    El=np.array(list(key.keys())); Ll=np.array(list(key.values()))
    deg=np.bincount(El.ravel(),minlength=N)
    Ulink=-np.sum(M[El[:,0]]*M[El[:,1]]/(Ll*GPC))*G
    typ=(K[El[:,0]]+K[El[:,1]])
    return dict(A=A,L=Ll,deg=deg,El=El,U_link_per_node_J=float(Ulink/N),links_cc=int((typ==0).sum()),links_cv=int((typ==1).sum()),links_vv=int((typ==2).sum()))
def Uall(P):
    m=M/MC; u=0.0
    for i in range(N):
        d=P-P[i]; d-=L*np.round(d/L); r=np.hypot(d[:,0],d[:,1]); r[i]=np.inf; u+=m[i]*np.sum(m/r)
    return -0.5*u*JU
OUT={"a_Gpc":a,"n_nodes":N,"n_voids":int(K.sum()),"n_clusters":int((K==0).sum()),"clusters_per_void":float((K==0).sum()/K.sum())}
VIDX=np.where(K==1)[0]
for nm,P in (("placed (random angles)",Z["P0"]),("annealed (constrained ground state)",Z["P"])):
    g=geom(P); U=Uall(P)
    rr=[]; 
    for i in np.where(K==0)[0]:
        v=VIDX[OWN[i]]; d=P[i]-P[v]; d-=L*np.round(d/L); rr.append(np.hypot(*d)/RAD[v])
    rr=np.array(rr)
    OUT[nm]=dict(U_all_per_node_J=U/N,U_all_per_node_GMc2_over_a=U/N/(G*MC**2/(a*GPC)),U_link_per_node_J=g["U_link_per_node_J"],U_link_per_node_GMc2_over_a=g["U_link_per_node_J"]/(G*MC**2/(a*GPC)),
        tri_area_median_over_lattice=float(np.median(g["A"])/(np.sqrt(3)/4*a*a)),tri_area_mean_over_lattice=float(np.mean(g["A"])/(np.sqrt(3)/4*a*a)),tri_area_cv=float(g["A"].std()/g["A"].mean()),
        link_median_over_a=float(np.median(g["L"])/a),link_mean_over_a=float(np.mean(g["L"])/a),mean_degree=float(g["deg"].mean()),deg_void=float(g["deg"][K==1].mean()),deg_cluster=float(g["deg"][K==0].mean()),
        links_cc_cv_vv=[g["links_cc"],g["links_cv"],g["links_vv"]],wall_r_over_R_mean=float(rr.mean()),frac_wall_outer_quarter=float((rr>1.225).mean()))
    OUT[nm]["_g"]=g
# relation
Md=np.array([r["M_def"] for r in rows]); Nc=np.array([r["n_clusters"] for r in rows]); held=np.array([r["M_held"] for r in rows])/Md; Rv=np.array([r["R"] for r in rows]); mMc=np.array([r["mean_Mc"] for r in rows])
ok=Nc>0; sl,ic=np.polyfit(np.log10(Md[ok]/MC),np.log10(Nc[ok]),1); ls,li=np.polyfit(Md/MC,Nc,1); mcl=MC/ls
OUT["relation"]=dict(linear_fit="N = %.2f + M_def / %.2e kg"%(li,mcl),linear_slope_per_1p4e50=float(ls),linear_intercept=float(li),effective_Mc_kg=float(mcl),mean_Mc_untrimmed_expected_kg=2.21e50,fit="log10 N = %.3f log10(M_def/1.4e50 kg) + %.3f"%(sl,ic),slope=float(sl),intercept=float(ic),median_Mc_kg=float(np.median(M[K==0])),
   frac_voids_fully_compensated=float((held>0.999).mean()),min_held=float(held.min()),N_range=[int(Nc.min()),int(Nc.max())],R_range=[float(Rv.min()),float(Rv.max())],
   rule="N = ceil(M_def / M_c) with M_def=|delta_v| rho 4/3 pi R^3, M_c = delta_c rho 4/3 pi r_c^3 (excess mass), rho=%.3e kg/m^3"%C["RHO"])
fig,ax=plt.subplots(figsize=(8.5,7.2)); fig.subplots_adjust(bottom=0.22)
sc=ax.scatter(Md/MC,np.maximum(Nc,0.8),c=Rv,cmap="viridis",s=28,edgecolor="k",lw=0.3)
x=np.logspace(np.log10((Md/MC).min()),np.log10((Md/MC).max()),50); ax.plot(x,li+ls*x,"k--",lw=1.5,label=f"straight-line fit: N = {li:.1f} + M_def / ({mcl/1e50:.2f}e50 kg)")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("void mass deficit M_def (units of 1.40e50 kg, the web's median cluster)",fontsize=FS); ax.set_ylabel("clusters around the void",fontsize=FS)
cb=fig.colorbar(sc,ax=ax,pad=0.02); cb.set_label("void radius R (Gpc)",fontsize=FS); ax.grid(alpha=0.3,which="both"); ax.tick_params(labelsize=10.5)
ax.set_title("Clusters per void grow in proportion to the void's mass deficit",fontsize=FT)
fig.legend(loc="lower center",frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.02),title="each dot = one void; contrast −0.5 to −0.9, R 0.75–3 Gpc; clusters 100–250 Mpc, δ 30–300",title_fontsize=10)
fig.savefig("comp_clusters_per_void.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# layout before/after (40 Gpc window)
W0,W1=0,40
fig,ax=plt.subplots(2,1,figsize=(8.4,17.5)); fig.subplots_adjust(hspace=0.16,bottom=0.08,top=0.97)
for a_,nm in zip(ax,("placed (random angles)","annealed (constrained ground state)")):
    P=Z["P0"] if nm.startswith("placed") else Z["P"]; g=OUT[nm]["_g"]; El=g["El"]; d=P[El[:,1]]-P[El[:,0]]; s=(np.abs(d)<L/2).all(1)
    a_.add_collection(LineCollection(np.stack([P[El[s,0]],P[El[s,1]]],1),colors="0.7",lw=0.5,zorder=1))
    v=K==1; a_.add_collection(PatchCollection([Circle(p,r) for p,r in zip(P[v],RAD[v])],fc="#aed6f1",ec="#5dade2",lw=0.5,alpha=0.6,zorder=0.5))
    a_.add_collection(PatchCollection([Circle(p,1.3*r) for p,r in zip(P[v],RAD[v])],fc="none",ec="#5dade2",lw=0.5,ls=":",zorder=0.6))
    c=K==0; a_.scatter(P[c,0],P[c,1],s=np.clip(14*(M[c]/MC)**(2/3),3,60),c="#c0392b",lw=0,zorder=3)
    o=OUT[nm]; a_.set_xlim(W0,W1); a_.set_ylim(W0,W1); a_.set_aspect("equal"); a_.tick_params(labelsize=10.5); a_.set_xlabel("x (Gpc)",fontsize=FS); a_.set_ylabel("y (Gpc)",fontsize=FS)
    a_.set_title(f"{nm.upper()}\nenergy per node: {o['U_all_per_node_J']:+.2e} J (all pairs); clusters at {o['wall_r_over_R_mean']:.2f} R on average",fontsize=FT)
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=10,mfc="#aed6f1",mec="#5dade2",label="void, drawn to size (−mass)"),Line2D([],[],color="#5dade2",ls=":",label="outer edge of its wall, 1.3 R"),
    Line2D([],[],marker="o",ls="",ms=5,color="#c0392b",label="cluster (+mass), size ∝ mass"),Line2D([],[],color="0.7",lw=1,label="Delaunay link")],loc="lower center",ncol=2,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.0),
    title=f"40 × 40 Gpc corner of the {L:.0f} Gpc periodic box",title_fontsize=FL)
fig.savefig("comp_layout_before_after.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# geometry histograms vs lattice
fig,ax=plt.subplots(3,1,figsize=(9,12)); fig.subplots_adjust(hspace=0.42,bottom=0.12,top=0.93)
for nm,col,fill in (("placed (random angles)","0.5",True),("annealed (constrained ground state)","#1f3a93",False)):
    g=OUT[nm]["_g"]; kw=dict(color=col,alpha=0.35) if fill else dict(histtype="step",color=col,lw=1.8)
    A_=g["A"]/(np.sqrt(3)/4*a*a); ax[0].hist(A_,np.logspace(-3,1.5,60),weights=np.ones(len(A_))/len(A_),label=nm,**kw)
    ax[1].hist(g["deg"],np.arange(1.5,30.5,1),weights=np.ones(N)/N,**kw); ax[2].hist(g["L"]/a,np.linspace(0,5,60),weights=np.ones(len(g["L"]))/len(g["L"]),**kw)
for a_,x in zip(ax,(1,6,1)): a_.axvline(x,color="k",lw=1.8,ls="--")
ax[0].plot([],[],"k--",label="perfect triangular lattice (area √3/4 a², 6 links, length a)")
ax[0].set_xscale("log"); ax[0].set_xlabel("triangle area ÷ lattice triangle (√3/4 a²)",fontsize=FS); ax[1].set_xlabel("connections per node",fontsize=FS); ax[2].set_xlabel("link length ÷ a",fontsize=FS)
for a_,t in zip(ax,("(a) triangle areas","(b) connections per node","(c) link lengths")): a_.set_title(t,fontsize=FT,loc="left"); a_.set_ylabel("fraction",fontsize=FS); a_.grid(alpha=0.25); a_.tick_params(labelsize=10.5)
fig.suptitle(f"Compensated void + cluster layout vs perfect lattice (a = {a:.2f} Gpc, same node density)",fontsize=FT+1,fontweight="bold")
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.0),ncol=1)
fig.savefig("comp_geometry_hist.png",dpi=140,bbox_inches="tight"); plt.close(fig)
for nm in ("placed (random angles)","annealed (constrained ground state)"): OUT[nm].pop("_g")
json.dump(OUT,open("comp_results.json","w"),indent=1); print(json.dumps(OUT,indent=1))
