# fine psi scan 100..120 deg (1 deg) for K-1 with SURF=t and SURF=a, to characterise the hit-surface jump
cd /workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained
P=$(python3 -c "print(','.join(str(x) for x in range(100,121)))")
for S in t a; do SURF=$S PSI=$P U=$(python3 -c "print(68.64*27.737393*1000/299792.458/1.00016)") DRV=$(python3 -c "print(27.737393/14.30048287260589)") LDX=0.42296181187611975 ETA=3.068065 OMP_NUM_THREADS=1 python3 /workspace/coepdf/joint_solve/phase5/famP_mig.py mig 0 0.01847547274832157 /workspace/coepdf/joint_solve/phase5/surf/fine_K-1_S$S.npz > /workspace/coepdf/joint_solve/phase5/surf/fine_S$S.log 2>&1 & done; wait
