# Phase4b (NEW) pipeline: bg-only compensation -> field from scratch -> compensation with field (quadrupole-weighted) -> TT chi2.
set -x
BG=/workspace/coepdf/joint_solve/phase4b/p2_sm_D27.3385_H68.64_f0.0185_ns64.npy
[ -f /workspace/coepdf/joint_solve/phase4b/p5_bg.json ] || python3 /workspace/coepdf/joint_solve/phase4b/p5_compensate4b.py --bg $BG --R 21.74 --MMAX 1e22 --LC 24 --WH 1 --NB 3000 --tag bg
[ -f /workspace/coepdf/joint_solve/phase4b/p8_d4b.json ] || TAG=d4b DG=27.3385 R=21.74 DFIXF=/workspace/coepdf/joint_solve/phase4b/p5_bg.npz NP=6 python3 /workspace/coepdf/joint_solve/phase4b/p8_scratch4b.py
for WQ in 1 5; do [ -f /workspace/coepdf/joint_solve/phase4b/p5_fld_WQ$WQ.json ] || python3 /workspace/coepdf/joint_solve/phase4b/p5_compensate4b.py --bg $BG --R 21.74 --MMAX 1e22 --LC 24 --WH 1 --NB 3000 --field /workspace/coepdf/joint_solve/phase4b/p8_d4b.npz --f 0 --WQ $WQ --tag fld_WQ$WQ; done
python3 /workspace/coepdf/joint_solve/phase4b/p6_tt_compensated4b.py /workspace/coepdf/joint_solve/phase4b/p8_d4b.npz /workspace/coepdf/joint_solve/phase4b/p5_fld_WQ5.npz d4b
echo CHAIN_DONE
