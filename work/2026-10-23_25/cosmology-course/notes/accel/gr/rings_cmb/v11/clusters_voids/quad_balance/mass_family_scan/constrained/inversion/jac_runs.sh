cd /workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained/
xargs -P 7 -L1 sh -c '[ -f $PWD/inversion/runs/$0.npz ] || DRV=$6 C1=$1 C2=$2 C3=$3 C4=$4 timeout 7200 python3 famP.py poly 0 $5 $PWD/inversion/runs/$0.npz > inversion/runs/$0.log 2>&1' < inversion/${1:-jobs.txt}
