cd /workspace/cosmology-course/notes/age_refit
while ! grep -q ALL_DONE fwd/rerun_all.log; do sleep 20; done
mkdir -p fwd/after_newton; cp fwd/*.json fwd/*.log fwd/after_newton/
rm -f fwd/bracket_B.json fwd/bracket_sr_B.json
export GRZ=gr FROMCKPT=1
for l in exp soft eta pl; do VMAX=1 python3 fit_forward.py $l > fwd/${l}_vc.log 2>&1 & python3 fit_forward.py $l > fwd/$l.log 2>&1 & done; wait
for l in exp soft eta pl; do python3 sr_local.py $l > fwd/sr_$l.log 2>&1 & done
python3 bracket.py > fwd/bracket.log 2>&1 &
python3 bracket_sr.py > fwd/bracket_sr.log 2>&1 &
wait; mkdir -p fwd/after_gr; cp fwd/*.json fwd/*.log fwd/after_gr/; echo GR_DONE
