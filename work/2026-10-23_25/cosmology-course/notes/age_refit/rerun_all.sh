cd /workspace/cosmology-course/notes/age_refit
for l in exp soft eta pl; do VMAX=1 python3 fit_forward.py $l > fwd/${l}_vc.log 2>&1 & python3 fit_forward.py $l > fwd/$l.log 2>&1 & done; wait
echo FITS_DONE
for l in exp soft eta pl; do python3 sr_local.py $l > fwd/sr_$l.log 2>&1 & done
python3 fwd_rrec.py > fwd/rrec.log 2>&1 &
python3 bracket.py > fwd/bracket.log 2>&1 &
python3 bracket_sr.py > fwd/bracket_sr.log 2>&1 &
python3 step2_redo.py > fwd/step2_redo.log 2>&1 &
wait; echo ALL_DONE
