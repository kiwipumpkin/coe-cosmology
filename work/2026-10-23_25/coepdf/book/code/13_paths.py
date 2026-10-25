"""Path constants for CoE figure regeneration."""
from pathlib import Path

ROOT = Path("/workspace/coe-figures")
SRC = ROOT / "src"
OUT = Path("/workspace/figs_letter_exact_python")
FIGS_FROM_MAKERS = OUT
REF = ROOT / "ref"

COSMO = Path("/workspace/cosmology")
MIGKAS_CSV = COSMO / "migkas/fig6_digitized.csv"  # asymmetric A_err_low/high (letter gold)
MIGKAS_FIG7_CSV = COSMO / "migkas/digitised/migkas2020_fig7_OurSample_A_vs_longitude.csv"
BUBBLE = COSMO / "papers/spreadsheet_model/bubble"
LETTER_FIGS = (
    COSMO
    / "papers/A_equals_fourier/steps_v3/cmb_reconcile/letter_brief/figs"
)
CMB = COSMO / "papers/A_equals_fourier/steps_v3/cmb_reconcile"
CYCLE4 = CMB / "cycle4"
CYCLE6 = CMB / "cycle6_two_emitters"
CYCLE7 = CMB / "cycle7_recombination"

UNIVERSE_NEWT = BUBBLE / "universe_newt_omb.json"
UNIVERSE_REL = BUBBLE / "universe_rel_mcfree.json"
THIN_FREE = BUBBLE / "thin_b_free_hug.json"
NEWTON_CURVE_NPY = BUBBLE / "gr_ltb/refit_mass/newton_curve.npy"
LTB_CURVE_NPY = BUBBLE / "gr_ltb/refit_mass/ltb_mcfree_z03_curve.npy"
DTT_JSON = CYCLE4 / "dtt/dtt.json"
CANCEL_JSON = CYCLE6 / "simple_cmb_two_quad_cancel.json"
VS_KILL_JSON = CYCLE6 / "vs_kill_term_0_360.json"
TWO_EMITTERS_JSON = CYCLE6 / "two_emitters.json"
COE_ANALEMMA_JSON = LETTER_FIGS / "coe_analemma_texas.json"
COE_CARDINAL_JSON = LETTER_FIGS / "coe_cardinal_direction_texas.json"

