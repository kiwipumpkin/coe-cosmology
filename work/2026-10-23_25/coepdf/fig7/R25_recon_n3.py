# Generates R25_recon_n3.png. Real code: rerun/L25_run.py (run from this folder).
import runpy,os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
runpy.run_path("rerun/L25_run.py",run_name="__main__")
