# Generates R28_mapB.png. Real code: R26_29_sky_chain.py (run from this folder).
import runpy,os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
runpy.run_path("R26_29_sky_chain.py",run_name="__main__")
