import os,sys
from common import *
SH=f"{SEQ}/shallow"; sys.path.insert(0,SH); _c=os.getcwd(); os.chdir(SH); import shallow_core as sc; os.chdir(_c)
def pat(X): return sc.shallow_pattern(X,1.5)
