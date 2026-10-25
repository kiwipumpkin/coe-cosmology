# consistency: tabulated exp profile vs analytic bubble_core
from prof_core import run as prun
from bubble_core import run as brun
import numpy as np
a = dict(DR=2.83, LD=0.0877); r1 = prun("exp", a, 0.0389, 0.3); r2_ = brun(2.83, 0.0877, 0.0389, 0.3)
print(np.max(np.abs(r1["af"] - r2_["af"])), r1["contrast"], r2_["contrast"], r1["vesc2"], r2_["vesc2"], r1["c34"], r2_["c34"])
