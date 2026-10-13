"""Shared panel boxes (figure fractions) so both composites align when stacked."""
FIGSIZE = (10, 7.5)
L0, R0 = 0.068, 0.458   # left column x-extent
L1, R1 = 0.533, 0.923   # right column x-extent
B0, T0 = 0.062, 0.440   # bottom row
B1, T1 = 0.545, 0.923   # top row
SUP_Y = 0.992
CB_PAD, CB_W = 0.006, 0.010
def boxes():
    return [[L0, B1, R0-L0, T1-B1], [L1, B1, R1-L1, T1-B1],
            [L0, B0, R0-L0, T0-B0], [L1, B0, R1-L1, T0-B0]]
