"""Term-by-term per-slice solve with Chris's Newtonian integrator (ptshell + prof_core.run, as coe_steps_v3 Step 2 / cycle3b)."""
import sys, os, json, numpy as np
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB); _cwd = os.getcwd(); os.chdir(BUB)
import ptshell, prof_core as P
os.chdir(_cwd)
RV = 14.30048287260589; TH1 = np.arange(0, 360, 1.0); TH = np.radians(TH1)
G = 6.674e-11; MPC = 3.0857e19; GPCm = 3.0857e25
def rhoc(H): h = H / MPC; return 3 * h * h / (8 * np.pi * G)
OMB = json.load(open(BUB + "/universe_newt_omb.json"))["Rs14.3|sine"]
KIN, KOUT = (OMB["D"] - OMB["R_in"]) / RV, (OMB["R_out"] - OMB["D"]) / RV        # headline shell rule (1.028, 1.475)
def run(D, mcc, Rs, kin=KIN, kout=KOUT, Os=0.049):
    a = dict(DR=D / Rs, Os=Os, mcc=mcc, rin=1 - kin * RV / D, rout=1 + kout * RV / D, Lc=0.005, mc=0.0, Ob=0.0)
    cwd = os.getcwd(); os.chdir(BUB)
    try: r = P.run("ptshell", a, 0.0, 0.0)
    finally: os.chdir(cwd)
    return r
def shifted(af, l0): return np.interp((TH1 - l0) % 360, np.append(TH1, 360), np.append(af, af[0]))
def A_from(af, th, amp, phi): return (1 - amp * np.sin(np.radians(th) + phi)) * af
def harm(y, m): return 2 * np.mean(y * np.cos(m * TH)), 2 * np.mean(y * np.sin(m * TH))
def target(a1, p1, a2, p2):   # A = 1 + a1 cos(l - p1) + a2 cos 2(l - p2)
    return {1: (a1 * np.cos(np.radians(p1)), a1 * np.sin(np.radians(p1))), 2: (a2 * np.cos(np.radians(2 * p2)), a2 * np.sin(np.radians(2 * p2)))}
def dq(A):
    A = A / A.mean(); c1, s1 = harm(A, 1); c2, s2 = harm(A, 2)
    return dict(a1=np.hypot(c1, s1), p1=np.degrees(np.arctan2(s1, c1)) % 360, a2=np.hypot(c2, s2), p2=(np.degrees(np.arctan2(s2, c2)) / 2) % 180,
                a3=np.hypot(*harm(A, 3)))
