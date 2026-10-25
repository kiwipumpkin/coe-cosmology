import json, numpy as np, multiprocessing as mp
import rm_core as C
def one(r):
    try: return C.evaluate(r["Mc"], r["D"], r["R"])
    except Exception as e: return dict(fail=str(e))
if __name__ == "__main__":
    S = [r for r in json.load(open("flipsearch.json"))["sample"] if r["ok"] and r["Hr"] > r["Hp"]]
    with mp.Pool(7) as P: R = P.map(one, S, chunksize=4)
    V = [o for o in R if "R2_s" in o and o["twoE_min"] > -1 and o["along"]]
    json.dump(V, open("along_quick.json", "w"), default=float)
    print(len(S), len(V))
    for o in sorted(V, key=lambda o: -o["R2_s"])[:6]:
        print("Mc %.2e D %.3f R %.3f R2s %.3f R2ns %.3f l0s %.1f 2E %.3f age %.2f Hp %.1f Hr %.1f m2 %.3f corr %.2f" % (o["Mc_kg"], o["D"], o["R"], o["R2_s"], o["R2_ns"], o["fits"]["0.1"]["l0_s"], o["twoE"], o["age_Gyr"], o["H_perp"], o["H_r"], o["fits"]["0.1"]["m2"], o["fits"]["0.1"]["corr"]))
