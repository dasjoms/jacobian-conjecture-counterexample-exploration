"""
Note 20, stage D: the tour-wide fresh census + the corner-chases-shadow post-checks.
Locks:
  LC2 (census N=240k fresh, sigma=1.5, C=1):
      d=11 cone-0-real in [8.67, 8.84]% (published band; law center 8.7556)
      d=11 4-real in [7.2, 7.7]%  and > fresh d=9 4-real (row monotone)
      d=11 6-real < 0.01%
      fresh-vs-archived main buckets within 0.2pp for d = 5..10
  LC5 (coalescence, note-15 pre-reg windows):
      gap s_c - s*(11) in [1.3, 2.1]e-4   (series est 1.67e-4)
      |t1_c - t*| * 11^3 in [5.0, 7.5]e-3 (series {6.63,6.65,6.35}e-3)
      r_c - s_c in [-4.5, -2.0]e-4        (series {-1.49,-0.86,-0.53}e-3)
"""
import numpy as np, json, time
from collections import Counter
t0 = time.time()

def phi_coeffs(d):
    m = d*(d+1)
    co = np.zeros(d+2, float)
    co[0] = -1.0/(d+1); co[1] = 1.0/d
    co[d-2] = -(1 - 2/m); co[d-1] = (1 - 3/m)
    return co

def p_coeffs(d):
    m = d*(d+1)
    pc = np.zeros(d+1, float)
    pc[0] = -1.0; pc[1] = 1.0; pc[-2] = -(3 - 6/m); pc[-1] = (2 - 6/m)
    return pc

def census(d, N=240000, seed=11, chunk=20000, tol=1e-8):
    rng = np.random.default_rng(seed)
    cnt = Counter(); co = phi_coeffs(d); pc = p_coeffs(d)
    done = 0
    while done < N:
        n = min(chunk, N - done); done += n
        A = rng.normal(0, 1.5, n); B = rng.normal(0, 1.5, n)   # (s,r) = (B,A) at C=1
        coef = np.tile(co, (n, 1)); coef[:, -2] = -B; coef[:, -1] = A
        M = np.zeros((n, d+1, d+1))
        M[:, 0, :] = -coef[:, 1:] / coef[:, [0]]
        idx = np.arange(d)
        M[:, idx+1, idx] = 1.0
        Rt = np.linalg.eigvals(M)
        realmask = np.abs(Rt.imag) < tol
        pw = np.polyval(pc, Rt.real)
        keep = realmask & (np.abs(B[:, None] - pw) > 1e-9)
        cnt.update(keep.sum(1).tolist())
    return cnt

print("d | census (fresh 240k)                                    | top-bucket %")
fresh = {}
for d in range(4, 12):
    c = census(d)
    tot = sum(c.values())
    fresh[d] = {int(k): 100*v/tot for k, v in sorted(c.items())}
    print(f"{d:2d} | {dict(sorted(c.items()))}   [{time.time()-t0:.0f}s]")

# archived cross-check
arch = {5: ("atlas5_realcensus.json", (0, 2, 4)), 6: ("atlas6_realcensus.json", (1, 3)),
        7: ("atlas7_realcensus.json", (0, 2, 4)), 8: ("atlas8_realcensus.json", (1, 3)),
        9: ("atlas9_realcensus.json", (0, 2, 4)), 10: ("atlas10_realcensus.json", (1, 3))}
maxdiff = 0.0
for d, (fn, keys) in arch.items():
    a = json.load(open(fn)); tot = sum(a.values())
    for k in keys:
        dv = abs(100*a[str(k)]/tot - fresh[d].get(k, 0.0))
        maxdiff = max(maxdiff, dv)
print(f"\nfresh-vs-archived max bucket diff: {maxdiff:.4f} pp  [lock <= 0.2]")

c11 = fresh[11]
cone = c11.get(0, 0.0); four = c11.get(4, 0.0); six = c11.get(6, 0.0)
print(f"d=11: cone {cone:.4f}% [lock 8.67..8.84, center 8.7556] | 4-real {four:.4f}% [lock 7.2..7.7]"
      f" | 6-real {six:.5f}% [lock <0.01]")
print(f"      monotone 4-real row d=9 -> 11: {four:.4f} > {fresh[9].get(4,0.0):.4f}:", four > fresh[9].get(4,0.0))
ok2 = (8.67 <= cone <= 8.84) and (7.2 <= four <= 7.7) and six < 0.01 and (four > fresh[9].get(4,0.0)) and maxdiff <= 0.2
print("LC2 GREEN:", ok2)

# ---- LC5 coalescence ----
s_star = -0.9417878439676974; t_star = -1.0834586214993059682   # note-15 certified windows
bt = json.load(open("atlas11_bitangents.json"))
ND = [(complex(a.replace(' ','')), complex(b.replace(' ','')), complex(sv.replace(' ','')), complex(rv.replace(' ','')), k)
      for a, b, sv, rv, k in bt["nodes"]]
cr = next(q for q in ND if q[4] == "CRUNODE")
t1c, t2c, sc, rc = cr[0].real, cr[1].real, cr[2].real, cr[3].real
gap = sc - s_star; tgap = abs(t1c - t_star)*11**3; roff = rc - sc
print(f"\ncrunode exact (120-dps pairing): t1={t1c:.15f} t2={t2c:.15f}")
print(f"  s_c = {sc:.15f}  r_c = {rc:.15f}")
print(f"  gap s_c - s* = {gap:.6e}  [lock 1.3e-4..2.1e-4;   series est 1.67e-4]")
print(f"  |t1-t*|*11^3 = {tgap:.6e}  [lock 5.0e-3..7.5e-3;   series ~6.6e-3]")
print(f"  r_c - s_c    = {roff:.6e}  [lock -4.5e-4..-2.0e-4;  series -> ~-3.4e-4]")
deltas = [(5, 6.55e-4), (7, 3.98e-4), (9, 2.50e-4), (11, gap)]
alphas = [float(np.log(deltas[i][1]/deltas[i+1][1]) / np.log(deltas[i+1][0]/deltas[i][0])) for i in range(3)]
cvals = [d_*d_*dl for d_, dl in deltas]
print(f"  gap series: {[f'{v:.4g}' for _, v in deltas]}   interval exponents: {[f'{a:.3f}' for a in alphas]}"
      f"   c=d^2*gap: {[f'{v:.5f}' for v in cvals]}")
ok5 = (1.3e-4 <= gap <= 2.1e-4) and (5.0e-3 <= tgap <= 7.5e-3) and (-4.5e-4 <= roff <= -2.0e-4)
print("LC5 GREEN:", ok5)

json.dump({"fresh_census_pct": {str(d): fresh[d] for d in fresh}, "archived_maxdiff_pp": maxdiff,
           "cone11": cone, "four11": four, "six11": six, "gap": gap, "tgap_d3": tgap, "roff": roff,
           "alphas": alphas, "cvals": cvals, "crunode": {"t1": t1c, "t2": t2c, "s": sc, "r": rc},
           "LC2": bool(ok2), "LC5": bool(ok5)}, open("atlas11_stageD.json","w"), indent=1)
print(f"[stage D done {time.time()-t0:.0f}s]")
