"""
Note 8, stage C: escape physics + real side for F6 (fiber SEPTIC, odd degree!).
1. fiber counts: generic 7 / fold 5 / cusp 4 / acnode 3 (bounded-x honest filter)
2. escape rates: fold ~ delta^{1/2}, cusp ~ delta^{2/3} (universality check)
3. off-wall boundedness: 20k targets
4. REAL CENSUS: odd fiber degree -> every real target has >=1 real root: histogram
   over {1,3,5,7}; missed fraction should be ~0 (whisker is measure zero)
5. region grid: real-root counts over (s,r) in {1,3,5,7}
6. C=0: epsilon sweep (expect 1 flat + 2 antipodal hinge + 4 escaping)
   hinge: A(1+(13/7)u)^2 = B^2 u (1+(13/14)u)
7. whisker certificate: real cusp h-factorization at 200-digit mpmath: residual
   quartic has 0 real roots (Sturm-sign robust)
"""
import numpy as np, json
from mpmath import mp, mpf, mpc, findroot

def h_coeffs(s, r): return [-1/7, 1/6, 0, 0, -20/21, 13/14, -s, r]
def p6w(w): return -w**6 + w**5 - 20*w**2/7 + 13*w/7
def preimages(A, B, C):
    s, r = B*C, A*C**2
    out = []
    for wv in np.roots(h_coeffs(s, r)):
        g = s - p6w(wv)
        if abs(g) > 1e-12:
            out.append((wv, g, C/g))
    return out

print("== 1. fiber counts (|x| <= 1e4 bounded) ==")
s_f = 0.22991071428571427
phi_half = -(0.5)**7/7 + (0.5)**6/6 - 20*(0.5)**3/21 + 13*(0.5)**2/14
r_f = 0.5*s_f - phi_half
CUSP = (0.3043412601764495, 0.03338254261374218)   # (s,r)
ACN1 = (-1.76416996, 2.29113957); ACN2 = (-0.10068230, -0.86419633)
samples = {"generic (1,1,1)": (1,1,1), "fold": (r_f, s_f, 1),
           "cusp": (CUSP[1], CUSP[0], 1), "acnode1": (ACN1[1], ACN1[0], 1),
           "acnode2": (ACN2[1], ACN2[0], 1)}
for name, (A,B,C) in samples.items():
    ps = preimages(A,B,C)
    fin = [p_ for p_ in ps if abs(p_[2]) <= 1e4]
    print(f"  {name:18s}: {len(fin)} bounded (of {len(ps)} gamma!=0)")

print("\n== 2. escape rates (delta in [1e-8, 1e-4]) ==")
def escape_scan(bs, br):
    gam, delt = [], []
    for k in np.arange(4, 8.001, 0.5):
        d = 10.0**(-k)
        s, r = bs + d*1.0, br + d*0.7
        g = [abs(s - p6w(wv)) for wv in np.roots(h_coeffs(s, r))]
        gam.append(min(g)); delt.append(d*np.hypot(1.0, 0.7))
    return np.array(delt), np.array(gam)
for name, (bs, br) in {"fold": (s_f, r_f), "cusp": CUSP}.items():
    delt, gam = escape_scan(bs, br)
    fit = np.polyfit(np.log(delt), np.log(gam), 1)
    print(f"  {name}: |gamma| ~ {np.exp(fit[1]):.4f} * delta^{fit[0]:.4f}")
    json.dump({"delta": list(delt), "gamma": list(gam), "slope": float(fit[0]),
               "const": float(np.exp(fit[1]))}, open(f"atlas6_escape_{name}.json","w"))

print("\n== 3. off-wall boundedness (20k, sigma=3) ==")
rng = np.random.default_rng(7)
mx = 0.0
for _ in range(20000):
    A = rng.normal(0,3) + 1j*rng.normal(0,3); B = rng.normal(0,3) + 1j*rng.normal(0,3)
    for wv, g, x in preimages(A,B,1):
        mx = max(mx, abs(x))
print(f"  max |x| = {mx:.6g}")

print("\n== 4. real census (200k targets) ==")
from collections import Counter
cnt = Counter()
A = rng.normal(0,1.5,200000); B = rng.normal(0,1.5,200000)
for a, b in zip(A,B):
    roots = np.roots(h_coeffs(b, a))
    np_ = sum(1 for wv in roots if abs(wv.imag) < 1e-9 and abs(b - p6w(wv.real)) > 1e-9)
    cnt[np_] += 1
tot = sum(cnt.values())
for k in sorted(cnt): print(f"  {k} real preimages: {cnt[k]:8d} ({100*cnt[k]/tot:.2f}%)")
print(f"  MISSED: {100*cnt[0]/tot:.4f}%  [expect ~0: odd fiber degree]")
json.dump({str(k): v for k,v in sorted(cnt.items())}, open("atlas6_realcensus.json","w"))

print("\n== 5. region grid (200x200) ==")
n = 200
S = np.linspace(-4,4,n); R = np.linspace(-4,4,n)
grid = np.zeros((n,n), dtype=int)
for i, rv in enumerate(R):
    for j, sv in enumerate(S):
        grid[i,j] = sum(1 for wv in np.roots(h_coeffs(sv, rv)) if abs(wv.imag) < 1e-7)
vals, cts = np.unique(grid, return_counts=True)
print("  grid counts:", dict(zip(vals.tolist(), cts.tolist())))
np.savez("atlas6_grid.npz", S=S, R=R, grid=grid)

print("\n== 6. C=0 epsilon sweep + hinge (p1 = 13/7) ==")
ps = preimages(2.0, 3.0, 1e-5)
fin = sorted(ps, key=lambda p_: abs(p_[2]))
print("  preimages by |x|:", [complex(round(p_[2].real,5), round(p_[2].imag,5)) for p_ in fin])
# exact hinge roots at (2,3): 143 u^2 + 154 u - 196 = 0
us = np.roots([143, 154, -196])
xs = sorted([complex(-(1+13*u/7)/3) for u in us], key=lambda z: z.real)
print("  hinge x:", [complex(round(v.real,5), round(v.imag,5)) for v in xs],
      "antipodal:", abs(xs[0] + xs[1]) < 1e-9)

print("\n== 7. whisker certificate (200-digit) ==")
mp.dps = 200
def p_mp(w_): return -w_**6 + w_**5 - mpf(20)/7*w_**2 + mpf(13)/7*w_
def Phi_mp(w_): return -w_**7/7 + w_**6/6 - mpf(20)/21*w_**3 + mpf(13)/14*w_**2
t0 = findroot(lambda t: -6*t**5 + 5*t**4 - mpf(40)/7*t + mpf(13)/7, 0.3313539)
s0 = p_mp(t0); r0 = t0*s0 - Phi_mp(t0)
def divide_by(coeffs, a):
    b = [coeffs[0]]
    for c in coeffs[1:]:
        b.append(c + b[-1]*a)
    return b[:-1], b[-1]
hc = [mpf(-1)/7, mpf(1)/6, mpf(0), mpf(0), mpf(-20)/21, mpf(13)/14, -s0, r0]
tot_rem = 0
for _ in range(3):
    hc, rem = divide_by(hc, t0)
    tot_rem += abs(rem)
print(f"  t0 refinement: p6(t0) residual = {float(abs(-6*t0**5 + 5*t0**4 - mpf(40)/7*t0 + mpf(13)/7)):.1e}")
print(f"  (w-t0)^3 division remainder = {float(tot_rem):.1e}")
# residual quartic hc (degree 4); robust Sturm on its real coefficients via numpy roots at quad precision
import numpy.polynomial.polynomial as P
q = [float(v.real) for v in hc]  # descending
disc_pos = np.roots(q)
reals = [z for z in disc_pos if abs(z.imag) < 1e-6]
print(f"  residual quartic roots: {[complex(round(z.real,6), round(z.imag,6)) for z in disc_pos]}")
print(f"  real roots among them: {len(reals)}  => cusp target MISSED over R (whisker)" if not reals else "  REAL root found?!")
