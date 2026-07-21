"""
Note 10, stage D: escape physics + real side for F8 (fiber NONIC - odd chamber, whiskers).
Predictions: counts at generic/fold/cusp = 9/7/6; whisker tip: 0 real preimages;
census odd counts only, missed 0.0000%; slopes 1/2 and 2/3; hinge antipodal.
"""
import numpy as np, json
from collections import Counter
import sympy as sp

def h_coeffs(s, r): return [-1/9, 1/8, 0, 0, 0, 0, -35/36, 23/24, -s, r]
def p8w(w): return -w**8 + w**7 - 35*w**2/12 + 23*w/12
def preimages(A, B, C):
    s, r = B*C, A*C**2
    out = []
    for wv in np.roots(h_coeffs(s, r)):
        g = s - p8w(wv)
        if abs(g) > 1e-12:
            out.append((wv, g, C/g))
    return out

w_ = sp.symbols("w")
half = sp.Rational(1,2)
p8 = -w_**8 + w_**7 - sp.Rational(35,12)*w_**2 + sp.Rational(23,12)*w_
Phi9 = sp.expand(sp.integrate(p8, w_))
s_f = float(p8.subs(w_, half)); r_f = float(half*p8.subs(w_, half) - Phi9.subs(w_, half))
CUSP = None
for wv in sp.nroots(sp.diff(p8, w_), n=60):
    if abs(sp.im(wv)) < sp.Float("1e-50"):
        sv = complex(p8.subs(w_, wv)); rv = complex(wv*sv - Phi9.subs(w_, wv))
        CUSP = (sv.real, rv.real)
print(f"fold ({s_f:.6f},{r_f:.6f})  whisker-cusp {CUSP}")

print("== 1. fiber counts (bounded |x| <= 1e4) ==")
samples = {"generic (1,1,1)": (1,1,1), "fold": (r_f, s_f, 1), "cusp": (CUSP[1], CUSP[0], 1)}
for name, (A,B,C) in samples.items():
    ps = preimages(A,B,C)
    fin = [p_ for p_ in ps if abs(p_[2]) <= 1e4]
    nreal = sum(1 for p_ in fin if abs(p_[0].imag) < 1e-7)
    print(f"  {name:18s}: {len(fin)} bounded (of {len(ps)})  real: {nreal}   expect 9/7/6")

print("== 1b. WHISKER TIP over R: cusp target should have NO real preimage ==")
ps = preimages(CUSP[1], CUSP[0], 1)
real_pre = [p_ for p_ in ps if abs(p_[0].imag) < 1e-7 and abs(p_[2]) <= 1e10]
print(f"  real preimages at whisker tip: {len(real_pre)}   [predict 0 - the missed tip]")

print("\n== 2. escape rates (delta [1e-8,1e-4]) ==")
def escape_scan(bs, br):
    gam, delt = [], []
    for k in np.arange(4, 8.001, 0.5):
        d = 10.0**(-k)
        s, r = bs + d, br + 0.7*d
        g = [abs(s - p8w(wv)) for wv in np.roots(h_coeffs(s, r))]
        gam.append(min(g)); delt.append(d*np.hypot(1.0, 0.7))
    return np.array(delt), np.array(gam)
for name, (bs, br) in {"fold": (s_f, r_f), "cusp": CUSP}.items():
    delt, gam = escape_scan(bs, br)
    fit = np.polyfit(np.log(delt), np.log(gam), 1)
    print(f"  {name}: |gamma| ~ {np.exp(fit[1]):.4f} * delta^{fit[0]:.4f}   [predict 1/2, 2/3]")
    json.dump({"delta": list(delt), "gamma": list(gam), "slope": float(fit[0]),
               "const": float(np.exp(fit[1]))}, open(f"atlas8_escape_{name}.json","w"))

print("\n== 3. off-wall boundedness (20k, complex targets) ==")
rng = np.random.default_rng(11)
mx = 0.0
for _i in range(20000):
    A = rng.normal(0,3) + 1j*rng.normal(0,3); B = rng.normal(0,3) + 1j*rng.normal(0,3)
    for wv, g, x in preimages(A,B,1):
        mx = max(mx, abs(x))
print(f"  max |x| = {mx:.6g}")

print("\n== 4. real census (200k, normal(0,1.5^2), C=1) - odd chamber, expect {1,3,5,7,9} ==")
cnt = Counter()
A = rng.normal(0,1.5,200000); B = rng.normal(0,1.5,200000)
for a, b in zip(A,B):
    roots = np.roots(h_coeffs(b, a))
    np_ = sum(1 for wv in roots if abs(wv.imag) < 1e-9 and abs(b - p8w(wv.real)) > 1e-9)
    cnt[np_] += 1
tot = sum(cnt.values())
for k in sorted(cnt): print(f"  {k}: {cnt[k]:8d} ({100*cnt[k]/tot:.3f}%)")
print(f"  MISSED (0 preimages): {100*cnt.get(0,0)/tot:.4f}%   [predict 0.0000% - whisker tip is measure zero]")
json.dump({str(k): v for k,v in sorted(cnt.items())}, open("atlas8_realcensus.json","w"))

print("\n== 5. region grid (200x200, real-root counts) ==")
n = 200
S = np.linspace(-4,4,n); R = np.linspace(-4,4,n)
grid = np.zeros((n,n), dtype=int)
for i, rv in enumerate(R):
    for j, sv in enumerate(S):
        grid[i,j] = sum(1 for wv in np.roots(h_coeffs(sv, rv)) if abs(wv.imag) < 1e-7)
vals, cts = np.unique(grid, return_counts=True)
print("  grid counts:", dict(zip(vals.tolist(), cts.tolist())))
np.savez("atlas8_grid.npz", S=S, R=R, grid=grid)

print("\n== 6. C=0 sweep + hinge (p1 = 23/12, q2 = 23/24) ==")
ps = preimages(2.0, 3.0, 1e-6)
fin = sorted(ps, key=lambda p_: abs(p_[2]))
print("  by |x|:", [complex(round(p_[2].real,5), round(p_[2].imag,5)) for p_ in fin])
u_ = sp.symbols("u_")
hinge_eq = sp.expand(2*(1+sp.Rational(23,12)*u_)**2 - 9*u_*(1+sp.Rational(23,24)*u_))
print("  hinge polynomial:", sp.Poly(hinge_eq*72, u_).as_expr())
us = np.roots([float(c) for c in sp.Poly(hinge_eq, u_).all_coeffs()])
xs = sorted([complex(-(1+23*u/12)/3) for u in us], key=lambda z: z.real)
print("  hinge x:", [complex(round(v.real,6), round(v.imag,6)) for v in xs],
      "antipodal:", abs(xs[0]+xs[1]) < 1e-9)
