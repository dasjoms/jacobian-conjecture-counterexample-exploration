"""
Note 9, stage C: escape physics + real side for F7 (fiber OCTIC - even, cone returns).
"""
import numpy as np, json
from collections import Counter

def h_coeffs(s, r): return [-1/8, 1/7, 0, 0, 0, -27/28, 53/56, -s, r]
def p7w(w): return -w**7 + w**6 - 81*w**2/28 + 53*w/28
def preimages(A, B, C):
    s, r = B*C, A*C**2
    out = []
    for wv in np.roots(h_coeffs(s, r)):
        g = s - p7w(wv)
        if abs(g) > 1e-12:
            out.append((wv, g, C/g))
    return out

print("== 1. fiber counts (bounded |x| <= 1e4) ==")
s_f = 0.23102845535714287
import sympy as sp
w_ = sp.symbols("w")
half = sp.Rational(1,2)
p7 = -w_**7 + w_**6 - sp.Rational(81,28)*w_**2 + sp.Rational(53,28)*w_
Phi7 = sp.expand(sp.integrate(p7, w_))
s_f = float(p7.subs(w_, half)); r_f = float(half*p7.subs(w_, half) - Phi7.subs(w_, half))
CUSP = (0.3104763509023746, 0.034001356659445606)
CRU = (-0.90945711, -0.91031417)
samples = {"generic (1,1,1)": (1,1,1), "fold": (r_f, s_f, 1), "cusp": (CUSP[1], CUSP[0], 1),
           "crunode": (CRU[1], CRU[0], 1)}
for name, (A,B,C) in samples.items():
    ps = preimages(A,B,C)
    fin = [p_ for p_ in ps if abs(p_[2]) <= 1e4]
    print(f"  {name:18s}: {len(fin)} bounded (of {len(ps)})  expect 8/6/5/4")

print("\n== 2. escape rates (delta [1e-8,1e-4]) ==")
def escape_scan(bs, br):
    gam, delt = [], []
    for k in np.arange(4, 8.001, 0.5):
        d = 10.0**(-k)
        s, r = bs + d, br + 0.7*d
        g = [abs(s - p7w(wv)) for wv in np.roots(h_coeffs(s, r))]
        gam.append(min(g)); delt.append(d*np.hypot(1.0, 0.7))
    return np.array(delt), np.array(gam)
for name, (bs, br) in {"fold": (s_f, r_f), "cusp": CUSP}.items():
    delt, gam = escape_scan(bs, br)
    fit = np.polyfit(np.log(delt), np.log(gam), 1)
    print(f"  {name}: |gamma| ~ {np.exp(fit[1]):.4f} * delta^{fit[0]:.4f}")
    json.dump({"delta": list(delt), "gamma": list(gam), "slope": float(fit[0]),
               "const": float(np.exp(fit[1]))}, open(f"atlas7_escape_{name}.json","w"))

print("\n== 3. off-wall boundedness (20k) ==")
rng = np.random.default_rng(11)
mx = 0.0
for _ in range(20000):
    A = rng.normal(0,3) + 1j*rng.normal(0,3); B = rng.normal(0,3) + 1j*rng.normal(0,3)
    for wv, g, x in preimages(A,B,1):
        mx = max(mx, abs(x))
print(f"  max |x| = {mx:.6g}")

print("\n== 4. real census (200k) - even chamber, expect 0/2/4/6/8 ==")
cnt = Counter()
A = rng.normal(0,1.5,200000); B = rng.normal(0,1.5,200000)
for a, b in zip(A,B):
    roots = np.roots(h_coeffs(b, a))
    np_ = sum(1 for wv in roots if abs(wv.imag) < 1e-9 and abs(b - p7w(wv.real)) > 1e-9)
    cnt[np_] += 1
tot = sum(cnt.values())
for k in sorted(cnt): print(f"  {k}: {cnt[k]:8d} ({100*cnt[k]/tot:.3f}%)")
print(f"  MISSED: {100*cnt[0]/tot:.4f}%  [cone returns]")
json.dump({str(k): v for k,v in sorted(cnt.items())}, open("atlas7_realcensus.json","w"))

print("\n== 5. region grid (200x200, counts in {0,2,4,6,8}) ==")
n = 200
S = np.linspace(-4,4,n); R = np.linspace(-4,4,n)
grid = np.zeros((n,n), dtype=int)
for i, rv in enumerate(R):
    for j, sv in enumerate(S):
        grid[i,j] = sum(1 for wv in np.roots(h_coeffs(sv, rv)) if abs(wv.imag) < 1e-7)
vals, cts = np.unique(grid, return_counts=True)
print("  grid counts:", dict(zip(vals.tolist(), cts.tolist())))
np.savez("atlas7_grid.npz", S=S, R=R, grid=grid)

print("\n== 6. C=0 sweep + hinge (p1 = 53/28, q2 = 53/56) ==")
ps = preimages(2.0, 3.0, 1e-5)
fin = sorted(ps, key=lambda p_: abs(p_[2]))
print("  by |x|:", [complex(round(p_[2].real,5), round(p_[2].imag,5)) for p_ in fin])
# hinge: 2(1+53u/28)^2 = 9u(1+53u/56)  ->  multiply 784:  2*784(1+53u/28)^2 - 784*9u(1+53u/56)
#  1568 + 2*2*53*28 u + 2*53^2 u^2 - 7056u - 378*53/(2)*??  -> use sympy
u_ = sp.symbols("u_")
hinge_eq = sp.expand(2*(1+sp.Rational(53,28)*u_)**2 - 9*u_*(1+sp.Rational(53,56)*u_))
print("  hinge polynomial:", sp.Poly(hinge_eq*784, u_).as_expr())
us = np.roots([float(c) for c in sp.Poly(hinge_eq, u_).all_coeffs()])
xs = sorted([complex(-(1+53*u/28)/3) for u in us], key=lambda z: z.real)
print("  hinge x:", [complex(round(v.real,5), round(v.imag,5)) for v in xs],
      "antipodal:", abs(xs[0]+xs[1]) < 1e-9)
