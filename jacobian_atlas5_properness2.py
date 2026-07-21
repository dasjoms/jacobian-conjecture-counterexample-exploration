"""
Note 7, stage C addendum: honest recounts.
1. stratification with bounded-x filter (|x| <= 1e4)
2. escape exponent refit on delta in [1e-8, 1e-4]
3. C=0 certificate: quadratic hinge 81u^2+90u-100=0 vs epsilon-sweep, 3 random targets
4. node preimages: which pairs escape (both contacts?); acnode: complex contacts -> real target
   has no real roots among escaping pair? (count REAL roots at acnode)
"""
import numpy as np, json

def h_coeffs(s, r): return [-1/6, 1/5, 0, -14/15, 9/10, -s, r]
def p5w(w): return -w**5 + w**4 - 14*w**2/5 + 9*w/5
def preimages(A, B, C):
    s, r = B*C, A*C**2
    out = []
    for wv in np.roots(h_coeffs(s, r)):
        g = s - p5w(wv)
        if abs(g) > 1e-12:
            out.append((wv, g, C/g))
    return out

print("== 1. stratification, honest (|x| <= 1e4) ==")
s_f = 0.23125
phi_half = -(0.5)**6/6 + (0.5)**5/5 - 14*(0.5)**3/15 + 9*(0.5)**2/10
r_f = 0.5*s_f - phi_half
samples = {
    "generic (1,1,1)": (1,1,1),
    "fold exact": (r_f, s_f, 1),
    "cusp (3,1,1,1) -> 3": (0.03302630284803969, 0.2971583157757676, 1),
    "crunode (2,2,1,1) -> 2": (-0.88337174, -0.88188341, 1),
    "acnode (2,2,1,1) -> 2": (-0.46406377, 0.49202797, 1),
}
for name, (A,B,C) in samples.items():
    ps = preimages(A,B,C)
    fin = [p_ for p_ in ps if abs(p_[2]) <= 1e4]
    print(f"  {name:24s}: {len(fin)} bounded (of {len(ps)} with gamma!=0); |x| of excluded:"
          f" {sorted([round(abs(p_[2]),2) for p_ in ps if abs(p_[2])>1e4])[:4]}")

print("\n== 2. escape exponents refit, delta in [1e-8, 1e-4] ==")
def escape_scan(bs, br, ds, dr):
    gam, delt = [], []
    for k in np.arange(4, 8.001, 0.5):
        d = 10.0**(-k)
        s, r = bs + d*ds, br + d*dr
        g = [abs(s - p5w(wv)) for wv in np.roots(h_coeffs(s, r))]
        gam.append(min(g)); delt.append(d*np.hypot(ds, dr))
    return np.array(delt), np.array(gam)
for name, (bs, br) in {"fold": (s_f, r_f), "cusp": (0.2971583157757676, 0.03302630284803969)}.items():
    delt, gam = escape_scan(bs, br, 1.0, 0.7)
    fit = np.polyfit(np.log(delt), np.log(gam), 1)
    print(f"  {name}: |gamma| ~ {np.exp(fit[1]):.4f} * delta^{fit[0]:.4f}")

print("\n== 3. C=0 hinge certificate (p1 = 9/5, q2 = 9/10) ==")
rng = np.random.default_rng(7)
for (A, B) in [(2.0, 3.0)] + [(float(rng.normal(0,2)), float(rng.normal(0,2))) for _ in range(3)]:
    # quadratic: A(1+(9/5)u)^2 = B^2 u (1+(9/10)u)
    # expand: (A*81/25 - B^2*9/10) u^2 + (A*18/5 - B^2) u + A = 0
    qa = A*81/25 - B**2*9/10; qb = A*18/5 - B**2; qc = A
    us = np.roots([qa, qb, qc]) if abs(qa) > 1e-12 else np.roots([qb, qc])
    xs = [(-(1+1.8*u)/B if abs(B) > 1e-12 else np.nan) for u in us]
    # epsilon sweep: actual preimages for eps = 1e-5
    ps = preimages(A, B, 1e-5)
    fin = sorted([p_[2] for p_ in ps if abs(p_[2]) < 1e6], key=lambda z: abs(z))
    print(f"  (A,B)=({A:.4f},{B:.4f}): hinge x = {[complex(round(x.real,5), round(x.imag,5)) for x in xs]}"
          f"  |  eps-sweep 3 smallest |x|: {[complex(round(v.real,5), round(v.imag,5)) for v in fin[:3]]}")

print("\n== 4. real roots at strata (acnode vs crunode) ==")
for name, (A,B) in {"crunode": (-0.88337174, -0.88188341), "acnode": (-0.46406377, 0.49202797)}.items():
    roots = np.roots(h_coeffs(B, A))
    nreal = sum(1 for wv in roots if abs(wv.imag) < 1e-8)
    print(f"  {name}: {nreal} real roots of h (real fiber multiplicities visible: "
          f"{sorted([complex(round(z.real,5), round(z.imag,5)) for z in roots], key=lambda z: z.real)})")
