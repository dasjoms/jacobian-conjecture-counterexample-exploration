"""
Note 7, stage C: escape physics + real census for F5 (fiber sextic, C=1 => s=B, r=A).
1. fiber-count stratification: generic 6 / fold 4 / cusp 3 / node 2 preimages (gamma != 0 roots)
2. fold escape:  |gamma| ~ c * sqrt(delta)   => |x| ~ delta^{-1/2}
3. cusp escape:  |gamma| ~ c * delta^{2/3}   => |x| ~ delta^{-2/3}   [NEW exponent vs note 6]
4. off-wall boundedness: 20k random targets, max |x|
5. real census: real preimages of random real targets (A,B,1)
6. C=0: flat-sheet restriction + epsilon sweep of finite preimages
"""
import numpy as np, json
rng = np.random.default_rng(66)

def h_coeffs(s, r):
    return [-1/6, 1/5, 0, -14/15, 9/10, -s, r]
def p5w(w):
    return -w**5 + w**4 - 14*w**2/5 + 9*w/5

def preimages(A, B, C):
    """return list of (w, gamma, x) with gamma != 0"""
    s, r = B*C, A*C**2
    roots = np.roots(h_coeffs(s, r))
    out = []
    for wv in roots:
        g = s - p5w(wv)
        if abs(g) > 1e-9:
            out.append((wv, g, C/g))
    return out

print("== 1. fiber-count stratification (C=1 targets) ==")
samples = {
    "generic (A=1,B=1)": (1,1,1),
    "fold (s=0.23125,r=0.0036458+0.001i perturb?)": None,
    "cusp (0.0330263, 0.2971583)": (0.03302630284803969, 0.2971583157757676, 1),
    "crunode (-0.88337174, -0.88188341)": (-0.88337174, -0.88188341, 1),
    "acnode (r=-0.46406377,s=0.49202797)": (-0.46406377, 0.49202797, 1),
}
# fold target: on-wall exact parametrization t=1/2
s_f = 0.23125
phi_half = -(0.5)**6/6 + (0.5)**5/5 - 14*(0.5)**3/15 + 9*(0.5)**2/10
r_f = 0.5*s_f - phi_half
samples["fold (exact, A=r, B=s)"] = (r_f, s_f, 1)
for name, tgt in samples.items():
    if tgt is None: continue
    A,B,C = tgt
    ps = preimages(A,B,C)
    print(f"  {name:45s} -> {len(ps)} preimages")

print("\n== 2/3. escape rates (C=1; base wall points) ==")
def escape_scan(base_s, base_r, dir_s, dir_r, powers):
    gam = []
    for k in powers:
        d = 10.0**(-k)
        s, r = base_s + d*dir_s, base_r + d*dir_r
        roots = np.roots(h_coeffs(s, r))
        g = [abs(s - p5w(wv)) for wv in roots]
        gam.append(min(g))                    # escaping sheet: gamma -> 0
    return np.array(gam)

powers = np.arange(2, 8.001, 0.5)
for name, (bs, br) in {"fold": (s_f, r_f),
                       "cusp": (0.2971583157757676, 0.03302630284803969)}.items():
    gam = escape_scan(bs, br, 1.0, 0.7, powers)
    delt = 10.0**(-powers)*np.hypot(1.0, 0.7)
    fit = np.polyfit(np.log(delt), np.log(gam), 1)
    print(f"  {name}: gamma ~ delta^{fit[0]:.4f}  (const {np.exp(fit[1]):.4f})   => |x| ~ delta^{{-fit[0]:.4f}}")
    json_name = f"atlas5_escape_{name}.json"
    json.dump({"delta": list(delt), "gamma": list(gam), "slope": float(fit[0])}, open(json_name, "w"))

print("\n== 4. off-wall boundedness survey (C=1, 20000 targets) ==")
mx, rec = 0.0, None
for _ in range(20000):
    A = rng.normal(0, 3) + 1j*rng.normal(0, 3)
    B = rng.normal(0, 3) + 1j*rng.normal(0, 3)
    for wv, g, x in preimages(A, B, 1):
        if abs(x) > mx:
            mx, rec = abs(x), (A, B)
print(f"  max |x| over all preimages: {mx:.6g}")

print("\n== 5. real census (C=1, 200000 real targets (A,B)) ==")
from collections import Counter
cnt = Counter()
A = rng.normal(0, 1.5, 200000); B = rng.normal(0, 1.5, 200000)
for a, b in zip(A, B):
    roots = np.roots([ -1/6, 1/5, 0, -14/15, 9/10, -b, a ])
    npre = 0
    for wv in roots:
        if abs(wv.imag) < 1e-9 and abs(b - p5w(wv.real)) > 1e-9:
            npre += 1
    cnt[npre] += 1
tot = sum(cnt.values())
for k in sorted(cnt):
    print(f"  {k} real preimages: {cnt[k]:8d}  ({100*cnt[k]/tot:.2f}%)")
print(f"  MISSED (0 real preimages): {100*cnt[0]/tot:.3f}%   [note 5 reference: 451/4000 = 11.275%]")
json.dump({str(k): v for k, v in sorted(cnt.items())}, open("atlas5_realcensus.json","w"))

print("\n== 6. C=0 slice ==")
import sympy as sp
x, y, z, w = sp.symbols("x y z w")
p = 2*w - 3*w**2 + w*(1-w)*(w**3 - sp.Rational(1,5))
q = sp.expand(sp.integrate(w*sp.diff(p, w), w))
a = sp.Rational(-19, 14)
u = 1 + x*y; g_expr = 1 + a*x*y + x**2*z; ws = u*g_expr
f1 = sp.expand(sp.cancel((u + q.subs(w, ws)/g_expr**2)/x**2))
f2 = sp.expand(sp.cancel((1 + p.subs(w, ws)/g_expr)/x))
print("  flat sheet f1(0,y,z) degree:", sp.Poly(f1.subs(x,0), y, z).total_degree(),
      " f2(0,y,z) degree:", sp.Poly(f2.subs(x,0), y, z).total_degree())
J2 = sp.Matrix([f1.subs(x,0), f2.subs(x,0)]).jacobian([y, z]).det()
print("  2D Keller check det J(flat) =", sp.factor(sp.expand(J2)) if J2 != 0 else 0, flush=True)
print("  epsilon-sweep, target (A=2, B=3, eps): #finite preimages (|x|<1e4):")
for k in (2, 3, 4, 6, 8, 10):
    eps = 10.0**(-k)
    ps = preimages(2.0, 3.0, eps)
    fin = [px for px in ps if abs(px[2]) < 1e4]
    print(f"    eps=1e-{k}: {len(fin)} finite of {len(ps)}   x-values {[complex(round(p_[2].real,6), round(p_[2].imag,6)) for p_ in fin[:6]]}")
