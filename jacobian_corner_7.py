"""
NOTE 15, stage 5 (series engine + remaining audits).
FF := F*eps^2 = (1+2eps)(U+2)^2 (X-2) + eps^2 (4U^2+17U+19)
                - eps X (U+1)(U+2) - eps^2 X (U+1)(1+2(U+2)),  X = (1+U)^(1/eps).
Solve U(eps) = u1 eps + u2 eps^2 + u3 eps^3 + u4 eps^4 order by order EXACTLY.
Then: LSQ vs exact roots (d=9..45); s*-series; ghost(43)-finish; C3; M3(1e6); F1/F2 publish.
"""
import json, time
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, diff, integrate, Rational as R, Poly, ln, series, nsimplify

mp.mp.dps = 90
w = symbols("w")
t0 = time.time()
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

# ------------ load/recompute exact roots quickly ------------
roots = {}
def diag_root(d):
    p = seed(d); Phi = integrate(p, w)
    G = expand(Phi - (w-1)*p)
    fc = [float(c) for c in Poly(sp.expand(G*(d*(d+1))), w).all_coeffs()]
    Gf = sp.lambdify(w, G, "mpmath")
    best = None
    for z in np.roots(fc):
        if abs(z.imag) < 1e-7:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(Gf, xr, tol=mp.mpf("1e-80"), verify=False)
            except Exception: pass
            if abs(Gf(xr)) < mp.mpf("1e-40") and xr < -1.000001 and (best is None or xr > best):
                best = xr
    return best, sp.lambdify(w, p, "mpmath")(best)
for d in [3] + list(range(5, 46, 2)):
    roots[d] = diag_root(d)
print(f"  roots ok (d=3..45 odd) [{time.time()-t0:.0f}s]", flush=True)

# ------------ exact eps-series solve ------------
print("="*84); print("exact order-by-order solve")
eps = symbols("eps", positive=True)
u1, u2, u3, u4 = symbols("u1 u2 u3 u4")
U_ser = u1*eps + u2*eps**2 + u3*eps**3 + u4*eps**4
# log(1+U)/eps to eps^3:
z = sp.Symbol("z")
E = sp.expand((z - z**2/2 + z**3/3 - z**4/4 + z**5/5).subs(z, U_ser))
E = sp.expand(E/eps)
# explicit vi exponents
v1 = sp.expand(E).coeff(eps, 1)
v2 = sp.expand(E).coeff(eps, 2)
v3 = sp.expand(E).coeff(eps, 3)
DmT = v1*eps + v2*eps**2 + v3*eps**3
X_ser = sp.expand(sp.exp(u1)*(1 + DmT + DmT**2/2 + DmT**3/6 + DmT**4/24))
def trunc3(expr):
    expr = sp.expand(expr)
    P = sp.Poly(expr, eps)
    return sum(c*eps**k for (k,), c in P.terms() if k <= 3)
X_ser = trunc3(X_ser)
FF = trunc3((1+5*eps+6*eps**2)*(U_ser+2)**2*(X_ser-2) + eps**2*(4*U_ser**2+17*U_ser+19)
            - eps*X_ser*(U_ser+1)*(U_ser+2) - eps**2*X_ser*(U_ser+1)*(1+2*(U_ser+2)))
eqs = []
cur = FF
for i in range(4):
    c0v = cur.subs(eps, 0)
    cur = sp.expand(P if False else cur)
    if c0v != 0:
        c0v = sp.expand(c0v)
    eqs.append(c0v)
    cur = sp.expand(cur - c0v)
    cur = sp.expand(cur/eps)
E1 = sp.exp(u1)
print("  O(eps^0):", sp.collect(eqs[0], E1))
print("  O(eps^1):", sp.collect(eqs[1], E1))
print("  O(eps^2):", sp.collect(eqs[2], E1))
print("  O(eps^3):", sp.collect(eqs[3], E1)[:600] if len(str(eqs[3]))>600 else sp.collect(eqs[3], E1))
sol = {}
e0 = eqs[0].subs(sp.exp(u1), 2)
if sp.expand(e0) != 0:
    print("  !! eps^0 not killed by E1=2:", e0)
else:
    print("  eps^0: forces exp(u1)=2 ✓ (u1 = ln 2)")
s2 = sp.solve(sp.Eq(eqs[1].subs(sp.exp(u1), 2), 0), u2)[0]
sol[u2] = sp.expand(s2); print("  u2 =", sol[u2])
s3 = sp.solve(sp.Eq(sp.expand(eqs[2].subs({sp.exp(u1): 2, u2: sol[u2]})), 0), u3)[0]
sol[u3] = sp.nsimplify(sp.expand(s3), [ln(2)]) if not s3.has(u1) else sp.expand(s3)
sol[u3] = sp.expand(s3).subs(u1, ln(2)) if s3.has(u1) else s3
print("  u3 =", sp.expand(s3).subs(u1, ln(2)))
s4e = sp.expand(eqs[3].subs({sp.exp(u1): 2, u2: sol[u2], u3: sol[u3]}))
s4 = sp.solve(sp.Eq(s4e, 0), u4)[0]
print("  u4 =", sp.expand(s4).subs(u1, ln(2)))
out["u_series"] = {"u1": "ln(2)", "u2": str(sol[u2]), "u3": str(sp.expand(s3).subs(u1, ln(2))),
                   "u4": str(sp.expand(s4).subs(u1, ln(2)))}
N = lambda x: mp.mpf(str(sp.N(x.subs(u1, ln(2)) if hasattr(x, "subs") else x, 90)))
u2n, u3n, u4n = N(sol[u2]), N(s3), N(s4)
print(f"  numeric: u2 = {mp.nstr(u2n, 15)}, u3 = {mp.nstr(u3n, 15)}, u4 = {mp.nstr(u4n, 15)}")

# LSQ validation on exact data: U(d) * dn^2 - ln2*dn -> u2 + u3/dn + u4/dn^2
Ln2 = mp.log(2)
res_tab = []
for d, (tv, sv) in sorted(roots.items()):
    if d < 5: continue
    dn = d - 2
    U = -tv - 1
    pred = Ln2/dn + u2n/dn**2 + u3n/dn**3 + u4n/dn**4
    res = abs(U - pred)
    rel = res / (mp.mpf(1)/dn**2)
    res_tab.append((d, U, res, rel*dn**2))
    if d <= 15 or d % 4 == 1:
        print(f"  d={d:2d}: |U - series4| = {mp.nstr(res, 6)}   (x dn^4 scale {mp.nstr(res*dn**4, 6)})", flush=True)
maxrel = max(r for _, _, r, _ in res_tab if _ is not None or True)
out["U_series_maxres"] = mp.nstr(max(res for _, _, res, _ in [(0,0,0,0)] + [(d,U,res,rel) for d,U,res,rel in res_tab]), 6) if res_tab else None
maxr = max(res for _, _, res, rel in res_tab)
print(f"  max |U-series4| over d=5..45: {mp.nstr(maxr, 6)}")
out["U_series_ok"] = bool(maxr < mp.mpf("3e-5"))

# U-series at dn^4 scale should approach a constant
print("  dn^4-scaled residuals (-> u5?):", [(d, mp.nstr(rel, 5)) for d, _, _, rel in res_tab[-6:]])

# s*-series: sigma*(d) = s*+1 = -4 -8U -3U^2 + (2+3U+U^2)(X + c0), c0 = 6/(d(d+1))
print("="*84); print("s*+1 series")
c0_ser = 6*eps**2/((1+eps)*(1+2*eps))     # careful: d(d+1) = (1/eps+2)(1/eps+3)
c0_ser = sp.expand(sp.series(6/((1+2*eps)*(1+3*eps)), eps, 0, 2).removeO())*eps**2
Sser = sp.expand(-4 - 8*U_ser - 3*U_ser**2 + (2+3*U_ser+U_ser**2)*(X_ser + c0_ser))
Ss = sp.expand(series(Sser/eps*0 + Sser, eps, 0, 3).removeO())
# sig* = O(eps): divide by eps, series coefficients of sig*/eps
sig_over_eps = sp.expand(Sser)
cur = sig_over_eps
_sig = []
for i in range(3):
    c0v = sp.expand(cur.subs(eps, 0))
    _sig.append(c0v)
    cur = sp.expand((cur - c0v)/eps)
    cur = trunc3(cur)
sig_coeffs = _sig
print("  sig*+1 = (s1 eps + s2 eps^2 + s3 eps^3) with:")
for i, c in enumerate(sig_coeffs):
    cc = sp.expand(c.subs({sp.exp(u1): 2, u1: ln(2), u2: sol[u2], u3: sol[u3]}))
    ccv = sp.N(cc, 60)
    print(f"   s{i+1} = {cc}   ~ {ccv}")
    out.setdefault("s_coeffs", {})[i+1] = str(cc)
json.dump(out, open("jcorner_stage5_partial.json", "w"))
print(f"checkpoint [{time.time()-t0:.0f}s]")
