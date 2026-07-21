"""
NOTE 15, stage 3: THE ASYMPTOTIC LAWS (certified) + validation-set locks.
==========================================================================
Structure: fit ONLY on d <= 11 (and d <= 41 cusp data already audited);
predict; THEN compute validation d = 13, 15 (shadow), d = 47 (ghost),
d = 13,15 (left cusp) and test the locks.
LOCKS (pre-computation):
 (A1) |U(13)| := |t*(13)+1| in (0.060, 0.068);  |U(15)| in (0.052, 0.058);
      s*(13)+1 in (-0.056, -0.050).
 (A2) ghost: 47^2 (1/3 - t_g(47)) in (0.3247, 0.3268).
 (A3) left-cusp fixed-point (recurrence u -> 1 - R(u)^{1/(d-2)})
      predicts exact u(13) to within 3e-4; also first-order law
      u*(d-2)/ln((d-2)/4) in (0.9, 1.35) for d = 13, 15.
 (F1) chamber-12 locks: crunode s(11) in (-0.9424, -0.9409),
      r(11) in (-0.9428, -1.9413+1.0->) = (-0.9428, -0.9413);
      whisker acnode s(12) in (-0.945, -0.895), r(12) in (-1.16, -1.10).
 (M3) missed mass census d = 3: m(3) in [8.45%, 8.85%].
Certification duties: ghost series a1 = -1/3 c0-normalized, a2 = 2/3 [from
v = -1/(3d(d+1))(1 + 2/(d(d+1)) + ...) match at d=41: 0.325775 predicted];
corner-shadow A1 = ln 2 and A2 derived + residuals audit; contact-gap series
delta(d) = s_c - s* ~ c/d^alpha, alpha -> 2 check; tau-slopes.
==========================================================================
"""
import json, time
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, diff, integrate, Rational as R, Poly

mp.mp.dps = 80
w = symbols("w")
t0 = time.time()
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

def diag_root(d):
    """exact near-corner diagonal root (t*,s*) via np + polish on G = Phi-(w-1)p."""
    p = seed(d); Phi = integrate(p, w)
    G = expand(Phi - (w-1)*p)
    P = Poly(sp.expand(G*(d*(d+1))), w)
    fc = [float(c) for c in P.all_coeffs()]
    Gf = sp.lambdify(w, G, "mpmath")
    Gfp = sp.lambdify(w, diff(G, w), "mpmath")
    best = None
    for z in np.roots(fc):
        if abs(z.imag) < 1e-6 and z.real < -1.0001:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(Gf, xr, tol=mp.mpf("1e-70"), verify=False)
            except Exception: pass
            if abs(Gf(xr)) < mp.mpf("1e-40") and (best is None or xr > best):
                best = xr
    sv = sp.lambdify(w, p, "mpmath")(best)
    return best, sv

# ---------------- data d = 3..11 (fit set) ----------------
fitset = {}
for d in (3, 5, 7, 9, 11):
    tv, sv = diag_root(d)
    fitset[d] = (tv, sv)
    print(f"  d={d:2d}: t* = {mp.nstr(tv, 20)}   s* = {mp.nstr(sv, 20)}", flush=True)
    # audit this is the same as stage-1's (for d<=11)
U = {d: abs(tv + 1) for d, (tv, _) in fitset.items()}
print("-"*84)
for d in fitset:
    print(f"  |U({d})|={mp.nstr(U[d],12)}  |U|*(d-2)={mp.nstr(U[d]*(d-2),10)}  "
          f"s*+1 = {mp.nstr(fitset[d][1]+1,12)}  (s*+1)(d-2) = {mp.nstr((fitset[d][1]+1)*(d-2),10)}", flush=True)

# ---------------- (B) corner-shadow certified series ----------------
print("="*84); print("(B) corner-shadow asymptotics: cluster extraction"); print("-"*84)
# G(t) = Phi - (t-1)p ; evaluate at t = -1-U, separate the (1+U)^(d-2) cluster.
d_sym = sp.symbols("d", positive=True, integer=True)
U_ = sp.symbols("U")
# p_d and Phi_d as explicit families:
c0 = 6/(d_sym*(d_sym+1))
wU = -1 - U_
# (w)^(d-2) for odd d: = (-1-U)^(d-2) = -(1+U)^(d-2); (w)^d = -(1+U)^d; (w)^(d+1) = +(1+U)^(d+1)
pwU = 2*wU - 3*wU**2 + wU*(1-wU)*(-(1+U_)**(d_sym-2) - c0)          # odd d
PhiU = wU**2 - wU**3 - (1+U_)**d_sym/d_sym - (1+U_)**(d_sym+1)/(d_sym+1) \
       - c0*(wU**2/2 - wU**3/3)
GU = sp.expand(PhiU - (wU - 1)*pwU)
# isolate the cluster: collect powers of (1+U)^(d-2)-family on a set-spectrum
X = sp.symbols("X")
GUx = GU.subs((1+U_)**(d_sym-2), X).subs((1+U_)**d_sym, X*(1+U_)**2).subs(
      (1+U_)**(d_sym+1), X*(1+U_)**3)
GUx = sp.expand(GUx)
E = sp.Poly(GUx, X).coeff_monomial(X)          # cluster coeff (poly in U, d)
Pp = sp.Poly(GUx, X).coeff_monomial(1)         # the rest
E = sp.expand(E); Pp = sp.expand(Pp)
print(f"  polynomial part  P(U) = {Pp}")
print(f"  cluster coeff    E(U) = {E}")
# numeric sanity at d=11, U=0.0834..
GUf = sp.lambdify((U_, d_sym), GU, "mpmath")
EPf = sp.lambdify((U_, d_sym), E*(1+U_)**(d_sym-2) + Pp, "mpmath")
for dd, uu in ((11, mp.mpf("0.083458621")), (7, mp.mpf("0.155950961"))):
    lhs = GUf(uu, dd)
    rhs = EPf(uu, dd)
    print(f"  d={dd}: G(-1-U) = {mp.nstr(lhs,6)}   P + E(1+U)^(d-2) = {mp.nstr(rhs,6)}  match {abs(lhs-rhs)<1e-40}")
# asymptotic balance: solve P(U,d) + E(U,d)*(1+U)^(d-2) = 0 with U = sum A_k L^k (d-2)^-k
# Ansatz: U = a1/dn + a2/dn^2, dn := d-2, a1 = ln2 expected.
dn = sp.symbols("dn", positive=True)
eps = sp.symbols("eps")   # 1/dn
a1, a2 = sp.symbols("a1 a2")
Use = a1*eps + a2*eps**2
ln1pU = sp.series(sp.log(1+sp.symbols('z')), sp.symbols('z'), 0, 3).removeO().subs(sp.symbols('z'), Use)
ln1pU = sp.expand(ln1pU)
Xser = sp.exp(ln1pU/eps)                      # (1+U)^(d-2) = exp(ln(1+U)/eps)
Xser = sp.expand(Xser)
# print series of exponent
print(f"  ln(1+U)/eps = {sp.expand(ln1pU/eps)}")
print(f"  (1+U)^(d-2) ~ {Xser}")
Pser = sp.series(Pp.subs({U_: Use, d_sym: dn+2}), eps, 0, 3).removeO()
Eser = sp.series(E.subs({U_: Use, d_sym: dn+2}) * Xser, eps, 0, 3).removeO()
tot = sp.expand(Pser.subs(dn, 1/eps) + Eser.subs(dn, 1/eps))
tot = sp.expand(tot)
print(f"  G series (orders): const: {sp.expand(tot.subs(eps,0))}")
c1 = sp.expand(tot.coeff(eps, 0)); c2 = sp.expand(tot.coeff(eps, 1)); c3 = sp.expand(tot.coeff(eps, 2))
print(f"   O(1):   {c1}")
print(f"   O(eps): {c2}")
print(f"   O(eps^2): {c3}")
s1 = sp.solve(sp.Eq(c1, 0) if c1 != 0 else sp.Eq(c2, 0), a1)
print(f"   balance1: a1 = {s1}")
out["shadow_cluster"] = {"P": str(Pp), "E": str(E), "c_orders": [str(c1), str(c2), str(c3)]}

# s*-series: s* = p(t*) ; use the same machinery
s_use = pwU.subs(d_sym, dn+2).subs(U_, Use).subs((1+Use)**(dn), Xser)  # (1+U)^(d-2)->Xser
sser = sp.expand(sp.series(s_use.subs(dn, 1/eps), eps, 0, 2).removeO())
print(f"  s* series: {sser}")
json.dump(out, open("jcorner_stage3.json", "w"))
print(f"checkpoint B [{time.time()-t0:.0f} s]")
