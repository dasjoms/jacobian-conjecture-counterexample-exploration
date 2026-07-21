"""
Note 20, stage C: modality, census of the whole tour, escape physics, det, corner polish.
================================================================
LOCKED PREDICTIONS:
  LC1 counts (bounded |x|<=1e4): generic 12 / fold 10 / cusp 9 / crunode 8 / pocket-acnode 8
  LC2 census (fresh 400k, N(0,1.5^2), C=1): d=11 cone-0 in [8.67, 8.84]% (law center 8.7556);
      4-real in [7.2, 7.7]%; 6-real < 0.01%; fresh-vs-archived |Delta| <= 0.2pp for d=5..10
      (archived: 8.693/8.737/8.726, 3-root 17.46/17.40/17.35)
  LC3 escape slopes at fold/cusp: 1/2 and 2/3 (+-0.03)
  LC4 det JF_12 = 1 exactly at 30 rational points + hinge p1 = 2 q2 exact
  LC5 coalescence post-checks (note 15 pre-reg): gap d11 in [1.3,2.1]e-4 (est 1.67e-4);
      contact |t1-t*|*11^3 in [5.0,7.5]e-3; r-offset (r_c-s_c) in [-4.5,-2.0]e-4
  LC6 no escape above the noise everywhere off-wall (max|x| over 20k targets finite, report)
================================================================
"""
import numpy as np, json, time
from collections import Counter
import sympy as sp
import mpmath as mp
mp.mp.dps = 120
t0 = time.time()

def p_d_np(w, d):
    m = d*(d+1)
    return (-w**d + w**(d-1) - (3-6/m)*w**2  + (2-6/m)*w)*1.0 if False else \
           (-w**d + w**(d-1) - (3 - 6/m)*w**2 + (2 - 6/m)*w)
# sanity: p_d(w) = 2w -3w^2 + w^{d-1} - w^{d} -(6/m)w + (6/m)w^2 = -w^d + w^{d-1} -(3-6/m)w^2 + (2-6/m)w
def phi_coeffs(d):
    m = d*(d+1)
    co = np.zeros(d+2, float)
    co[0] = -1.0/(d+1); co[1] = 1.0/d
    co[d-2] = -(1 - 2/m); co[d-1] = (1 - 3/m)   # positions: index k <-> power d+1-k
    # descending: idx 0 -> w^{d+1}, idx 1 -> w^d, ..., idx d-1 -> w^2, idx d -> w^1, idx d+1 -> w^0
    return co
def p_d_exact(d):
    w_ = sp.symbols("w"); m = d*(d+1)
    return -w_**d + w_**(d-1) - sp.Rational(3*m-6, m)*w_**2 + sp.Rational(2*m-6, m)*w_
# cross-check templates at d=11 against stage A's exact forms
w_ = sp.symbols("w")
assert sp.expand(p_d_exact(11) - (-w_**11 + w_**10 - sp.Rational(65,22)*w_**2 + sp.Rational(43,22)*w_)) == 0
co11 = phi_coeffs(11)
assert abs(co11[0] + 1/12) < 1e-15 and abs(co11[1] - 1/11) < 1e-15 and abs(co11[9] + 65/66) < 1e-14 and abs(co11[10] - 43/44) < 1e-14

def h_coeffs(d, s, r):
    co = phi_coeffs(d).astype(complex if isinstance(s+r, complex) else float)
    co[-2] = -s; co[-1] = r; return co

# ---- sample points ----
p11 = p_d_exact(11); Phi12 = sp.expand(sp.integrate(p11, w_))
half = sp.Rational(1,2)
s_f = float(p11.subs(w_, half)); r_f = float(half*p11.subs(w_, half) - Phi12.subs(w_, half))
_c = [complex(v) for v in sp.nroots(sp.diff(p11, w_), n=60, maxsteps=8000)]
_c = sorted([z.real for z in _c if abs(z.imag) < 1e-40])
CU = []
for tc in _c:
    CU.append((float(p11.subs(w_, sp.Float(repr(tc), 40))),
               float(sp.re(sp.Float(repr(tc),40)*p11.subs(w_, sp.Float(repr(tc),40)) - Phi12.subs(w_, sp.Float(repr(tc),40))))))
bt = json.load(open("atlas11_bitangents.json"))
ND = []
for a, b, sv, rv, k in bt["nodes"]:
    sv_ = complex(sv.replace(" ","")); rv_ = complex(rv.replace(" ",""))
    ND.append((sv_.real, rv_.real, k))
CR = next((sv, rv) for sv, rv, k in ND if k == "CRUNODE")
ACp = min(((sv, rv) for sv, rv, k in ND if k == "ACNODE"), key=lambda q: abs(q[0]+0.8937))
print(f"fold ({s_f:.6f},{r_f:.6f})  cusps {[tuple(round(v,6) for v in c) for c in CU]}")
print(f"crunode ({CR[0]:.10f},{CR[1]:.10f})  pocket-acnode ({ACp[0]:.10f},{ACp[1]:.10f})")

# ---- LC1 fiber counts ----
def preimages(d, A, B, C):
    s, r = B*C, A*C**2
    out = []
    for wv in np.roots(h_coeffs(d, s, r)):
        g = s - p_d_np(wv, d)
        if abs(g) > 1e-12:
            out.append((wv, g, C/g))
    return out
print("\n== LC1 fiber counts (bounded |x| <= 1e4) - predict 12/10/9/8/8 ==")
samples = {"generic (1,1,1)": (1,1,1), "fold t=1/2": (r_f, s_f, 1), "cusp (small)": (CU[0][1], CU[0][0], 1),
           "cusp (neg)": (CU[1][1], CU[1][0], 1), "crunode": (CR[1], CR[0], 1), "pocket acnode": (ACp[1], ACp[0], 1)}
fc = {}
for name, (A,B,C) in samples.items():
    ps = preimages(11, A, B, C)
    fin = [p_ for p_ in ps if abs(p_[2]) <= 1e4]
    nreal = sum(1 for p_ in fin if abs(p_[0].imag) < 1e-7)
    fc[name] = (len(fin), nreal)
    print(f"  {name:16s}: bounded {len(fin):2d} (of {len(ps):2d})   real {nreal}")
ok1 = fc["generic (1,1,1)"][0]==12 and fc["fold t=1/2"][0]==10 and fc["cusp (small)"][0]==9 \
      and fc["cusp (neg)"][0]==9 and fc["crunode"][0]==8 and fc["pocket acnode"][0]==8
print("LC1 GREEN:", ok1)

# ---- LC3 escape rates ----
print("\n== LC3 escape rates ==  [predict 1/2, 2/3]")
def escape_scan(d, bs, br):
    gam, delt = [], []
    for k in np.arange(4, 8.001, 0.5):
        dd = 10.0**(-k)
        s, r = bs + dd, br + 0.7*dd
        g = [abs(s - p_d_np(wv, d)) for wv in np.roots(h_coeffs(d, s, r))]
        gam.append(min(g)); delt.append(dd*np.hypot(1.0, 0.7))
    return np.array(delt), np.array(gam)
es = {}
for name, (bs, br) in {"fold": (s_f, r_f), "cusp": CU[0]}.items():
    delt, gam = escape_scan(11, bs, br)
    fit = np.polyfit(np.log(delt), np.log(gam), 1)
    es[name] = float(fit[0])
    print(f"  {name}: slope {fit[0]:.4f}")
ok3 = abs(es["fold"]-0.5) <= 0.03 and abs(es["cusp"]-2/3) <= 0.03
print("LC3 GREEN:", ok3)

# ---- LC4 det + hinge ----
print("\n== LC4 det JF_12 at 30 exact rational points + hinge ==")
x, y, z = sp.symbols("x y z")
q_ = sp.expand(sp.integrate(w_*sp.diff(p11, w_), w_))
kap = sp.diff(p11, w_).subs(w_, 1); a_ = sp.Rational(-(1+kap)/(2+kap))
u_ = 1 + x*y; g_ = 1 + a_*x*y + x**2*z; ws_ = u_*g_
al_ = u_ + q_.subs(w_, ws_)/g_**2; be_ = 1 + p11.subs(w_, ws_)/g_
f1 = sp.expand(sp.cancel(al_/x**2)); f2 = sp.expand(sp.cancel(be_/x)); f3 = sp.expand(x*g_)
print(f"  F_12 term counts: {[len(sp.Poly(f, x, y, z).terms()) for f in (f1,f2,f3)]}"
      f"  degrees: {[sp.Poly(f, x, y, z).total_degree() for f in (f1,f2,f3)]}")
print(f"  hinge p1 = [w]p11 = {sp.Poly(p11,w_).coeff_monomial(w_)}  ==  2*q2 = 2*[w^2]Phi12 ="
      f" {2*sp.Poly(Phi12,w_).coeff_monomial(w_**2)}:",
      sp.expand(sp.Poly(p11,w_).coeff_monomial(w_) - 2*sp.Poly(Phi12,w_).coeff_monomial(w_**2)) == 0)
import random
random.seed(7)
pts = [(sp.Rational(1,1),1,1), (2,1,3), (sp.Rational(1,3),2,1), (3,-1,2), (sp.Rational(-2,5),1,4)]
for _ in range(25):
    pts.append((sp.Rational(random.randint(-9,9), random.randint(1,9)),
                sp.Rational(random.randint(-9,9), random.randint(1,9)),
                sp.Rational(random.randint(-9,9), random.randint(1,9))))
J = [[sp.diff(f, v) for v in (x, y, z)] for f in (f1, f2, f3)]
mx = 0
for px_, py_, pz_ in pts:
    if sp.expand(g_.subs({x:px_,y:py_,z:pz_})) == 0: continue
    M = [[ji.subs({x:px_,y:py_,z:pz_}) for ji in row] for row in J]
    dv = sp.det(sp.Matrix(M))
    mx = max(mx, abs(int(sp.expand(dv - 1))))
print(f"  max |det-1| over {len(pts)} exact rational points: {mx}  (0 => det==1 certificate)")
ok4 = (mx == 0)
print("LC4 GREEN:", ok4)

# ---- LC6 off-wall boundedness ----
rng = np.random.default_rng(13)
mxb = 0.0; esc = 0
for _i in range(20000):
    A = rng.normal(0,3) + 1j*rng.normal(0,3); B = rng.normal(0,3) + 1j*rng.normal(0,3)
    for wv, g, xv in preimages(11, A, B, 1):
        mxb = max(mxb, abs(xv))
    if not preimages(11, A, B, 1): esc += 1
print(f"\n== LC6 boundedness 20k targets: max |x| = {mxb:.6g} (empty fibers: {esc}) ==")
ok6 = np.isfinite(mxb) and esc == 0

json.dump({"fiber_counts": fc, "escape_slopes": es, "det_max_resid": int(mx),
           "bounded_max_x": float(mxb), "Cr": CR, "ACp": ACp,
           "LC1": bool(ok1), "LC3": bool(ok3), "LC4": bool(ok4), "LC6": bool(ok6)},
          open("atlas11_stageC.json","w"), indent=1, default=str)
print(f"[stage C part 1 done {time.time()-t0:.0f}s]")
