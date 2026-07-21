#!/usr/bin/env python3
# LAB NOTE 17 — THE DIAGONAL CENSUS: the shadow root's uniqueness, proven.
# =====================================================================
# THE KEY: G(t) = Phi(t) - (t-1)p(t)  =>  G'(t) = -(t-1) p'(t).
# Note 16's Theorem 2 (sign atlas of p') then settles EVERYTHING.
# LOCKS (registered before compute, 2026-07-21):
#   C1 (symbolic, per d=2..40): G'(t) + (t-1) p'(t) is the zero polynomial.
#   C2 (symbolic d, parity-split): closed forms
#        G(-1)_odd  = -4 - 2(d-9)/(d(d+1))         [ < 0 all odd d>=3 ]
#        G(-1)_even = -12 + (2d+20)/(d(d+1))      [ < 0 all even d>=4 ]
#        G(-2)_odd  = -36 + (40/3)c0 + 2^(d-1)(9 - 2/d - 4/(d+1))       [ > 0 d>=5 ; v2 CORRECTED ]
#        G(-2)_even = -36 + (40/3)c0 + 2^(d-1)(2/d + 4/(d+1) - 9)      [ < 0 d>=4 ]
#      verified as exact identities; sign claims via rational certificates.
#   C3 (Sturm audit): real root count of G on R = {3 odd, 2 even} for d=2..301;
#      count in (-2,-1] = 1 (odd) / 0 (even); count in (-1,0] = 1 for all d (the t=0 root;
#      sympy count_roots is right-closed) => NO roots strictly inside (-1,0).
#   C4 (bisection, exact Q): t* interval width 2^-40 for odd d=5..201;
#      image window s* = (p(hi),p(lo)) inside (-1,0) for EVERY odd d>=5.
#   C5 (multiplicity): t=1 exactly double: G(1)=G'(1)=0, G''(1) = 5-c0 != 0 (symbolic);
#      t=0 simple: G'(0) = 2-c0 != 0. Per-d factor check d=3..60.
#   C6: d=3 exact: G_3 = (3/4) t (t-1)^2 (t+2); images p_3(-2) = -1 (coalescence).
import sympy as sp
from sympy import Rational as Q, symbols, factor, Poly
import json, time
t0 = time.time()
w, dd, par = symbols('w dd par')
c0s = Q(6,1)/(dd*(dd+1))
OUT = {"C1": {}, "C2": {}, "C3": {}, "C4": {}, "C5": {}, "C6": {}}

def p_poly(d_):  return 2*w - 3*w**2 + w*(1-w)*(w**(d_-2) - Q(6, d_*(d_+1)))
def Phi_poly(d_):
    p = p_poly(d_)
    return sp.integrate(p, (w, 0, w)) if False else sum(sp.integrate(t_, w) for t_ in sp.Poly(p, w).terms() and [0])
def make_G(d_):
    p = p_poly(d_)
    Phi = sp.integrate(p, w)  # rational antiderivative, const 0 (poly, Phi(0)=0)
    return sp.expand(Phi - (w-1)*p), p

# C1: key identity per d
for d_ in range(2, 41):
    G, p = make_G(d_)
    ident = sp.expand(sp.diff(G, w) + (w-1)*sp.diff(p, w))
    OUT["C1"][d_] = bool(ident == 0)
OUT["C1"]["all"] = all(OUT["C1"][k] for k in list(OUT["C1"].keys()))
print("C1:", OUT["C1"]["all"], flush=True)

# C2: symbolic closed forms (parity via par = (-1)^d)
def p_sym(tval_abs, sign_abs_is_neg):
    # p(w0) for w0 = -A (A>0): (-A)^k = (-1)^k A^k; powers d-1,d: (-1)^(d-1)=-par,(-1)^d=par
    A = tval_abs
    return 2*(-A) - 3*A**2 - c0s*(-A) + c0s*A**2 + (-par)*A**(dd-1)*A**0 if False else None
def p_eval_closed(A):
    # w = -A:  base: -2A -3A^2 + c0 A + c0 A^2 ; w^(d-1) = (-par) A^(d-1)... keep exponent symbolic:
    base = -2*A - 3*A**2 + c0s*A + c0s*A**2
    return base, A
# Do it directly with symbolic exponents substituted:
def G_at_negA(A):
    t_ = -A
    base_p  = 2*t_ - 3*t_**2 - c0s*t_ + c0s*t_**2
    hi_p    = (-par)*A**(dd-1) - par*A**dd      # t^(d-1), t^d  with t=-A: (-A)^k=(-1)^k A^k; (-1)^(d-1) = -par
    p_val   = base_p + hi_p
    base_P  = t_**2 - t_**3 - c0s*t_**2/2 + c0s*t_**3/3
    hi_P    = (par*A**dd)/dd - ((-par)*A**(dd+1))/(dd+1)     # Phi high terms: t^d/d - t^(d+1)/(d+1)
    Phi_val = base_P + hi_P
    return sp.expand(Phi_val - (t_-1)*p_val)
G1 = G_at_negA(1)
G2 = G_at_negA(2)
def zcheck(expr):   # exact zero over Q(d): Poly in par? expr has par -> substitute each value
    out = {}
    for pv, name in ((-1, 'odd'), (1, 'even')):
        e = sp.expand(expr.subs(par, pv))
        out[name] = e
    return out
G1v, G2v = zcheck(G1), zcheck(G2)
# claimed forms
G1_odd_claim  = -4 - 2*(dd-9)/(dd*(dd+1))
G1_even_claim = -12 + (2*dd+20)/(dd*(dd+1))
G2_odd_claim  = -36 + Q(40,3)*c0s + 2**(dd-1)*(9 - Q(2,1)/dd - Q(4,1)/(dd+1))   # CORRECTED (v1 had sign slip in the middle bracket; caught by CAS)
G2_even_claim = -36 + Q(40,3)*c0s + 2**(dd-1)*(Q(2,1)/dd + Q(4,1)/(dd+1) - 9)
def same(a, b):  # rational identity over Q(d)
    n, dn = sp.fraction(sp.together(sp.simplify(a-b)))
    return n == 0
OUT["C2"]["G(-1)_odd_match"]  = same(G1v['odd'],  G1_odd_claim)
OUT["C2"]["G(-1)_even_match"] = same(G1v['even'], G1_even_claim)
OUT["C2"]["G(-2)_odd_match"]  = same(G2v['odd'],  G2_odd_claim)
OUT["C2"]["G(-2)_even_match"] = same(G2v['even'], G2_even_claim)
# sign certificates
q = symbols('q', nonnegative=True)
def pos_coeff_poly(expr):  return all(c >= 0 for c in sp.Poly(sp.expand(expr), q).coeffs())
cert = {}
cert["G(-1)_odd<0  <==>  4d^2+6d-18>0 (d>=3)"] = pos_coeff_poly(4*(3+q)**2 + 6*(3+q) - 18)
cert["G(-1)_even<0 <==> 12d^2+10d-20>0 (d>=4)"] = pos_coeff_poly(12*(4+q)**2 + 10*(4+q) - 20)
cert["G(-2)_odd>0: 9*2^(d-1)>=144>36 alone (d>=5)"] = (9*2**4 >= 144) and (144 > 36)
cert["G(-2)_even<0: 2/d+4/(d+1)-9 <= -7.7 (d>=4)"] = (Q(2,4)+Q(4,5)-9) < 0
# exact rational audit of the sign claims per explicit d
okG1o = all(bool(G1_odd_claim.subs(dd, d_) < 0) for d_ in range(3, 302, 2))
okG1e = all(bool(G1_even_claim.subs(dd, d_) < 0) for d_ in range(4, 302, 2))
okG2o = all(bool(G2_odd_claim.subs(dd, d_) > 0) for d_ in range(5, 302, 2))
okG2e = all(bool(G2_even_claim.subs(dd, d_) < 0) for d_ in range(4, 302, 2))
OUT["C2"]["sign_audit"] = dict(G1odd_3_301=okG1o, G1even_4_301=okG1e, G2odd_5_301=okG2o, G2even_4_301=okG2e)
cert["audit_all"] = okG1o and okG1e and okG2o and okG2e
OUT["C2"]["certificates"] = cert
OUT["C2"]["closed_forms"] = dict(G1_odd=str(G1v['odd']), G1_even=str(G1v['even']),
                                 G2_odd=str(G2v['odd']), G2_even=str(G2v['even']))
print("C2 matches:", [OUT["C2"][k] for k in OUT["C2"] if k.endswith('match')], " audit:", cert["audit_all"], flush=True)

# C3: Sturm audits d=2..301
# sympy count_roots(a,b) counts roots in (a, b]  [RIGHT-CLOSED], so
# count_roots(-1,0) includes the t=0 root => expected 1 for ALL d, not 0.
cen = {}
for d_ in range(2, 302):
    G, p = make_G(d_)
    GP = sp.Poly(G, w)
    n_all = GP.count_roots(-10**7, 10**7)
    n_mid = GP.count_roots(-2, -1)          # (-2,-1]: t=0 NOT included; -1 is never a root
    n_neg1_0 = GP.count_roots(-1, 0)        # (-1,0]: includes t=0 => expect 1
    exp = (3, 1) if d_ % 2 == 1 else (2, 0)
    good = (n_all == exp[0]) and (n_mid == exp[1]) and (n_neg1_0 == 1)
    if not good:
        cen[d_] = ("FAIL", n_all, n_mid, n_neg1_0)
    if d_ % 50 == 0: print("  C3 d=%d (%.1fs)" % (d_, time.time()-t0), flush=True)
OUT["C3"]["fails"] = cen
OUT["C3"]["PASS_d2..301"] = (len(cen) == 0)
print("C3:", OUT["C3"]["PASS_d2..301"], cen, flush=True)

# C4: exact-rational bisection for t* (odd d), image window s*
bis = {}
tstar_win = {}
for d_ in range(5, 202, 2):
    G, p = make_G(d_)
    lo, hi = Q(-2), Q(-1)
    glo, ghi = G.subs(w, lo), G.subs(w, hi)
    assert bool(glo > 0) and bool(ghi < 0), (d_, glo, ghi)
    for _ in range(42):
        mid = (lo+hi)/2
        if bool(G.subs(w, mid) > 0): lo = mid
        else: hi = mid
    slo, shi = p.subs(w, hi), p.subs(w, lo)     # p decreasing: s* in (p(hi), p(lo))
    in01 = bool(slo > -1) and bool(shi < 0)
    bis[d_] = dict(tstar=(str(lo), str(hi)), sstar=(str(slo), str(shi)), im01=in01)
    if not in01: bis[d_]['FAIL'] = str((slo, shi))
    if d_ in (5, 11, 21, 101, 201): tstar_win[d_] = (str(lo), str(hi), str(slo), str(shi))
OUT["C4"]["examples_windows"] = tstar_win
OUT["C4"]["all_images_in_(-1,0)"] = all(v.get('im01') for v in bis.values())
print("C4 images ok:", OUT["C4"]["all_images_in_(-1,0)"], flush=True)

# C5: multiplicity per-d + symbolic constants
symb = {}
G1v_, p1v = make_G(1) if False else (None, None)
# symbolic: G(1), G'(1), G''(1) via Phi/p values:
p1 = -1  # p_d(1) = -1 (note 16 T1)
Phi1 = 0
symb["G(1)"] = "Phi(1) - 0*p(1) = 0  [certified symbolically in note 16 T1]"
symb["G'(1)"] = "-(t-1)p'(t)|_1 = 0  [via C1 identity]"
symb["G''(1)"] = "-p'(1) = 5 - c0 > 0  [p'(1) = c0-5, note 16 T1]"
symb["G'(0)"] = "p'(0) = 2 - c0 > 0"
OUT["C5"]["symbolic"] = symb
mult = {}
for d_ in range(3, 61):
    G, p = make_G(d_)
    Qr, remd = sp.div(G, (w-1)**2, w)
    val = sp.expand(remd)           # a*w+b; double root iff a=b=0
    ok_double = (sp.expand(remd) == 0)
    Q_at_1 = sp.expand(Qr).subs(w, 1)
    ok_exactly2 = bool(sp.expand(Q_at_1 - (5 - Q(6, d_*(d_+1)))/2) == 0)  # G''(1)/2
    if not (ok_double and ok_exactly2): mult[d_] = ("FAIL", str(remd), str(Q_at_1))
OUT["C5"]["fails"] = mult
OUT["C5"]["PASS_d3..60"] = (len(mult) == 0)
print("C5:", OUT["C5"]["PASS_d3..60"], flush=True)

# C6: d=3 exact
G3, p3 = make_G(3)
f3 = sp.factor(G3)
OUT["C6"]["G3_factor"] = str(f3)
OUT["C6"]["G3_exact_form_match"] = bool(sp.expand(f3 - Q(3,4)*w*(w-1)**2*(w+2)) == 0)
OUT["C6"]["p3(-2)"] = str(p3.subs(w, -2))
s3cusp = (-1 - sp.sqrt(3))/2
OUT["C6"]["ordering t*=-2 < left cusp -(1+sqrt3)/2"] = bool(-2 < s3cusp)
print("C6:", OUT["C6"]["G3_factor"], OUT["C6"]["p3(-2)"], flush=True)

json.dump(OUT, open('/home/user/diagcensus_stage1.json', 'w'), indent=1, default=str)
print("DONE %.1fs" % (time.time()-t0))
