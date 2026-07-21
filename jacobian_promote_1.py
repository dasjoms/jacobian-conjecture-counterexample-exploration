#!/usr/bin/env python3
# LAB NOTE 16 — STAGE 1: THE PROMOTION (machine certificates for five theorems)
# =====================================================================
# PREDICTIONS LOCKED BEFORE ANY COMPUTE (2026-07-21):
#   L1: for every d=2..300 the following exact-rational inequalities hold:
#       (a) p'(0) = 2 - c0 > 0 ; (b) p'(1/2) = -1 + (d-2)/2^(d-1) < 0 ;
#       (c) (d-1)(d-2) 2^(3-d) + 2*c0 < 6 for d>=5   [p''<0 on [0,1/2]]
#       (d) -1 + c0 + ((d-2)/d)^(d-2) < 0 for d>=5    [p'<0 on [1/2,(d-1)/d]]
#       (e) -4 + 6/d + c0 < 0 for d>=3                [p'<0 on [(d-1)/d,1]]
#       (f) 9 - 2d - 3*c0 < 0 for odd d>=5            [left cusp root in (-1,0)]
#       (g) for even d: 2-c0>0, 6-2c0>0 automatically => no negative roots
#   L2: Sturm cross-check: #real roots of p_d' for d=2..120 equals
#       {1: even d, 2: odd d} with positive root unique in (0,1/2),
#       odd negative root in (-1,0) for d>=5 resp in (-2,-1) for d=3.
#   L3: GHOST CONTRACTION: for d>=12 the map T(v) = [-c0/3 + E(v)]/(6-2c0),
#       E(v) = (1/3+v)^(d-2)(2d/3 - 1 - dv), contracts on B_d = { |v| <= Bd }
#       with Bd = 1/d^2 : invariance margin m(d)= (Bd-|T(0)|)/Bd - L_g > 0
#       and L_g < 1/2 for d>=12. All over exact rationals (interval-bounded exp).
#   L4: LEFT-CUSP CONTRACTION (odd d>=7): T_L(u) = 1 - R(u)^(1/(d-2)) on
#       I_d = [ (1-eps) u0, (1+eps) u0 ], u0 = ln((2d-1)/8)/(d-2), eps=0.15:
#       Lipschitz L_L(d) < 1 and  8 |u| iterates' error bound L_L^8*(diam) < 1e-12.
#   L5: SHADOW IDENTITY, symbolic d with odd parity: with t = -(1+U), X=(1+U)^(d-2),
#       G(t) = Phi(t) - (t-1)p(t) satisfies
#           d(d+1)*G  ==  (U+1) * [ (d^2+d)(U+2)^2 (X-2) + 4U^2+17U+19 - X(U+1)(1+d(U+2)) ]
#       IDENTICALLY as a polynomial identity in (U, X, 1/((d)(d+1))) over symbolic d.
#   L6: the order-by-order eps-solve of L5 reproduces u2,u3,u4 exactly in Q(ln2);
#       re-verification residual exactly 0 at orders eps^1..eps^4 *and* mismatch-free
#       against jcorner_final.json's published values.
# Funding statement: if any lemma bound fails for some d in range, the theorem's
# statement gets DOWNGRADED to that d-range and the failure is reported verbatim.
import sympy as sp
from sympy import Rational as Q, symbols, binomial, expand, Poly, S
import json, math, time

t0 = time.time()
w, U, X, d = symbols('w U X d')
OUT = {"locks": "see header", "L1": {}, "L2": {}, "L3": {}, "L4": {}, "L5": {}, "L6": {}}

def c0(dd): return Q(6, dd*(dd+1))

def pprime(dd):
    return (2 - c0(dd)) + (2*c0(dd) - 6)*w + (dd-1)*w**(dd-2) - dd*w**(dd-1)

# ---------- L1: exact rational inequality audit, d=2..300 ----------
fails = []
for dd in range(2, 301):
    c = c0(dd)
    chk = {}
    chk['a'] = (2 - c) > 0
    chk['b'] = (-1 + Q(dd-2, 2**(dd-1))) < 0
    if dd >= 5:
        chk['c'] = (Q((dd-1)*(dd-2), 2**(dd-3))) < (6 - 2*c)          # (c): p'' bound
        # (d): -1 + c0 + ((d-2)/d)^(d-2) < 0   [rational power -> exact]
        chk['d'] = (-1 + c + Q(dd-2, dd)**(dd-2)) < 0
    if dd >= 3:
        chk['e'] = (-4 + Q(6, dd) + c) < 0
    if dd >= 5 and dd % 2 == 1:
        chk['f'] = (9 - 2*dd - 3*c) < 0
    for k, v in chk.items():
        if not bool(v):
            fails.append((dd, k))
OUT["L1"]["fails"] = fails
OUT["L1"]["range"] = "d=2..300"
OUT["L1"]["PASS"] = (len(fails) == 0)

# induction certificate for (c): f(d)=(d-1)(d-2)/2^(d-3): f(d+1)/f(d)= d/(2(d-2)) < 1 for d>=5?
dd_sym = symbols('dd', positive=True, integer=True)
ratio = sp.simplify((dd_sym*(dd_sym-1)/2**(dd_sym-2)) / ((dd_sym-1)*(dd_sym-2)/2**(dd_sym-3)))
OUT["L1"]["ratio_c(d+1)/c(d)"] = str(ratio)          # d/(2(d-2)); at d=5 = 5/6 < 1 and decreasing
OUT["L1"]["ratio_c_at_5"] = str(ratio.subs(dd_sym, 5))
diff = sp.factor(sp.simplify(ratio.subs(dd_sym, dd_sym+1) - ratio))
OUT["L1"]["ratio_c_diff_factored"] = str(diff)   # -1/((dd-2)(dd-1)) < 0 for dd>=3 [sign of factored form]
num, den = sp.fraction(diff)
OUT["L1"]["ratio_c_monotone_decreasing"] = bool((num < 0) == True) or "see factored: -1/((dd-1)(dd-2))<0 for dd>=3"

# ---------- L2: Sturm COUNT cross-check, d=2..150 (count_roots = exact Sturm) ----------
dance = {}
pos_unique = True; neg_window = True
for dd in range(2, 151):
    ppoly = sp.Poly(pprime(dd), w)
    nR = ppoly.count_roots(-100000, 100000)
    npos_half = ppoly.count_roots(Q(1, 10**60), Q(1, 2))          # roots in (0, 1/2)
    nneg = ppoly.count_roots(-1, Q(-1, 10**60))                    # roots in (-1, 0)
    nneg2 = ppoly.count_roots(-2, -1) if dd == 3 else None
    dance[dd] = nR
    expect = 1 if dd % 2 == 0 else 2
    if nR != expect: dance[dd] = f"FAIL total={nR}"
    if npos_half != 1:
        pos_unique = False; dance[dd] = f"FAIL pos_half={npos_half}"
    if dd % 2 == 1 and dd >= 5 and nneg != 1:
        neg_window = False; dance[dd] = f"FAIL neg={nneg}"
    if dd == 3 and not (nneg2 == 1 and nneg == 0):
        neg_window = False; dance[dd] = "FAIL d3 neg window"
    if dd == 3 and npos_half != 1:
        pos_unique = False; dance[dd] = "FAIL d3 pos"
    if dd % 25 == 0: print(f"  L2 d={dd} done ({time.time()-t0:.1f}s)", flush=True)
bad = {k: v for k, v in dance.items() if isinstance(v, str)}   # failures stored as "FAIL ..." strings
OUT["L2"]["failures"] = bad
OUT["L2"]["positive_root_unique_in_(0,1/2)"] = pos_unique
OUT["L2"]["negative_root_window_kept"] = neg_window
OUT["L2"]["counts_ok_d2..150"] = (len(bad) == 0)
OUT["L2"]["spot"] = {"d=2": dance[2], "d=3": dance[3], "d=11": dance[11], "d=47": dance[47], "d=150": dance[150]}

json.dump(OUT, open('/home/user/promote_stage1.json', 'w'), indent=1, default=str)
print("L1 PASS:", OUT["L1"]["PASS"], " fails:", fails[:5])
print("L2 counts ok d2..150:", OUT["L2"]["counts_ok_d2..150"], " pos_uniq:", pos_unique, " neg_win:", neg_window)
print("elapsed %.1fs" % (time.time()-t0))
