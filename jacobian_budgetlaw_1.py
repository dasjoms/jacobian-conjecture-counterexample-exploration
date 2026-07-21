"""
Note 21, stage 1: THE BUDGET LAW - the universal-algebra certificates.
=====================================================================
THEOREM (local + algebraic, all d):
  E_d(t2) = Res_{t1}(eq1, eq2) = p'_d(t2) * Res_{t1}(eq1, te2),  where
  eq2 = (t2 - t1) * te2  EXACTLY (the diagonal is a component of {eq2 = 0}),
  so  E_d = p'_d^2 * Cof_d  with  Cof_d = Res_{t1}(eq1, te2) / p'_d in Q[t2].
Locks (registered BEFORE computation):
  L-A1: (w2 - w1) | eq2 exactly, for d = 3..14 (sympy div, remainder 0)
  L-A2: eq1(t, t) = p'_d(t) exactly, d = 3..14
  L-A3: te2(t, t) = p'_d(t)/2 exactly, d = 3..14
  L-B : Res_{t1}(eq1, t2 - t1) = eq1(t2, t2) = p'_d(t2) exactly, d = 3..12
        (sign EXACTLY +1, not -1)
  L-C : Res_{t1}(eq1, eq2) == p'_d(t2) * Res_{t1}(eq1, te2) exactly, d = 3..10
  L-D (germ certificate, universal - not tower-specific):
        with Phi(tau + T) = c0 + c1 T + c3 T^3 + c4 T^4 + ... (c2 = p'/2 = 0):
        eq1 = 3 c3 (u + x) + O(2)        -> tangent slope -1
        eq2 = (u - x) * te2 EXACTLY (through the germ's order)
        te2 = c3 (u + 2 x) + O(2)        -> tangent slope -2,  distinct from -1, +1
        all as exact symbolic identities in the c_k (order-6 germ)
  L-E : p'_d squarefree over Q (gcd(p', p'') = 1) for d = 3..30 (dance machinery)
=====================================================================
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, cancel, resultant, Rational as R, div, gcd

t0 = time.time()
w, w1, w2, u = symbols("w w1 w2 u")
out = {"LA1": {}, "LA2": {}, "LA3": {}, "LB": {}, "LC": {}, "LE": []}

def tower(d):
    m = d*(d+1)
    p = 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - R(6, m))
    return sp.expand(p), sp.expand(integrate(p, w))

eqs = {}
for d in range(3, 15):
    p, Ph = tower(d)
    eq1 = sp.expand(cancel((p.subs(w, w2) - p.subs(w, w1))/(w2 - w1)))
    eq2 = sp.expand(cancel((Ph.subs(w, w2) - Ph.subs(w, w1))/(w2 - w1) - p.subs(w, w1)))
    quo, rem = sp.div(sp.Poly(eq2, w1, w2), sp.Poly(w2 - w1, w1, w2))
    assert expand(rem.as_expr()) == 0
    te2 = expand(quo.as_expr())
    eqs[d] = (p, Ph, eq1, eq2, te2)
    out["LA1"][d] = True
    out["LA2"][d] = (sp.expand(eq1.subs(w2, w1) - diff(p, w).subs(w, w1)) == 0)
    out["LA3"][d] = (sp.expand(te2.subs(w2, w1) - diff(p, w).subs(w, w1)/2) == 0)
    if d == 3 or d == 14 or not (out["LA2"][d] and out["LA3"][d]):
        print(f"  d={d}: LA1=True LA2={out['LA2'][d]} LA3={out['LA3'][d]}", flush=True)
print(f"L-A all true d=3..14: {all(out[k][d] for k in ('LA1','LA2','LA3') for d in range(3,15))}  [{time.time()-t0:.0f}s]", flush=True)

for d in range(3, 13):
    p, Ph, eq1, eq2, te2 = eqs[d]
    RB = resultant(eq1, w2 - w1, w1)
    pd2 = sp.expand(diff(p, w).subs(w, w2))
    exact_plus = sp.expand(RB - pd2) == 0
    print(f"  L-B d={d}: Res(eq1, t2-t1) == + p'_d(t2): {exact_plus}", flush=True)
    out["LB"][d] = exact_plus
    if d <= 10:
        Efull = resultant(eq1, eq2, w1)
        Eprod = sp.expand(pd2 * resultant(eq1, te2, w1))
        ok = sp.expand(Efull - Eprod) == 0
        out["LC"][d] = ok
        print(f"  L-C d={d}: E == p' * Res(eq1, te2) EXACT: {ok}   [deg E {sp.degree(Efull, w2)} = d(d-1) {d*(d-1)}: {sp.degree(Efull, w2)==d*(d-1)}]", flush=True)

# --------- L-D: the universal germ certificate ----------
print("\nL-D: (1/2)  N := Phi(t2)-Phi(t1)-(t2-t1)p(t1) divisible by (t2-t1)^2 EXACTLY, d=3..14:")
div2 = {}
for d in range(3, 15):
    p, Ph, eq1, eq2, te2 = eqs[d]
    N = sp.expand(Ph.subs(w,w2) - Ph.subs(w,w1) - (w2-w1)*p.subs(w,w1))
    quo2, rem2 = sp.div(sp.Poly(N, w1, w2), sp.Poly((w2-w1)**2, w1, w2))
    ok = (expand(rem2.as_expr()) == 0) and (expand(quo2.as_expr() - te2) == 0)
    div2[d] = ok
print(f"      N/(t2-t1)^2 == te2 exactly: {all(div2.values())}  (te2 = second divided difference)")
out["LA1b"] = all(div2.values())

c0, c1, c3, c4, c5, c6 = symbols("c0 c1 c3 c4 c5 c6")
x = symbols("x")
Tsym = symbols("T")
germ = c0 + c1*Tsym + c3*Tsym**3 + c4*Tsym**4 + c5*Tsym**5 + c6*Tsym**6   # c2 = p'/2 = 0 (flex)
pgerm = diff(germ, Tsym)
eq1g = sp.expand(cancel((pgerm.subs(Tsym, u) - pgerm.subs(Tsym, x))/(u - x)))
rem1 = sp.Poly(sp.expand(eq1g - 3*c3*(u + x)), u, x)
mindeg1 = min(a + b for (a, b), c in rem1.terms() if c != 0)
print(f"      eq1-cone: 3*c3*(u+x), residual min degree {mindeg1} (>=2) -> tangent slope -1")

Ng = sp.expand(germ.subs(Tsym, u) - germ.subs(Tsym, x) - (u - x)*pgerm.subs(Tsym, x))
_q1, _r1 = sp.div(sp.Poly(Ng, u, x), sp.Poly(u - x, u, x))
assert sp.expand(_r1.as_expr()) == 0          # (u-x) | N
_q2, _r2 = sp.div(_q1, sp.Poly(u - x, u, x))
assert sp.expand(_r2.as_expr()) == 0          # (u-x)^2 | N
te2g = sp.expand(_q2.as_expr())
expect_N = c3*(u-x)**2*(u+2*x)
remN = sp.Poly(sp.expand(Ng - expect_N), u, x)
mindegN = min(a + b for (a, b), c in remN.terms() if c != 0)
rem2c = sp.Poly(sp.expand(te2g - c3*(u + 2*x)), u, x)
mindeg2 = min(a + b for (a, b), c in rem2c.terms() if c != 0)
print(f"      N_germ = c3*(u-x)^2*(u+2x) + O(4): residual min degree {mindegN} (>=4 predict)")
print(f"      te2-cone: c3*(u+2x), residual min degree {mindeg2} (>=2) -> tangent slope -2")
dg = sp.expand(te2g.subs(u, x))
print(f"      te2(x,x) = {sp.factor(dg)} ; p'(x)/2 = {sp.expand(sp.diff(pgerm,Tsym).subs(Tsym,x)/2)} ; equal:",
      sp.expand(dg - sp.diff(pgerm, Tsym).subs(Tsym, x)/2) == 0)
out["LD"] = {"cones": ["3*c3*(u+x)", "(u-x) [diag]", "c3*(u+2x)"],
             "mindeg1": mindeg1, "mindegN": mindegN, "mindeg2": mindeg2,
             "diag_id": sp.expand(dg - sp.diff(pgerm, Tsym).subs(Tsym, x)/2) == 0}

# --------- L-E: p' squarefree (dance) ----------
bad = []
for d in range(3, 31):
    p, Ph = tower(d)
    g = sp.gcd(diff(p, w), diff(p, w, 2))
    if sp.degree(g, w) > 0: bad.append((d, str(g)))
    out["LE"].append((d, sp.degree(g, w) == 0))
print(f"\nL-E p'_d squarefree d=3..30: {all(v for _, v in out['LE'])}  bad: {bad}", flush=True)

green = (all(out[k][d] for k in ('LA1','LA2','LA3') for d in range(3,15)) and out["LA1b"]
         and all(out['LB'].values()) and all(out['LC'].values())
         and mindeg1 >= 2 and mindegN >= 4 and mindeg2 >= 2 and not bad)
json.dump(out, open("budgetlaw_stage1.json", "w"), indent=1, default=str)
print(f"[stage 1 done {time.time()-t0:.0f}s]  UNIVERSAL GREEN: {green}", flush=True)
