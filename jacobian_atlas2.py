"""
Stage A part 2: finish bitangent nodes (backsolve lex GB), print Disc_4 cleanly,
and singular-point classification of the wall curve.
"""
import sympy as sp
from sympy import symbols, diff, integrate, expand, factor, resultant, Rational as R

w, s, r = symbols("w s r")
w1, w2 = symbols("w1 w2")

p4 = -w**4 + w**3 - R(27,10)*w**2 + R(17,10)*w
Phi4 = expand(integrate(p4, w))
h4 = Phi4 - s*w + r
D4 = sp.expand(resultant(h4, diff(h4, w), w))
D4 = sp.expand(D4 / sp.gcd_list(list(sp.Poly(D4, s, r).coeffs())))
print("Disc_4(s,r): degree", sp.Poly(D4, s, r).total_degree(), ",", len(sp.Poly(D4,s,r).terms()), "terms")
print(D4)

# wall curve singularities: D4 = 0, dD4/ds = 0, dD4/dr = 0
Ds, Dr = diff(D4, s), diff(D4, r)
# eliminate r: resultant(Ds?? ) - use resultant of Ds, Dr over r? The singular locus is 0-dim:
print("\neliminating r from {D4=0, D4_r=0} ...", flush=True)
Sing_r = sp.factor(resultant(D4, Dr, r))
print("res(D, D_r ; r) =", Sing_r if len(str(Sing_r)) < 400 else f"<deg {sp.degree(sp.Poly(Sing_r,s)) if S else '?'}>")
Sr = sp.Poly(Sing_r, s)
print("degree in s:", Sr.degree())
s_vals = sp.nroots(Sr, n=25, maxsteps=300)
print("roots s:", [complex(v) for v in s_vals])
nodes_cusps = []
seen_sr = set()
for sv in s_vals:
    # find r with D=0 and D_r=0 (double root condition): gcd of D(sv, r) and Dr(sv, r)
    g = sp.gcd(sp.Poly(D4.subs(s, sv), r), sp.Poly(Dr.subs(s, sv), r))
    if g.degree() >= 1:
        for rv in sp.nroots(sp.Poly(g.as_expr(), r), n=25, maxsteps=400):
            dd = float(abs(Ds.subs({s: sv, r: rv})))
            key = (round(complex(sv).real,4), round(complex(sv).imag,4), round(complex(rv).real,4), round(complex(rv).imag,4))
            if key in seen_sr: continue
            seen_sr.add(key)
            nodes_cusps.append((complex(sv), complex(rv), dd))
for sv, rv, dd in nodes_cusps:
    tru = "SINGULAR" if dd < 1e-3 else "smooth (vertical tangent)"
    print(f"   ({sv:.8f},{rv:.8f})  |D_s| = {dd:.4g}  {tru}")

# --- bitangents via eliminant ---
eq1 = sp.expand(sp.cancel((p4.subs(w, w2) - p4.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(sp.cancel((p4.subs(w, w2) - p4.subs(w, w1))/(w2 - w1) - p4.subs(w, w1)))
elim_w2 = sp.factor(resultant(eq1, eq2, w1))
print("\neliminant in w2:", sp.factor(elim_w2))
# strip known diagonal factor (w2 - w1 limit) is absent; factors carry multiplicities
roots_all = []
for fac, mul in sp.factor_list(sp.factor(elim_w2))[1]:
    if w2 in fac.free_symbols:
        rr = sp.nroots(sp.Poly(fac, w2), n=30, maxsteps=500)
        roots_all += [(mul, complex(v)) for v in rr]
        print(f"  factor {fac} (x{mul}): roots {[complex(v) for v in rr]}")
w2_roots = [v for _, v in roots_all]
pairs = []
import itertools
for wv2 in w2_roots:
    # find w1 roots of eq1,eq2 numerically via nsolve random starts
    found = []
    for g in (eq1, eq2):
        p1sub = sp.lambdify(w1, g.subs(w2, wv2), "mpmath")
        # solve univariate numerically via nroots of the poly in w1
        rr = sp.nroots(sp.Poly(sp.expand(g.subs(w2, wv2)), w1), n=25, maxsteps=300)
        found.append(set((round(complex(rv).real,6), round(complex(rv).imag,6)) for rv in rr))
    common = found[0] & found[1]
    for wtup in common:
        pairs.append((complex(wtup[0], wtup[1]), complex(wv2)))

def ckey(cc, nd=5):
    return (round(cc.real, nd), round(cc.imag, nd))
seen = set(); final = []
for a, b in pairs:
    key = tuple(sorted((ckey(a), ckey(b))))
    if key not in seen:
        seen.add(key); final.append((a,b))
print(f"\nbitangent/cusp solution pairs: {len(final)} (unordered)")
for a, b in final:
    sv = p4.subs(w, a); rv = a*sv - Phi4.subs(w, a)
    sv2 = p4.subs(w, b)
    real_sr = abs(complex(sv).imag) < 1e-8 and abs(complex(rv).imag) < 1e-8
    real_contacts = abs(complex(a).imag) < 1e-8 and abs(complex(b).imag) < 1e-8
    kind = "cusp-diag" if abs(a-b) < 1e-6 else ("REAL node" if real_sr and real_contacts else ("real-line node (complex contacts)" if real_sr else "complex"))
    print(f"  w1={a:.6f} w2={b:.6f} -> (s,r)=({complex(sv):.6f},{complex(rv):.6f}) [{kind}]  check s equal: {abs(complex(sv-sv2))<1e-6}")
