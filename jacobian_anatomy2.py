"""
Anatomy part 2: for each map compute the wall curve (multiple-root locus of the
fiber equation), cusps (triple escape), and bitangents (double-double escape),
which determine the missing-from-image target curves. Then numeric spot checks.

Map      fiber equation                     escape gamma
F        w^2 - w^3      = Pw - Q (P=BC/4,  Q=AC^2/4)   gamma = P - (2w - 3w^2)
G        (w^2 - w^4)/2  = Pw - Q (P=BC,    Q=AC^2/2)   gamma = P - (w - 2w^3)
H        3/2 w^4 - 4w^3 + 5/2 w^2 = Pw - Q (P=BC, Q=AC^2)  gamma = P - (5w-12w^2+6w^3)
"""
import sympy as sp
import numpy as np

w, P, Q, PP, QQ = sp.symbols("w P Q PP QQ")

def analyze(name, Phi, gamma_minus):
    p = sp.diff(Phi, w)
    print(f"==== {name}: Phi(w) = {sp.expand(Phi)}, gamma = P - ({p})")
    # wall curve: resultant of (PP - p(w)) and (QQ - (w*p - Phi))
    Qw = sp.expand(w * p - Phi)
    R = sp.factor(sp.resultant(PP - p, QQ - Qw, w))
    print("  wall curve R(PP,QQ) =", R)
    # cusps: Phi'' = 0 (triple root with line; 3 sheets escape)
    cusps = sp.solve(sp.diff(Phi, w, 2), w)
    for wc in cusps:
        Pc = sp.simplify(p.subs(w, wc)); Qc = sp.simplify(Qw.subs(w, wc))
        print(f"  cusp: w0 = {wc}  ->  (P*, Q*) = ({Pc}, {Qc})   [order-{sp.degree(Phi)+0} tangency]")
    # bitangents (quartic case): (w-w1)^2 (w-w2)^2 match, both roots escape
    deg = sp.degree(Phi)
    if deg == 4:
        s, t = sp.symbols("s t")
        lead = sp.Rational(sp.LC(sp.Poly(Phi, w)))
        prod = sp.expand(lead * (w**2 - s * w + t)**2)
        diff = sp.Poly(sp.expand(Phi - P * w + Q - prod), w)
        eqs = [e for m, e in diff.terms() if m[0] in (3, 2)]
        solst = sp.solve(eqs, [s, t], dict=True)
        print("  bitangent (s,t) solutions:", solst)
        for sol in solst:
            w1, w2 = sp.solve(sp.Eq(w**2 - sol[s] * w + sol[t], 0), w)
            Pm = sp.simplify(sp.expand(Phi - P * w + Q - prod).subs({**sol}).coeff(w, 1) * -1)
            Qm = sp.simplify(sp.expand(prod).subs(sol).subs(w, 0))
            esc1 = sp.simplify(gamma_minus.subs({P: Pm, w: w1}))
            esc2 = sp.simplify(gamma_minus.subs({P: Pm, w: w2}))
            print(f"    bitangent: roots {w1}, {w2}; P = {Pm}, Q = {Qm}; gamma at roots: {esc1}, {esc2}")
            if esc1 == 0 and esc2 == 0:
                print("    ==> BOTH escape: fiber EMPTY. Missing curve present.")
    print()

analyze("F", w**2 - w**3, P - (2 * w - 3 * w**2))
analyze("G", (w**2 - w**4) / 2, P - (w - 2 * w**3))
analyze("H", sp.Rational(3, 2) * w**4 - 4 * w**3 + sp.Rational(5, 2) * w**2,
        P - (5 * w - 12 * w**2 + 6 * w**3))

# --- numeric spot checks of emptiness / wall fiber counts using reconstructions ---
def census_generic(Phi_coeffs, pfun, backfun, T, Pfun, Qfun):
    A, B, C = [complex(v) for v in T]
    Pv, Qv = Pfun(A, B, C), Qfun(A, B, C)
    cfs = np.array(Phi_coeffs, dtype=complex)
    cfs[-2] -= Pv; cfs[-1] += Qv
    roots = np.roots(cfs)
    n = 0
    for wr in roots:
        g = Pv - pfun(wr)
        if abs(g) < 1e-9:
            continue
        n += 1
    return n, len(roots)

print("numeric: F at missing point (4/27, 4/3, 1):",
      census_generic([-1, 1, 0, 0], lambda wv: 2*wv-3*wv**2, None,
                     (sp.Rational(4,27), sp.Rational(4,3), 1),
                     lambda A,B,C: B*C/4, lambda A,B,C: A*C**2/4))
print("numeric: F at wall non-cusp, e.g. (0,0,1):",
      census_generic([-1, 1, 0, 0], lambda wv: 2*wv-3*wv**2, None, (0,0,1),
                     lambda A,B,C: B*C/4, lambda A,B,C: A*C**2/4))
print("numeric: G at bitangent target (-1/4, 0, 1):",
      census_generic([-1/2, 0, 1/2, 0, 0], lambda wv: wv - 2*wv**3, None,
                     (-sp.Rational(1,4), 0, 1), lambda A,B,C: B*C, lambda A,B,C: A*C**2/2))
print("numeric: G near-miss generic (1,0,1):",
      census_generic([-1/2, 0, 1/2, 0, 0], lambda wv: wv - 2*wv**3, None, (1,0,1),
                     lambda A,B,C: B*C, lambda A,B,C: A*C**2/2))
print("numeric: H at bitangent target (1/216, -2/9, 1):",
      census_generic([3/2, -4, 5/2, 0, 0], lambda wv: 5*wv - 12*wv**2 + 6*wv**3, None,
                     (sp.Rational(1,216), -sp.Rational(2,9), 1), lambda A,B,C: B*C, lambda A,B,C: A*C**2))
