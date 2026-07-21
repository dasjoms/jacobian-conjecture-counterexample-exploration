"""
Fiber census V2 -- via the exact univariate fiber equations (deterministic,
high precision), plus mod-p bijection tests.

Each map has a fiber equation  Phi(w) = P*w - c*Q  with P, Q read off the
target (A,B,C), and an exact reconstruction of each preimage from a root w.

  F (Alpoge):   Phi = w^2 - w^3,      P = B*C/4,   Q = A*C^2/4
                gamma = P - p(w), p=2w-3w^2; x = C/(2 gamma), u = w/gamma,
                y = (u-1)/x, z = (2(1 - 1.5(u-1) - gamma))/x^2
  G (explainer):Phi = (w^2 - w^4)/2,  P = B*C,     Q = A*C^2/2
                gamma = P - p(w), p = w - 2w^3; x = C/gamma, u = w/gamma,
                y = (u-1)/(3x), z = (1 - (4/3)(u-1) - gamma)/x^2
  H (homemade): Phi = (9/10)w^4 - 4w^3 + (5/2)w^2,  P = B*C,  Q = A*C^2
                gamma = P - p(w), p = 5w-12w^2+6w^3; x = C/gamma, u = w/gamma,
                y = (u-1)/x, z = (gamma - 1)/x^2
"""
import sympy as sp
import numpy as np

x, y, z, w = sp.symbols("x y z w")

uG = 1 + 3 * x * y; gG = 1 - 4 * x * y - x**2 * z
G = [sp.expand(sp.cancel(sp.expand(2 * uG + uG**2 - 3 * uG**4 * gG**2) / x**2)),
     sp.expand(sp.cancel(sp.expand(1 + uG - 2 * uG**3 * gG**2) / x)),
     sp.expand(x * gG)]

pH = 5 * w - 12 * w**2 + 6 * w**3
qH = sp.expand(sp.integrate(w * sp.diff(pH, w), w))
uH = 1 + x * y; gH = 1 + x**2 * z; ws = uH * gH
M1 = sp.cancel(sp.expand(uH * gH**2 + qH.subs(w, ws)) / gH**2)
H1 = sp.expand(sp.cancel(sum(ci * x**m[0] * y**m[1] * z**m[2]
                             for m, ci in sp.Poly(M1, x, y, z).terms()) / x**2))
betaH = sp.cancel(sp.expand(gH + pH.subs(w, ws)) / gH)
H = [H1, sp.expand(sp.cancel(betaH / x)), sp.expand(x * gH)]

MAPS = {}
def make_map(name, Fpolys, Phi_coeffs, Pfun, Qfun, gamma_fun, back):
    Ff = sp.lambdify([x, y, z], Fpolys, "numpy")
    MAPS[name] = dict(Ff=Ff, Phi=np.array(Phi_coeffs, float), P=Pfun, Q=Qfun,
                      gamma=gamma_fun, back=back)
w0 = 1 + x * y; A0 = w0**2 * z + y**2 * (4 + 3 * x * y)
Fpolys = [sp.expand(w0 * A0), sp.expand(y + 3 * x * A0), sp.expand(2 * x - x**2 * (3 * y + x * z))]
make_map("F", Fpolys, [-1, 1, 0, 0],           # w^2 - w^3  -> poly coeffs descending
         lambda A, B, C: B * C / 4, lambda A, B, C: A * C**2 / 4,
         lambda wv, P: P - (2 * wv - 3 * wv**2),
         lambda wv, gm, C: (lambda xx: (xx, (wv / gm - 1) / xx,
                                        2 * (1 - 1.5 * (wv / gm - 1) - gm) / xx**2))(C / (2 * gm)))
make_map("G", G, [-1/2, 0, 1/2, 0, 0],         # (w^2 - w^4)/2
         lambda A, B, C: B * C, lambda A, B, C: A * C**2 / 2,
         lambda wv, P: P - (wv - 2 * wv**3),
         lambda wv, gm, C: (lambda xx: (xx, (wv / gm - 1) / (3 * xx),
                                        (1 - 4 / 3 * (wv / gm - 1) - gm) / xx**2))(C / gm))
make_map("H", H, [9/10, -4, 5/2, 0, 0],        # quartic
         lambda A, B, C: B * C, lambda A, B, C: A * C**2,
         lambda wv, P: P - (5 * wv - 12 * wv**2 + 6 * wv**3),
         lambda wv, gm, C: (lambda xx: (xx, (wv / gm - 1) / xx,
                                        (gm - 1) / xx**2))(C / gm))

# ---- census through the fiber equation ----
def census(name, T):
    A, B, C = [complex(v) for v in T]
    P, Q = MAPS[name]["P"](A, B, C), MAPS[name]["Q"](A, B, C)
    cfs = MAPS[name]["Phi"].copy()
    cfs[-2] -= P      # Phi(w) - P*w + Q = 0
    cfs[-1] += Q
    wroots = np.roots(cfs)
    pts = []
    for wr in wroots:
        gm = MAPS[name]["gamma"](wr, P)
        if abs(gm) < 1e-9 or abs(C) < 1e-12:
            pts.append((wr, None)); continue
        pt = MAPS[name]["back"](wr, gm, C)
        Fv = np.array(MAPS[name]["Ff"](*pt), dtype=complex).reshape(3)
        res = np.linalg.norm(Fv - np.array([A, B, C]))
        pts.append((wr, (pt, res)))
    return pts

tests = {
    "F": [(0, 0, 1), (-sp.Rational(1, 4), 0, 1), (-sp.Rational(8, 27), 0, 1), (1, 2, 3)],
    "G": [(0, 0, 1), (-sp.Rational(1, 8), 0, 1), (1, 2, 3), (1, 0, 1)],
    "H": [(sp.Rational(441, 32), 0, 1), (1, 2, 3), (0, 0, 1)],
}
for name, ts in tests.items():
    for T in ts:
        pts = census(name, T)
        valid = [p for p in pts if p[1] is not None]
        ok = [p for p in valid if p[1][1] < 1e-8]
        nreal = sum(all(abs(im) < 1e-9 for im in np.imag(p[1][0])) for p in ok)
        print(f"{name} over {tuple(str(t) for t in T)}: {len(ok)} reconstructed preimages "
              f"({nreal} real), max residual {max((p[1][1] for p in ok), default=0):.3g}")
        for wr, dat in ok:
            pt, res = dat
            print("     w =", np.round(wr, 6), "-> pt =", tuple(np.round(c, 5) for c in pt))

# ---- discriminate: generic fiber has 3 distinct roots for F -----------------
P_, Q_ = sp.symbols("P_ Q_")
cub = w**3 - w**2 + P_ * w - Q_
disc = sp.factor(sp.discriminant(cub, w))
print("discriminant of F's fiber cubic:", disc, " (nonzero poly -> generic fiber = 3 points)")

# ---- mod-p bijection tests ---------------------------------------------------
def modp_bijection(Fpolys, p=101):
    idx = np.arange(p**3, dtype=np.int64)
    xa, ya, za = idx % p, (idx // p) % p, idx // (p * p)
    out = np.zeros((3, len(idx)), dtype=np.int64)
    for k, f in enumerate(Fpolys):
        acc = np.zeros(len(idx), dtype=np.int64)
        for mon, coef in sp.Poly(f, x, y, z).terms():
            rc = sp.Rational(coef)
            c = int(rc.p % p) * pow(int(rc.q % p), -1, p) % p
            term = np.full(len(idx), c, dtype=np.int64)
            for var, e in zip((xa, ya, za), mon):
                for _ in range(e):
                    term = (term * var) % p
            acc = (acc + term) % p
        out[k] = acc
    code = out[0] * p * p + out[1] * p + out[2]
    return bool(np.unique(code).size == p**3)

for name, Mp in (("F", Fpolys), ("G", G), ("H", H)):
    det = sp.factor(sp.Matrix(Mp).jacobian([x, y, z]).det())
    print(f"{name}: det = {det}; bijective on F_101^3? -> {modp_bijection(Mp, 101)}")
Fx = [x**2 + y**3, y**2, z**2]
print("contrast map (det =", sp.factor(sp.Matrix(Fx).jacobian([x, y, z]).det()), "):",
      modp_bijection(Fx, 101), "(expect False)")
