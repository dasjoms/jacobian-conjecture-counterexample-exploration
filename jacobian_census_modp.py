"""
Numeric fiber census (batched complex Newton) for the three Keller maps,
plus finite-field bijection tests (JC is TRUE over finite fields -> each map
must be a permutation of F_p^3; a great falsification test of our algebra).
"""
import sympy as sp
import numpy as np

x, y, z, w = sp.symbols("x y z w")

# ---- rebuild the three maps -------------------------------------------------
w0 = 1 + x * y
A0 = w0**2 * z + y**2 * (4 + 3 * x * y)
F = [sp.expand(w0 * A0), sp.expand(y + 3 * x * A0), sp.expand(2 * x - x**2 * (3 * y + x * z))]

uG = 1 + 3 * x * y; gG = 1 - 4 * x * y - x**2 * z
N1 = sp.expand(2 * uG + uG**2 - 3 * uG**4 * gG**2)
N2 = sp.expand(1 + uG - 2 * uG**3 * gG**2)
G = [sp.expand(sp.cancel(N1 / x**2)), sp.expand(sp.cancel(N2 / x)), sp.expand(x * gG)]

p = 5 * w - 12 * w**2 + 6 * w**3
q = sp.expand(sp.integrate(w * sp.diff(p, w), w))
uH = 1 + x * y; gH = 1 + x**2 * z; ws = uH * gH
M1 = sp.expand(uH * gH**2 + q.subs(w, ws))
alpha = sp.cancel(M1 / gH**2)
H1 = sp.expand(sum(ci * x**m[0] * y**m[1] * z**m[2]
                   for m, ci in sp.Poly(alpha, x, y, z).terms()))
H1 = sp.expand(sp.cancel(H1 / x**2))
beta = sp.cancel(sp.expand(gH + p.subs(w, ws)) / gH)
H2 = sp.expand(sp.cancel(beta / x))
H = [H1, H2, sp.expand(x * gH)]

MAPS = {"F (Alpoge's)": F, "G (explainer's)": G, "H (homemade)": H}
for name, M in MAPS.items():
    det = sp.factor(sp.Matrix(M).jacobian([x, y, z]).det())
    print(f"{name}: det = {det}")

# ---- batched complex Newton census ------------------------------------------
def census(Fpolys, target, n=5000, iters=100, seed=1):
    J = sp.Matrix(Fpolys).jacobian([x, y, z])
    Ff = sp.lambdify([x, y, z], Fpolys, "numpy")
    Jf = sp.lambdify([x, y, z], J, "numpy")
    rng = np.random.default_rng(seed)
    rad = 10.0 ** rng.uniform(-2, 2.5, n)
    ang = rng.uniform(0, 2 * np.pi, (n, 3))
    X = rad[:, None] * np.exp(1j * ang) * np.abs(rng.standard_normal((n, 3)))[:, :]
    X = X.astype(complex)
    T = np.asarray(target, dtype=complex)
    for _ in range(iters):
        xa, ya, za = X[:, 0], X[:, 1], X[:, 2]
        Fv = np.array(Ff(xa, ya, za), dtype=complex).T - T
        Jm = np.array(Jf(xa, ya, za), dtype=complex).transpose(2, 0, 1)
        try:
            d = np.linalg.solve(Jm, Fv[..., None])[..., 0]
        except np.linalg.LinAlgError:
            d = np.zeros_like(X)  # det is a nonzero constant; singular solves shouldn't occur
        X = X - d
    Fv = np.array(Ff(X[:, 0], X[:, 1], X[:, 2]), dtype=complex).T - T
    conv = X[np.linalg.norm(Fv, axis=1) < 1e-9]
    keys = {}
    for v in conv:
        k = (round(v[0].real, 5), round(v[0].imag, 5), round(v[1].real, 5),
             round(v[1].imag, 5), round(v[2].real, 5), round(v[2].imag, 5))
        keys[k] = v
    return list(keys.values()), len(conv)

targets = {
    "F (Alpoge's)": [(0, 0, 0), (-sp.Rational(1, 4), 0, 0), (sp.Rational(1), 2, 3),
                     (-sp.Rational(8, 27), 0, 1)],
    "G (explainer's)": [(0, 0, 1), (-sp.Rational(1, 8), 0, 1), (1, 2, 3)],
    "H (homemade)": [(sp.Rational(441, 32), 0, 1), (1, 2, 3), (0, 0, 1)],
}
for name, ts in targets.items():
    for T in ts:
        roots, nconv = census(MAPS[name], [complex(t) for t in T])
        nreal = sum(all(abs(v[i].imag) < 1e-7 for i in range(3)) for v in roots)
        print(f"{name} over {tuple(str(t) for t in T)}: {len(roots)} distinct roots "
              f"({nreal} real), from {nconv} converged starts")
        for v in sorted(roots, key=lambda v: (v[0].real, v[0].imag)):
            print("    ", tuple(np.round(v, 6)))

# ---- finite-field bijection tests -------------------------------------------
def modp_bijection(Fpolys, p=101):
    Ps = [sp.Poly(f, x, y, z) for f in Fpolys]
    idx = np.arange(p**3, dtype=np.int64)
    xa, ya, za = idx % p, (idx // p) % p, idx // (p * p)
    out = np.zeros((3, len(idx)), dtype=np.int64)
    for k, P in enumerate(Ps):
        acc = np.zeros(len(idx), dtype=np.int64)
        for mon, coef in P.terms():
            c = (sp.Rational(coef).p * pow(sp.Rational(coef).q % p, -1, p)) % p
            term = np.full(len(idx), int(c), dtype=np.int64)
            for var, e in zip((xa, ya, za), mon):
                for _ in range(e):
                    term = (term * var) % p
            acc = (acc + term) % p
        out[k] = acc
    code = out[0] * p * p + out[1] * p + out[2]
    return np.unique(code).size == p**3

for name, M in MAPS.items():
    print(f"{name}: bijective on F_101^3? -> {modp_bijection(M, 101)}")

# contrast: a map whose Jacobian determinant is NOT constant
Fx = [x**2 + y**3, y**2, z**2]
detx = sp.factor(sp.Matrix(Fx).jacobian([x, y, z]).det())
print("contrast map det =", detx, "-> bijective?", modp_bijection(Fx, 101), "(expected False)")
