"""
Complete image enumeration over F_101^3 for F, G, H:
compute the exact complement of the image and compare with the missing curve.

Predictions:
  F: complement = {( 4/(27 C^2),   4/(3 C),   C) : C in F_101*}
  G: complement = {(-1/(4  C^2),   0,         C) : C in F_101*}
  H: complement = {( 1/(216 C^2), -2/(9 C),   C) : C in F_101*}
"""
import sympy as sp
import numpy as np

x, y, z, w = sp.symbols("x y z w")
P = 101

# ---------- rebuild maps ----------
w0 = 1 + x * y
A0 = w0**2 * z + y**2 * (4 + 3 * x * y)
F = [sp.expand(w0 * A0), sp.expand(y + 3 * x * A0), sp.expand(2 * x - x**2 * (3 * y + x * z))]

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

# ---------- generic mod-p evaluator ----------
def image_codes(Fpolys, p=101):
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
    return out[0] * p * p + out[1] * p + out[2]

def curve_codes(param, p=101):
    """param: function C(int) -> (a,b,c) mod p rational params via Fractions"""
    from fractions import Fraction as Fr
    codes = set()
    for C in range(1, p):
        vals = param(Fr(C))
        trip = tuple((v.numerator % p) * pow(v.denominator % p, -1, p) % p for v in vals)
        codes.add(trip[0] * p * p + trip[1] * p + trip[2])
    return codes

from fractions import Fraction as Fr
curves = {
    "F": lambda C: (Fr(4, 27) / C**2, Fr(4, 3) / C, C),
    "G": lambda C: (Fr(-1, 4) / C**2, Fr(0), C),
    "H": lambda C: (Fr(1, 216) / C**2, Fr(-2, 9) / C, C),
}

for name, Mp in (("F", F), ("G", G), ("H", H)):
    codes = image_codes(Mp)
    hit = np.zeros(101**3, dtype=bool)
    hit[codes] = True
    complement = set(np.nonzero(~hit)[0].tolist())
    M = curve_codes(curves[name])
    print(f"{name}: |image| = {hit.sum()}, |complement| = {len(complement)} "
          f"(curve has {len(M)} pts);  complement == M mod 101 ? -> {complement == M}")
    if complement != M:
        extra = complement - M
        missing = M - complement
        print("   complement ⊄ M:", len(extra), " excess;", len(missing), " curve pts missing from complement")
