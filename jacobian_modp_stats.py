"""
The finite-field analog of the anatomy theorem: for every target in F_p^3,
#preimages is governed EXACTLY by the fiber equation:
  #preimages = #{ F_p-roots w of Phi(w) - P w + Q with gamma(w) != 0 }
Full-enumeration check: 40k random sampled targets per map, predicted vs actual.

Over C, every cubic/quartic has roots (FTA), so escape is the ONLY failure mode
(missing curve M). Over F_p the dominant phenomenon is different: the fiber
polynomial may have NO F_p-root at all (~1/3 of cubics are irreducible!).
"""
import sympy as sp
import numpy as np

x, y, z, w = sp.symbols("x y z w")
p = 101

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

def image_counts(Fpolys):
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
    return np.bincount(code, minlength=p**3)

inv4 = pow(4, -1, p)
W = np.arange(p, dtype=np.int64)[None, :]          # (1, 101)

def check(name, Fpolys, phi_desc, pfun, PQfun, n_targets=40000, seed=7):
    counts = image_counts(Fpolys)
    rng = np.random.default_rng(seed)
    sample = rng.choice(p**3, size=n_targets, replace=False)
    A = sample // (p * p); B = (sample // p) % p; C = sample % p
    Pv, Qv = PQfun(A, B, C)
    # fiber polynomial mod p: coefficient arrays per target (low->high order, len 5)
    d = sp.degree(phi_desc, w)
    cd = {}
    for (me,), coef in sp.Poly(phi_desc, w).terms():
        rc = sp.Rational(coef)
        red = (rc.p % p) * pow(int(rc.q % p), -1, p) % p
        cd[int(me)] = np.full(A.shape, np.int64(red))
    cd.setdefault(1, np.zeros_like(A)); cd.setdefault(0, np.zeros_like(A))
    cd[1] = (cd[1] - Pv) % p
    cd[0] = (cd[0] + Qv) % p
    val = np.zeros((n_targets, p), dtype=np.int64)
    for e in range(d, -1, -1):
        val = (val * W + cd.get(e, np.zeros_like(A))[:, None]) % p
    roots_is = (val == 0)
    # gamma(w) = P - pfun(w)
    pw = pfun(W)
    gam = (Pv[:, None] - pw) % p
    pred = (roots_is & (gam != 0)).sum(axis=1)
    actual = counts[sample]
    nonzeroC = C != 0
    mism = (pred != actual) & nonzeroC
    note = (pred != actual) & ~nonzeroC
    print(f"{name}: targets checked {nonzeroC.sum()}; mismatches (C!=0): {mism.sum()}; "
          f"(C=0 targets, handled by sheet analysis: {note.sum()} of {(~nonzeroC).sum()})")
    print(f"   histogram of actual preimage counts: "
          f"{dict(zip(*np.unique(counts[sample[nonzeroC]], return_counts=True)))}")
    # theoretical: distribution of #Fp-roots (gamma!=0) -- maximum model detail
    print(f"   histogram predicted (gamma-escape-adjusted): "
          f"{dict(zip(*np.unique(pred[nonzeroC], return_counts=True)))}")

check("F", F, -w**3 + w**2, lambda W_: (2 * W_ - 3 * W_ * W_) % p,
      lambda A, B, C: ((B * C % p) * inv4 % p, (A * C % p * C % p) * inv4 % p))
check("G", G, -w**4 / 2 + w**2 / 2, lambda W_: (W_ - 2 * W_ * W_ * W_) % p,
      lambda A, B, C: (B * C % p, (A * C % p * C % p) * pow(2, -1, p) % p))
check("H", H, sp.Rational(3, 2) * w**4 - 4 * w**3 + sp.Rational(5, 2) * w**2,
      lambda W_: (5 * W_ - 12 * W_ * W_ + 6 * W_ * W_ * W_) % p,
      lambda A, B, C: (B * C % p, A * C % p * C % p))
