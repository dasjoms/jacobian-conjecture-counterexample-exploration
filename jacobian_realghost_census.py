"""
Stage 4: F_101 census for F4 — actual fiber sizes from full enumeration of
101^3 source points vs the fiber-equation model:
  C != 0: count = #{w in F_p : Phi(w) - BC w + AC^2 = 0,  BC - p(w) != 0}
  C = 0 : count = 1 (flat sheet, S affine-invertible) + gamma-sheet roots
Everything exact over F_101. Reports mismatches and the count histogram.
"""
import sympy as sp
import numpy as np
from sympy import symbols, diff, integrate, expand, cancel, Rational as R

x, y, z, w = symbols("x y z w")
P = 101

def build(d):
    p = 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - sp.Rational(6, d*(d+1)))
    q = sp.expand(sp.integrate(w*sp.diff(p, w), w))
    kap = sp.diff(p, w).subs(w, 1)
    a = sp.Rational(-(1 + kap) / (2 + kap))
    u = 1 + x*y
    g = 1 + a*x*y + x**2*z
    ws = u*g
    alpha = u + q.subs(w, ws)/g**2
    beta = 1 + p.subs(w, ws)/g
    f1 = sp.expand(sp.cancel(alpha/x**2))
    f2 = sp.expand(sp.cancel(beta/x))
    f3 = sp.expand(x*g)
    return sp.expand(p), sp.expand(sp.integrate(p, w)), f1, f2, f3

print("building F4 ...")
p, Phi, f1, f2, f3 = build(4)
print("  terms:", [len(sp.Poly(f, x, y, z).terms()) for f in (f1, f2, f3)])

inv = lambda a, m: pow(int(a) % m, -1, m)
def terms_mod(f):
    out = []
    for (i, j, k), c in sp.Poly(f, x, y, z).terms():
        num, den = sp.fraction(c)
        out.append((i, j, k, int(num) % P * inv(int(den), P) % P))
    return out

X = np.arange(P, dtype=np.int64)[:, None, None]
Y = np.arange(P, dtype=np.int64)[None, :, None]
Z = np.arange(P, dtype=np.int64)[None, None, :]

def evaluate(terms):
    acc = np.zeros((P, P, P), dtype=np.int64)
    Xp = [ (X**i) % P for i in range(18) ]
    Yp = [ (Y**j) % P for j in range(18) ]
    Zp = [ (Z**k) % P for k in range(8) ]
    for i, j, k, c in terms:
        acc = (acc + c * ((Xp[i] * Yp[j]) % P) * Zp[k]) % P
    return acc

print("evaluating F4 on F_101^3 ...")
T1 = evaluate(terms_mod(f1)); print("  f1 done")
T2 = evaluate(terms_mod(f2)); print("  f2 done")
T3 = evaluate(terms_mod(f3)); print("  f3 done")

idx = (T1 * (P*P) + T2 * P + T3).ravel()
actual = np.bincount(idx, minlength=P**3).astype(np.int64)
print("enumeration done; distinct targets hit:", int((actual > 0).sum()))

# ---------- model ----------
pco = [sp.expand(p).coeff(w, k) for k in range(5)]
Phico = [sp.expand(Phi).coeff(w, k) for k in range(6)]
pc = np.array([int(sp.fraction(c)[0]) % P * inv(sp.fraction(c)[1], P) % P for c in pco], dtype=np.int64)
Phc = np.array([int(sp.fraction(c)[0]) % P * inv(sp.fraction(c)[1], P) % P for c in Phico], dtype=np.int64)
wv = np.arange(P, dtype=np.int64)
pw = np.zeros(P, dtype=np.int64)
for k in range(5): pw = (pw + pc[k] * wv**k) % P
Phw = np.zeros(P, dtype=np.int64)
for k in range(6): Phw = (Phw + Phc[k] * wv**k) % P

# table over (s, r): roots of h = Phi - s w + r with gamma = s - p(w) != 0
model = np.zeros(P**3, dtype=np.int64)
gam_nz = np.zeros((P, P), dtype=bool)   # gam_nz[s, w] = (s - p(w)) % P != 0
for s in range(P):
    gam_nz[s] = ((s - pw) % P) != 0
h_base = Phw[None, :].repeat(P, axis=0) - (np.arange(P)[:, None] * wv[None, :])  # h_s(w) = Phi(w) - s w
h_base %= P
print("building (s,r) table ...")
for s in range(P):
    hs = h_base[s]
    for r in range(P):
        roots = ((hs + r) % P == 0) & gam_nz[s]
        n = int(roots.sum())
        # all (A,B,C) with C != 0, BC = s, AC^2 = r
        for C in range(1, P):
            B = s * inv(C, P) % P
            A = r * inv(C*C % P, P) % P
            model[(A*P + B)*P + C] = n

# C = 0: flat sheet always exactly 1 preimage (S has det 1 and is affine/quadratic:
# y = 27 B / 10, z from first coord). Gamma sheet quadratic: A*(1+p1 u)^2 - B^2 u (1+q2 u).
p1v = 17 * inv(10, P) % P
q2v = 17 * inv(20, P) % P
def gamma_count(A, B):
    if B != 0:
        a2 = (A*p1v*p1v - B*B*q2v) % P
        a1 = (2*A*p1v - B*B) % P
        a0 = A % P
        cnt = 0
        for uu in range(P):
            if (a2*uu*uu + a1*uu + a0) % P == 0:
                if (1 + p1v*uu) % P != 0:  # x = (1+p1 u)/B must be nonzero
                    cnt += 1
        return cnt
    else:
        if A == 0:
            return 0
        # B = 0: u = -1/p1 forced; x^2 = u(1+q2 u)/A
        uu = (-inv(p1v, P)) % P
        val = uu * (1 + q2v*uu) % P * inv(A, P) % P
        if val == 0: return 0
        return 2 if pow(val, (P-1)//2, P) == 1 else 0
print("adding C=0 model ...")
for A in range(P):
    for B in range(P):
        model[(A*P + B)*P + 0] = 1 + gamma_count(A, B)

mism = actual != model
print("MISMATCHES:", int(mism.sum()))
if mism.sum():
    bad = np.nonzero(mism)[0][:10]
    for t in bad:
        print(f"   target ({t//P//P},{t//P%P},{t%P}): actual {actual[t]} model {model[t]}")

# histogram of actual preimage counts
mx = int(actual.max())
print("\npreimage-count histogram (actual, all 101^3 targets):")
for k in range(mx+1):
    print(f"   {k}: {int((actual == k).sum())}")
print("\nmax fiber size:", mx)
