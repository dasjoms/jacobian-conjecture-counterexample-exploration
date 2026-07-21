"""
Stage 6: d=5 (fiber sextic) over R -- hunt for real targets with NO real roots
(open gaps in R-image), and make figures (cusp geometry for F4 vs the rescued
constructed seed; census histogram).
"""
import sympy as sp
from sympy import symbols, diff, integrate, expand, Rational as R
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import base64, io

w = symbols("w")

def seed_family(d):
    p = 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - R(6, d*(d+1)))
    Phi = expand(integrate(p, w))
    return sp.expand(p), Phi

p5, Phi5 = seed_family(5)
print("d=5:  Phi6(w) =", Phi5)

rng = np.random.default_rng(7)
nogap = 0; tot = 0; examples = []
f = sp.lambdify(w, Phi5, "numpy")
for _ in range(4000):
    A0, B0 = float(rng.normal(0, 2)), float(rng.normal(0, 2))
    C0 = 1.0
    co = np.polynomial.polynomial.polyfit(
        np.linspace(-3, 3, 400), f(np.linspace(-3, 3, 400)) - B0*np.linspace(-3, 3, 400) + A0, 30)
    # better: direct coefficients
    coeffs = [Phi5.expand().coeff(w, k) for k in range(7)]
    coeffs = [float(c) for c in coeffs]
    coeffs[1] -= B0
    coeffs[0] += A0
    roots6 = np.roots(list(reversed(coeffs)))
    nreal = sum(1 for r in roots6 if abs(r.imag) < 1e-8)
    tot += 1
    if nreal == 0:
        nogap += 1
        if len(examples) < 3:
            examples.append((round(A0,4), round(B0,4), 1))
print(f"  random real targets with NO real fiber root: {nogap}/{tot}")
print("  examples (missed over R, open neighborhoods):", examples)

# ---------- figures ----------
fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

# (a) F4 cusp: Phi with inflection tangent NOT recrossing
p4, Phi4 = seed_family(4)
crits = sp.nroots(diff(p4, w), n=30)
w0 = float([r for r in crits if abs(sp.im(r)) < 1e-20][0])
s0 = float(p4.subs(w, w0)); r0 = w0*s0 - float(Phi4.subs(w, w0))
ww = np.linspace(-1.6, 2.2, 800)
Phi4n = sp.lambdify(w, Phi4, "numpy")
ax = axes[0]
ax.plot(ww, Phi4n(ww), lw=2.2, color="#1f77b4", label=r"$\Phi_4(w)$")
ax.plot(ww, s0*ww - r0, "--", color="crimson", lw=1.8, label="cusp tangent (slope $BC$, order-3 contact)")
ax.plot([w0], [float(Phi4.subs(w, w0))], "o", color="crimson", ms=9)
ax.annotate(f"cusp $w_0$={w0:.4f}\nouter roots complex\n(no recrossing => missed over R)",
            xy=(w0, float(Phi4.subs(w, w0))), xytext=(0.55, -0.8),
            arrowprops=dict(arrowstyle="->"), fontsize=10, color="crimson")
ax.set_title(r"(a) $F_4$: the cusp of $\Phi_4$ — real target $(%.5f, %.5f, 1)$ missed over $\mathbb{R}^3$" % (r0, s0),
             fontsize=11)
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# (b) constructed rescued seed: p~(w) = 5w^4 - 3w^2  ->  Phi = w^5 - w^3 ; cusp at 0
ax = axes[1]
Phic = lambda t: t**5 - t**3
ax.plot(ww, Phic(ww), lw=2.2, color="#1f77b4", label=r"$\tilde\Phi(w) = w^5 - w^3$")
ax.plot(ww, 0.0*ww, "--", color="seagreen", lw=1.8, label="cusp tangent at $w_0=0$ (the $w$-axis)")
ax.plot([0], [0], "o", color="seagreen", ms=9)
for rr in (-1, 1):
    ax.plot([rr], [0], "s", color="k", ms=8)
ax.annotate("rescued: two real\nrecrossings at w = ±1", xy=(1, 0), xytext=(0.2, 0.45),
            arrowprops=dict(arrowstyle="->"), fontsize=10, color="seagreen")
ax.set_title("(b) but rescued cusps DO exist off the normalized locus\n(seed $5w^4-3w^2$: violates $p(1)=-1$)", fontsize=11)
ax.set_ylim(-1.6, 1.6); ax.legend(fontsize=9); ax.grid(alpha=0.3)

# (c) census histogram
ax = axes[2]
hist = {0: 377700, 1: 386651, 2: 167800, 3: 91350, 4: 0, 5: 6800}
ax.bar(list(hist.keys()), list(hist.values()), color="#555588")
for k, v in hist.items():
    ax.text(k, v*1.01 + 3000, f"{v:,}", ha="center", fontsize=9)
ax.set_xlabel("preimages over $\\mathbb{F}_{101}$"); ax.set_ylabel("# targets")
ax.set_title("(c) $F_4$ over $\\mathbb{F}_{101}^3$: model = reality,\n0 mismatches on $101^3$ targets (note: no 4-fibers!)", fontsize=11)
ax.grid(alpha=0.3, axis="y")

plt.tight_layout()
plt.savefig("/home/user/realghost.png", dpi=150)
print("figure saved")
