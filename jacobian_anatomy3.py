"""
Anatomy part 3: exact (proof-grade) verification of the emptiness claims
for the three missing curves, plus figures.

Missing curves (targets that have NO preimage):
  F:  M_F = { (4/(27 C^2), 4/(3 C), C) : C != 0 }        (cusp, triple escape)
  G:  M_G = { (-1/(4 C^2), 0, C)        : C != 0 }        (bitangent, 2x2 escape)
  H:  M_H = { (1/(216 C^2), -2/(9 C), C): C != 0 }        (bitangent, 2x2 escape)
Mechanism: preimages over C != 0 targets correspond bijectively to roots w of
Phi(w) = Pw - Q with gamma(w) = P - Phi'(w) != 0 (gamma = 0 puts the preimage
at infinity). Multiple roots automatically have gamma = 0 (tangency), and at
the special points below ALL roots are multiple -> nothing finite remains.
"""
import sympy as sp

w, P, Q = sp.symbols("w P Q")

def prove_empty(name, Phi, P0, Q0, note):
    print(f"==== {name}: fiber polynomial at (P,Q) = ({P0}, {Q0})   {note}")
    poly = sp.factor(Phi - P0 * w + Q0)
    print("   Phi - P0 w + Q0 =", poly)
    fac = sp.factor_list(poly)[1]
    all_mult = all(m >= 2 for _, m in fac)
    gam = P - sp.diff(Phi, w)
    for rpt, m in fac:
        rts = sp.solve(rpt, w)
        for r in rts:
            g = sp.simplify(gam.subs({P: P0, w: r}))
            print(f"     root {r} (mult {m}): gamma = {g}")
    print("   every root multiple:", all_mult,
          " => every root escapes => fiber EMPTY (x=0 impossible since C != 0)")

PhiF = w**2 - w**3
PhiG = (w**2 - w**4) / 2
PhiH = sp.Rational(3, 2) * w**4 - 4 * w**3 + sp.Rational(5, 2) * w**2

prove_empty("F", PhiF, sp.Rational(1, 3), sp.Rational(1, 27), "(cusp)")
prove_empty("G", PhiG, 0, -sp.Rational(1, 8), "(bitangent)")
prove_empty("H", PhiH, -sp.Rational(2, 9), sp.Rational(1, 216), "(bitangent)")

# extra consistency: F at cusp-lift target (4/27, 4/3, 1):
print("\ncheck P,Q from F-target (4/27,4/3,1): P = BC/4 =", sp.Rational(4, 3) / 4,
      ", Q = AC^2/4 =", sp.Rational(4, 27) / 4)
print("check P,Q from G-target (-1/4,0,1): P = BC =", 0, ", Q = AC^2/2 =", -sp.Rational(1, 4) / 2)
print("check P,Q from H-target (1/216,-2/9,1): P = BC =", -sp.Rational(2, 9),
      ", Q = AC^2 =", sp.Rational(1, 216))

# ---------------- figures ----------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    # panel 1: Phi_F with three tangent/intersection regimes
    ww = np.linspace(-1.2, 1.6, 800)
    phi = ww**2 - ww**3
    axes[0].plot(ww, phi, "k-", lw=2, label="Phi(w) = w^2 - w^3")
    axes[0].plot(ww, 0 * ww - 1 / 27 + 1 / 3 * ww, "r--", lw=1.5,
                 label="inflection tangent (P=1/3, Q=1/27): triple root -> MISSING")
    axes[0].plot(ww, 0 * ww - (-0.074), "g--", lw=1.5, label="secant: 3 roots -> 3 preimages")
    axes[0].plot(ww, 0 * ww - (1 / 27) + 0.0 * ww, "b--", lw=1.5, label="wall tangent at w=0: 1 preimage")
    axes[0].set_title("Fiber cubic of F: line meets cubic")
    axes[0].legend(fontsize=8); axes[0].set_xlabel("w"); axes[0].set_ylabel("y")

    # panel 2: wall curve Gamma_F in (P,Q) space with cusp
    wt = np.linspace(-0.9, 1.7, 1200)
    Pc = 2 * wt - 3 * wt**2
    Qc = wt**2 - 2 * wt**3
    axes[1].plot(Pc, Qc, "purple", lw=2, label="wall curve Gamma = (p(w), q(w))")
    axes[1].plot([1 / 3], [1 / 27], "r*", ms=16, label="cusp (1/3, 1/27) => missing targets")
    axes[1].annotate("fiber = 3 off-curve\nfiber = 1 on curve\nfiber = 0 at cusp lifts",
                     xy=(0.05, -0.4), fontsize=9)
    axes[1].set_title("Wall curve of F in (P,Q)-space")
    axes[1].set_xlabel("P = BC/4"); axes[1].set_ylabel("Q = AC^2/4")
    axes[1].legend(fontsize=8, loc="lower left")
    plt.tight_layout()
    plt.savefig("/home/user/anatomy_F.png", dpi=130)

    # panel 3 (separate fig): G's bitangent
    fig2, axq = plt.subplots(figsize=(6, 4.2))
    phi4 = (ww**2 - ww**4) / 2
    axq.plot(ww, phi4, "k-", lw=2, label="Phi_G(w) = (w^2 - w^4)/2")
    axq.plot(ww, np.full_like(ww, 0.125), "r--", lw=1.5,
             label="bitangent y = 1/8: both contact roots escape -> MISSING")
    axq.plot([-1 / np.sqrt(2)], [0.125], "ro"); axq.plot([1 / np.sqrt(2)], [0.125], "ro")
    axq.set_title("Fiber quartic of G: the horizontal bitangent")
    axq.legend(fontsize=8); axq.set_xlabel("w")
    plt.tight_layout()
    plt.savefig("/home/user/anatomy_G.png", dpi=130)
    print("\nfigures written: anatomy_F.png, anatomy_G.png")
except ImportError as e:
    print("matplotlib unavailable, skip figures:", e)
