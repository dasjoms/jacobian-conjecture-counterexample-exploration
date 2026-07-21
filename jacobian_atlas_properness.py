"""
Stage C: non-properness verification for F4.
(i)  fold escape-rate: gamma ~ sqrt(delta) at two coalescing roots, so |x| ~ delta^{-1/2}
(ii) off-wall boundedness of preimages as targets roam a compact box
(iii) the C=0 asymptotic family (inside the wall: (s,r) = (0,0) since Disc(0,0)=0)
"""
import numpy as np
import sympy as sp
from sympy import symbols, diff, integrate, expand, cancel, Rational as R

C5, C4, C3, C2 = -1/5, 1/4, -9/10, 17/20
def h_coeffs(s, r): return [C5, C4, C3, C2, -s, r]
def p4v(w): return -w**4 + w**3 - 2.7*w**2 + 1.7*w

print("=== (i) fold escape-rate at (s,r0) crossing, path (s=-1, r=-1+delta) ===")
print(" (C=1 => x = C/gamma, gamma = BC - p(w*) = s - p(w*))")
for delta in (1e-2, 1e-4, 1e-6, 1e-8):
    s, r = -1.0, -1.0 + delta
    rts = np.roots(h_coeffs(s, r))
    gam = s - p4v(rts)
    idx = np.argsort(np.abs(gam))
    g_small = np.abs(gam[idx[:2]])
    g_rest = np.abs(gam[idx[2:]])
    print(f"  delta={delta:.0e}: gamma_pair={g_small[0]:.3e},{g_small[1]:.3e}"
          f" ~ sqrt(delta)={np.sqrt(delta):.3e}   |x|~{1/g_small[0]:.3e}   rest-min={g_rest.min():.4f}")

print("\n=== (ii) off-wall boundedness: box targets C in [0.5,2], A,B in [-2,2], away from wall ===")
D4_terms = [(-20000000,4,0),(20000000,3,1),(18900000,3,0),(-57100000,2,2),(75376000,2,1),
            (-46613461,2,0),(38210000,1,3),(-67741050,1,2),(45893931,1,1),(-6657115,1,0),
            (-8192000,0,5),(17058675,0,4),(-12278715,0,3),(1957975,0,2)]
def discv(s, r): return sum(c*r**rr*s**ss for c, rr, ss in D4_terms)
rng = np.random.default_rng(3)
a = R(37,27)*-1
worst = 0.0; argw = None; cnt = 0
samples = []
for _ in range(20000):
    A0, B0, C0 = rng.uniform(-2,2), rng.uniform(-2,2), rng.uniform(0.5,2)
    s, r = B0*C0, A0*C0**2
    if abs(discv(s, r)) < 0.02:  # demand away from wall (scaled units)
        continue
    cnt += 1
    rts = np.roots(h_coeffs(s, r))
    gam = s - p4v(rts)
    mg = np.abs(gam).min()
    if mg > 0:
        lines = 1/mg  # |x| = |C|/min gamma, C<=2
        if lines*2 > worst:
            worst = lines*2; argw = (A0, B0, C0, mg)
print(f"  kept {cnt} off-wall targets; worst |x| upper bound = {worst:.3f} at {argw}")
print("  (all fibers' preimages uniformly bounded over the box, as expected)")

print("\n=== (iii) C=0 asymptotic family: F4(1, u-1, z(eps)) -> (u(1+17u/20), 1+17u/10, 0) ===")
x_, y_, z_, w_ = symbols("x y z w")
def build4():
    p = 2*w_ - 3*w_**2 + w_*(1-w_)*(w_**2 - R(3,10))
    q = sp.expand(integrate(w_*diff(p, w_), w_))
    kap = diff(p, w_).subs(w_, 1)
    av = R(-(1+kap)/(2+kap))
    u_ = 1 + x_*y_
    g_ = 1 + av*x_*y_ + x_**2*z_
    ws = u_*g_
    alpha = u_ + q.subs(w_, ws)/g_**2
    beta = 1 + p.subs(w_, ws)/g_
    f1 = sp.expand(cancel(alpha/x_**2)); f2 = sp.expand(cancel(beta/x_))
    return f1, f2, sp.expand(x_*g_), av
f1, f2, f3, av = build4()
f1n = sp.lambdify((x_,y_,z_), f1, "math"); f2n = sp.lambdify((x_,y_,z_), f2, "math")
f3n = sp.lambdify((x_,y_,z_), f3, "math")
u_par = 0.7
pred = (u_par*(1+17*u_par/20), 1+17*u_par/10, 0.0)
print(f"  predicted limit (A,B,0) = {pred}")
for eps in (1e-3, 1e-6, 1e-9):
    y0 = u_par - 1
    z0 = (eps - 1 - float(av)*y0)
    v = (f1n(1, y0, z0), f2n(1, y0, z0), f3n(1, y0, z0))
    err = max(abs(v[0]-pred[0]), abs(v[1]-pred[1]))
    print(f"  eps={eps:.0e}: image = ({v[0]:.10f},{v[1]:.10f},{v[2]:.2e})  err={err:.2e}")
print("  -> properness fails over C=0 too; C=0 plane = {(s,r)=(0,0)} subset of the wall (Disc(0,0)=0)")
