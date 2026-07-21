#!/usr/bin/env python3
# LAB NOTE 18 — THE 2-D CHAMBER / MOH'S SHADOW, stage 1: the compression algebra.
# =====================================================================
# The explainer recipe in 3D: F = (alpha/x^2, beta/x, x*gamma) with
# alpha = u + q(w)/gamma^2, beta = c + p(w)/gamma, u = x^2 y, w = u*gamma,
# gamma = 1 + a v + b t  (two warp knobs),  q' = w p'/c,  det = b c.
# THE COMPRESSION (one room fewer): F_2D(x,y) = (alpha/x^2, beta/x),
#   u = x^2 y,  v = x y,  gamma = gamma(v) univariate,  w = u gamma = x^2 y gamma(v).
# LOCKS (registered 2026-07-21 BEFORE compute):
#   K1: after the x-powers close via x = w/(v gamma), det JF_2D is EXACTLY a
#       function H(w, v) (no free x, no free y) for ARBITRARY functions p,q,gamma.
#   K2 (polynomial collapse): with p,q,gamma polynomial, gamma(0)=1, deg p<=4,
#       deg q<=5, deg gamma in {1,2}, the identity H(w,v) = B0 has NO solution
#       with deg p >= 2 (i.e. only the Moh-trivial classes: p affine or q const
#       or gamma const). Machine: exact Q-Groebner certified ideal equations.
#   K3 (entire rescue): gamma = exp(lambda v): there EXIST constants (lambda,b,c)
#       and a q-recipe making H(w,v) = b c for EVERY polynomial seed p with q
#       given by the explainer's own rule q' = w p'/c (q(0)=0). Predict: forced.
#   K4: the 2-D fiber equations are EXACTLY
#         (I)  Y w = c v gamma(v) + v p(w)
#         (II) X w^2 = v^2 (gamma(v) w + q(w))
#       for every preimage with x != 0 (verified symbolically).
#   K5: with seed p = 2w - 3w^2 (fiber-3), c = 4/3, lambda = b = -1/4 (explainer-ish
#       normalizations TBD by K3's answer), the map is NON-INJECTIVE and the
#       injectivity defect is visible in target-space walls = images of the
#       seed's wall curve (numeric, 50 digits, multi-start audits).
import sympy as sp
from sympy import Rational as Q, symbols, Function, exp, lambdify, groebner, solve
import json, time
t0 = time.time()
x, y, V, W, lam, c, b, B0 = symbols('x y V W lam c b B0')
p = Function('p'); q = Function('q'); ga = Function('ga')

v_ = symbols('v_')                                   # independent warp symbol  (stands for x*y)
w_expr = x**2 * y * ga(v_)

alpha = x**2*y + q(w_expr)/ga(v_)**2
beta  = c + p(w_expr)/ga(v_)
X_ = alpha/x**2
Y_ = beta/x
dX_x = sp.diff(X_, x); dX_y = sp.diff(X_, y)
dY_x = sp.diff(Y_, x); dY_y = sp.diff(Y_, y)
det = sp.expand(dX_x*dY_y - dX_y*dY_x)
# structural rewrites: express det in (W, V):   x*y = v_ ;  x = W/(v_*gamma)
ga_, dg_, p_, dp_, q_, dq_ = symbols('ga_ dg_ p_ dp_ q_ dq_')
det1 = det
# generic derivative-atom rewrite (head-driven)
for at in sorted(det1.atoms(sp.Derivative), key=str, reverse=True):
    fn = at.expr.func if hasattr(at.expr, 'func') else None
    if fn is ga:
        det1 = det1.subs(at, dg_*ga_)
    elif fn is p:
        det1 = det1.subs(at, dp_)
    elif fn is q:
        det1 = det1.subs(at, dq_)
det1 = det1.subs(ga(v_), ga_)
det1 = det1.subs(p(w_expr), p_).subs(q(w_expr), q_)
# close BOTH coordinates: x = W/(V ga_), then y = V^2 ga_/W  (from w=x^2 y ga, v=xy)
det2 = det1.subs(x, W/(V*ga_))
det2 = det2.subs(y, V**2*ga_/W).subs(v_, V)
det2 = sp.together(sp.expand(det2))
free_syms = set(det2.free_symbols) - {p_, dp_, q_, dq_, ga_, dg_, W, V, c, b, lam}
OUT = {"K1": {}, "log": []}
OUT["K1"]["free_after_closure"] = [str(s) for s in free_syms]
OUT["K1"]["x_free"] = (x not in det2.free_symbols) and (y not in det2.free_symbols) and (v_ not in det2.free_symbols)
OUT["K1"]["det_expr_before_gamma_choice"] = str(det2)[:3000]
print("K1 x_free:", OUT["K1"]["x_free"], " residual syms:", free_syms, flush=True)

json.dump(OUT, open('/home/user/keller2d_stage1a.json', 'w'), indent=1, default=str)
# save the raw det (function-form) for stage B/C reuse
open('/home/user/keller2d_H_raw.txt', 'w').write(sp.sstr(det2))
print("stage1a done %.1fs" % (time.time()-t0))
