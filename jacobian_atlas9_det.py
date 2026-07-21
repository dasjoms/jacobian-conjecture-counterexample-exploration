"""
Note 11, stage F: pointwise-exact determinants. Tonight: d = 9 (F9) = 1 at 5/5,
PLUS the queue debt from note 10: retro det checks for d = 6 and d = 7.
"""
import sympy as sp, random
x, y, z, w = sp.symbols("x y z w")

def det_at_points(d, npts=5, seed=99):
    p = sp.expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - sp.Rational(6, d*(d+1))))
    c = sp.Integer(1); b = sp.Integer(1)
    q = sp.expand(sp.integrate(w*sp.diff(p, w)/c, w))
    kap = sp.Rational(sp.diff(p, w).subs(w, 1))/c
    a = sp.Rational(-(1+kap)/(2+kap))
    u = 1 + x*y
    g = 1 + a*x*y + b*x**2*z
    ws = u*g
    alpha = u + q.subs(w, ws)/g**2
    beta = c + p.subs(w, ws)/g
    M = sp.Matrix([alpha/x**2, beta/x, x*g]).jacobian([x, y, z])
    random.seed(seed + d)
    ok = 0
    for _ in range(npts):
        while True:
            pt = {x: sp.Rational(random.randint(-4,4), random.randint(1,5)),
                  y: sp.Rational(random.randint(-4,4), random.randint(1,5)),
                  z: sp.Rational(random.randint(-4,4), random.randint(1,5))}
            if pt[x] != 0 and sp.simplify(g.subs(pt)) != 0: break
        vals = [entry.subs(pt) for entry in M]
        det = sp.Matrix(3, 3, lambda ii, jj: vals[3*ii+jj]).det()
        if sp.simplify(det - 1) == 0: ok += 1
    return ok

for d in (6, 7, 9):
    ok = det_at_points(d)
    print(f"d = {d} (fiber {d+1}): det JF = 1 at {ok}/5 exact rational points", flush=True)
