import sympy as sp
import numpy as np

x, y, z = sp.symbols("x y z")
w0 = 1 + x * y
A0 = w0**2 * z + y**2 * (4 + 3 * x * y)
F = [sp.expand(w0 * A0), sp.expand(y + 3 * x * A0), sp.expand(2 * x - x**2 * (3 * y + x * z))]
J = sp.Matrix(F).jacobian([x, y, z])

Ff = sp.lambdify([x, y, z], F, "numpy")
Jf = sp.lambdify([x, y, z], J, "numpy")

n = 5000
rng = np.random.default_rng(0)
X = (rng.uniform(-20, 20, (n, 3)) + 1j * rng.uniform(-20, 20, (n, 3))).astype(complex)
T = np.array([-0.25, 0, 0], dtype=complex)

for it in range(150):
    xa, ya, za = X[:, 0], X[:, 1], X[:, 2]
    Fv = np.array(Ff(xa, ya, za), dtype=complex).reshape(3, -1).T - T
    Jm = np.array(Jf(xa, ya, za), dtype=complex)  # (3,3,n)
    Jm = np.moveaxis(Jm, -1, 0)                    # (n,3,3)
    d = np.linalg.solve(Jm, Fv[..., None])[..., 0]
    X = X - d
    if it % 30 == 29:
        res = np.linalg.norm(Fv, axis=1)
        print(it + 1, "residual quantiles:", np.percentile(res, [10, 50, 90, 99]))

res = np.linalg.norm(np.array(Ff(X[:, 0], X[:, 1], X[:, 2]), dtype=complex).reshape(3, -1).T - T, axis=1)
print("converged:", np.sum(res < 1e-9), "of", n)
conv = X[res < 1e-9]
seen = {}
for v in conv:
    k = tuple(np.round(v, 4))
    seen.setdefault(k, v)
print("distinct roots:", len(seen))
for k in sorted(seen):
    print(k)
