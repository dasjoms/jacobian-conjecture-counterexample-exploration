"""Generate & batch-run Singular cases for all-multiple contact patterns."""
import subprocess
import sympy as sp

w = sp.symbols("w")

def family_seed(d):
    return 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - sp.Rational(6, d*(d+1)))

FIBERS = {}
for d in range(2, 7):
    FIBERS[f"family_d{d}"] = sp.expand(sp.integrate(family_seed(d), w))
FIBERS["G_seed"] = (w**2 - w**4) / 2                       # sanity: has bitangent
FIBERS["H_seed"] = sp.Rational(3,2)*w**4 - 4*w**3 + sp.Rational(5,2)*w**2

PATTERNS = {3: [(3,)], 4: [(2,2), (4,)], 5: [(3,2), (5,)],
            6: [(2,2,2), (3,3), (4,2), (6,)]}

def polystr(e, v="w"):
    frees = [str(s) for s in e.free_symbols]
    if v not in frees and len(frees) == 1:
        v = frees[0]
    e = sp.expand(e)
    if e == 0:
        return "0"
    terms = sp.Poly(e, sp.Symbol(v) if isinstance(v, str) else v).terms() if e.has(sp.Symbol(v) if isinstance(v, str) else v) else [((0,), e)]
    out = []
    for mon, coef in terms:
        if mon[0] == 0:
            out.append(f"({coef})")
        else:
            out.append(f"({coef})*{v}^{mon[0]}")
    return " + ".join(out)

cases = []
for name, Phi in FIBERS.items():
    n = sp.degree(Phi)
    for patt in PATTERNS.get(n, []):
        cases.append((name, Phi, patt))

lines = ["option(redSB);"]
for name, Phi, patt in cases:
    n = sp.degree(Phi)
    m = len(patt)
    Ws = [f"w{i+1}" for i in range(m)]
    K = m * (m - 1) // 2
    Ss = [f"s{i+1}" for i in range(K)]
    varstr = ",".join(Ws + ["P", "Q"] + Ss)
    label = f"{name} (n={n}) pattern {patt}"
    blk = [f'ring R{len(lines)} = 0,({varstr}),dp;']
    eqs = []
    for i, wi in enumerate(Ws):
        k = patt[i]
        Pi = sp.expand(Phi.subs(w, sp.Symbol(wi)))
        eqs.append(polystr(Pi, wi) + f" - P*{wi} - Q")
        eqs.append(polystr(sp.diff(Pi, sp.Symbol(wi)), wi) + " - P")
        for j in range(2, k):
            eqs.append(polystr(sp.diff(Pi, sp.Symbol(wi), j), wi))
    aux = 0
    for i in range(m):
        for j in range(i+1, m):
            aux += 1
            eqs.append(f"s{aux}*({Ws[i]}-{Ws[j]}) - 1")
    blk.append("ideal I = " + ",\n".join(eqs) + ";")
    blk.append(f'"=== CASE: {label} ===";')
    blk.append("I = groebner(I);")
    blk.append('"emptiness flag (0 => EMPTY / no such line):", size(reduce(1, I));')
    blk.append('if (size(reduce(1, I)) != 0) { "basis:"; I; }')
    lines.extend(blk)

open("/home/user/surjectivity.sing", "w").write("\n".join(lines))
print("generated", len(cases), "cases")
r = subprocess.run(["bash", "-c",
    "export PATH=/home/user/tools/bin:$PATH; timeout 1500 Singular -q -t /home/user/surjectivity.sing"],
    capture_output=True, text=True, timeout=1600)
print(r.stdout[-6000:])
print("STDERR tail:", r.stderr[-500:])
