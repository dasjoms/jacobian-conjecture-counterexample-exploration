# The Jacobian counterexample — sandbox lab notes
*2026-07-20. All symbolic checks below were run exactly (rational/algebraic arithmetic) unless labeled "numeric". Scripts: `verify_jacobian.py`, `jacobian_structure.py`, `jacobian_family_and_fibers.py`, `jacobian_homemade.py`, `jacobian_census2.py` (+ inline experiments).*

## 0 · TL;DR

| # | Experiment | Result |
|---|---|---|
| 1 | Alpöge's announced map F: det J constant? | ✅ `det JF = −2` exactly |
| 2 | Three distinct points collide? | ✅ all → `(−1/4, 0, 0)` |
| 3 | Hidden structure of F | ✅ factored form, odd symmetry, weighted scaling |
| 4 | Mechanism of failure | ✅ explicit "escape curve": properness fails |
| 5 | Fiber census of F | ✅ generic fiber = 3 points; special fibers 1 |
| 6 | Explainer's claimed degree-4 map G | ✅ fully verified: polynomial, `det JG = −6`, integer collision |
| 7 | **Homegrown NEW counterexample H** | ✅ minted: `det JH = 1`, exact algebraic collision |
| 8 | Float64 pitfall | ⚠️ numeric naivety destroys the det identity at moderate scales |
| 9 | Reduction mod p | ✅ collision descends to F₁₀₁: det-constant ≠ injective over F_p |

## 1 · The announced map, in its hidden normal form

With `w = 1+xy`, `A = w²z + y²(4+3xy)`, `B = 3y + xz`, the posted map factors as

$$F = (\,wA,\;\; y + 3xA,\;\; 2x - x^2 B\,)$$

and SymPy confirms `det JF = −2` identically, plus an unexpected **odd symmetry**

$$F(-x,-y,z) = (f_1, -f_2, -f_3),$$

which explains why the three collision points come as a symmetric pair `(±1, ∓3/2, 13/2)` plus a fixed point `(0,0,−1/4)`.

## 2 · Why a constant-Jacobian map can fail — the escape hatch

The inverse function theorem gives a local inverse *near every finite point*. Global failure can only occur by a preimage escaping to infinity while its image stays bounded (that's exactly what "non-proper" means). The map advertises its own escape route:

$$F\big(t,\; -1/t,\; 5/t^2\big) = (0,\; 2/t,\; 0) \;\longrightarrow\; (0,0,0) \quad\text{as } t\to\infty,$$

and more generally the whole divisor `{xy = −1}` maps **bijectively** onto `{0}×C*×C`: `F(x, −1/x, z) = (0, 2/x, 5x − x³z)`. (An early internet objection — "the map isn't proper, so no contradiction" — had it backwards: non-properness is not a rescue of the conjecture, it *is* the loophole the counterexample exploits.)

## 3 · Fiber census (exact fiber-equation method)

Fibers of F over `(A,B,C)` correspond bijectively (off `C=0`) to roots of one fixed cubic:

$$w^2 - w^3 = Pw - Q, \qquad P = BC/4,\; Q = AC^2/4,$$

each root reconstructing a unique preimage. The cubic's discriminant is the nonzero polynomial `−4P³ + P² + 18PQ − 27Q² − 4Q`, so **generic fibers have exactly 3 points**; reconstructions verified against F with residuals ≤ 1e-13:

| Target | Preimages found |
|---|---|
| `(0,0,1)` | **1** (two sheets escaped to ∞; boundary fiber) |
| `(−1/4,0,1)` | 3 real |
| `(−8/27,0,1)` | 3 real, incl. the explainer's exact radical triple |
| `(1,2,3)` | 3 (one real + conjugate pair) |

The fiber over the origin is provably just the origin (short case analysis: on `w=0` the map has no zeros; on `A=0` fibers force `y=z=x=0`).

## 4 · The explainer's degree-4 map G — independently confirmed

Claimed on jacobianfun.org, unverified until run here:
`G = ((2u+u²−3u⁴γ²)/x², (1+u−2u³γ²)/x, xγ)` with `u = 1+3xy`, `γ = 1−4xy−x²z`.

| Claim | Check |
|---|---|
| G is polynomial (divisions cancel) | ✅ degrees (12, 11, 4) |
| det J constant | ✅ `−6` exactly |
| Collision certificate | ✅ `G(1,0,0) = G(−1,0,2) = (0,0,1)` |
| Four real preimages over `(−1/8,0,1)` | ✅ exact; roots are `±cos(π/8), ±sin(π/8)` |
| Generic fiber degree 4 | ✅ fiber quartic `(w²−w⁴)/2 = Pw − Q` |

## 5 · Minting a counterexample of my own (map H)

To test whether the construction is a *machine* (not two lucky artifacts), I rolled a different seed satisfying the recipe's three conditions — `p(w) = 5w − 12w² + 6w³` with `p(0)=0`, `p(1)=−1`, `∫₀¹p = 0` — and ran the lift with `b = c = 1, a = 0`:

$$H = \big(\alpha/x^2,\; \beta/x,\; x\gamma\big),\quad \alpha = u + \tfrac{q(w)}{\gamma^2},\ \beta = 1 + \tfrac{p(w)}{\gamma},\ u = 1+xy,\ \gamma = 1 + x^2z,\ w = u\gamma .$$

Results, all exact: H is polynomial (degrees 12,11,4); `det JH = 1` identically; and designing the fiber quartic `Φ = (3/2)w⁴ − 4w³ + (5/2)w²` to have `Φ(2+d) = Φ(2−d)` gives `d² = −5/4`, i.e. two exact conjugate algebraic points which SymPy verifies map to the same target `(441/32, 0, 1)`. Census: generic fibers = 4 points; the fiber over `(0,0,1)` has 2 (routes through the integer point `(1,0,0)`).

**So the sandbox now holds three distinct explicit counterexamples — F (fiber degree 3), G (4), H (4, different det and coefficients).**

## 6 · Two cautionary tales from the numeric trenches

- **Float64 destroys the determinant.** Far from the origin, evaluating det JF numerically involves canceling monomials of size ~10¹⁵ down to the constant −2 — machine epsilon swallows it entirely (measured det values ranged from 0 to 892,258 at |coordinates| ~ 28). Naive Newton-based root finding therefore fails; exact arithmetic or the fiber-equation structure is mandatory. A fitting epistemic echo of the whole story: *check symbolically or don't check at all.*
- **Mod-p reality check.** The collision points are rational with denominators 2 and 4, so they descend to F_p for every odd prime. Exact check mod 101: all three points → `(25, 0, 0)`, and `det ≡ −2 ≢ 0`. So the constant-Jacobian hypothesis fails to give injectivity over F₁₀₁ as well, and none of F, G, H is a permutation of F₁₀₁³. (My initial expectation of the opposite taught me something; the reduction is a good extra certificate that the collision is a *rational identity*, not a fluke at special reals.)

## 7 · What I'd try next

1. **Formalize the certificates in Lean** — det identity and the collisions are short, routine, and would be the highest-trust version of today's news.
2. **The repair program**: implement the proposed patched conjecture (constant Jacobian + proper *on the nose*, or injective-on-the-axis variants) and probe it on F/G/H.
3. **Census the boundary**: the map of missing/escaped targets (fiber size < degree) is still only partially charted — §10 of the explainer flags it as open.
4. Hunt for a degree-3 seed giving fiber degree 4 with *integer* coefficients smaller than H's, and generally minimize "coefficient height" over the seed family.
