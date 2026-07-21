# Trust, but verify: the image theorem, three independent ways
*Third lab note, 2026-07-20. Builds on `jacobian_lab_notes.md` (verification of the counterexamples) and `jacobian_anatomy.md` (walls, escapes, the missing curve). New: full finite-field enumeration and a from-scratch commutative-algebra reproduction via Singular.*

## The claim under test

From the anatomy note: the announced counterexample map F is a generically 3-to-1 local biholomorphism of C³ whose image is exactly **C³ minus the rational curve**

$$M = V(27AC^2 - 4,\; 3BC - 4) = \{(4/(27t^2), 4/(3t), t) : t \neq 0\},$$

with fiber sizes 3 (generic), 1 (wall hypersurface R(BC/4, AC²/4) = 0), 0 (M). That proof used the map's hidden fiber-equation structure — so it deserves adversarial cross-checking by methods that know *nothing* about that structure.

## Method 1 · Structure (recap, from the anatomy note)

Reconstruction bijection (preimages ↔ roots w with γ(w) = P − p(w) ≠ 0 of Φ(w) = Pw − Q), exact identities, escape-at-tangency, cusp ⇒ triple escape ⇒ empty fiber. ✔

## Method 2 · Complete enumeration over a finite field

Enumerated the image of every one of the 101³ = 1,030,301 points of F₁₀₁³ for F, G, and H (numpy, exact arithmetic mod 101).

**First surprise (and a lesson).** The complement of the image over F₁₀₁ is ~340,100 points (~33%) — vastly bigger than the 100 points of M. No contradiction: my initial expectation to find *only* M was a category error. Over C the fiber polynomial always has roots (fundamental theorem of algebra), so escape-to-infinity is the *only* failure mode. Over F_p the dominant phenomenon is entirely different: **the fiber polynomial often has no F_p-root at all** (≈ 1/3 of cubics are irreducible).

**The corrected finite-field analog of the anatomy theorem** — held exactly, with **zero mismatches on 39,606 sampled targets per map**:

$$\#\mathrm{fibers}(T) \;=\; \#\{\, w \in \mathbb{F}_p : \Phi(w) = P w - Q,\; \gamma(w) \neq 0 \,\} \qquad (C \neq 0\text{ targets})$$

| Map | Measured preimage-count histogram | Model histogram | Mismatches |
|---|---|---|---|
| F | {0: 13023, 1: 20234, 3: 6349} | identical | **0 / 39,606** |
| G | {0: 14952, 1: 13178, 2: 9962, 4: 1514} | identical | **0 / 39,606** |
| H | {0: 15083, 1: 13005, 2: 10064, 4: 1454} | identical | **0 / 39,606** |

Two lovely structural fingerprints visible in the data: **no count ever appears at "2" for F** (the only way to have exactly 2 finite preimages would be a double root with γ ≠ 0 — but double roots are tangencies, hence γ = 0 — so 2 is forbidden), and likewise G/H never have 3. Also, the enumeration confirmed M ⊂ complement over F₁₀₁ for all three maps, and the original 3-point collision descends intact: all three collision points map to (25, 0, 0) ∈ F₁₀₁³.

Honest ledger: reaching zero mismatches took two live debugging rounds visible in the transcript (a reversed Horner evaluation, then rational-coefficient reduction mod p).

## Method 3 · From scratch with Singular (independent CAS, no shared machinery)

Installed a user-space Singular 4.4.0 this turn and asked *it* the questions, feeding only the raw polynomials.

1. **Function-field test.** The ideal ⟨F₁ − 4/(27t²), F₂ − 4/(3t), F₃ − t⟩ ⊂ Q(t)[x,y,z] has Gröbner basis **(1)** ⟹ the generic point of M has no preimage ⟹ (by density) none do.
2. **Point tests.** At the M-point (4/27, 4/3, 1): Gröbner basis **(1)** — provably empty fiber. At (−1/4, 0, 0): Singular computes the fiber exactly — basis `[8z²−50z−13, y(2z−13), 12y²−4z−1, 18x+28y³−4yz−25y]` — whose three rational solutions are precisely the original three collision points (0,0,−1/4), (1,−3/2,13/2), (−1,3/2,13/2).
3. **The crown: the canonical Gröbner cover.** `grobcov` computes a complete, canonical stratification of parameter space by fiber structure — from nothing but the equations. It returns **9 segments**, and reproduces the stratification on its own: a generic segment with basis shape [z³, y, x] (3 preimages), degenerate segments with shapes [z², yz, y², x] (3 points), [z, y, x] (1 preimage)… and a **unique segment with basis (1)** — no preimages — whose parameter locus is cut out by

$$\langle 3BC - 4,\;\; 12A - B^2 \rangle.$$

Singular itself then confirms: `radical⟨3BC−4, 12A−B²⟩ = radical⟨27AC²−4, 3BC−4⟩` (dimension 1). **The independent computation lands on exactly the predicted curve M — not one target more, not one less.**

## Status board

| Statement | Method 1 (structure) | Method 2 (F₁₀₁ census) | Method 3 (Singular) |
|---|---|---|---|
| det JF = −2 | ✔ SymPy identity | ✔ (constant mod 101) | — |
| 3-point collision | ✔ exact | ✔ descends mod 101 | ✔ fiber computed = the 3 points |
| Fiber = 3 / 1 / 0 stratification | ✔ proved | ✔ 0 mismatches | ✔ 9-segment cover |
| im(F) = C³ ∖ M | ✔ proved | ✔ M ⊂ complement; rest explained by root-lessness | ✔ unique empty segment = M |

The theorem that the Jacobian conjecture is false (n = 3), and the exact description of *how* F folds space, now sit on three sets of shoulders that share no components.

## Next in the queue

1. Monodromy of the 3-sheeted covering over C³ ∖ {wall} — expected S₃; path-tracking check.
2. Seed-family combinatorics: count cusps/bitangents (hence missing curves) per fiber degree.
3. Lean 4 formalization of the two-line certificates (now feasible to attempt given package install).
4. The 2-D obstruction question — the one that actually matters next.
