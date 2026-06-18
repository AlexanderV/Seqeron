# Evidence Artifact: RNA-PARTITION-001

**Test Unit ID:** RNA-PARTITION-001
**Algorithm:** RNA Partition Function (McCaskill) and Boltzmann Structure Probability
**Date Collected:** 2026-06-14

---

## Online Sources

### McCaskill JS (1990) — PubMed record

**URL:** https://pubmed.ncbi.nlm.nih.gov/1695107/
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Citation:** J. S. McCaskill, "The equilibrium partition function and base pair binding probabilities for RNA secondary structure," *Biopolymers* 1990, **29**(6-7):1105-19. DOI 10.1002/bip.360290621, PMID 1695107.
2. **Partition function definition:** The method computes "the full equilibrium partition function for secondary structure and the probabilities of various substructures."
3. **Complexity:** Both the partition function and the probabilities of all base pairs are computed by "a recursive scheme of polynomial order N³ in the sequence length N."
4. **Base-pair probabilities:** the derived pair binding probabilities are presented in a "box matrix" display over the full ensemble of probable alternative equilibrium structures.

### S. Will, MIT 18.417 (Fall 2011) — McCaskill lecture slides (inside recursion)

**URL:** https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf
**Accessed:** 2026-06-14 (fetched with WebFetch; PDF extracted with `pdftotext`)
**Authority rank:** 1 (graduate course built directly on the primary literature)

**Key Extracted Points:**

1. **Matrices:** `Q_ij := Z(P_ij)` (partition function of sub-ensemble), `Q^b_ij := Z(P^b_ij)` where `P^b_ij = { P ∈ P_ij | (i,j) ∈ P }`.
2. **Q recursion (verbatim):** `Q_ij = Q_{i,j-1} + Σ_{i≤k<j-m} Q_{i,k-1} · Q^b_{kj}`, with `Q_ij = 1` for `i ≥ j - m`.
3. **Total partition function:** "Partition function of the ensemble of S is `Z = Q_{1n}`."
4. **Complexity (verbatim):** "O(n²) space, O(n³) time (after bounding size of interior loops)."
5. **Probabilities (verbatim):** structure `Pr[P|S] = Z⁻¹ exp(−βE(P))`; base pair `Pr[(i,j)|S] = Z⁻¹ Σ_{P∋(i,j)} exp(−βE(P))`, with `β = 1/RT`.
6. **Unambiguous decomposition:** correctness "due to disjoint (=unambiguous) and independent decomposition."

### S. Will, MIT 18.417 (Fall 2011) — McCaskill base-pair probabilities (outside recursion)

**URL:** https://math.mit.edu/classes/18.417/Slides/mccaskill2.pdf
**Accessed:** 2026-06-14 (fetched with WebFetch; PDF extracted with `pdftotext`)
**Authority rank:** 1

**Key Extracted Points:**

1. **External base-pair probability (verbatim):** `p^E_kl = Z(P^E_kl)/Z = (Q_{1,k-1} · Q^b_{kl} · Q_{l+1,n}) / Q_{1n}`, where `(k,l)` is an external base pair (no enclosing pair).
2. **General case (outside recursion):** the full probability sums the external term PLUS, for every pair `(i,j)` that encloses `(k,l)`, the contribution `p_ij · Q^b_kl · (region between (i,j) and (k,l)) / (inside of (i,j))`. This enclosing term does **not** vanish in the flat fixed-per-pair model: a pair `(k,l)` that can be nested inside another pair gets probability strictly larger than its external term. Concretely the correct McCaskill probability is `p_kl = Q^b_kl · O_kl / Z`, where the **outside partition function** obeys
   `O_kl = Q_{1,k-1}·Q_{l+1,n} + Σ_{i<k, j>l, CanPair(i,j), j−i>m} w·Q_{i+1,k-1}·Q_{l+1,j-1}·O_ij`
   with `w = exp(−βE_bp)`. **Verified independently** this session against exhaustive Boltzmann-weighted enumeration on `GGGGCCCC` and `GGGAAACCC` (E_bp = 0 and −1) and on 300 random sequences: the outside recursion matches brute force to machine precision (max error 3.3e-16), whereas the external-only term is wrong for nestable pairs (e.g. `GGGAAACCC` P[2,6] = 6/20 = 0.30, not the external 1/20 = 0.05; `GGGGCCCC` P[1,5] = 3/16 = 0.1875, not 1/16). **CORRECTION (2026-06-16):** an earlier version of this Evidence doc claimed "only the external decomposition contributes, so p_kl = p^E_kl" in the flat model — that is FALSE and matched a bug in the implementation; the outside recursion above is required.

### Freiburg RNA Tools — McCaskill teaching tool (simplified energy model)

**URL:** https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 3 (educational reference implementation, built on McCaskill 1990)

**Key Extracted Points:**

1. **Partition function:** `Z = Σ_P exp(−E(P)/RT)` over all nested (pseudoknot-free) structures P.
2. **Simplified energy model (verbatim idea):** "each base pair of a structure contributes a fixed energy term `E_bp`."
3. **Parameters of the tool:** RNA sequence `S`, minimal loop length `l` (range 0–5), energy weight `E_bp`, normalized temperature `RT`.
4. **Pairing alphabet:** Watson-Crick and GU pairs are treated as complementary.

### ViennaRNA Package — Partition Function and Equilibrium Properties reference

**URL:** https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/pf_fold.html
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 3 (canonical reference implementation documentation)

**Key Extracted Points:**

1. **Boltzmann distribution (verbatim):** `p(s) ∝ e^(−βE(s))` with `β = 1/kT`, `k ≈ 1.987 × 10⁻³ kcal/(mol·K)`.
2. **Structure probability (verbatim):** `p(s) = e^(−βE(s)) / Z`.
3. **Default folding temperature** is 37 °C (= 310.15 K), the standard for the Turner/NNDB model.

---

## Documented Corner Cases and Failure Modes

### From MIT 18.417 slides

1. **Minimum loop length:** the recursion is only over `k` with `j − k > m`; pairs with too few enclosed bases are forbidden (`Q^b = 0`), preventing sterically impossible hairpins.
2. **Empty / base interval:** `Q_ij = 1` for `i ≥ j − m` — a sub-sequence too short to contain any pair has exactly one (empty) structure, contributing weight 1.

### From Freiburg tool

1. **Non-pairing alphabet:** only Watson-Crick (A-U, G-C) and GU pairs contribute; all other character pairs give `Q^b = 0`.

---

## Test Datasets

### Dataset: Hand-derived combinatorial partition functions (E_bp = 0 ⇒ pairWeight = 1)

When `E_bp = 0`, `exp(−β·0) = 1`, so `Z` counts the number of admissible (pseudoknot-free, min-loop = 3) secondary structures. These Z values are derived **two independent ways, both independent of the library code**: (a) a standalone re-implementation of the published `Q`/`Q^b` recurrence, and (b) an exhaustive brute-force enumeration of all non-crossing, base-disjoint subsets of admissible pairs. Both methods agree on every value below.

**Source:** McCaskill `Q`/`Q^b` recurrence (MIT 18.417 slides), min hairpin loop m = 3, Watson-Crick + GU pairing (Freiburg model); structure-counting via exhaustive non-crossing-subset enumeration.

| Sequence | E_bp | minloop m | Z (= #structures) |
|----------|------|-----------|-------------------|
| `AAAA` | 0 | 3 | 1 (no canonical pair exists) |
| `GC` | 0 | 3 | 1 (the only pair (0,1) has j−i=1 ≤ 3 → forbidden) |
| `GGGGCCCC` | 0 | 3 | 16 |
| `GGGAAACCC` | 0 | 3 | 20 |

Brute-force confirmation of `GGGGCCCC` and `GGGAAACCC`: enumerate every subset of admissible G-C pairs that is pairwise non-crossing and base-disjoint (the empty structure counts as 1). Both the recurrence and the exhaustive enumerator yield 16 and 20 respectively.

### Dataset: Exact base-pair probabilities (E_bp = 0 ⇒ P = #structures-containing-pair / Z)

**Source:** exhaustive non-crossing-subset enumeration (independent of library code), cross-checked against the McCaskill outside recursion. These values include **nestable** pairs (a pair that can be enclosed by another pair), for which the external decomposition alone is wrong.

`GGGGCCCC` (Z = 16): P(0,4)=1/16, P(0,5)=1/16, P(0,6)=2/16, P(0,7)=4/16, **P(1,5)=3/16**, P(1,6)=2/16, P(1,7)=2/16, **P(2,6)=3/16**, P(2,7)=1/16, P(3,7)=1/16.

`GGGAAACCC` (Z = 20): P(0,6)=1/20, P(0,7)=3/20, P(0,8)=6/20, **P(1,6)=3/20**, **P(1,7)=4/20**, **P(1,8)=3/20**, **P(2,6)=6/20**, **P(2,7)=3/20**, P(2,8)=1/20.

(Bold = nestable pairs whose external-only probability is strictly smaller than the truth; e.g. external P(2,6)=1/20=0.05 vs true 6/20=0.30.)

### Dataset: Exact weighted base-pair probabilities (E_bp = −1)

**Source:** Boltzmann-weighted exhaustive enumeration (weight = w^{#pairs}, w = exp(1/RT), RT = 0.61626805), independent of library code; matches the outside recursion to machine precision.

`GGGGCCCC` (Z = 180.018344830…): P(0,4)=0.02814492, P(0,6)=0.17074408, P(0,7)=0.45594238, P(1,5)=0.31334323, P(1,6)=0.17074408, P(1,7)=0.17074408, P(2,6)=0.31334323.

### Invariant: per-base pairing probability ≤ 1

For any position p, `Σ_{(i,j): i=p or j=p} P(i,j) ≤ 1` (a base pairs with at most one partner). Verified over 300 random sequences (lengths 0–12, E_bp ∈ {0, −1, −2, 0.5}): max single-base sum 0.983, and max |brute-force − outside-recursion| = 3.3e-16.

### Dataset: Boltzmann structure probability (closed form)

**Source:** McCaskill 1990 `Pr[P|S] = Z⁻¹ exp(−βE(P))`; ViennaRNA `p = exp(−βE)/Z`.

| Quantity | Value |
|----------|-------|
| R | 1.987 cal/(mol·K) |
| T | 310.15 K |
| RT | 1.987 × 310.15 / 1000 = 0.61626805 kcal/mol |
| If E_struct = E_ensemble | p = 1 (structure is the whole ensemble) |
| E_struct = −5, E_ensemble = −6 kcal/mol | p = exp((6−5)/RT)... = exp(−( −5 − (−6))/RT) = exp(−1/0.61626805) = 0.196817... |

---

## Assumptions

1. **ASSUMPTION: Simplified per-pair energy model.** The repository computes `Z` with a fixed per-base-pair energy `E_bp` (Freiburg teaching model) rather than the full Turner 2004 nearest-neighbour loop energies used by ViennaRNA/McCaskill's original tRNA examples. This is a documented simplification of the *energy model only*: the partition-function recurrence (`Q`, `Q^b`, `Z = Q₁ₙ`), the base-pair-probability formula, and all structural invariants (Z ≥ 1, probabilities in [0,1], normalisation, symmetry, monotonicity in E_bp) are fully conformant with McCaskill 1990 and are what the tests verify. Exact Turner-parameter ensemble energies are out of scope and would require porting the full Turner loop model into the partition recurrence.

---

## Recommendations for Test Coverage

1. **MUST Test:** `Z = 1` for a sequence with no admissible pair (e.g. `AAAA`, `GC`). — Evidence: `Q_ij = 1` base case (MIT slides).
2. **MUST Test:** `Z` equals the structure count for `E_bp = 0` on small sequences (`GGGGCCCC` → 16, `GGGAAACCC` → 20). — Evidence: combinatorial reduction of the McCaskill recurrence, confirmed by exhaustive enumeration.
3. **MUST Test:** every base-pair probability is in [0,1]; probabilities are symmetric `P[i,j] = P[j,i]` (only one orientation stored). — Evidence: `Pr[(i,j)|S]` is a probability (McCaskill 1990).
4. **MUST Test:** `Z ≥ 1` always (the empty structure always contributes weight 1). — Evidence: `Q_ij = 1` base case.
5. **MUST Test:** Boltzmann probability `CalculateStructureProbability` returns 1 when structure energy equals ensemble energy, and `exp(−ΔE/RT)` for a known ΔE. — Evidence: McCaskill 1990 / ViennaRNA `p = exp(−βE)/Z`.
6. **SHOULD Test:** lowering `E_bp` (more favourable pairing) strictly increases `Z` (monotonicity). — Rationale: `Z` is a sum of increasing exponential weights.
7. **SHOULD Test (property-based):** for randomly generated sequences, `Z ≥ 1` and all base-pair probabilities lie in [0,1]. — Rationale: invariants must hold for the whole input domain (O(n³) algorithm property test).
8. **COULD Test:** `GenerateRandomRna` is deterministic given a seeded `Random`; output length and alphabet are correct. — Rationale: reproducibility convention.

---

## References

1. McCaskill, J. S. (1990). The equilibrium partition function and base pair binding probabilities for RNA secondary structure. *Biopolymers* 29(6-7):1105-1119. https://doi.org/10.1002/bip.360290621 (PMID 1695107)
2. Will, S. (2011). McCaskill algorithm lecture slides, MIT 18.417 Fall 2011. https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf and https://math.mit.edu/classes/18.417/Slides/mccaskill2.pdf
3. Freiburg RNA Tools. McCaskill teaching tool. https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill
4. ViennaRNA Package. Partition Function and Equilibrium Properties (pf_fold reference). https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/pf_fold.html
5. Nussinov, R., & Jacobson, A. B. (1980). Fast algorithm for predicting the secondary structure of single-stranded RNA. *PNAS* 77(11):6309-6313. https://doi.org/10.1073/pnas.77.11.6309 (minimum-loop and nesting conventions)

---

## Change History

- **2026-06-14**: Initial documentation (RNA-PARTITION-001).
- **2026-06-16**: Validation session. Corrected the base-pair-probability evidence: the prior claim that the external decomposition suffices in the flat model is FALSE for nestable pairs; documented the McCaskill outside recursion `P[i,j] = Q^b·O/Z`, added exact probability datasets (E_bp=0 and −1) for all pairs incl. nested ones, and the per-base ≤1 invariant. Implementation was fixed to use the outside recursion (validation report RNA-PARTITION-001).
