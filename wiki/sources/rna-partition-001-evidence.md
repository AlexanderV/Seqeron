---
type: source
title: "Evidence: RNA-PARTITION-001 (RNA partition function — McCaskill algorithm, base-pair probabilities)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-PARTITION-001-Evidence.md
sources:
  - docs/Evidence/RNA-PARTITION-001-Evidence.md
source_commit: c74e2076b6e891a02c917198a54544896c4dbafa
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-PARTITION-001

The validation-evidence artifact for test unit **RNA-PARTITION-001** — the **RNA partition function
(McCaskill algorithm)**: the Boltzmann-weighted ensemble computation of the equilibrium partition
function `Z`, per-base-pair binding probabilities, and the Boltzmann probability of a given structure.
It is one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the synthesizing concept is [[rna-partition-function-mccaskill]]. [[test-unit-registry]]
tracks the unit.

This is the **probabilistic / ensemble** layer of the RNA secondary-structure family — the
Boltzmann-weighted counterpart of [[rna-minimum-free-energy-folding|MFE folding]] (which returns a
*single* optimal structure). Instead of one minimum-energy fold it computes the full equilibrium
distribution over all pseudoknot-free structures and, from it, the probability that any given pair
`(i,j)` forms.

## What this file records

- **Algorithm:** the **McCaskill (1990)** recursion. Both `Z` and all base-pair probabilities are
  computed by a recursive scheme of polynomial order **O(n³) time, O(n²) space** (after bounding
  interior-loop size), over all nested (pseudoknot-free) structures.

- **Authoritative sources:**
  - **McCaskill, J. S. (1990)** *Biopolymers* 29(6-7):1105–1119 (rank 1, primary; PMID 1695107) —
    "the full equilibrium partition function for secondary structure and the probabilities of various
    substructures"; base-pair binding probabilities in a "box matrix" over the ensemble.
  - **S. Will, MIT 18.417 (2011)** McCaskill lecture slides (rank 1) — the inside `Q`/`Q^b` recursion
    and the outside base-pair-probability recursion, verbatim.
  - **Freiburg RNA Tools — McCaskill teaching tool** (rank 3) — the simplified **fixed-per-pair energy
    model** (`E_bp`), min-loop length `l`, normalized temperature `RT`, Watson-Crick + GU pairing.
  - **ViennaRNA Package — pf_fold reference** (rank 3) — the Boltzmann distribution `p(s) = e^(−βE(s))/Z`,
    `β = 1/kT`, `k ≈ 1.987×10⁻³ kcal/(mol·K)`, default folding temperature **37 °C = 310.15 K**.
  - **Nussinov & Jacobson (1980)** *PNAS* 77(11):6309 — minimum-loop and nesting conventions.

- **Inside recursion (partition of a sub-ensemble):** with `Q_ij := Z(P_ij)` and
  `Q^b_ij := Z(P^b_ij)` (sub-ensembles requiring pair `(i,j)`):
  `Q_ij = Q_{i,j-1} + Σ_{i≤k<j-m} Q_{i,k-1} · Q^b_{kj}`, base case `Q_ij = 1` for `i ≥ j − m`.
  The **total partition function is `Z = Q_{1n}`**. Correctness follows from a disjoint (unambiguous)
  and independent decomposition.

- **Probabilities:** structure `Pr[P|S] = Z⁻¹ exp(−βE(P))`; base pair
  `Pr[(i,j)|S] = Z⁻¹ Σ_{P∋(i,j)} exp(−βE(P))`, with `β = 1/RT`.

- **Base-pair probability requires the OUTSIDE recursion (2026-06-16 correction):** the correct
  McCaskill probability is `p_kl = Q^b_kl · O_kl / Z`, where the **outside partition function** obeys
  `O_kl = Q_{1,k-1}·Q_{l+1,n} + Σ_{i<k, j>l, CanPair(i,j), j−i>m} w·Q_{i+1,k-1}·Q_{l+1,j-1}·O_ij`
  with `w = exp(−βE_bp)`. The external-only term `p^E_kl = (Q_{1,k-1}·Q^b_kl·Q_{l+1,n})/Q_{1n}` is
  **wrong for nestable pairs** (a pair that can be enclosed by another). An earlier version of this
  Evidence doc claimed the external term suffices in the flat model — that claim was FALSE and matched
  a since-fixed implementation bug.

## Oracles and datasets (all independent of the library code)

- **Structure counting (`E_bp = 0` ⇒ pairWeight = 1, so `Z` = # admissible structures, min-loop m = 3):**
  `AAAA` → **1** (no canonical pair), `GC` → **1** (only pair has j−i ≤ 3, forbidden),
  `GGGGCCCC` → **16**, `GGGAAACCC` → **20**. Derived two independent ways — a standalone re-implementation
  of the `Q`/`Q^b` recurrence and an exhaustive brute-force enumeration of all non-crossing, base-disjoint
  admissible-pair subsets — which agree on every value.
- **Exact base-pair probabilities (`E_bp = 0`, `P = #structures-containing-pair / Z`), incl. nestable pairs:**
  `GGGGCCCC` (Z=16): P(0,4)=1/16, P(0,7)=4/16, **P(1,5)=3/16**, **P(2,6)=3/16**, …;
  `GGGAAACCC` (Z=20): **P(2,6)=6/20** (external-only wrongly gives 1/20=0.05), **P(1,7)=4/20**, ….
  Bold = nestable pairs whose external-only probability is strictly too small.
- **Weighted probabilities (`E_bp = −1`, `RT = 0.61626805`, w = exp(1/RT)):** e.g. `GGGGCCCC`
  Z = 180.018344830…, P(1,5) = 0.31334323, P(0,7) = 0.45594238. Matches Boltzmann-weighted exhaustive
  enumeration and the outside recursion to machine precision (max |brute − outside| = **3.3e-16**).
- **Boltzmann structure probability (closed form):** `R = 1.987 cal/(mol·K)`, `T = 310.15 K`,
  `RT = 0.61626805 kcal/mol`; `p = 1` when a structure's energy equals the ensemble energy; for
  `E_struct = −5`, `E_ensemble = −6` → `p = exp(−1/RT) = 0.196817…`.

## Invariants

- **`Z ≥ 1` always** — the empty structure always contributes weight 1 (`Q_ij = 1` base case).
- **`Z = 1`** for a sequence with no admissible pair (`AAAA`, `GC`).
- **Every base-pair probability ∈ [0,1]** and symmetric `P[i,j] = P[j,i]` (one orientation stored).
- **Per-base pairing sum ≤ 1** — `Σ_{(i,j): i=p or j=p} P(i,j) ≤ 1` (a base pairs with ≤ 1 partner);
  verified over 300 random sequences (max single-base sum 0.983).
- **Monotonicity** — lowering `E_bp` (more favourable pairing) strictly increases `Z`.
- **Minimum loop / base case** — pairs with `j − k ≤ m` are forbidden (`Q^b = 0`); a subsequence too
  short to contain any pair has one (empty) structure of weight 1. Only Watson-Crick (A-U, G-C) and GU
  pairs contribute; all other character pairs give `Q^b = 0`.

## Deviations and assumptions

1. **Simplified per-pair energy model.** The repository computes `Z` with a fixed per-base-pair energy
   `E_bp` (the Freiburg teaching model) rather than the full **Turner 2004 nearest-neighbour loop
   energies** used by ViennaRNA and McCaskill's original tRNA examples. This is a documented
   simplification of the **energy model only** — the partition-function recurrence (`Q`, `Q^b`,
   `Z = Q₁ₙ`), the base-pair-probability formula, and all structural invariants (Z ≥ 1, probabilities
   in [0,1], normalisation, symmetry, monotonicity in `E_bp`) are fully conformant with McCaskill 1990.
   Exact Turner-parameter ensemble energies are out of scope.

**No source contradictions** — McCaskill 1990, the MIT 18.417 slides, the Freiburg tool, and ViennaRNA
agree on the recurrence, the Boltzmann probability form, and the O(n³)/O(n²) complexity. The only
recorded historical item is the corrected external-vs-outside base-pair-probability claim (fixed
2026-06-16).
