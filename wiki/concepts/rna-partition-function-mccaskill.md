---
type: concept
title: "RNA partition function (McCaskill algorithm — base-pair probabilities & ensemble free energy)"
tags: [rna, algorithm]
sources:
  - docs/Evidence/RNA-PARTITION-001-Evidence.md
source_commit: c74e2076b6e891a02c917198a54544896c4dbafa
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-partition-001-evidence
      evidence: "Test Unit ID: RNA-PARTITION-001 ... Algorithm: RNA Partition Function (McCaskill) and Boltzmann Structure Probability"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:rna-minimum-free-energy-folding
      source: rna-partition-001-evidence
      evidence: "McCaskill 1990 computes 'the full equilibrium partition function for secondary structure and the probabilities of various substructures' over the whole ensemble of probable alternative structures; ViennaRNA pf_fold p(s)=e^(−βE(s))/Z is the Boltzmann-weighted ensemble counterpart of MFE folding (a single optimal structure) — same pseudoknot-free structure space, probabilistic vs single-optimum output"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-free-energy-turner-model
      source: rna-partition-001-evidence
      evidence: "Z = Σ_P exp(−E(P)/RT) Boltzmann-weights structure free energies; ViennaRNA's default folding temperature is 37°C = 310.15 K, the Turner/NNDB standard (RT = 0.61626805 kcal/mol). The canonical McCaskill partition function consumes the Turner nearest-neighbour terms — though Seqeron's implementation uses a documented simplified fixed-per-pair E_bp energy model"
      confidence: medium
      status: current
---

# RNA partition function (McCaskill algorithm)

The **probabilistic / ensemble layer** of the RNA secondary-structure family (test unit
**RNA-PARTITION-001**): the **McCaskill (1990)** algorithm computing the equilibrium **partition
function `Z`**, the **base-pair binding probability** `P(i,j)` for every pair, and the **Boltzmann
probability of a given structure**. Where [[rna-minimum-free-energy-folding|MFE folding]] returns the
*single* lowest-energy fold, the partition function Boltzmann-weights and sums over the *whole*
ensemble of pseudoknot-free structures — it is the probabilistic counterpart of MFE. The record is
[[rna-partition-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

## 1. What the algorithm does

Given an RNA sequence, McCaskill's recursion computes both the total partition function and every
base-pair probability in **O(n³) time and O(n²) space** (after bounding interior-loop size), over all
nested (pseudoknot-free) structures. The partition function is the Boltzmann-weighted sum

```
Z = Σ_P exp(−E(P)/RT)
```

over all admissible structures `P`. When the per-pair energy `E_bp = 0`, every weight is 1 and `Z`
simply **counts** the admissible structures.

## 2. Inside recursion (`Q` / `Q^b`)

With `Q_ij := Z(P_ij)` the partition function of the sub-ensemble on `[i,j]`, and
`Q^b_ij := Z(P^b_ij)` the sub-ensemble of structures in which `(i,j)` **is** paired:

```
Q_ij = Q_{i,j-1} + Σ_{i≤k<j-m} Q_{i,k-1} · Q^b_{kj},   with Q_ij = 1 for i ≥ j − m
```

The **total partition function is `Z = Q_{1n}`**. Correctness follows from a **disjoint (unambiguous)
and independent** decomposition (MIT 18.417 slides): every structure is counted exactly once.

## 3. Base-pair probabilities need the OUTSIDE recursion

The Boltzmann probability of a structure is `Pr[P|S] = Z⁻¹ exp(−βE(P))` and of a pair is
`Pr[(i,j)|S] = Z⁻¹ Σ_{P∋(i,j)} exp(−βE(P))` (β = 1/RT). The correct McCaskill pair probability is

```
p_kl = Q^b_kl · O_kl / Z
```

where the **outside partition function** `O_kl` accumulates the ensemble *outside* the pair:

```
O_kl = Q_{1,k-1}·Q_{l+1,n} + Σ_{i<k, j>l, CanPair(i,j), j−i>m} w·Q_{i+1,k-1}·Q_{l+1,j-1}·O_ij,  w = exp(−βE_bp)
```

**Correction (crucial):** the external-only term `p^E_kl = (Q_{1,k-1}·Q^b_kl·Q_{l+1,n})/Z` is **wrong
for nestable pairs** — a pair that can be enclosed by another gets probability strictly larger than its
external term. E.g. `GGGAAACCC` P(2,6) = 6/20 = 0.30, not the external 1/20 = 0.05; `GGGGCCCC` P(1,5) =
3/16, not 1/16. The outside recursion matches Boltzmann-weighted brute-force enumeration to machine
precision (max error 3.3e-16); a prior "external term suffices" claim was a since-fixed bug.

## 4. Worked oracles

**Structure counting (`E_bp = 0`, min-loop m = 3, WC + GU pairing):**

| Sequence | `Z` (= # admissible structures) |
|----------|--------------------------------|
| `AAAA` | 1 (no canonical pair exists) |
| `GC` | 1 (only pair has j−i ≤ 3 → forbidden) |
| `GGGGCCCC` | 16 |
| `GGGAAACCC` | 20 |

Each derived two independent ways — a standalone `Q`/`Q^b` re-implementation **and** exhaustive
non-crossing base-disjoint subset enumeration.

**Base-pair probabilities** (`E_bp = 0`): `GGGGCCCC` P(0,7)=4/16, **P(1,5)=P(2,6)=3/16**; `GGGAAACCC`
**P(2,6)=6/20**, **P(1,7)=4/20**. **Weighted** (`E_bp = −1`, RT = 0.61626805): `GGGGCCCC` Z ≈
180.0183, P(1,5) = 0.31334323, P(0,7) = 0.45594238.

**Boltzmann structure probability** (closed form, R = 1.987 cal/(mol·K), T = 310.15 K, RT = 0.61626805
kcal/mol): `p = 1` when structure energy = ensemble energy; `E_struct = −5`, `E_ensemble = −6` →
`p = exp(−1/RT) = 0.196817…`.

## 5. Invariants and edge cases

- **`Z ≥ 1` always** — the empty structure always contributes weight 1 (`Q_ij = 1` base case); `Z = 1`
  for any sequence with no admissible pair.
- **Every base-pair probability ∈ [0,1]** and symmetric `P[i,j] = P[j,i]` (one orientation stored).
- **Per-base pairing sum ≤ 1** — `Σ_{(i,j): i=p or j=p} P(i,j) ≤ 1` (a base pairs with ≤ 1 partner);
  verified over 300 random sequences (max single-base sum 0.983).
- **Monotonicity** — lowering `E_bp` (more favourable pairing) strictly increases `Z`.
- **Minimum loop** — pairs with `j − k ≤ m` are forbidden (`Q^b = 0`), preventing sterically
  impossible hairpins; only Watson-Crick (A-U, G-C) and GU pairs contribute, all others `Q^b = 0`.

## 6. Scope, relationships, and assumptions

A [[scientific-rigor|research-grade]] ensemble method. It is the **Boltzmann-weighted counterpart** of
[[rna-minimum-free-energy-folding]] — same pseudoknot-free structure space, but a *probability
distribution over all structures* (partition function + pair probabilities) rather than one optimal
fold. It Boltzmann-weights the folding free energies of [[rna-free-energy-turner-model]] (the canonical
McCaskill/ViennaRNA partition function consumes the Turner nearest-neighbour terms; ViennaRNA's default
37 °C = 310.15 K is the Turner/NNDB standard). It builds on the [[rna-base-pairing]] Watson-Crick + G-U
wobble pairing rule and shares the pseudoknot-free nesting convention with the whole family.

**One documented assumption — simplified per-pair energy model.** Seqeron computes `Z` with a fixed
per-base-pair energy `E_bp` (the Freiburg teaching model) rather than the full Turner 2004
nearest-neighbour loop energies used by ViennaRNA. This is a simplification of the **energy model
only** — the recurrence (`Q`, `Q^b`, `Z = Q₁ₙ`), the base-pair-probability formula, and all structural
invariants (Z ≥ 1, probabilities in [0,1], normalisation, symmetry, monotonicity in `E_bp`) are fully
conformant with McCaskill 1990. Exact Turner-parameter ensemble energies are out of scope.

**No source contradictions** — McCaskill 1990, the MIT 18.417 slides, the Freiburg teaching tool, and
ViennaRNA agree on the recurrence, the Boltzmann probability form, and the O(n³)/O(n²) complexity. The
only historical item is the corrected external-vs-outside base-pair-probability claim (fixed 2026-06-16).
