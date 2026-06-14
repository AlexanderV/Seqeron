# Disorder Propensity (TOP-IDP scale & Dunker amino-acid classification)

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinPred |
| Test Unit ID | DISORDER-PROPENSITY-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

This unit exposes the per-amino-acid building blocks used by the intrinsic-disorder predictor: a single-residue **disorder propensity lookup** based on the TOP-IDP scale [1] and a **categorical order/disorder classification** of the 20 standard amino acids based on Dunker et al. (2001) [2]. `GetDisorderPropensity(char)` returns the TOP-IDP value of one residue; `IsDisorderPromoting(char)` answers whether a residue belongs to the disorder-promoting class; and `DisorderPromotingAminoAcids` / `OrderPromotingAminoAcids` / `AmbiguousAminoAcids` expose the three classification sets. All operations are exact table lookups (O(1)); no statistics or windowing are involved.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Intrinsically disordered proteins (IDPs) lack a stable folded structure under physiological conditions. Their amino-acid composition is biased: disordered regions are enriched in hydrophilic/charged residues and depleted in bulky hydrophobic residues [2]. Two complementary descriptions of this bias are used here — a continuous per-residue scale (TOP-IDP) and a discrete three-way classification (Dunker).

### 2.2 Core Model

**TOP-IDP scale.** Campen et al. (2008) optimised, over 517 surveyed amino-acid scales, a per-residue scale that maximally separates ordered from disordered residues. The published Table 2 values (the lookup returned by `GetDisorderPropensity`) are [1]:

W = -0.884, F = -0.697, Y = -0.510, I = -0.486, M = -0.397, L = -0.326, V = -0.121, N = 0.007, C = 0.020, T = 0.059, A = 0.060, G = 0.166, R = 0.180, D = 0.192, H = 0.303, Q = 0.318, S = 0.341, K = 0.586, E = 0.736, P = 0.987.

Negative values are order-promoting, positive values disorder-promoting; W is the global minimum and P the global maximum [1].

**Dunker classification.** Dunker et al. (2001) partition the 20 standard residues into [2]:

- disorder-promoting: {A, R, G, Q, S, P, E, K} (8);
- order-promoting: {W, C, F, I, Y, V, L, N} (8);
- ambiguous: {H, M, T, D} (4).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `GetDisorderPropensity` returns the exact Table 2 value for each of the 20 standard residues | Lookup table copied verbatim from Campen et al. (2008) Table 2 [1] |
| INV-02 | -0.884 ≤ propensity ≤ 0.987 over the 20 residues; min at W, max at P | Scale extrema in Table 2 [1] |
| INV-03 | `IsDisorderPromoting(c)` ⇔ c ∈ `DisorderPromotingAminoAcids` | Predicate and property read the same backing set [2] |
| INV-04 | The disorder (8) / order (8) / ambiguous (4) sets are pairwise disjoint and cover all 20 standard residues | Direct partition from Dunker et al. (2001) [2] |
| INV-05 | `GetDisorderPropensity` and `IsDisorderPromoting` are case-insensitive | Input is upper-cased before lookup |

### 2.5 Comparison with Related Methods

| Aspect | TOP-IDP scale | Dunker classification |
|--------|---------------|------------------------|
| Output type | continuous real value | discrete class (3-way) |
| Granularity | per residue, ordered/ranked | set membership |
| Primary use | smoothed disorder score in a sliding window | quick compositional check |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| aminoAcid | char | required | Single one-letter amino-acid code | Case-insensitive; 20 standard codes defined |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `GetDisorderPropensity` | double | TOP-IDP value; 0.0 for residues outside the 20-residue scale |
| `IsDisorderPromoting` | bool | true iff the residue is in the disorder-promoting set |
| `DisorderPromotingAminoAcids` | IReadOnlyList&lt;char&gt; | {A, E, G, K, P, Q, R, S} (sorted) |
| `OrderPromotingAminoAcids` | IReadOnlyList&lt;char&gt; | {C, F, I, L, N, V, W, Y} (sorted) |
| `AmbiguousAminoAcids` | IReadOnlyList&lt;char&gt; | {D, H, M, T} (sorted) |

### 3.3 Preconditions and Validation

Input is upper-cased (lowercase 'p' equals 'P'). Any character not among the 20 standard residues returns propensity 0.0 and `IsDisorderPromoting` false — see §6.1. No exceptions are thrown for unknown characters.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case the input residue.
2. For `GetDisorderPropensity`, look up the residue in the TOP-IDP table; return its value or 0.0 if absent.
3. For `IsDisorderPromoting`, test membership in the disorder-promoting set.
4. The properties return pre-sorted cached copies of the three Dunker sets.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The two reference tables — the 20-residue TOP-IDP scale [1] and the three Dunker classification sets [2] — fully determine behavior. Both are stored as static dictionaries / hash sets initialised verbatim from the cited sources; the property lists are sorted once and cached.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GetDisorderPropensity` | O(1) | O(1) | single dictionary lookup |
| `IsDisorderPromoting` | O(1) | O(1) | single hash-set lookup |
| classification properties | O(1) | O(1) | return cached sorted list |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)

- `DisorderPredictor.GetDisorderPropensity(char)`: TOP-IDP value lookup, default 0.0.
- `DisorderPredictor.IsDisorderPromoting(char)`: membership in the disorder-promoting set.
- `DisorderPredictor.DisorderPromotingAminoAcids` / `OrderPromotingAminoAcids` / `AmbiguousAminoAcids`: cached sorted classification sets.

### 5.2 Current Behavior

The TOP-IDP scale is held in a `Dictionary<char,double>` and the three Dunker classes in `HashSet<char>`; the public list properties expose `OrderBy`-sorted cached copies. No search/matching is performed, so the repository suffix tree is not applicable (N/A) — these are O(1) table lookups.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All 20 TOP-IDP per-residue values, exactly as Campen et al. (2008) Table 2 [1].
- The disorder/order/ambiguous classification, exactly as Dunker et al. (2001) [2].

**Intentionally simplified:**

- (none)

**Not implemented:**

- The TOP-IDP windowed prediction / 0.542 cutoff is out of scope here; **users should rely on:** [Disorder_Prediction.md](Disorder_Prediction.md) (`PredictDisorder`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Unknown residue → 0.0 | Assumption | Out-of-scale chars score as neutral, not undefined | accepted | Implementation `GetValueOrDefault`; not source-defined |
| 2 | S/K rank-string order | Assumption | None on the scope methods | accepted | Table 2 numeric values are authoritative over the rendered rank string |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Unknown residue (X, Z, B, '*') | propensity 0.0, `IsDisorderPromoting` false | Scale undefined outside 20 residues; lookup default [1] |
| Lowercase input ('p') | same as uppercase ('P') | Input upper-cased (INV-05) |
| W / P | -0.884 / 0.987 (scale extrema) | Anchor residues [1] |

### 6.2 Limitations

Defined only for the 20 standard amino acids. Single-residue propensity in isolation is a weak disorder predictor; the windowed `PredictDisorder` aggregates it. The Dunker classification is a coarse 3-way partition, not a probability.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double pPro = DisorderPredictor.GetDisorderPropensity('P'); // 0.987 (max)
double pTrp = DisorderPredictor.GetDisorderPropensity('W'); // -0.884 (min)
bool eDisorder = DisorderPredictor.IsDisorderPromoting('E'); // true
bool wDisorder = DisorderPredictor.IsDisorderPromoting('W'); // false
var promoting = DisorderPredictor.DisorderPromotingAminoAcids; // [A,E,G,K,P,Q,R,S]
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [DisorderPredictor_GetDisorderPropensity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_GetDisorderPropensity_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [DISORDER-PROPENSITY-001-Evidence.md](../../../docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md)
- Related algorithms: [Disorder_Prediction.md](Disorder_Prediction.md)

## 8. References

1. Campen A, Williams RM, Brown CJ, Meng J, Uversky VN, Dunker AK. 2008. TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder. Protein Pept Lett 15(9):956-963. https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/
2. Dunker AK, Lawson JD, Brown CJ, et al. 2001. Intrinsically disordered protein. J Mol Graph Model 19(1):26-59. https://pubmed.ncbi.nlm.nih.gov/11381529/
