# Primer3 Weighted Penalty Objective (Per-Primer)

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PRIMER-TM-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-06-23 |

## 1. Overview

The Primer3 penalty objective scores a single PCR primer by how far its measured
properties deviate from user-defined optima, as a weighted sum of one-sided deviations.
It is the de-facto field-standard primer-scoring objective: Primer3 evaluates every
candidate primer with this function and selects the candidate with the **lowest** value
[1][3][4]. This implementation faithfully reproduces the left/right-primer branch of
Primer3's reference `p_obj_fn` with its documented default weights and optima, so the
score is no longer a heuristic without a ground truth — it computes the Primer3 objective
and is cross-validated against the published formula. It is exact (a deterministic
weighted sum), not heuristic or probabilistic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In PCR primer design the "best" primer is operationally defined as the one minimizing a
penalty function over primer properties (melting temperature, length, GC content,
self-complementarity, ambiguous bases) [1]. Primer3 is the canonical open-source tool
implementing this approach [1][2].

### 2.2 Core Model

For a left or right primer, the penalty is the weighted sum of one-sided deviations from
the optima (Primer3 `p_obj_fn`, manual §19) [3][4]:

```
penalty =
    (Tm  > OPT_TM)  ? WT_TM_GT  * (Tm  - OPT_TM)  : 0
  + (Tm  < OPT_TM)  ? WT_TM_LT  * (OPT_TM  - Tm)  : 0
  + (GC% > OPT_GC)  ? WT_GC_GT  * (GC% - OPT_GC)  : 0
  + (GC% < OPT_GC)  ? WT_GC_LT  * (OPT_GC - GC%)  : 0
  + (len > OPT_SIZE)? WT_SIZE_GT* (len - OPT_SIZE): 0
  + (len < OPT_SIZE)? WT_SIZE_LT* (OPT_SIZE - len): 0
  + WT_SELF_ANY * SELF_ANY
  + WT_SELF_END * SELF_END
  + WT_NUM_NS   * (number of N bases)
```

Each term is added only when its weight is non-zero and the deviation has the matching
sign; a parameter exactly at its optimum contributes nothing [3]. GC content is a
percentage in [0, 100] (`gc_content = 100·num_gc/num_gcat` in Primer3 source) [3].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | `SELF_ANY` / `SELF_END` are local-alignment self-complementarity scores supplied by the caller (Primer3 computes them with its `dpal` engine, scoring +1.00 complementary / −1.00 mismatch / −2.00 gap [4]). | If a caller passes a value on a different scale, the self-complementarity term is misweighted. Under the default weights (0) these terms are inert. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | penalty ≥ 0 for all inputs | every term is weight·(non-negative deviation) [3] |
| INV-02 | penalty = 0 ⇔ Tm=OPT_TM, len=OPT_SIZE, GC=OPT_GC and SELF_ANY=SELF_END=N=0 | sign-gated terms; all contributions vanish at the optimum [3] |
| INV-03 | a parameter at its optimum contributes 0 to its term | strict `>` / `<` gates exclude the equality case [3] |
| INV-04 | each term scales linearly with its weight | term = weight·deviation [3][4] |
| INV-05 | default weights TM/SIZE = 1, GC/SELF/NUM_NS = 0; optima OPT_TM=60 °C, OPT_SIZE=20 bases, OPT_GC=50 % | Primer3 `pr_set_default_global_args_2` (TM/SIZE/GC) and manual (OPT_GC 50.0) [3][4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `inputs.Tm` | double | required | Primer melting temperature | °C |
| `inputs.Length` | int | required | Primer length | bases |
| `inputs.GcPercent` | double | required | GC content | percent, [0, 100] |
| `inputs.SelfAny` | double | 0 | Self-complementarity local-alignment score (PRIMER_SELF_ANY) | ≥ 0 |
| `inputs.SelfEnd` | double | 0 | 3'-self-complementarity score (PRIMER_SELF_END) | ≥ 0 |
| `inputs.NumNs` | int | 0 | Number of N bases | ≥ 0 |
| `weights` | Primer3PenaltyWeights? | `DefaultPrimer3Weights` | `PRIMER_WT_*` weights | one-sided _gt/_lt |
| `optima` | Primer3Optima? | `DefaultPrimer3Optima` | `PRIMER_OPT_*` optima | OPT_TM/SIZE/GC |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | double | Primer3 objective-function value (penalty); lower is better; 0 = all terms at optimum |

### 3.3 Preconditions and Validation

The method is a pure weighted sum and does not throw for in-range numeric inputs. GC is
interpreted as a percentage (0–100), consistent with Primer3 [3]; passing a fraction
(0–1) would silently mis-scale the GC term. Length is in bases; Tm in °C. No sequence is
parsed here — the caller supplies measured properties (the repository computes Tm via the
SEQ-THERMO-001-validated routines and GC via `CalculateGcContent`).

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve weights and optima (defaults if not supplied).
2. Add the Tm term: `WT_TM_GT·(Tm−OPT_TM)` if `Tm>OPT_TM`, else `WT_TM_LT·(OPT_TM−Tm)` if `Tm<OPT_TM`.
3. Add the GC% term symmetrically using `WT_GC_GT` / `WT_GC_LT` about `OPT_GC`.
4. Add the size term symmetrically using `WT_SIZE_GT` / `WT_SIZE_LT` about `OPT_SIZE`.
5. Add `WT_SELF_ANY·SELF_ANY`, `WT_SELF_END·SELF_END`, `WT_NUM_NS·N`.
6. Return the sum.

### 4.2 Decision Rules, Scoring, Reference Tables

Default weights and optima (Primer3 source / manual) [3][4]:

| Parameter | Default | Source token |
|-----------|---------|--------------|
| WT_TM_GT, WT_TM_LT | 1.0 | `weights.temp_gt`, `weights.temp_lt` |
| WT_SIZE_GT, WT_SIZE_LT | 1.0 | `weights.length_gt`, `weights.length_lt` |
| WT_GC_PERCENT_GT, WT_GC_PERCENT_LT | 0.0 | `weights.gc_content_gt`, `weights.gc_content_lt` |
| WT_SELF_ANY, WT_SELF_END | 0.0 | `weights.compl_any`, `weights.compl_end` |
| WT_NUM_NS | 0.0 | `weights.num_ns` |
| OPT_TM | 60.0 °C | `opt_tm` |
| OPT_SIZE | 20 bases | `opt_size` |
| OPT_GC_PERCENT | 50.0 % | manual `PRIMER_OPT_GC_PERCENT` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculatePrimer3Penalty` | O(1) | O(1) | fixed number of arithmetic terms |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.CalculatePrimer3Penalty(Primer3PenaltyInputs, Primer3PenaltyWeights?, Primer3Optima?)`: computes the per-primer penalty (this objective).
- `PrimerDesigner.DefaultPrimer3Weights`: sourced default weights.
- `PrimerDesigner.DefaultPrimer3Optima`: sourced default optima.

### 5.2 Current Behavior

Reproduces the left/right-primer branch of Primer3 `p_obj_fn` for the core terms
(Tm, size, GC%, self_any, self_end, num_ns). The legacy convenience `Score`
(`CalculatePrimerScore`, used by `EvaluatePrimer`/`DesignPrimers`) is **kept unchanged**
and available for backward compatibility; the new method is the validated, Primer3-anchored
objective. No search/matching is involved, so the repository suffix tree is N/A here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- One-sided Tm term with separate `WT_TM_GT` / `WT_TM_LT` weights about `OPT_TM` [3][4].
- One-sided GC% term (`WT_GC_PERCENT_GT/LT` about `OPT_GC_PERCENT`), GC as percent [3][4].
- One-sided size term (`WT_SIZE_GT/LT` about `OPT_SIZE`) [3][4].
- Linear self_any, self_end and num_ns terms [3][4].
- Default weights and optima taken verbatim from Primer3 source / manual [3][4].

**Intentionally simplified:**

- self_any/self_end alignment scores are caller-supplied: the penalty arithmetic on them is
  exact, but Primer3's `dpal` local-alignment computation of those scores is not reproduced;
  **consequence:** with the default weights (0) the term is inert; with non-zero weights the
  caller must supply a Primer3-scale alignment score (+1.00/−1.00/−2.00 [4]).

**Not implemented:**

- The thermodynamic-alignment penalty branch (`*_TH` terms, `temp_cutoff`), `pos_penalty`,
  `end_stability`, `seq_quality`, `repeat_sim`, `template_mispriming`; **users should rely on:**
  Primer3 itself for those terms (all default to weight 0, so they do not affect the default objective).
- The pair-level objective (`PRIMER_PAIR_*`, Tm-difference, product size); **users should rely on:**
  a future pair-penalty unit or Primer3.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | self_any/self_end scores caller-supplied | Assumption | only when WT_SELF_* ≠ 0 (default 0) | accepted | ASM-01 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| All parameters at optimum | penalty = 0 | INV-02 |
| Parameter exactly at optimum | that term contributes 0 | INV-03 (strict gates) |
| Weight = 0 for a term | term never contributes | weight-gated [3] |
| Tm above vs below optimum | uses `WT_TM_GT` vs `WT_TM_LT` respectively | one-sided weights [3][4] |

### 6.2 Limitations

Per-primer only (no pair penalty); the `*_TH` thermodynamic-alignment, position, end-stability,
sequence-quality, repeat and template-mispriming terms are not implemented (they default to
weight 0 in Primer3, so the default objective is unaffected). self_any/self_end alignment scores
are caller-supplied (§5.3).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Default weights/optima: penalty = |Tm-60| + |len-20|.
var p = PrimerDesigner.CalculatePrimer3Penalty(
    new Primer3PenaltyInputs(Tm: 62.5, Length: 22, GcPercent: 50));
// 1.0*(62.5-60) + 1.0*(22-20) = 4.5
```

**Numerical walk-through:** For Tm=57, len=18, GC=50, default weights: Tm term =
`WT_TM_LT·(60−57)=3`, size term = `WT_SIZE_LT·(20−18)=2`, GC term = 0 (weight 0), total = 5.0.

### 7.2 Applications and Use Cases

- **PCR primer ranking:** score and sort candidate primers by their Primer3 objective value,
  matching Primer3's selection of the lowest-penalty primer [1][3].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PrimerDesigner_Primer3Penalty_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_Primer3Penalty_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [PRIMER-TM-001-Evidence.md](../../../docs/Evidence/PRIMER-TM-001-Evidence.md)
- TestSpec: [PRIMER-TM-001-Penalty.md](../../../tests/TestSpecs/PRIMER-TM-001-Penalty.md)
- Related algorithms: [Melting_Temperature](../Molecular_Tools/Melting_Temperature.md), [Primer_Design](./Primer_Design.md)

## 8. References

1. Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. 2012. Primer3—new capabilities and interfaces. Nucleic Acids Research 40(15):e115. https://doi.org/10.1093/nar/gks596
2. Koressaar T, Remm M. 2007. Enhancements and modifications of primer design program Primer3. Bioinformatics 23(10):1289–1291. https://doi.org/10.1093/bioinformatics/btm091
3. Primer3 source code, `src/libprimer3.cc` (`p_obj_fn`, `pr_set_default_global_args_2`) and `src/libprimer3.h` (`DEFAULT_OPT_GC_PERCENT`), branch main. https://github.com/primer3-org/primer3
4. Primer3 manual, §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE" and global input tags. https://primer3.org/manual.html
