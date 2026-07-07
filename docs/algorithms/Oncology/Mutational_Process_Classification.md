# Mutational Process Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Mutational Signatures |
| Test Unit ID | ONCO-SIG-004 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Given the per-signature exposures (activities) returned by signature fitting, this algorithm classifies which
**mutational processes** (mutagenic aetiologies) are active in a tumour. Each COSMIC SBS signature has a
proposed aetiology (e.g. SBS4 = tobacco smoking, SBS2/SBS13 = APOBEC) [3]. Exposures are converted to
normalized relative contributions; signatures contributing below a fixed cutoff are discarded; the survivors
are mapped to their processes and aggregated, and the process with the largest aggregated contribution is the
dominant process. It is a deterministic, threshold-driven classification built on the deconstructSigs
presence rule [1] and the COSMIC SBS→aetiology map [3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Mutational signatures are characteristic patterns of single-base substitutions imprinted on a genome by
distinct mutagenic processes. The COSMIC catalogue enumerates SBS signatures and assigns each a proposed
aetiology [3][4]. Signature *fitting* (ONCO-SIG-002) decomposes an observed catalogue into a non-negative
combination of reference signatures, yielding an exposure per signature. Translating those exposures into a
statement of which biological processes are active requires (a) a presence rule that separates real activity
from fitting noise, and (b) the signature→aetiology map.

### 2.2 Core Model

Let exposures be `{(sᵢ, eᵢ)}` with `eᵢ ≥ 0`. Define the total `T = Σ eᵢ` and the normalized relative
contribution `wᵢ = eᵢ / T` (deconstructSigs: "the weights W are normalized between 0 and 1") [1]. A signature
is **present/active** iff `wᵢ ≥ τ`, where `τ = 0.06` is the deconstructSigs cutoff — "any signature with
Wᵢ < 6% is excluded" [1]; the reference code applies `weights[weights < signature.cutoff] <- 0` with default
`signature.cutoff = 0.06` [2]. The comparison is strict less-than, so `wᵢ = 0.06` is retained [2].

Each surviving signature is mapped to a process `p(sᵢ)` via the COSMIC aetiology map [3]. The aggregated
contribution of a process P is `C(P) = Σ_{i: p(sᵢ)=P, wᵢ≥τ} wᵢ` (additive weights) [1]. The **active set** is
`{P : C(P) > 0}` and the **dominant process** is `argmax_P C(P)`.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A process's contribution is the **sum** of its member signatures' surviving contributions (additive weights) | Non-additive aggregation would change per-process totals and the dominant call |
| ASM-02 | The 6% cutoff is applied per signature, then aggregated by process (matching deconstructSigs, which operates per signature) | Applying the cutoff to per-process totals would change which signatures survive |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each surviving `wᵢ ∈ [0,1]`; Σ surviving contributions ≤ 1 | `wᵢ = eᵢ/T` with `eᵢ ≥ 0, T = Σeᵢ`; sub-cutoff dropped [1] |
| INV-02 | A signature is active iff `wᵢ ≥ 0.06` (strict `<` excludes; 0.06 retained) | deconstructSigs `weights < signature.cutoff`, cutoff 0.06 [1][2] |
| INV-03 | `C(P) = Σ` surviving member contributions | additive weights [1] + ASM-01 |
| INV-04 | Dominant = active process with max `C(P)`; `Unknown` when none active | derivation |
| INV-05 | `T = 0` ⇒ no active processes, dominant = `Unknown` | normalization undefined for zero total |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `exposures` | `IReadOnlyList<(string Signature, double Exposure)>` | required | Per-signature COSMIC label + non-negative activity | non-null; labels non-null; exposures ≥ 0 |
| `contributionCutoff` | `double` | `0.06` | Minimum normalized contribution to be active | `[0, 1)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ActiveProcesses` | `IReadOnlyList<ProcessActivity>` | Active processes (descending contribution, then process enum); each carries the summed contribution ∈ [0,1] |
| `DominantProcess` | `MutationalProcess` | Active process with largest contribution, or `Unknown` when none |

### 3.3 Preconditions and Validation

`exposures` null → `ArgumentNullException`; a null label → `ArgumentNullException`; a negative or NaN exposure
→ `ArgumentException`; `contributionCutoff` NaN or outside `[0, 1)` → `ArgumentOutOfRangeException`. An empty
list, or one whose total is 0, yields an empty active set and `Unknown` dominant. Signature-label matching is
case-insensitive; labels outside the COSMIC map resolve to `Unknown` and contribute to no recognized process.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; compute total `T = Σ eᵢ`.
2. If `T ≤ 0`, return empty active set and `Unknown` dominant (INV-05).
3. For each signature: `wᵢ = eᵢ / T`; skip if `wᵢ < τ` (cutoff) or if its process is `Unknown`.
4. Accumulate surviving `wᵢ` into a per-process sum.
5. Order active processes by descending contribution (then by process enum); dominant = first.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

**Presence cutoff:** `τ = 0.06` (deconstructSigs `signature.cutoff`) [1][2].

**COSMIC SBS → process map** (proposed aetiologies, verbatim) [3]:

| Process | SBS signatures | COSMIC aetiology string |
|---------|----------------|--------------------------|
| Aging (clock-like) | SBS1, SBS5 | "Spontaneous deamination of 5-methylcytosine (clock-like)"; "Unknown (clock-like)" |
| APOBEC | SBS2, SBS13 | "Activity of APOBEC family of cytidine deaminases" |
| Tobacco smoking | SBS4 | "Tobacco smoking" |
| Ultraviolet light | SBS7a, SBS7b, SBS7c, SBS7d | "Ultraviolet light exposure" |
| MMR deficiency | SBS6, SBS15, SBS20, SBS26 | "Defective DNA mismatch repair"; SBS20 "Concurrent POLD1 mutations and defective DNA mismatch repair" |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ClassifyMutationalProcess` | O(k log k) | O(k) | k signatures; the log k is the final ordering of ≤ 5 processes (effectively O(k)) |
| `GetMutationalProcess` | O(1) | O(1) | dictionary lookup |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ClassifyMutationalProcess(exposures, contributionCutoff)`: normalizes exposures, applies the 6% cutoff, aggregates surviving contributions by COSMIC process, returns active set + dominant.
- `OncologyAnalyzer.GetMutationalProcess(signatureLabel)`: case-insensitive COSMIC SBS label → `MutationalProcess`.

### 5.2 Current Behavior

The cutoff is applied to per-signature normalized contributions (matching deconstructSigs, which operates per
signature), then survivors are grouped by process. Surviving contributions can sum to less than 1 because
sub-cutoff mass is dropped (deconstructSigs attributes the remainder to "unknown") [1]. Signatures mapping to
`Unknown` (unmapped aetiology) are excluded from the active set. This is not a search/matching unit, so the
repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Normalized relative contributions `wᵢ = eᵢ / Σe` [1].
- Presence cutoff `τ = 0.06` with strict `<` (a 0.06 contribution is retained) [1][2].
- COSMIC SBS→aetiology map for Aging, APOBEC, Tobacco, UV, MMR deficiency [3].

**Intentionally simplified:**

- Only the five aetiologies named for ONCO-SIG-004 are mapped; other COSMIC signatures resolve to `Unknown`; **consequence:** processes such as platinum chemotherapy or AID are not reported. Extend `SbsToProcess` to add them.

**Not implemented:**

- Confidence-based presence (bootstrap-CI lower bound > 0); **users should rely on:** `OncologyAnalyzer.BootstrapExposures` (ONCO-SIG-003) for interval estimates, combined with this cutoff rule.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| All-zero exposures | empty active set, `Unknown` dominant | T = 0, normalization undefined (INV-05) |
| Empty exposure list | empty active set, `Unknown` dominant | T = 0 |
| Contribution exactly 0.06 | retained | strict `<` cutoff [2] |
| Contribution just below 0.06 | excluded | `weights < signature.cutoff` [2] |
| Unmapped SBS label | contributes to no process | COSMIC unknown aetiology |
| Negative / NaN exposure | `ArgumentException` | invalid input |

### 6.2 Limitations

Aetiology assignments follow the caller-independent COSMIC catalogue; the method does not validate that the
supplied labels correspond to genuine COSMIC profiles (profiles are caller-supplied during fitting). The 6%
cutoff is the deconstructSigs default; thresholds tuned for other panels can be passed explicitly.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var exposures = new (string, double)[]
{
    ("SBS2", 50), ("SBS13", 30), ("SBS1", 15), ("SBS4", 5)
};
var result = OncologyAnalyzer.ClassifyMutationalProcess(exposures);
// normalized: SBS2 0.50, SBS13 0.30, SBS1 0.15, SBS4 0.05 (< 0.06 → dropped)
// APOBEC = 0.80, Aging = 0.15; DominantProcess == MutationalProcess.Apobec
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ClassifyMutationalProcess_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ClassifyMutationalProcess_Tests.cs) — covers INV-01..INV-05
- Evidence: [ONCO-SIG-004-Evidence.md](../../../docs/Evidence/ONCO-SIG-004-Evidence.md)
- Related algorithms: [Mutational_Signature_Fitting](Mutational_Signature_Fitting.md), [Mutational_Signature_Exposure_Bootstrap](Mutational_Signature_Exposure_Bootstrap.md)

## 8. References

1. Rosenthal R, McGranahan N, Herrero J, Taylor BS, Swanton C. 2016. deconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. *Genome Biology* 17:31. https://doi.org/10.1186/s13059-016-0893-4
2. deconstructSigs reference implementation, `whichSignatures.R` (`signature.cutoff = 0.06`; `weights[weights < signature.cutoff] <- 0`). https://github.com/raerose01/deconstructSigs/blob/master/R/whichSignatures.R
3. COSMIC Mutational Signatures — SBS (proposed aetiologies). Wellcome Sanger Institute. https://cancer.sanger.ac.uk/signatures/sbs/
4. Alexandrov LB, Kim J, Haradhvala NJ, et al. 2020. The repertoire of mutational signatures in human cancer. *Nature* 578(7793):94–101. https://doi.org/10.1038/s41586-020-1943-3
