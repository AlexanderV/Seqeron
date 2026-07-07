# Driver Mutation Detection (20/20 Rule)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-DRIVER-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Driver mutation detection separates the few mutations that "drive" tumorigenesis from the many neutral passenger mutations. This unit implements the **20/20 rule** of Vogelstein et al. (2013) [1]: a heuristic that classifies a gene as an **oncogene** when more than 20% of its mutations are missense at recurrent positions, or as a **tumor suppressor gene** when more than 20% of its mutations are inactivating (truncating) [1][2]. A mutation is then called a driver if it lies in a so-classified gene or at a caller-supplied known hotspot. The method is a deterministic heuristic over a gene's mutation spectrum, not a trained pathogenicity model; curated gene/hotspot catalogs are caller-supplied inputs.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Across cancer genomes a small number of genes are altered in many tumors ("mountains") while most are altered rarely ("hills") [1]. Oncogenes are activated by gain-of-function mutations that recur at specific amino-acid positions (hotspots); tumor suppressors are inactivated by loss-of-function (truncating) mutations dispersed along the gene [1][3]. The 20/20 rule formalizes this asymmetry into a position/consequence test.

### 2.2 Core Model

For a gene with `N` recorded coding mutations, define two fractions:

- **Recurrent-missense fraction** `f_OG = (# missense mutations at recurrent positions) / N`, where a position is *recurrent* if it carries at least 2 missense mutations [4].
- **Truncating fraction** `f_TSG = (# inactivating mutations) / N`, where inactivating = nonsense (stop gain/loss), frameshift indel, or splice donor/acceptor mutation [2][3][4].

The 20/20 rule [1][2]:

- Gene is an **oncogene** if `f_OG > 0.20`.
- Gene is a **tumor suppressor** if `f_TSG > 0.20`.

Tokheim & Karchin (2020) restate it verbatim: "OGs have >20% mutations causing missense changes at recurrent positions and TSGs have >20% mutations causing inactivating changes" [2]. The threshold comparison is strict (`>`) per the "more than 20%" wording [1][2].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The submitted mutations are a representative somatic spectrum of the gene | Skewed sampling inflates/deflates either fraction and can misclassify the gene |
| ASM-02 | Truncating mutations indicate loss of function | Passenger truncations from drift can falsely raise `f_TSG` in oncogenes [2] |
| ASM-03 | Recurrence at a position reflects positive selection | Mutational-process hotspots (not selection) can falsely raise `f_OG` |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Returned driver mutations ⊆ input mutations | `IdentifyDriverMutations` only filters the input list, never synthesizes entries |
| INV-02 | Oncogene ⟺ `f_OG > 0.20` (and, on a dual pass, `f_OG ≥ f_TSG`) | Direct encoding of [1][2] with the dominant-fraction tie-break |
| INV-03 | TumorSuppressor ⟺ `f_TSG > 0.20` (and, on a dual pass, `f_TSG > f_OG`) | Direct encoding of [1][2] |
| INV-04 | `0 ≤ ScoreDriverPotential ≤ 1` | It is `max` of two fractions of counts over `N` |
| INV-05 | A recurrent missense position requires ≥ 2 missense at the same position | Miller et al. (2017) hotspot definition [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `mutations` | `IEnumerable<GeneMutation>` | required | Coding mutations (gene, 1-based protein position, consequence) | non-null; for `ClassifyGene` should be one gene's mutations |
| `knownHotspots` | `IReadOnlySet<(string,int)>?` | null (→ empty) | Caller-supplied (gene, position) hotspot catalog | null treated as empty |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `DriverGeneClassification.Role` | `DriverGeneRole` | Oncogene / TumorSuppressor / Ambiguous |
| `DriverGeneClassification.TruncatingFraction` | `double` | `f_TSG` in [0, 1] |
| `DriverGeneClassification.RecurrentMissenseFraction` | `double` | `f_OG` in [0, 1] |
| `ScoreDriverPotential` | `double` | `max(f_TSG, f_OG)` in [0, 1] |
| `IdentifyDriverMutations` | `IReadOnlyList<GeneMutation>` | Input subset in driver genes or at known hotspots, in input order |

### 3.3 Preconditions and Validation

Null `mutations`/`knownHotspots` throw `ArgumentNullException`. Empty mutation list → `Ambiguous` with both fractions 0. Protein position is 1-based and used only for recurrence/hotspot equality. Gene symbols are matched with ordinal (case-sensitive) comparison.

## 4. Algorithm

### 4.1 High-Level Steps

1. Group mutations by gene.
2. For each gene compute `f_TSG` (truncating count / N) and `f_OG` (recurrent-missense count / N).
3. Apply the 20/20 rule: `f_OG > 0.20` → Oncogene, `f_TSG > 0.20` → TumorSuppressor, else Ambiguous; on a dual pass take the dominant fraction (exact tie → Ambiguous).
4. A mutation is a driver if its gene is Oncogene/TumorSuppressor or its (gene, position) is in the hotspot set.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Constant | Value | Source |
|----------|-------|--------|
| Driver fraction threshold (strict `>`) | 0.20 | Vogelstein 2013 [1]; Tokheim 2020 [2] |
| Recurrent-position min count | 2 | Miller 2017 [4] |
| Truncating consequences | Nonsense, Frameshift, SpliceSite | Schroeder 2014 [3]; Miller 2017 [4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ClassifyGene` | O(N) | O(P) | N mutations, P distinct missense positions |
| `IdentifyDriverMutations` | O(N) | O(G + H) | G genes, H hotspots; group + hash lookups |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ClassifyGene(mutations)`: 20/20-rule role + criterion fractions for one gene.
- `OncologyAnalyzer.ScoreDriverPotential(mutations)`: `max(f_TSG, f_OG)` driver-signal score.
- `OncologyAnalyzer.MatchCancerHotspots(mutation, knownHotspots)`: caller-supplied hotspot membership.
- `OncologyAnalyzer.IdentifyDriverMutations(mutations, knownHotspots)`: per-gene classification + hotspot match → driver subset.

### 5.2 Current Behavior

Recurrence is detected by counting missense mutations per protein position in a dictionary; positions with ≥ 2 missense contribute their full count to `f_OG`. Hotspot catalogs are not embedded — `MatchCancerHotspots`/`IdentifyDriverMutations` consult a caller-supplied set so no unverifiable curated list is hardcoded. **Search reuse:** the repository suffix tree was evaluated and is not used — there is no substring/sequence search here; the "lookup" is set membership on (gene, position) keys (O(1) hash), for which a suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Oncogene if `f_OG > 0.20` (recurrent missense) and TumorSuppressor if `f_TSG > 0.20` (truncating), strict threshold [1][2].
- Truncating = nonsense + frameshift + splice site [3][4]; recurrent = same position ≥ 2 times [4].

**Intentionally simplified:**

- `ScoreDriverPotential`: returns the 20/20 driver-signal fraction `max(f_TSG, f_OG)` instead of a trained CADD/SIFT/PolyPhen score; **consequence:** the score is the rule's own evidence strength, not an external pathogenicity probability.
- Dual-criterion genes: resolved by the dominant fraction (tie → Ambiguous); **consequence:** the source does not prescribe a single label for a gene passing both, so an atypical dual pass may differ from a curated catalog.

**Not implemented:**

- Curated driver-gene lists and hotspot catalogs (COSMIC, OncoKB, Cancer Hotspots, ClinVar); **users should rely on:** supplying `knownHotspots` (and pre-filtered mutations) from those databases as inputs.
- Statistical driver tests (MutSigCV, 20/20+ random forest, OncodriveFML); **users should rely on:** dedicated external tools — the 20/20 rule is a deterministic heuristic only.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Strict `>` vs `≥` 0.20 | Assumption | Boundary genes at exactly 20% | accepted | Vogelstein/Tokheim "more than 20%" [1][2]; OncodriveROLE writes "≥" [3] — see ASM/Assumption #2 |
| 2 | Dual-pass tie-break by dominant fraction | Assumption | Atypical genes passing both criteria | accepted | Source silent; deterministic resolution |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty mutation list | Ambiguous, fractions 0 | No evidence to classify |
| Truncating fraction exactly 0.20 | Not a TumorSuppressor | Strict `>` per "more than 20%" [1] |
| Single missense at a position | Not recurrent (`f_OG` unaffected) | Recurrence needs ≥ 2 [4] |
| Null inputs | `ArgumentNullException` | Contract |

### 6.2 Limitations

The 20/20 rule misses lowly recurrent drivers [3] and can be misled by drift-induced truncations in oncogenes [2]; it requires a reasonable number of mutations to be stable. It is a classification heuristic, not a significance test, and does not model copy-number, fusion, or non-coding drivers.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

// IDH1: ten missense, all at codon 132 → recurrent-missense fraction 1.0 → Oncogene.
var idh1 = Enumerable.Range(0, 10)
    .Select(_ => new GeneMutation("IDH1", 132, MutationConsequence.Missense));
DriverGeneClassification c = ClassifyGene(idh1);
// c.Role == DriverGeneRole.Oncogene; c.RecurrentMissenseFraction == 1.0
```

**Numerical walk-through:** IDH1 with 10 missense at codon 132: position 132 carries 10 ≥ 2 missense, so recurrent-missense count = 10, `f_OG = 10/10 = 1.0 > 0.20` → Oncogene; `f_TSG = 0/10 = 0` [1][4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_IdentifyDriverMutations_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_IdentifyDriverMutations_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ONCO-DRIVER-001-Evidence.md](../../../docs/Evidence/ONCO-DRIVER-001-Evidence.md)

## 8. References

1. Vogelstein B, Papadopoulos N, Velculescu VE, Zhou S, Diaz LA Jr, Kinzler KW. 2013. Cancer Genome Landscapes. Science 339(6127):1546–1558. https://doi.org/10.1126/science.1235122
2. Tokheim C, Karchin R. 2020. Somatic selection distinguishes oncogenes and tumor suppressor genes. Bioinformatics 36(6):1712–1719. https://doi.org/10.1093/bioinformatics/btz759
3. Schroeder MP, Rubio-Perez C, Tamborero D, Gonzalez-Perez A, Lopez-Bigas N. 2014. OncodriveROLE classifies cancer driver genes in loss of function and activating mode of action. Bioinformatics 30(17):i549–i555. https://doi.org/10.1093/bioinformatics/btu467
4. Miller ML, Reznik E, Gauthier NP, et al. 2017. Identification and analysis of mutational hotspots in oncogenes and tumour suppressors. Oncotarget 8(20):33321–33333. https://doi.org/10.18632/oncotarget.15514
