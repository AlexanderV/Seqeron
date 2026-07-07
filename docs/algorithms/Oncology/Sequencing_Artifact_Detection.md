# Sequencing Artifact Detection (OxoG / FFPE Deamination / Strand Bias)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-ARTIFACT-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Detects and filters sequencing artifacts that masquerade as low-frequency somatic variants. Two substitution-class artifacts are recognized: FFPE cytosine-deamination (C>T / G>A) and OxoG 8-oxoguanine oxidation (G>T / C>A). Oxidation is confirmed by the read-orientation imbalance captured by the GIV (Global Imbalance Value) score, and a Phred-scaled Fisher strand-bias score (FS) is reported per variant. The algorithm is a deterministic, rule-based classifier driven by published substitution signatures and thresholds — not a probabilistic damage model. It is used to clean somatic call sets before downstream interpretation [1][4][6].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Library preparation and tissue fixation introduce base-level damage that is read out as recurrent low-VAF substitutions. 8-oxoguanine pairs with adenine, so an oxidized guanine on the template produces an apparent G>T on read 1 (and its reverse complement C>A on read 2) [1]. Formalin fixation deaminates cytosine to uracil, which pairs with adenine, producing C>T (and G>A on the antisense strand) [4]. Strand bias — an allele appearing predominantly on one strand — is a generic artifact indicator [6].

### 2.2 Core Model

**Substitution classes** (disjoint) [1][4]:

- FFPE deamination: ref→alt ∈ {C>T, G>A}.
- OxoG oxidation: ref→alt ∈ {G>T, C>A}.

**GIV (Global Imbalance Value)** for a substitution type [1][2][3]:

```
GIV = (alt-supporting reads in read 1) / (alt-supporting reads in read 2)
```

GIV = 1 for an undamaged, balanced library; GIV > 1.5 is defined as damaged [3]. For the canonical OxoG substitution G>T, an excess in read 1 over read 2 indicates oxidative damage [1].

**FisherStrand (FS)** [6]: build the 2×2 strand contingency table with cell order [refForward, refReverse, altForward, altReverse], compute the two-sided Fisher exact p-value, and Phred-scale it:

```
FS = -10 · log10( max(p, MIN_PVALUE) )      MIN_PVALUE = 1E-320
```

FS = 0 when there is no strand bias (p = 1); FS grows as the alleles segregate by strand [6].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `FilterArtifacts` output ⊆ input, preserving order | output is built by skipping flagged inputs in order |
| INV-02 | GIV ≥ 0; GIV = 1 for balanced read1 = read2 counts | ratio of non-negative counts; equal counts ⇒ 1 [1][3] |
| INV-03 | FS ≥ 0; FS = 0 for a balanced table (p = 1) | −10·log10(1) = 0; p ≤ 1 ⇒ FS ≥ 0 [6] |
| INV-04 | FFPE class {C>T, G>A} and OxoG class {G>T, C>A} are disjoint | the four ordered pairs are distinct [1][4] |
| INV-05 | Greater strand segregation ⇒ FS non-decreasing | Fisher p decreases as the table becomes more extreme [6] |

### 2.5 Comparison with Related Methods

| Aspect | This (substitution-class + GIV + FS) | Probabilistic orientation-bias model (GATK F1R2 / OrientationBiasFilter) |
|--------|--------------------------------------|--------------------------------------------------------------------------|
| Model | rule-based substitution class + ratio threshold | Bayesian artifact-vs-real posterior per read-orientation |
| Inputs | per-strand and per-read-mate counts | full F1R2/F2R1 read tensors and priors |
| Output | discrete artifact class + GIV + FS | per-variant artifact probability |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| observation | `ArtifactObservation` | required | ref/alt bases + per-strand and per-read-mate counts | bases A/C/G/T; counts ≥ 0 |
| r1Count, r2Count | `int` | required | read-1 / read-2 alt counts for GIV | ≥ 0 |
| refForward..altReverse | `int` | required | 2×2 strand table cells | ≥ 0 |
| variants | `IEnumerable<ArtifactObservation>` | required | candidate variant set | not null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ArtifactCall.Type` | `ArtifactType` | None / FfpeDeamination / OxoG |
| `ArtifactCall.GivScore` | `double` | GIV ratio; 1.0 balanced; +∞ when only read-2 count is 0 |
| `ArtifactCall.StrandBiasPhred` | `double` | Phred-scaled FisherStrand FS (≥ 0) |
| `ArtifactCall.IsArtifact` | `bool` | flagged as artifact |
| `FilterArtifacts` | `IReadOnlyList<ArtifactObservation>` | non-artifact subset, input order |
| `DetectOxoGArtifacts` | `IReadOnlyList<ArtifactCall>` | OxoG calls with GIV > 1.5 |

### 3.3 Preconditions and Validation

Bases are upper-cased before classification (case-insensitive). Negative counts throw `ArgumentOutOfRangeException`; null sequences throw `ArgumentNullException`. GIV with a zero read-2 count returns 1.0 when read-1 is also 0 (no imbalance) and `double.PositiveInfinity` otherwise. An all-zero strand table yields p = 1 (FS = 0).

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case ref/alt; map the ordered (ref, alt) pair to an `ArtifactType` (C>T/G>A → FFPE; G>T/C>A → OxoG; else None).
2. Compute GIV = read1Alt / read2Alt.
3. Compute FS = −10·log10(two-sided Fisher exact p) on the strand table.
4. Flag as artifact: FFPE always; OxoG iff GIV > 1.5.
5. `FilterArtifacts` keeps non-flagged; `DetectOxoGArtifacts` returns flagged OxoG calls.

### 4.2 Decision Rules, Scoring, Reference Tables

| ref→alt | Class | Flag rule |
|---------|-------|-----------|
| C>T, G>A | FfpeDeamination | always |
| G>T, C>A | OxoG | GIV > 1.5 |
| other | None | never |

Constants: `DamagedGivThreshold = 1.5` [3]; `UndamagedGivScore = 1.0` [3]; `MinFisherPValue = 1E-320` [6].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ClassifyArtifact | O(N) | O(1) | N = column sum of the strand table (Fisher sum over feasible cells) |
| FilterArtifacts / DetectOxoGArtifacts | O(n·N) | O(n) | n = variant count |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ClassifyArtifact(ArtifactObservation)`: classifies one variant.
- `OncologyAnalyzer.CalculateGivScore(int, int)`: GIV ratio.
- `OncologyAnalyzer.CalculateStrandBias(int, int, int, int)`: Phred-scaled FisherStrand FS.
- `OncologyAnalyzer.DetectOxoGArtifacts(IEnumerable<ArtifactObservation>)`: OxoG calls.
- `OncologyAnalyzer.FilterArtifacts(IEnumerable<ArtifactObservation>)`: removes artifacts.

### 5.2 Current Behavior

The two-sided Fisher exact p-value is computed by summing hypergeometric probabilities of all same-margin tables whose probability is ≤ the observed table's, using a Lanczos log-gamma for numerically stable binomial coefficients. No substring search / matching is involved, so the repository suffix tree is **not used** (N/A — this is a substitution-class and statistical-test unit, not an occurrence-enumeration unit).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- FFPE C>T / G>A and OxoG G>T / C>A substitution classes [1][4].
- GIV = read1/read2 ratio, threshold 1.5 = damaged, 1.0 = undamaged [1][3].
- FisherStrand FS = −10·log10(two-sided Fisher p) on the [refFwd, refRev, altFwd, altRev] table with MIN_PVALUE floor [6].

**Intentionally simplified:**

- OxoG confirmation uses the GIV ratio threshold directly rather than a per-context (tri-nucleotide) GIV; **consequence:** the call is per-substitution-type, not sequence-context-specific.

**Not implemented:**

- Bayesian read-orientation artifact posterior (GATK F1R2/OrientationBiasFilter); **users should rely on:** GATK Mutect2 + `LearnReadOrientationModel` / `FilterMutectCalls` for a probabilistic model.
- BAM parsing of strand/read-mate counts; **users should rely on:** an upstream caller to populate `ArtifactObservation` counts.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | strand/read-mate counts supplied on the record, not parsed from a BAM | Assumption | API shape only; rules unchanged | accepted | repo has no BAM reader (ASM-01) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| GIV with read2 = 0, read1 = 0 | GIV = 1.0 | no imbalance evidence [1] |
| GIV with read2 = 0, read1 > 0 | GIV = +∞ (damaged) | maximal one-sided imbalance [1] |
| balanced strand table | FS = 0 | p = 1 [6] |
| substitution not in the four artifact pairs | Type = None, not flagged | classes are specific [1][4] |
| empty variant list | empty result | nothing to filter |
| null variant list | `ArgumentNullException` | API contract |

### 6.2 Limitations

Rule-based: it does not estimate a per-variant artifact probability, does not model tri-nucleotide context, and does not distinguish OxoG from other G>T sources beyond the read-orientation imbalance. FFPE flagging is by substitution class alone, so a genuine somatic C>T at sufficient VAF would also be flagged — intended as a conservative pre-filter, not a final caller.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var v = new OncologyAnalyzer.ArtifactObservation(
    ReferenceAllele: 'G', AlternateAllele: 'T',
    RefForward: 20, RefReverse: 18, AltForward: 12, AltReverse: 10,
    AltReadsR1: 200, AltReadsR2: 100);          // GIV = 2.0 > 1.5
var call = OncologyAnalyzer.ClassifyArtifact(v); // Type = OxoG, IsArtifact = true
```

**Numerical walk-through:** strand table [20,0,0,20] → two-sided Fisher p = 1.4508889×10⁻¹¹ → FS = −10·log10(p) = 108.384 [6].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_FilterArtifacts_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_FilterArtifacts_Tests.cs) — covers INV-01..INV-05
- Evidence: [ONCO-ARTIFACT-001-Evidence.md](../../../docs/Evidence/ONCO-ARTIFACT-001-Evidence.md)

## 8. References

1. Chen L., Liu P., Evans T.C. Jr., Ettwiller L.M. 2017. DNA damage is a pervasive cause of sequencing errors, directly confounding variant identification. Science 355(6326):752–756. https://www.science.org/doi/10.1126/science.aai8690
2. Ettwiller L. Damage-estimator (reference implementation; GIV score). https://github.com/Ettwiller/Damage-estimator
3. Nature Methods. 2017. DNA variants or DNA damage? Nat Methods 14:330. https://www.nature.com/articles/nmeth.4254
4. Do H., Dobrovic A. 2015. Sequence Artifacts in DNA from Formalin-Fixed Tissues / Deamination Effects in FFPE Tissue Samples. Clin Chem / ScienceDirect. https://www.sciencedirect.com/science/article/pii/S152515781630188X
5. Comment on "DNA damage is a pervasive cause of sequencing errors" (GIV_G_T interpretation). PMC7350422. https://pmc.ncbi.nlm.nih.gov/articles/PMC7350422/
6. Broad Institute. GATK FisherStrand (FS) — Fisher's Exact Test for strand bias; StrandBiasTest.java. https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandBiasTest.java
