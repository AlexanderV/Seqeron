# Centromere Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Chromosome Analysis |
| Test Unit ID | CHROM-CENT-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Centromere analysis identifies a candidate centromeric region in a chromosome-scale sequence and classifies the chromosome by centromere position. In this repository, `AnalyzeCentromere` uses a sliding-window heuristic that favors repeat-rich regions with low GC-content variability, then derives a cytogenetic class from the centromere midpoint. The method is intended for computational sequence analysis rather than clinical or reference-grade centromere annotation. It is useful when an assembled sequence is available and approximate centromere localization is sufficient.[1][2][3]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The centromere is a specialized chromosome region that joins sister chromatids and provides the kinetochore attachment site during cell division.[1][4] It separates the short `p` arm from the long `q` arm of the chromosome.[1]

The existing repository documentation also records the following centromere-associated sequence features in human chromosomes.[1][4]

| Feature | Description |
|---------|-------------|
| Alpha-satellite DNA | Approximately 171 bp tandem repeat units |
| Repeat content | High repetitive DNA content, often arranged in higher-order repeat units |
| Chromatin state | Constitutive heterochromatin |
| GC variability | Lower variability than many gene-rich regions |

Centromere-position nomenclature follows Levan et al. (1964) and is based on the arm-length ratio `q/p`.[2]

| Classification | Arm Ratio (`q/p`) | Description |
|----------------|-------------------|-------------|
| Metacentric | `1.0` to `1.7` | Arms approximately equal in length |
| Submetacentric | `> 1.7` to `3.0` | Arms unequal but both substantial |
| Subtelocentric | `> 3.0` and `< 7.0` | One arm clearly shorter |
| Acrocentric | `>= 7.0` | One arm very short |
| Telocentric | `p = 0` | Centromere at an end |

### 2.2 Core Model

The biological classification model is the centromere-arm-ratio system of Levan et al. (1964): if the centromere midpoint splits the chromosome into arms of lengths `p` and `q`, the chromosome class depends on the ratio `q/p`.[2] The sequence-localization part of the repository implementation is a heuristic model rather than a published centromere finder: it searches for windows that are simultaneously repeat-rich and relatively uniform in GC content.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The supplied sequence is long enough and contiguous enough for the centromere to appear as a repeat-rich interval. | The algorithm can return `Unknown` or localize a non-centromeric repetitive region instead. |
| ASM-02 | Repeat-rich, low-GC-variability windows are a useful proxy for the target centromeric region in the input sequence. | Non-centromeric repeats or unusual centromeres can produce false positives or missed calls. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | When a centromere is found, `Start <= End` and `Length = End - Start`. | The result length is computed directly from the two stored boundaries. |
| INV-02 | `CentromereType` is one of `Metacentric`, `Submetacentric`, `Subtelocentric`, `Acrocentric`, `Telocentric`, or `Unknown`. | `DetermineCentromereType` returns only those values. |
| INV-03 | `IsAcrocentric` is true if and only if `CentromereType == "Acrocentric"`. | The implementation sets the flag from the final type string. |
| INV-04 | `AlphaSatelliteContent >= 0`. | The stored value is derived from a non-negative score accumulator. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `chromosomeName` | `string` | required | Identifier copied into the result. | Preserved exactly in `CentromereResult.Chromosome`. |
| `sequence` | `string` | required | DNA sequence to scan for a centromere-like interval. | `null` or empty input returns `Unknown` with null boundaries. Sequences shorter than the scan window also return `Unknown`. |
| `windowSize` | `int` | `100000` | Sliding-window size used during the scan. | Windows are advanced by `windowSize / 4` during the initial scan. |
| `minAlphaSatelliteContent` | `double` | `0.3` | Minimum score threshold required before a candidate region is accepted. | Boundary extension uses `70%` of this threshold. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CentromereResult.Chromosome` | `string` | Chromosome identifier provided by the caller. |
| `CentromereResult.Start` | `int?` | Start coordinate of the detected region, or `null` when no centromere-like region is found. |
| `CentromereResult.End` | `int?` | End coordinate of the detected region, or `null` when no centromere-like region is found. |
| `CentromereResult.Length` | `int` | `End - Start` when a region is found, otherwise `0`. |
| `CentromereResult.CentromereType` | `string` | Cytogenetic class derived from the centromere midpoint and chromosome length. |
| `CentromereResult.AlphaSatelliteContent` | `double` | Stored score for the best candidate region. |
| `CentromereResult.IsAcrocentric` | `bool` | Convenience flag for the `Acrocentric` class only. |

### 3.3 Preconditions and Validation

`AnalyzeCentromere` returns `Unknown` with null boundaries when the input sequence is `null`, empty, or effectively too short for the initial scan loop to execute. The method uppercases the sequence before analysis. It does not require a reference genome or alpha-satellite database; instead it derives the result entirely from the supplied sequence with a repeat-content and GC-variability heuristic.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return `Unknown` immediately when the input sequence is `null` or empty.
2. Uppercase the sequence and scan overlapping windows of size `windowSize` with a step of `windowSize / 4`.
3. For each window, estimate repeat content with 15-mer counting and GC variability with 1 kb sub-windows.
4. Compute the candidate score as `repeatContent * (1 - gcVariability)` and retain the highest-scoring window above `minAlphaSatelliteContent`.
5. Extend the chosen region left and right while neighboring half-windows maintain at least `0.7 * minAlphaSatelliteContent` repeat content.
6. Compute the centromere midpoint, derive the `q/p` arm ratio, classify the chromosome, and return a `CentromereResult`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository-specific scoring rule for the candidate centromeric region is:

```text
score = repeatContent * (1 - gcVariability)
```

Repeat content is estimated from repeated 15-mers, and GC variability is the standard deviation of GC fractions over 1 kb sub-windows. Classification then follows the Levan arm-ratio system shown in Section 2.1.[2]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `AnalyzeCentromere` | `O(n)` | `O(k)` | Existing repository documentation characterizes the scan as linear in sequence length, with `k` representing the temporary k-mer dictionary used per window. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ChromosomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs)

- `ChromosomeAnalyzer.AnalyzeCentromere(...)`: scans the input sequence, selects and extends the best-scoring region, and returns a `CentromereResult`.
- `ChromosomeAnalyzer.EstimateRepeatContent(...)`: estimates repetitiveness from repeated 15-mers.
- `ChromosomeAnalyzer.CalculateGcVariability(...)`: measures GC variability over fixed sub-windows.
- `ChromosomeAnalyzer.DetermineCentromereType(...)`: assigns the Levan-based centromere class from the final midpoint.
- `ChromosomeAnalyzer.DetectAlphaSatellite(...)` / `FindCenpBBoxes(...)` / `DetectHigherOrderRepeat(...)`: opt-in alpha-satellite-specific detection (171-bp tandem + AT + CENP-B box; HOR structure).
- `ChromosomeAnalyzer.AssignSuprachromosomalFamily(...)` / `LoadBundledAlphaSatelliteReference()`: **opt-in** suprachromosomal-family (SF) assignment against a bundled **CC0** Dfam alpha-satellite reference (additive; does not change the detectors above).

### 5.2 Current Behavior

The initial scan uses overlapping windows and a 15-mer repeat heuristic, while GC variability is measured with 1 kb sub-windows. Boundary extension uses a relaxed repeat-content threshold of `0.7 * minAlphaSatelliteContent`. The classification logic uses the Levan `q/p` arm-ratio thresholds and returns `Subtelocentric` for ratios between `3.0` and `7.0`. The method returns only the single best-scoring region and is explicitly intended for computational analysis rather than clinical diagnostics.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Centromere classes are assigned with the Levan et al. (1964) arm-ratio nomenclature, including `Metacentric`, `Submetacentric`, `Subtelocentric`, `Acrocentric`, and `Telocentric`.[2]
- The centromere is treated as the sequence feature that divides the chromosome into `p` and `q` arms.[1][2]

**Intentionally simplified:**

- The algorithm estimates centromere likelihood from repeated 15-mers and GC variability instead of aligning to alpha-satellite reference libraries; **consequence:** results are approximate and depend on sequence composition rather than explicit repeat-family matching.
- Only the best-scoring candidate region is returned; **consequence:** multiple centromere-like repetitive regions are collapsed to one output.

**Not implemented:**

- **SF1-vs-SF2 separation and chromosome-specific HOR identity.** `AssignSuprachromosomalFamily` assigns SF3 (pentameric), SF4 (monomeric A-type) and SF5 (irregular A/B) and narrows dimeric arrays to {SF1, SF2} from the bundled CC0 reference (Dfam ALR/ALRa/ALRb) + HOR period + A/B-box composition (Shepelev 2009; McNulty & Sullivan 2018[5][6]). Separating SF1 from SF2 (both dimeric, identical A→B pattern), and tagging SF3 arrays whose period is not a multiple of 5 (e.g. the dodecameric DXZ1), need the **SF-resolved consensus monomer library** (J1/J2/D1/D2/W1–W5/M1/R1/R2); those sequences are only in unlicensed third-party HMM repos and are not CC0/redistributable. **users should rely on:** a caller-supplied `reference` of SF-resolved consensus monomers.
- Multi-candidate or uncertainty-ranked centromere reporting; **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `AlphaSatelliteContent` stores the best composite score, not a direct alpha-satellite fraction. | Deviation | Users should interpret the field as a heuristic centromere score rather than a literal repeat fraction. | accepted | Directly confirmed from the source: the returned value is `maxScore = repeatContent * (1 - gcVariability)`. |
| 2 | `Telocentric` is supported in the classifier but effectively unreachable through `AnalyzeCentromere` for ordinary detected windows. | Assumption | The public method will ordinarily report a non-zero centromere midpoint when it finds a candidate region. | accepted | Documented in [CHROM-CENT-001.md](../../../tests/TestSpecs/CHROM-CENT-001.md). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null sequence | Returns `Unknown` with null boundaries and length `0`. | The method special-cases empty input. |
| Sequence shorter than the analysis window | Returns `Unknown`. | The initial scan loop does not execute, leaving the candidate region unset. |
| Non-repetitive sequence | Returns `Unknown`. | No window exceeds the minimum candidate score threshold. |
| Strongly repetitive sequence | Returns a detected region with non-null boundaries. | Repeat-rich windows can exceed the acceptance threshold. |
| Lowercase input | Produces the same result as uppercase input. | The method normalizes to uppercase before scoring. |

### 6.2 Limitations

This implementation is heuristic. It does not use reference alpha-satellite libraries, experimentally defined centromeres, or multi-signal centromere models. Artificially repetitive sequences can be scored as centromeric, and only the best-scoring region is returned even when several repetitive regions are present.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ChromosomeAnalyzer_Centromere_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Centromere_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Tests (SF assignment): [ChromosomeAnalyzer_SuprachromosomalFamily_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_SuprachromosomalFamily_Tests.cs)
- Bundled CC0 reference + provenance: [Resources/README.md](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/Resources/README.md)
- Test specification: [CHROM-CENT-001.md](../../../tests/TestSpecs/CHROM-CENT-001.md)
- Evidence: [CHROM-CENT-001-Evidence.md](../../../docs/Evidence/CHROM-CENT-001-Evidence.md)
- Related algorithms: [Karyotype_Analysis.md](Karyotype_Analysis.md)

## 8. References

1. Wikipedia contributors. 2026. Centromere. Wikipedia. https://en.wikipedia.org/wiki/Centromere
2. Levan A, Fredga K, Sandberg AA. 1964. Nomenclature for centromeric position on chromosomes. Hereditas. N/A
3. Mehta GD, Agarwal MP, Ghosh SK. 2010. Centromere identity: a challenge to be faced. Molecular Genetics and Genomics. N/A
4. Wikipedia contributors. 2026. Karyotype. Wikipedia. https://en.wikipedia.org/wiki/Karyotype
5. McNulty SM, Sullivan BA. 2018. Alpha satellite DNA biology: finding function in the recesses of the genome. Chromosome Research 26:115–138. https://pmc.ncbi.nlm.nih.gov/articles/PMC6121732/
6. Shepelev VA, Uralsky LI, Alexandrov AA, Yurov YB, Rogaev EI, Alexandrov IA. 2009. The evolutionary origin of man can be traced in the layers of defunct ancestral alpha satellites … PLOS Genetics 5(9):e1000641. https://doi.org/10.1371/journal.pgen.1000641
7. Storer J, Hubley R, Rosen J, Wheeler TJ, Smit AF. 2021. The Dfam community resource of transposable element families, sequence models, and genome annotations. Mobile DNA 12:2 (Dfam = CC0). https://doi.org/10.1186/s13100-020-00230-y
