# Probe Validation

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-VALID-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-24 |

## 1. Overview

Probe validation assesses whether a hybridization probe is likely to bind specifically to its intended targets with limited cross-hybridization. In this repository, validation combines substitution-tolerant fixed-length window matching against reference sequences with self-complementarity and secondary-structure screening. The result is a simple validation record with a specificity score, issue list, and probe-quality flags.

An opt-in **gapped** off-target scan (`ScanOffTargetsGapped`) supplements the default ungapped scan: it reuses the library's validated Smith–Waterman local aligner [1] to find off-target sites reachable through insertions or deletions — the indel-aware "BLAST-grade" improvement over a pure ungapped Hamming scan [2] — and it separates the single intended on-target match from genuine off-target hits (correcting the on/off-target pooling of `ValidateProbe`'s `OffTargetHits`). The default `ValidateProbe`/`CheckSpecificity` behaviour is unchanged.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Cross-hybridization occurs when a probe binds sequences other than its intended target, and this is a central design concern for FISH, DNA microarrays, qPCR, and related assays. The original document also notes that mismatch tolerance, assay stringency, and probe length affect specificity, and that self-complementarity and low-complexity content can increase the risk of non-specific behavior. Sources: Wikipedia (Hybridization probe, DNA microarray, BLAST), Altschul et al. (1990), Amann & Ludwig (2000).

Off-target detection by *local alignment* follows the Smith–Waterman recurrence [1]: `H(i,j) = max{ H(i-1,j-1) + s(a_i,b_j), H(i-1,j) - W, H(i,j-1) - W, 0 }`, whose zero floor returns the best-scoring local subsequence match and whose gap terms admit insertions/deletions. The gapped local alignment is the recognized improvement over ungapped matching for similarity/homology search [2]. The identity threshold above which a non-target sequence is treated as a cross-hybridizing off-target is taken from Kane et al. (2000) [7]: non-target transcripts **>75% similar** over the probe length may cross-hybridize; this is the default `minIdentity = 0.75`, exposed as a caller parameter.

### 2.2 Core Model

The validation workflow normalizes the probe to uppercase, counts approximate matches across all supplied reference sequences, and then maps hit counts to specificity:

$$
specificity =
\begin{cases}
0.0, & offTargetHits = 0 \\
1.0, & offTargetHits = 1 \\
1.0 / offTargetHits, & offTargetHits > 1
\end{cases}
$$

It also computes self-complementarity as a fraction of aligned matches against the reverse complement and checks for secondary-structure potential with a sequence-level hairpin screen.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0.0 <= SpecificityScore <= 1.0` | The score is explicitly mapped from hit counts to `0`, `1`, or `1/hits` |
| INV-02 | `SelfComplementarity` is non-negative and at most `1.0` | It is computed as a fraction of aligned positions |
| INV-03 | `OffTargetHits >= 0` | Hit counts are accumulated from match enumeration |
| INV-04 | `OffTargetHits == 1` implies `SpecificityScore == 1.0` | That mapping is explicit in source |
| INV-05 | Each gapped hit has `0.0 <= Identity <= 1.0` and `0.0 <= Coverage <= 1.0` | Both are counts of aligned columns divided by probe length |
| INV-06 | The intended on-target (first perfect ungapped full-coverage exact match) is excluded from `OffTargetHits` | On/off separation in `ScanOffTargetsGapped` |
| INV-07 | An indel-only off-target has `HasGaps == true` and is found by `ScanOffTargetsGapped` but not by the ungapped `ValidateProbe` scan | Gapped local alignment admits indels [1][2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `probeSequence` | `string` | required | Probe sequence to validate | Null input throws `ArgumentNullException`; empty string yields a structured invalid result |
| `referenceSequences` | `IEnumerable<string>` | required | Reference sequences scanned for approximate matches | Null input throws `ArgumentNullException` |
| `maxMismatches` | `int` | `3` | Maximum mismatch tolerance for approximate matching | Passed through to the internal approximate-match search |
| `selfComplementarityThreshold` | `double` | `0.3` | Threshold above which self-complementarity is recorded as an issue | Default documented in source |
| `genomeIndex` | `ISuffixTree` | required for `CheckSpecificity(...)` | Pre-built suffix tree for exact hit counting | Used only by the suffix-tree specificity helper |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `IsValid` | `bool` | Validation outcome based on issue count, hit count, and self-complementarity |
| `SpecificityScore` | `double` | Specificity value in the range `0.0-1.0` |
| `OffTargetHits` | `int` | Total approximate hits across all reference sequences |
| `SelfComplementarity` | `double` | Fraction of self-complementary positions |
| `HasSecondaryStructure` | `bool` | Hairpin-potential flag |
| `Issues` | `IReadOnlyList<string>` | Recorded validation issues |

### 3.3 Preconditions and Validation

`ValidateProbe(...)` uppercases the input probe before analysis. Null probe or reference collections raise `ArgumentNullException`. An empty probe sequence returns a structured invalid result with `SpecificityScore = 0.0`, `OffTargetHits = 0`, and an `"Empty probe sequence"` issue. `CheckSpecificity(...)` uppercases the probe before querying the suffix tree.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the probe sequence to uppercase.
2. Search every reference sequence for approximate matches within the mismatch tolerance.
3. Accumulate the total number of hits across all references.
4. Compute probe self-complementarity and secondary-structure potential.
5. Map hit count to the specificity score and derive `IsValid` from the issue set plus the `offTargetHits <= 1 && selfComp <= 0.4` fallback rule.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Validation defaults preserved from the original document and source:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `maxMismatches` | `3` | Approximate-match tolerance |
| `selfComplementarityThreshold` | `0.3` | Threshold for recording a self-complementarity issue |
| Secondary-structure check | enabled | `ValidateProbe(...)` always checks for hairpin potential |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ValidateProbe` | `O(n × g × m)` | `O(1)` auxiliary | The original document describes dependence on probe length `n`, reference count `g`, and reference lengths `m` |
| `CheckSpecificity` | `O(m)` | `O(1)` | Suffix-tree exact hit counting for probe length `m` |
| `ScanOffTargetsGapped` | `O(g × n × m)` | `O(n × m)` per window | Sliding Smith–Waterman: for each of `g` reference start positions, a local alignment of probe length `n` against a window of width `m = n + 2`. Reuses `SequenceAligner.LocalAlign`. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProbeDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs)

- `ProbeDesigner.ValidateProbe(string, IEnumerable<string>, int, double)`: Performs ungapped approximate-match, self-complementarity, and secondary-structure validation (default behaviour, unchanged).
- `ProbeDesigner.CheckSpecificity(string, ISuffixTree)`: Computes exact-hit specificity from suffix-tree occurrence counts.
- `ProbeDesigner.ScanOffTargetsGapped(string, IEnumerable<string>, double, ScoringMatrix?)`: Opt-in gapped (Smith–Waterman) off-target scan. Returns a `GappedSpecificityResult` separating `OnTargetHits` (the perfect ungapped full-coverage exact match) from `OffTargetHits` (imperfect/indel hits ≥ `minIdentity`, default 0.75). Reuses `SequenceAligner.LocalAlign` for the indel-aware alignment.

### 5.2 Current Behavior

The current validator treats an empty probe as invalid rather than throwing. It records a cross-hybridization issue when more than one hit is found across all references and records a self-complementarity issue when the computed fraction exceeds the configured threshold. `ValidateProbe`'s approximate matching is implemented as an ungapped fixed-length sliding scan with mismatch tolerance rather than a gapped local-alignment search (its `OffTargetHits` pools the on-target with off-targets). `IsValid` becomes `true` either when no issues are found or when the probe has at most one hit and self-complementarity no greater than `0.4`. `CheckSpecificity(...)` uses exact suffix-tree hits rather than approximate matching.

`ScanOffTargetsGapped(...)` is the opt-in gapped alternative. For each reference it slides a window of length `probeLen + 2` and runs `SequenceAligner.LocalAlign` (Smith–Waterman, `BlastDna` scoring by default), computes identity = identical aligned columns / probe length and coverage = ungapped columns / probe length, keeps every site whose identity ≥ `minIdentity`, and collapses overlapping window detections to one best (highest-identity, then highest-coverage, then leftmost) hit per site via greedy non-overlapping selection. The first perfect ungapped full-coverage exact match (identity = coverage = 1.0, no gaps) is classified as the intended on-target and excluded from the off-target count; all imperfect or indel-containing hits — and any additional perfect repeats — are off-targets.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Ungapped (Hamming) approximate-match cross-hybridization screening (`ValidateProbe`).
- Gapped (Smith–Waterman) local-alignment off-target scan with indel handling [1][2] and on/off-target separation (`ScanOffTargetsGapped`).
- Self-complementarity and secondary-structure checks.
- Specificity scoring as `0`, `1`, or `1 / hits` depending on hit count (`ValidateProbe`).
- Off-target identity threshold (default 0.75 over the probe length) per Kane et al. (2000) [7].

**Intentionally simplified:**

- `ValidateProbe`'s specificity collapses all multi-hit outcomes to `1 / hits`; **consequence:** it distinguishes hit multiplicity but not mismatch severity, thermodynamics, or genomic context.
- `ValidateProbe`'s approximate matching is substitution-only and fixed-length, and its `OffTargetHits` pools the on-target match with off-targets; **consequence:** for indel-aware detection and on/off separation use `ScanOffTargetsGapped` instead.
- Suffix-tree specificity uses exact hits only; **consequence:** approximate off-targets are only modeled through `ValidateProbe(...)`/`ScanOffTargetsGapped(...)`, not through `CheckSpecificity(...)`.
- On/off-target labelling: the first perfect ungapped full-coverage exact match is taken as the intended on-target; **consequence:** when several identical perfect sites exist, the first is on-target and the rest are off-targets.

**Not implemented:**

- A seeded BLAST k-mer index over a whole genome; `ScanOffTargetsGapped` is an exhaustive sliding Smith–Waterman scan (O(g · n·m)), not a genome-scale seed-and-extend index [2]; **users should rely on:** an external seeded aligner for genome-scale off-target search.
- Thermodynamic hybridization (duplex-Tm) modeling and E-value ranking; **users should rely on:** external experimental validation or richer sequence-search tools when required.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty probe sequence | Returns a structured invalid result | Explicit special case in source |
| Unique probe hit | `SpecificityScore = 1.0` | Exact mapping in implementation |
| No probe hits | `SpecificityScore = 0.0` | Probe does not match the references |
| Multiple hits | `SpecificityScore = 1.0 / hits` | Cross-hybridization penalty |

### 6.2 Limitations

The implementation is a screening tool. The opt-in `ScanOffTargetsGapped` adds indel-aware (gapped) off-target detection and on/off-target separation, but it remains an exhaustive sliding Smith–Waterman scan, not a seeded BLAST k-mer index over a whole genome; it does not incorporate thermodynamic binding (duplex-Tm) models, mismatch-position weighting, assay stringency, or E-value ranking, and the suffix-tree helper only captures exact-hit uniqueness.

## 8. References

1. Smith TF, Waterman MS (1981) - Identification of common molecular subsequences, J. Mol. Biol. 147(1):195–197. https://doi.org/10.1016/0022-2836(81)90087-5 (recurrence via https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
2. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) - Basic local alignment search tool, J. Mol. Biol. 215(3):403–410. https://doi.org/10.1016/S0022-2836(05)80360-2 (via https://en.wikipedia.org/wiki/BLAST_(biotechnology))
3. Wikipedia: Hybridization probe - https://en.wikipedia.org/wiki/Hybridization_probe
4. Wikipedia: DNA microarray - https://en.wikipedia.org/wiki/DNA_microarray
5. Wikipedia: Off-target genome editing - https://en.wikipedia.org/wiki/Off-target_genome_editing
6. Amann R, Ludwig W (2000) - Ribosomal RNA-targeted nucleic acid probes for studies in microbial ecology, FEMS Microbiology Reviews.
7. Kane MD, Jatkoe TA, Stumpf CR, Lu J, Thomas JD, Madore SJ (2000) - Assessment of the sensitivity and specificity of oligonucleotide (50mer) microarrays, Nucleic Acids Research 28(22):4552–4557. https://pmc.ncbi.nlm.nih.gov/articles/PMC113865/
