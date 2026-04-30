# Probe Validation

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-VALID-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Probe validation assesses whether a hybridization probe is likely to bind specifically to its intended targets with limited cross-hybridization. In this repository, validation combines substitution-tolerant fixed-length window matching against reference sequences with self-complementarity and secondary-structure screening. The result is a simple validation record with a specificity score, issue list, and probe-quality flags.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Cross-hybridization occurs when a probe binds sequences other than its intended target, and this is a central design concern for FISH, DNA microarrays, qPCR, and related assays. The original document also notes that mismatch tolerance, assay stringency, and probe length affect specificity, and that self-complementarity and low-complexity content can increase the risk of non-specific behavior. Sources: Wikipedia (Hybridization probe, DNA microarray, BLAST), Altschul et al. (1990), Amann & Ludwig (2000).

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

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProbeDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs)

- `ProbeDesigner.ValidateProbe(string, IEnumerable<string>, int, double)`: Performs approximate-match, self-complementarity, and secondary-structure validation.
- `ProbeDesigner.CheckSpecificity(string, ISuffixTree)`: Computes exact-hit specificity from suffix-tree occurrence counts.

### 5.2 Current Behavior

The current validator treats an empty probe as invalid rather than throwing. It records a cross-hybridization issue when more than one hit is found across all references and records a self-complementarity issue when the computed fraction exceeds the configured threshold. Approximate matching is implemented as an ungapped fixed-length sliding scan with mismatch tolerance rather than a gapped local-alignment search. `IsValid` becomes `true` either when no issues are found or when the probe has at most one hit and self-complementarity no greater than `0.4`. `CheckSpecificity(...)` uses exact suffix-tree hits rather than approximate matching.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Approximate-match-based cross-hybridization screening.
- Self-complementarity and secondary-structure checks.
- Specificity scoring as `0`, `1`, or `1 / hits` depending on hit count.

**Intentionally simplified:**

- Specificity collapses all multi-hit outcomes to `1 / hits`; **consequence:** it distinguishes hit multiplicity but not mismatch severity, thermodynamics, or genomic context.
- Approximate probe matching is substitution-only and fixed-length; **consequence:** gaps, bulges, and local realignment effects are not modeled in `ValidateProbe(...)`.
- Suffix-tree specificity uses exact hits only; **consequence:** approximate off-targets are only modeled through `ValidateProbe(...)`, not through `CheckSpecificity(...)`.

**Not implemented:**

- Thermodynamic hybridization modeling and BLAST-like alignment ranking; **users should rely on:** external experimental validation or richer sequence-search tools when required.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty probe sequence | Returns a structured invalid result | Explicit special case in source |
| Unique probe hit | `SpecificityScore = 1.0` | Exact mapping in implementation |
| No probe hits | `SpecificityScore = 0.0` | Probe does not match the references |
| Multiple hits | `SpecificityScore = 1.0 / hits` | Cross-hybridization penalty |

### 6.2 Limitations

The current implementation is a lightweight screening tool. It does not incorporate thermodynamic binding models, mismatch-position weighting, assay stringency, or database-scale alignment heuristics, and the suffix-tree helper only captures exact-hit uniqueness.

## 8. References

1. Wikipedia: Hybridization probe - https://en.wikipedia.org/wiki/Hybridization_probe
2. Wikipedia: DNA microarray - https://en.wikipedia.org/wiki/DNA_microarray
3. Wikipedia: Off-target genome editing - https://en.wikipedia.org/wiki/Off-target_genome_editing
4. Wikipedia: BLAST (biotechnology) - https://en.wikipedia.org/wiki/BLAST_(biotechnology)
5. Altschul et al. (1990) - Basic local alignment search tool, J. Mol. Biol.
6. Amann R, Ludwig W (2000) - Ribosomal RNA-targeted nucleic acid probes for studies in microbial ecology, FEMS Microbiology Reviews.
