# Off-Target Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | CRISPR-OFF-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Off-target analysis estimates potential unintended cleavage sites for CRISPR guide RNAs by scanning a genome or reference sequence for near matches adjacent to valid PAMs. In this repository, off-target sites are defined as PAM-supported guide-length targets with at least one mismatch and no more than a configured mismatch limit, and they are scored by a simple position-dependent mismatch penalty. The resulting scores are intended for heuristic ranking rather than experimental prediction.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

CRISPR off-target activity arises when a guide RNA binds and cleaves genomic loci that are similar, but not identical, to the intended target. The original document highlights the importance of the PAM-proximal seed region, cites observed off-target activity with 3-5 mismatches, and notes that PAM-proximal mismatches are especially important for specificity. Sources: Hsu et al. (2013), Fu et al. (2013), Wikipedia (Off-target genome editing).

### 2.2 Core Model

The repository's off-target search applies three criteria:

1. The candidate site must have a valid PAM for the selected CRISPR system.
2. The PAM-adjacent target must have the correct guide length.
3. The guide-to-target mismatch count must satisfy `0 < mismatches <= maxMismatches`.

The off-target score is computed by summing mismatch penalties by position:

| Region | Weight |
|--------|--------|
| Seed region mismatch | `5` |
| Non-seed mismatch | `2` |

The seed region is the last 12 bp for PAM-after-target systems such as Cas9 and the first 12 bp for PAM-before-target systems such as Cas12a. `CalculateSpecificityScore(...)` finds off-targets with up to 4 mismatches, sums their penalties, and returns `max(0, 100 - totalPenalty)`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Exact matches are not returned as off-targets | `FindOffTargets(...)` yields only when `mismatches > 0` |
| INV-02 | Returned mismatch counts are bounded by the requested `maxMismatches` | The source filters on `mismatches <= maxMismatches` |
| INV-03 | Specificity scores are clamped to the range `[0, 100]` | `CalculateSpecificityScore(...)` applies `Math.Max(0, 100 - totalPenalty)` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `guideSequence` | `string` | required | Guide RNA sequence to analyze | Must match the guide length of the selected CRISPR system |
| `genome` | `DnaSequence` | required | Genome or reference sequence to search | Null input throws `ArgumentNullException` |
| `maxMismatches` | `int` | `3` | Maximum number of mismatches allowed in `FindOffTargets(...)` | Must be between `0` and `5` inclusive |
| `systemType` | `CrisprSystemType` | `SpCas9` | CRISPR system definition used for PAM and guide length | Must resolve through `GetSystem(...)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Position of the off-target PAM site |
| `Sequence` | `string` | Target sequence adjacent to the PAM |
| `Mismatches` | `int` | Guide-to-target mismatch count |
| `MismatchPositions` | `IReadOnlyList<int>` | Zero-based mismatch positions |
| `IsForwardStrand` | `bool` | Strand orientation of the off-target site |
| `OffTargetScore` | `double` | Position-weighted mismatch penalty |
| `specificity` | `double` | Overall specificity score in the range `0-100` |

### 3.3 Preconditions and Validation

`FindOffTargets(...)` throws `ArgumentNullException` for null guide or genome input and `ArgumentOutOfRangeException` when `maxMismatches` is outside `0-5`. It also throws `ArgumentException` when the guide length does not match the expected guide length for the selected CRISPR system. `CalculateSpecificityScore(...)` always evaluates off-targets using a fixed mismatch cap of 4.

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve the CRISPR system into its PAM pattern, guide length, and orientation.
2. Find all PAM sites in the genome on both strands.
3. Extract the guide-length target sequence for each PAM site.
4. Count mismatches between the guide and each target.
5. Yield sites whose mismatch count is positive and within the configured threshold.
6. Score each site by assigning a higher penalty to mismatches inside the 12-base seed region.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

CRISPR system variations documented for off-target analysis:

| System | PAM | Guide Length | PAM Position | Seed Region |
|--------|-----|--------------|--------------|-------------|
| SpCas9 | NGG | 20 bp | After target | Last 12 bp |
| SaCas9 | NNGRRT | 21 bp | After target | Last 12 bp |
| Cas12a | TTTV | 23 bp | Before target | First 12 bp |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindOffTargets` | `O(n × m)` | `O(k)` | `n` is genome length, `m` is guide length, `k` is number of off-targets |
| `CalculateSpecificityScore` | `O(n × m)` | `O(k)` | Reuses `FindOffTargets(...)` with a mismatch cap of 4 |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CrisprDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs)

- `CrisprDesigner.FindOffTargets(string, DnaSequence, int, CrisprSystemType)`: Finds PAM-constrained near matches for a guide sequence.
- `CrisprDesigner.CalculateSpecificityScore(string, DnaSequence, CrisprSystemType)`: Converts off-target penalties into a `0-100` specificity score.

### 5.2 Current Behavior

The current implementation enumerates PAM sites first and then compares the guide against each PAM-adjacent target sequence. Mismatch positions are collected explicitly, and exact matches are excluded from the off-target results. `CalculateSpecificityScore(...)` always uses up to 4 mismatches, regardless of the `FindOffTargets(...)` default of 3. Seed mismatches add `5` points and non-seed mismatches add `2` points to the off-target score.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PAM-constrained candidate enumeration.
- A PAM-proximal seed region with stronger mismatch penalties.
- A specificity score that decreases as off-target penalties accumulate.

**Intentionally simplified:**

- Off-target scoring uses only two positional weights (`5` in the seed region and `2` outside it); **consequence:** base-specific and context-specific models such as CFD or MIT-style scoring are not represented.
- The implementation uses strict PAM matching for the chosen system; **consequence:** off-targets at other PAMs are ignored unless that PAM is modeled as its own system (for example `SpCas9_NAG`).
- Only base mismatches are considered; **consequence:** bulges and gap-containing alignments are out of scope.

**Not implemented:**

- Chromatin accessibility, bulge-aware alignment, and richer off-target scoring models; **users should rely on:** external experimental or specialized computational off-target tools when those factors matter.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No off-targets found | Specificity score is `100` | Maximum specificity per implementation |
| `maxMismatches < 0` or `> 5` | Throws `ArgumentOutOfRangeException` | Explicit source guard |
| Guide length mismatches the system | Throws `ArgumentException` | The method enforces system-specific guide length |
| Exact on-target match | Not returned as an off-target | Exact matches are excluded by the `mismatches > 0` filter |

### 6.2 Limitations

The current implementation is a simplified position-weighted sequence model. It does not include chromatin state, indel/bulge alignments, or richer empirical scoring functions, and it relies on the CRISPR system's explicit PAM definition for site discovery.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
string guide = "ACGTACGTACGTACGTACGT";
var genome = new DnaSequence("ACGT...long genome...ACGT");

var offTargets = CrisprDesigner.FindOffTargets(
    guide,
    genome,
    maxMismatches: 3,
    CrisprSystemType.SpCas9);

double specificity = CrisprDesigner.CalculateSpecificityScore(
    guide, genome, CrisprSystemType.SpCas9);
```

## 8. References

1. Hsu PD, Scott DA, Weinstein JA, et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. *Nature Biotechnology*, 31(9):827-832. doi:10.1038/nbt.2647
2. Fu Y, Foden JA, Khayter C, et al. (2013). High-frequency off-target mutagenesis induced by CRISPR-Cas nucleases in human cells. *Nature Biotechnology*, 31(9):822-826. doi:10.1038/nbt.2623
3. Wikipedia. Off-target genome editing. https://en.wikipedia.org/wiki/Off-target_genome_editing
