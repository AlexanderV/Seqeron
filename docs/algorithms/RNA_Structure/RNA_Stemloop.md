# RNA Stem-Loop Detection

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-STEMLOOP-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Stem-loop detection scans a single RNA sequence for hairpin-like motifs formed by a paired stem flanking an unpaired loop. These motifs are a fundamental secondary-structure pattern and provide the candidate structures used elsewhere in the repository's RNA prediction code. The implementation performs an exhaustive local scan over loop starts and loop sizes, extends stems while pairing remains valid, and returns every candidate meeting the configured minimum stem length. It is therefore a motif detector, not a global structure optimizer.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Stem-loops consist of a paired stem and a terminal loop. In RNA, standard admissible pairs include Watson-Crick `A-U` and `G-C` pairs plus wobble `G-U` pairs. Hairpin loops shorter than three nucleotides are sterically disfavored and are treated as impossible in standard RNA secondary-structure models.

### 2.2 Core Model

The repository's detector treats a stem-loop as a candidate loop interval surrounded by antiparallel complementary sequence. For each candidate loop, the algorithm extends outward from the loop boundaries while base pairing remains valid. Once extension stops, the candidate is retained only if the resulting stem length is at least the configured minimum.

The resulting `TotalFreeEnergy` is:

$$
\Delta G_{total} = \Delta G_{stem} + \Delta G_{hairpin}
$$

where `ΔG_stem` is computed from stacked base pairs and `ΔG_hairpin` from the hairpin-loop energy helper.

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Stem-loop motifs can be identified by local complementary extension around a loop candidate | More global structural context is ignored during motif detection |
| ASM-02 | Wobble pairing is acceptable when `allowWobble = true` | Candidate sets change when wobble pairs are biologically inappropriate for the use case |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Returned stem-loops have `Loop.Size >= 3` | `minLoopSize` is clamped to `3` before scanning |
| INV-02 | Returned stem-loops satisfy `Stem.Length >= minStemLength` | Candidates shorter than the threshold are discarded |
| INV-03 | Returned loop type is always `Hairpin` | The constructor uses `LoopType.Hairpin` for this API |
| INV-04 | `DotBracketNotation.Length = End - Start + 1` | The dot-bracket helper allocates exactly the detected span length |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `rnaSequence` | `string` | required | RNA sequence to analyze | Empty or too-short input yields no results |
| `minStemLength` | `int` | `3` | Minimum number of paired bases in a returned stem | Candidates below this threshold are rejected |
| `minLoopSize` | `int` | `3` | Minimum hairpin loop size | Values below `3` are clamped to `3` |
| `maxLoopSize` | `int` | `10` | Maximum scanned loop size | Restricts candidate enumeration |
| `allowWobble` | `bool` | `true` | Whether `G-U` / `U-G` wobble pairs are allowed during stem extension | When `false`, wobble pairs stop extension |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | Inclusive 5' start of the detected stem-loop span |
| `End` | `int` | Inclusive 3' end of the detected stem-loop span |
| `Stem` | `Stem` | Stem coordinates, base pairs, and stem free energy |
| `Loop` | `Loop` | Hairpin-loop coordinates, size, and sequence |
| `TotalFreeEnergy` | `double` | Sum of stem and loop free energies |
| `DotBracketNotation` | `string` | Local dot-bracket string covering only the detected stem-loop span |

### 3.3 Preconditions and Validation

The input sequence is uppercased internally. If `minLoopSize < 3`, it is clamped to `3`. When the sequence is empty or shorter than `minStemLength * 2 + minLoopSize`, the method yields no candidates. Stem extension stops at the first invalid base pair, or at the first wobble pair when `allowWobble` is `false`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Clamp `minLoopSize` to `3` and reject too-short sequences.
2. Enumerate loop-start positions across the sequence.
3. Enumerate loop sizes between `minLoopSize` and `maxLoopSize` for each start.
4. Extend the stem outward from the candidate loop boundaries while pairing remains valid.
5. Discard candidates whose stem length is below `minStemLength`.
6. Build `Stem`, `Loop`, `TotalFreeEnergy`, and local dot-bracket notation for each accepted candidate.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Base-pair extension uses `GetBasePairType`, so canonical pairs and optionally wobble pairs are accepted. Accepted base pairs are reversed before the result is built so that they are ordered from 5' to 3'. Stem energy comes from `CalculateStemEnergy`, and loop energy comes from `CalculateHairpinLoopEnergy`, including the special `G-U` closure detection path used by the source.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindStemLoops` | `O(n^2 * L)` | `O(k)` | `n` = sequence length, `L` = scanned loop-size range, `k` = number of returned stem-loops |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.FindStemLoops(string, int, int, int, bool)`: Enumerates candidate stem-loop motifs.
- `RnaSecondaryStructure.DetectPseudoknots(IReadOnlyList<BasePair>)`: Post hoc crossing-pair detector for already assembled base-pair sets.

### 5.2 Current Behavior

The detector returns all candidates that satisfy the local constraints; it does not suppress overlaps or attempt to choose a single best motif. Each result includes both structural coordinates and a locally scored energy. Pseudoknots are not predicted directly from sequence in this method, but the same class exposes `DetectPseudoknots` for already assembled base-pair sets.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Minimum hairpin-loop size of three nucleotides.
- Watson-Crick and optional wobble pairing for stem extension.
- Hairpin free-energy scoring through Turner 2004-based helpers.

**Intentionally simplified:**

- The method returns every qualifying candidate instead of choosing a globally optimal or non-overlapping subset; **consequence:** callers may need a downstream selection step.
- Only hairpin stem-loops are modeled directly; **consequence:** internal-loop and whole-structure reasoning is left to other APIs.

**Not implemented:**

- Sequence-level pseudoknot prediction in this method; **users should rely on:** `RnaSecondaryStructure.DetectPseudoknots(...)` for post hoc detection or `RnaSecondaryStructure.PredictStructure(...)` for higher-level assembly.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or too-short sequence | No results | No valid loop-plus-stem span can be formed |
| `allowWobble = false` | Wobble pairs terminate stem extension | Explicit branch in the extension loop |
| Overlapping candidates | All qualifying candidates are yielded | The method performs no overlap filtering |

### 6.2 Limitations

This API is restricted to local hairpin-like motifs. It does not optimize over the whole sequence, does not incorporate pseudoknots during sequence scanning, and still depends on simplified downstream interpretation when many overlapping candidates are returned.

## 7. Examples and Related Material

- [RNA-STEMLOOP-001](../../../tests/TestSpecs/RNA-STEMLOOP-001.md) documents the repository's stem-loop test specification.
- [RNA_Secondary_Structure.md](./RNA_Secondary_Structure.md) documents the higher-level structure-prediction API that consumes these candidates.

## 8. References

1. Wikipedia contributors. Stem-loop. Wikipedia. https://en.wikipedia.org/wiki/Stem-loop
2. Wikipedia contributors. Tetraloop. Wikipedia. https://en.wikipedia.org/wiki/Tetraloop
3. Wikipedia contributors. Pseudoknot. Wikipedia. https://en.wikipedia.org/wiki/Pseudoknot
4. Woese, C. R., S. Winker, and R. R. Gutell. 1990. Architecture of ribosomal RNA: constraints on the sequence of tetra-loops. Proceedings of the National Academy of Sciences 87(21):8467-8471.
5. Heus, H. A., and A. Pardi. 1991. Structural features that give rise to the unusual stability of RNA hairpins containing GNRA loops. Science 253(5016):191-194.
6. Svoboda, P., and A. Cara. 2006. Hairpin RNA: A secondary structure of primary importance. Cellular and Molecular Life Sciences 63(7):901-908.
