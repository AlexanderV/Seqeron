# Synteny Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Chromosome Analysis |
| Test Unit ID | CHROM-SYNT-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Synteny analysis identifies conserved blocks of gene order between two genomes and uses those blocks to summarize candidate chromosomal rearrangements. In this repository, `FindSyntenyBlocks` groups ortholog pairs by chromosome pair, looks for collinear runs with consistent orientation and bounded gaps, and emits `SyntenyBlock` values. `DetectRearrangements` then inspects those blocks for inversions, translocations, deletions, and duplications. The implementation is coordinate-based and heuristic, so it is suited to comparative block finding rather than full whole-genome alignment or sequence-identity analysis.[1][3]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Synteny originally referred to loci located on the same chromosome; in comparative genomics it is commonly used in the stronger sense of conserved gene order, or collinearity, between species.[4][5] Conserved synteny is useful for evolutionary analysis, orthology support, and cross-species genome comparison.[4][5]

The current repository documentation identifies four rearrangement types that can disrupt synteny.[6]

| Type | Description | Detection Signature |
|------|-------------|---------------------|
| Inversion | Segment reversed in orientation | Strand change within the same chromosome pair |
| Translocation | Segment moved to a different chromosome | Change in the target chromosome between adjacent blocks |
| Deletion | Segment missing in one genome relative to the other | Asymmetric gap between adjacent blocks |
| Duplication | Segment copied | Overlapping source blocks that map to different target locations |

### 2.2 Core Model

The repository model is a coordinate-and-order view of synteny: ortholog pairs are first grouped by chromosome pair and ordered by the first genome, then split into blocks when orientation changes or the inter-gene gap exceeds a threshold. Rearrangement inference is then performed from the relationships among the resulting blocks rather than from base-level sequence comparison.[1][3]

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The input ortholog pairs are correct and represent comparable loci between the two genomes. | Collinear blocks and inferred rearrangements can reflect orthology errors rather than real genome structure. |
| ASM-02 | A megabase-scale `maxGap` threshold is appropriate for deciding whether adjacent orthologs belong to the same block. | Blocks can split too aggressively or merge across biologically distinct regions. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `SyntenyBlock.Strand` is always `'+'` or `'-'`. | The implementation derives orientation from monotonicity in the second genome and emits only those two values. |
| INV-02 | `Species1Start <= Species1End` and `Species2Start <= Species2End`. | Block boundaries are built from the minimum/maximum coordinates of the first and last genes in the run. |
| INV-03 | `SequenceIdentity` is `NaN` for coordinate-only input. | The implementation explicitly emits `double.NaN` because it does not compute sequence identity. |
| INV-04 | Rearrangement `Type` is one of `Inversion`, `Translocation`, `Deletion`, or `Duplication`. | `DetectRearrangements` emits only those labels. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[FindSyntenyBlocks] orthologPairs` | `IEnumerable<(string Chr1, int Start1, int End1, string Gene1, string Chr2, int Start2, int End2, string Gene2)>` | required | Ortholog gene pairs with coordinates in both genomes. | Empty input or fewer than `minGenes` total pairs produce no blocks. |
| `[FindSyntenyBlocks] minGenes` | `int` | `3` | Minimum number of genes required for a block. | Runs shorter than this threshold are discarded. |
| `[FindSyntenyBlocks] maxGap` | `int` | `10` | Maximum allowed inter-gene gap in megabases. | The implementation multiplies the value by `1,000,000` before comparing gaps. |
| `[DetectRearrangements] syntenyBlocks` | `IEnumerable<SyntenyBlock>` | required | Synteny blocks to inspect for rearrangements. | Empty or single-block input produces no rearrangements. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `SyntenyBlock.Species1Chromosome` | `string` | Chromosome name from the first genome. |
| `SyntenyBlock.Species1Start` | `int` | Start coordinate of the block in the first genome. |
| `SyntenyBlock.Species1End` | `int` | End coordinate of the block in the first genome. |
| `SyntenyBlock.Species2Chromosome` | `string` | Chromosome name from the second genome. |
| `SyntenyBlock.Species2Start` | `int` | Start coordinate of the block in the second genome. |
| `SyntenyBlock.Species2End` | `int` | End coordinate of the block in the second genome. |
| `SyntenyBlock.Strand` | `char` | Relative orientation of the block (`'+'` or `'-'`). |
| `SyntenyBlock.GeneCount` | `int` | Number of ortholog pairs retained in the block. |
| `SyntenyBlock.SequenceIdentity` | `double` | `NaN` placeholder because sequence identity is not computed from coordinate-only input. |
| `ChromosomalRearrangement.Type` | `string` | Rearrangement class. |
| `ChromosomalRearrangement.Chromosome1` | `string` | Chromosome in the first genome where the event is reported. |
| `ChromosomalRearrangement.Position1` | `int` | Primary breakpoint or anchor position. |
| `ChromosomalRearrangement.Chromosome2` | `string?` | Secondary chromosome when relevant, such as translocations or duplications. |
| `ChromosomalRearrangement.Position2` | `int?` | Secondary breakpoint or target position when relevant. |
| `ChromosomalRearrangement.Size` | `int?` | Event size when computed by the heuristic. |
| `ChromosomalRearrangement.Description` | `string?` | Human-readable summary assembled by the implementation. |

### 3.3 Preconditions and Validation

`FindSyntenyBlocks` materializes the ortholog-pair sequence and returns no output when the total pair count is below `minGenes`. The method groups input by `(Chr1, Chr2)` and interprets `maxGap` in megabases by multiplying it by `1,000,000`. `DetectRearrangements` sorts blocks by first-genome chromosome and start coordinate before applying its event-specific heuristics, and it does not require sequence data beyond the block coordinates and strand.

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize the ortholog pairs and group them by chromosome pair.
2. Sort each chromosome-pair group by the first-genome coordinate.
3. Track collinear runs by comparing orientation and gap size between consecutive pairs.
4. Emit a `SyntenyBlock` when a run ends and contains at least `minGenes` orthologs.
5. Sort emitted blocks by first-genome chromosome and coordinate.
6. Compare adjacent blocks for inversions, translocations, and deletions.
7. Compare block pairs for overlapping first-genome intervals that map to different targets, and emit duplications when found.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The rearrangement signatures preserved from the original file are:

| Type | Heuristic Signature |
|------|---------------------|
| Inversion | Same first-genome chromosome, same second-genome chromosome, opposite strand |
| Translocation | Same first-genome chromosome, different second-genome chromosome |
| Deletion | Same first-genome chromosome and strand, but a much larger gap in genome 1 than in genome 2 |
| Duplication | Overlapping genome-1 blocks that map to different genome-2 locations |

`FindSyntenyBlocks` emits `SequenceIdentity = NaN` because identity is not computable from coordinate-only input.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindSyntenyBlocks` | `O(n log n)` | `O(n)` | Sorting chromosome-pair groups dominates the scan. |
| `DetectRearrangements` | `O(n^2)` | `O(n)` | The adjacent-block pass is linear after sorting, but duplication detection uses an all-pairs overlap scan. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ChromosomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs)

- `ChromosomeAnalyzer.FindSyntenyBlocks(...)`: detects collinear forward or reverse blocks from ortholog coordinates.
- `ChromosomeAnalyzer.DetectRearrangements(...)`: infers inversions, translocations, deletions, and duplications from `SyntenyBlock` values.

### 5.2 Current Behavior

The implementation multiplies `maxGap` by `1,000,000`, so the public parameter is interpreted in megabases. It determines block orientation from the direction of movement in the second genome and sets `SequenceIdentity` to `double.NaN`. Deletions are called only when the gap in genome 1 is more than twice the corresponding gap in genome 2, and duplications are called when overlapping source blocks map to different target coordinates.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Conserved forward and reverse collinearity are represented as synteny blocks with `+` and `-` strand labels.[1][4]
- The repository reports inversions, translocations, deletions, and duplications as rearrangement classes.[6]

**Intentionally simplified:**

- Block detection uses only coordinate order, chromosome pairing, and gap thresholds; **consequence:** sequence-level conservation and gene-family scoring are not part of the block definition.
- Rearrangement detection is heuristic and based on pairwise block relationships; **consequence:** complex rearrangements are summarized through a limited set of signatures.

**Not implemented:**

- Sequence-identity calculation for `SyntenyBlock`; **users should rely on:** no current alternative.
- Whole-genome alignment or multi-genome synteny chaining; **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `SequenceIdentity` is always `NaN`. | Deviation | Downstream consumers cannot treat the field as a real similarity estimate. | accepted | Explicitly tested in [CHROM-SYNT-001.md](../../../tests/TestSpecs/CHROM-SYNT-001.md). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty ortholog input | No synteny blocks are emitted. | The implementation returns immediately when the input is empty. |
| Fewer than `minGenes` ortholog pairs | No blocks are emitted. | The total pair count is checked before scanning. |
| Single synteny block | No rearrangements are emitted. | Rearrangement heuristics require comparisons between blocks. |
| Collinear blocks without signature changes | No rearrangements are emitted. | None of the inversion/translocation/deletion/duplication rules are triggered. |
| Gap larger than `maxGap` | The run is split into separate blocks. | Gap thresholds terminate the current collinear block. |

### 6.2 Limitations

This implementation is coordinate-only and heuristic. It does not compute sequence identity, does not evaluate more elaborate genome-rearrangement models, and uses a bounded set of signatures for event calling. As a result, the output is best interpreted as a synteny-and-rearrangement summary rather than a complete comparative-genomics analysis.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
{
    ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
    ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
    ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
};

var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(orthologPairs, minGenes: 3, maxGap: 10);
var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks);
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ChromosomeAnalyzer_Synteny_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Synteny_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [CHROM-SYNT-001.md](../../../tests/TestSpecs/CHROM-SYNT-001.md)

## 8. References

1. Wang Y, Tang H, Debarry JD, Tan X, Li J, Wang X, Lee TH, Jin H, Marler B, Guo H, Kissinger JC, Paterson AH. 2012. MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. Nucleic Acids Research. N/A
2. Goel M, Sun H, Jiao WB, Schneeberger K. 2019. SyRI: finding genomic rearrangements and local sequence differences from whole-genome assemblies. Genome Biology. N/A
3. Liu D, Hunt M, Tsai IJ. 2018. Inferring synteny between genome assemblies: a systematic evaluation. BMC Bioinformatics. N/A
4. Wikipedia contributors. 2026. Synteny. Wikipedia. N/A
5. Wikipedia contributors. 2026. Comparative genomics. Wikipedia. N/A
6. Wikipedia contributors. 2026. Chromosomal rearrangement. Wikipedia. N/A
