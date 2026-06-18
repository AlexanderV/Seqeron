# Synteny / Collinearity Block Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-SYNTENY-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Synteny block detection finds chromosomal regions whose orthologous genes occur in the
same order (collinearity) between two genomes. Given a set of orthologous anchor pairs, the
algorithm chains adjacent anchors into collinear runs, scores each run with the MCScanX
dynamic-programming scheme (reward per anchored pair, penalty per intervening gene), and
reports the non-overlapping runs that meet a minimum score / anchor count [1]. It is a
heuristic, deterministic procedure used to identify conserved gene-order blocks (forward or
inverted) for comparative and evolutionary analyses.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Synteny originally denoted homologous genes co-located on a chromosome; in modern comparative
genomics it denotes preservation of gene order across genomes descended from a common
ancestor [3]. *Collinearity* is the stricter form: "Collinearity, a more specific form of
synteny, requires conserved gene order." [1]. Anchors are homologous (orthologous or
paralogous) gene pairs; "anchor genes are more likely to be homologs" [1].

### 2.2 Core Model

MCScanX detects collinear blocks by dynamic programming over anchor pairs sorted by position,
rewarding adjacent collinear pairs and penalizing the distance between them [1]:

```
Score(v) = max( MatchScore(v),
                max_u ( Score(u) + MatchScore(v) + GapPenalty × NumberofGaps(u,v) ) )
```

with defaults `MatchScore = 50`, `GapPenalty = −1`, and `NumberofGaps(u,v) < 25`
(maximum intervening genes between anchors `u` and `v`) [1]. Non-overlapping chains with
score over 250 — i.e. at least 5 collinear gene pairs — are reported [1]. Anchors are
generated from BLASTP matches (default E-value cutoff `10^−5`); near-duplicate matches that
share a gene whose partners are within five genes are collapsed to the lowest-E-value pair [1].
Anchors are sorted "in both transcriptional directions", which yields both forward and
inverted (reverse) collinear blocks [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Supplied anchors are true orthologs (homologs). | Spurious anchors produce false collinear blocks. |
| ASM-02 | Genes are provided in chromosomal order per genome. | Out-of-order input misrepresents gene adjacency, corrupting chains. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported block has `GeneCount ≥ minAnchors` (≥ 5 by default). | Report rule requires ≥ 5 collinear pairs [1]. |
| INV-02 | Every reported block has MCScanX score ≥ 250 (default). | Only chains scoring over 250 are reported [1]. |
| INV-03 | Within a block, consecutive anchors keep one genome-2 direction and `NumberofGaps < maxGap`. | Chaining rejects direction changes and oversized gaps [1]. |
| INV-04 | `IsInverted = true` iff the genome-2 order decreases along the block. | Reverse direction = inverted block [1]. |
| INV-05 | `Start1 ≤ End1`, `Start2 ≤ End2`; coordinates lie within parent gene spans. | Coordinates derived from min/max parent gene positions. |
| INV-06 | Reported chains are non-overlapping (each anchor in ≤ 1 block). | Greedy single-pass chaining assigns each anchor once [1]. |

### 2.5 Comparison with Related Methods

| Aspect | This implementation | MCScanX reference [1] |
|--------|---------------------|-----------------------|
| Anchor generation | external (`orthologMap` input) | BLASTP + E-value cutoff + collapsing |
| Chaining | greedy single pass, MCScanX scoring/cutoffs | full DP over all anchor predecessors |
| Defaults | MatchScore 50, GapPenalty −1, maxGap 25, min 5 pairs | same defaults |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `genome1Genes` | `IReadOnlyList<Gene>` | required | Genome-1 genes in chromosomal order | non-null; 0-based positions |
| `genome2Genes` | `IReadOnlyList<Gene>` | required | Genome-2 genes in chromosomal order | non-null |
| `orthologMap` | `IReadOnlyDictionary<string,string>` | required | genome-1 gene id → genome-2 gene id anchors | non-null |
| `minAnchors` | `int` | 5 | Minimum collinear pairs per block | ≥ 1; MCScanX default 5 |
| `maxGap` | `int` | 25 | Max intervening genes between anchors | MCScanX `NumberofGaps < 25` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `SyntenicBlock.Start1/End1` | `int` | Genome-1 coordinate span (min..max parent gene positions) |
| `SyntenicBlock.Start2/End2` | `int` | Genome-2 coordinate span |
| `SyntenicBlock.IsInverted` | `bool` | True when genome-2 order is reversed |
| `SyntenicBlock.GeneCount` | `int` | Number of anchored gene pairs in the block |
| `SyntenicBlock.Identity` | `double` | Fixed 1.0 (anchor-level identity not recomputed here) |

### 3.3 Preconditions and Validation

Null `genome1Genes`, `genome2Genes`, or `orthologMap` throws `ArgumentNullException`. An empty
genome or an anchor set smaller than `minAnchors` returns an empty result (no exception).
Positions are 0-based gene indices; coordinates are taken from the `Gene.Start`/`Gene.End`
fields as supplied. Ortholog entries whose target id is absent from genome 2 are silently
skipped.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; return empty for empty genomes.
2. Index genome-2 gene positions; collect anchor pairs `(pos1, pos2)` from `orthologMap`, ordered by genome-1 position.
3. Walk anchors in genome-1 order, extending the current chain while genome-2 direction stays consistent and `NumberofGaps < maxGap`; otherwise flush and start a new chain.
4. On flush, report the chain iff MCScanX `score ≥ 250` and `anchorCount ≥ minAnchors`.
5. Build `SyntenicBlock`s with coordinate spans and the inverted flag.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Constant | Value | Source |
|----------|-------|--------|
| `MatchScore` | 50 | MCScanX default [1] |
| `GapPenalty` | −1 | MCScanX default [1] |
| `DefaultMaxGaps` | 25 | `NumberofGaps < 25` [1] |
| `MinChainScore` | 250 | "scores over 250" [1] |
| `DefaultMinAnchors` | 5 | ≥ 5 collinear pairs [1] |

Chain score: `n × MatchScore + GapPenalty × Σ NumberofGaps`, where `NumberofGaps` between two
anchors is `|Δpos2| − 1`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindSyntenicBlocks` | O(n) anchors after O(n) index, sort O(a log a) | O(n) | n = genes; a = anchors. Single-pass greedy chaining; classified O(n²) in the Registry because the full MCScanX DP is O(a²). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.FindSyntenicBlocks(...)`: detects collinear blocks from anchors.
- `ComparativeGenomics.VisualizeSynteny(blocks)`: renders one text line per block.

### 5.2 Current Behavior

Anchor generation is delegated to the caller via `orthologMap` (orthologs are produced by
COMPGEN-ORTHO-001). The chaining is a deterministic greedy single pass rather than the full
predecessor-DP of MCScanX, but applies the same scoring constants and report cutoffs, so for
direction-consistent, gap-bounded inputs it produces the same blocks. `Identity` is reported as
1.0 (the algorithm scores gene-order collinearity, not per-block sequence identity).

**Search reuse:** the repository suffix tree was evaluated and **not used**. This unit chains
pre-supplied anchor *positions* by gene order and score; it performs no substring/occurrence
search, so the suffix tree (`Contains`/`FindAllOccurrences`) does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- MCScanX scoring constants `MatchScore = 50`, `GapPenalty = −1`, `NumberofGaps < 25` [1].
- Report rule: chain score ≥ 250 and ≥ 5 collinear gene pairs [1].
- Forward and inverted (reverse) collinear blocks via genome-2 direction [1].
- Non-overlapping chain reporting [1].

**Intentionally simplified:**

- Chaining uses a greedy single pass instead of full predecessor dynamic programming; **consequence:** for inputs with interleaved competing chains the greedy partition may segment anchors differently than MCScanX's optimal DP path.
- Anchor identity is fixed at 1.0; **consequence:** users do not get per-block sequence identity from this method.

**Not implemented:**

- BLASTP anchor generation, E-value cutoff (1e-5), and match collapsing; **users should rely on:** `ComparativeGenomics.FindOrthologs` / `FindReciprocalBestHits` (COMPGEN-ORTHO-001) to build the `orthologMap`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Greedy chaining vs full DP | Deviation | Different segmentation on interleaved anchors | accepted | ASM-02; same constants/cutoffs |
| 2 | Anchors supplied externally | Assumption | Quality of blocks depends on ortholog quality | accepted | ASM-01; delegated to COMPGEN-ORTHO-001 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty genome list | empty result | no anchors possible |
| Empty ortholog map | empty result | no anchors [1] |
| < 5 anchors | empty result | below report threshold [1] |
| Anchor gap ≥ maxGap | chain breaks | `NumberofGaps < 25` [1] |
| Reversed genome-2 order | one inverted block | both directions [1] |
| Null required argument | `ArgumentNullException` | input validation |
| Ortholog target absent in genome 2 | anchor skipped | robustness |

### 6.2 Limitations

Single chromosome pair per call; no whole-paranome / multi-genome multi-alignment; no tandem
collapsing; per-block sequence identity is not computed; greedy chaining is not guaranteed to
reproduce MCScanX's optimal DP segmentation on adversarial interleaved inputs.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var blocks = ComparativeGenomics.FindSyntenicBlocks(
    genome1Genes, genome2Genes, orthologMap).ToList();
// 5 adjacent forward anchors -> 1 block: GeneCount=5, IsInverted=false, score=5*50=250
string text = ComparativeGenomics.VisualizeSynteny(blocks);
```

**Numerical walk-through:** five adjacent anchors at genome-2 positions 0,1,2,3,4 have
`Σ NumberofGaps = 0`, so score `= 5 × 50 + (−1) × 0 = 250 ≥ 250` and `5 ≥ 5` → one forward
block reported.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_FindSyntenicBlocks_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindSyntenicBlocks_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [COMPGEN-SYNTENY-001-Evidence.md](../../../docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md)

## 8. References

1. Wang Y, Tang H, DeBarry JD, Tan X, Li J, Wang X, Lee T-H, Jin H, Marler B, Guo H, Kissinger JC, Paterson AH. 2012. MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research* 40(7):e49. https://doi.org/10.1093/nar/gkr1293 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336)
2. Wang Y, et al. 2012. MCScanX (HTML). *Nucleic Acids Research* 40(7):e49. https://academic.oup.com/nar/article/40/7/e49/1202057
3. Wikipedia contributors. Synteny. https://en.wikipedia.org/wiki/Synteny (accessed 2026-06-13).
