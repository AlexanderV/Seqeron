# Structural Variant Detection (Paired-End Mapping Signatures)

| Field | Value |
|-------|-------|
| Algorithm Group | StructuralVar |
| Test Unit ID | SV-DETECT-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Detects germline/somatic structural variants (SVs) from paired-end mappings by recognising
discordant read-pair *signatures*. A read pair whose mates map farther apart than the library
insert size signals a deletion; closer than expected signals an insertion; mates on the same
strand signal an inversion; mates on different chromosomes signal a translocation [1]. The
detector flags anomalous pairs using an insert-size cutoff in standard deviations [2], clusters
nearby supporting pairs, and reports an SV for each cluster meeting a minimum read-pair support
[2]. It is a heuristic, signature-based method (not exhaustive assembly or read-depth modelling).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A paired-end library sequences both ends of DNA fragments of a known approximate length (the
*insert size*). Concordant pairs map to the reference with the expected separation and the
forward-reverse (FR) orientation — one mate on the '+' strand, the mate on '−', pointing inward
(SAM proper-pair FLAG 0x02) [4]. Structural rearrangements in the donor genome perturb this
mapping, producing *discordant* pairs whose pattern (the *signature*) identifies the event type
[1].

### 2.2 Core Model

For a library with mean insert size μ and standard deviation σ, and a cutoff c (in units of σ),
a read pair with observed insert size s and mate chromosomes (chr₁, chr₂) / strands (st₁, st₂) is
classified as follows [1][2][4]:

- **Discordant by span:** `s < μ − c·σ` OR `s > μ + c·σ` (bounds = μ ± c·σ) [2].
- **Deletion (DEL):** same chromosome, span larger than the insert size (`s > μ + c·σ`) — "the mapped distance is greater than the insert size" [1].
- **Insertion (INS):** same chromosome, span smaller than the insert size (`s < μ − c·σ`) — "if the event is an insertion, then the distance is smaller" [1].
- **Inversion (INV):** same chromosome, mates on the same strand (`st₁ == st₂`) — "the orientation of the read, lying within the inversion, flipped" [1]; FF/RR is the abnormal orientation that supports an inversion [4].
- **Translocation (CTX):** mates on different chromosomes (`chr₁ ≠ chr₂`) — a linking signature that "can connect regions ... on different chromosomes" [1]; BreakDancer code CTX [2].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The library insert size is approximately normal with known μ and σ, so μ ± c·σ bounds the concordant span range [2]. | A skewed/multimodal insert distribution mis-flags concordant pairs or misses true discordant pairs. |
| ASM-02 | Inter-chromosomal mapping is classified as Translocation *before* orientation/span are examined (chromosome difference takes precedence) [1]. | A same-strand cross-chromosome pair would otherwise be mislabelled an inversion; inversion is defined only intra-chromosomally. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `chr₁ ≠ chr₂` ⇒ Translocation, irrespective of span/orientation. | Linking signature across chromosomes [1]; ASM-02 precedence. |
| INV-02 | Same chromosome and `st₁ == st₂` ⇒ Inversion. | Flipped-orientation signature [1][4]. |
| INV-03 | Same chromosome, opposite strands, `s > μ + c·σ` ⇒ Deletion. | Span larger than insert size [1]. |
| INV-04 | Same chromosome, opposite strands, `s < μ − c·σ` ⇒ Insertion. | Span smaller than insert size [1]. |
| INV-05 | A pair is discordant-by-span iff `s < μ − c·σ` or `s > μ + c·σ`. | Bounds = μ ± c·σ [2]. |
| INV-06 | `DetectSVs` emits an SV for a cluster iff its supporting-pair count ≥ minSupport. | Minimum-support gate (BreakDancer -r, default 2) [2]. |

### 2.5 Comparison with Related Methods

| Aspect | PEM-signature detection (this) | Read-depth (CNV) detection |
|--------|-------------------------------|----------------------------|
| Signal | Discordant pair span/orientation | Coverage depth per window |
| Copy-neutral events (inversion, balanced translocation) | Detectable [1] | Not detectable |
| Breakpoint resolution | Pair-level (approximate) | Window-level |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| readPairs | `IEnumerable<(string,string,int,char,string,int,char,int)>` | required | Mapped pairs: ReadId, Chr1, Pos1, Strand1, Chr2, Pos2, Strand2, InsertSize | Strand ∈ {'+','−'}; positions 0-based |
| expectedInsertSize | int | 400 | Library mean insert size μ | > 0 |
| insertSizeStdDev | int | 50 | Library insert-size σ | ≥ 0 |
| cutoffSd | double | 3.0 | Anomaly cutoff c in units of σ [2] | ≥ 0 |
| clusterDistance | int | 500 | Max coordinate gap to keep adjacent pairs in one cluster | ≥ 0 |
| minSupport | int | 2 | Minimum supporting pairs to call an SV [2] | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ClassifySV` | `SVType` | One of Deletion / Insertion / Inversion / Translocation / ComplexRearrangement |
| `DetectSVs` | `IEnumerable<StructuralVariant>` | One SV per qualifying cluster (Type, Start, End, Length, SupportingReads, Quality) |
| `FindDiscordantPairs` | `IEnumerable<ReadPairSignature>` | The pairs flagged anomalous |

### 3.3 Preconditions and Validation

Null `readPairs` / `discordantPairs` throws `ArgumentNullException`. Empty input yields an empty
result. Coordinates are treated as 0-based; strands use '+'/'−'. No normalization of chromosome
names is performed (string equality). The classifier assumes the input pair is already anomalous;
a concordant pair passed directly to `ClassifySV` is classified by the same rules but would not
normally be produced by `FindDiscordantPairs`.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Flag discordant pairs** (`FindDiscordantPairs`): mark a pair anomalous if interchromosomal, span outside μ ± c·σ, or non-FR/RF orientation [1][2][4].
2. **Cluster** (`ClusterDiscordantPairs`): sort by (chr, position) and group pairs whose coordinates lie within `clusterDistance` of the previous pair.
3. **Support gate:** discard clusters with fewer than `minSupport` pairs [2].
4. **Classify** (`ClassifySV`): assign the cluster's SV type from the representative pair's PEM signature [1].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Classification order (first match wins), per §2.2: (1) `chr₁ ≠ chr₂` → Translocation; (2)
`st₁ == st₂` → Inversion; (3) `s > μ + c·σ` → Deletion; (4) `s < μ − c·σ` → Insertion; (5)
otherwise → ComplexRearrangement. Numeric parameters: cutoff c default 3 σ [2]; min support
default 2 [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindDiscordantPairs | O(n) | O(1) streaming | Per-pair classification is O(1). |
| ClusterDiscordantPairs / DetectSVs | O(n log n) | O(n) | Dominated by the sort of discordant pairs. |
| ClassifySV | O(1) | O(1) | Constant-time signature test. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [StructuralVariantAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs)

- `StructuralVariantAnalyzer.DetectSVs(...)`: canonical entry point — find discordant → cluster → classify.
- `StructuralVariantAnalyzer.ClassifySV(...)`: maps a read-pair signature to an `SVType`.
- `StructuralVariantAnalyzer.FindDiscordantPairs(...)`: anomaly detection using the μ ± c·σ cutoff and FR/RF orientation test.

### 5.2 Current Behavior

The detector is single-pass and signature-based. `FindDiscordantPairs` additionally flags any pair
above a hard `maxInsertSize` guard (default 10000); such a pair with no matching basic signature is
classified `ComplexRearrangement`. Clustering is a simple linear sweep over the sorted pairs rather
than a windowing/CBS model. **Search reuse:** this unit performs no substring/pattern search over a
sequence (it operates on already-mapped coordinate records), so the repository suffix tree is not
applicable and was not used.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Deletion = larger span, Insertion = smaller span, Inversion = flipped (same-strand) orientation, Translocation = interchromosomal linking signature [1].
- Discordant-by-span cutoff bounds μ ± c·σ with default c = 3 [2].
- Minimum supporting read pairs default 2 [2].
- Concordant FR/RF (opposite-strand) orientation; FF/RR abnormal [4].

**Intentionally simplified:**

- Clustering: linear adjacency sweep; **consequence:** no statistical confidence model or windowing as in BreakDancer's connection scoring — cluster membership depends only on `clusterDistance`.
- Tandem duplication / everted-duplication and complex linking signatures are folded into `ComplexRearrangement` rather than separately resolved; **consequence:** duplications from PEM are not distinguished here.

**Not implemented:**

- Insertions larger than the fragment insert size (no span signature) [1]; **users should rely on:** split-read / assembly-based detection (out of scope for this unit).
- Inserted-sequence recovery and breakpoint base-resolution refinement; **users should rely on:** the split-read methods in this class (`FindSplitReads`, `ClusterSplitReads`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Inter-chromosomal precedence over orientation | Assumption | Determines Translocation vs Inversion for cross-chromosome same-strand pairs | accepted | ASM-02; inversion defined intra-chromosomally only [1] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty `readPairs` | Empty result | No pairs to classify |
| Null input | `ArgumentNullException` | Input validation contract |
| Span exactly μ ± c·σ | Concordant (not discordant) | Bound is inclusive; discordant iff strictly outside [2] |
| Cluster below minSupport | No SV emitted | Support gate [2] |
| Cross-chromosome, same strand | Translocation | ASM-02 precedence [1] |

### 6.2 Limitations

PEM signatures cannot see insertions longer than the fragment, do not recover inserted sequence,
and give only approximate (pair-level) breakpoints [1]. The normal-insert-size assumption (ASM-01)
degrades on skewed libraries. Tandem duplications are not separately classified.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var pairs = new[]
{
    // three concordant deletion-signature pairs (same chr, FR, span 5000 > 400 + 3*50)
    ("r1","chr1",1000,'+',"chr1",6000,'-',5000),
    ("r2","chr1",1010,'+',"chr1",6010,'-',5000),
    ("r3","chr1",1020,'+',"chr1",6020,'-',5000),
};

var svs = StructuralVariantAnalyzer.DetectSVs(
    pairs, expectedInsertSize: 400, insertSizeStdDev: 50,
    cutoffSd: 3.0, clusterDistance: 500, minSupport: 2).ToList();
// svs[0].Type == SVType.Deletion, svs[0].SupportingReads == 3
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [StructuralVariantAnalyzer_DetectSVs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_DetectSVs_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [SV-DETECT-001-Evidence.md](../../../docs/Evidence/SV-DETECT-001-Evidence.md)

## 8. References

1. Medvedev P, Stanciu M, Brudno M. 2009. Computational methods for discovering structural variation with next-generation sequencing. Nature Methods 6(11s):S13–S20. https://doi.org/10.1038/nmeth.1374
2. Chen K, Wallis JW, McLellan MD, et al. 2009. BreakDancer: an algorithm for high-resolution mapping of genomic structural variation. Nature Methods 6:677–681. https://doi.org/10.1038/nmeth.1363 (distribution README: https://raw.githubusercontent.com/genome/breakdancer/master/README)
3. Fan X, Abbott TE, Larson D, Chen K. 2014. BreakDancer: Identification of Genomic Structural Variation from Paired-End Read Mapping. Curr Protoc Bioinformatics 45:15.6.1–15.6.11. https://pmc.ncbi.nlm.nih.gov/articles/PMC3661775/
4. Kennedy A. 2012. Forward and reverse reads in paired-end sequencing (SAM proper-pair FLAG 0x02; BWA FR convention). https://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
