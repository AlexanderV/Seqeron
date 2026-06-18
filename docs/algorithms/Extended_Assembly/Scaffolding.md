# Scaffolding (Joining Ordered Contigs with N-Gaps)

| Field | Value |
|-------|-------|
| Algorithm Group | Extended Assembly |
| Test Unit ID | ASSEMBLY-SCAFFOLD-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Scaffolding orders and orients a set of contigs into longer sequences (scaffolds) using
paired-end / mate-pair link information, bridging the unknown stretches between contigs with runs
of the gap character `N` [1][4]. This unit consumes pre-computed links — each a triple
`(contig1, contig2, gapSize)` giving an estimated inter-contig distance — and emits the scaffolds:
contigs along a link path are concatenated, interspersed with a run of `N` whose length is the
distance estimate [1]. It is deterministic and specification-driven; the upstream maximum-likelihood
gap estimator [1][3] is out of scope. It should be used to assemble final scaffold strings once
link order and gap estimates are known.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

After contigs are assembled, scaffolding "link[s] together a non-contiguous series of genomic
sequences into a scaffold, consisting of sequences separated by gaps of known length" [4]. Links
are derived from read pairs that map to two different contigs; greedy scaffolders (e.g. Bambus)
join contigs with the most links first [4]. The gap between two contigs is the unknown distance `d`
separating them in the genome [3].

### 2.2 Core Model

For a path of ordered contigs `c0, c1, …, ck`, the scaffold is the concatenation of the contig
sequences, "interspersed with gaps represented by a run of the character N, whose length
corresponds to the estimate of the distance between those two contigs" [1]. Thus a link
`(i, j, g)` appends `contigs[j]` after `contigs[i]` preceded by `g` copies of the gap character
when `g > 0` [1].

A distance estimate may be negative, "indicating that the two contigs should in fact overlap" [1];
this case is frequent because a de Bruijn assembler splits contigs at a node leaving an overlap of
one k-mer length [3]. When the overlap is not resolved, no positive run of `N` represents it. The
governing file-format convention requires gap lengths to be positive — "Negative gaps and gap lines
with zero length are not valid" — and for "negative gaps, or gaps of unknown size, use … 100 as the
gap size, since 100 is the GenBank/EMBL/DDBJ standard for gaps of unknown size" [2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A `gapSize > 0` link emits a contiguous run of exactly `gapSize` gap characters between the two contigs | Run length = distance estimate [1] |
| INV-02 | A non-positive estimate (`gapSize ≤ 0`) emits exactly 100 gap characters | GenBank/EMBL/DDBJ unknown-gap length [2] |
| INV-03 | Scaffold length over a followed path = Σ(contig lengths) + Σ(emitted gap lengths) | Concatenation + gap runs [1] |
| INV-04 | Each contig appears in exactly one scaffold | Scaffold is a path of distinct contigs [4] |
| INV-05 | Contig substrings are preserved verbatim and in link order within a scaffold | "sequences … are concatenated" [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `contigs` | `IReadOnlyList<string>` | required | Contigs indexed by position | non-null; indices referenced by links |
| `links` | `IReadOnlyList<(int contig1, int contig2, int gapSize)>` | required | Ordering links with gap estimates (characters) | non-null; out-of-range / self indices ignored |
| `gapCharacter` | `char` | `'N'` | Character used to fill gaps | any char |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IReadOnlyList<string>` | One string per scaffold; contigs joined in link order with gap-character runs between them |

### 3.3 Preconditions and Validation

Null `contigs` or `links` raise `ArgumentNullException`. An empty `contigs` list returns an empty
result. Link indices are 0-based positions into `contigs`; a link whose `contig1` or `contig2` is
out of range, or whose endpoints are equal, is ignored. A contig is placed into at most one
scaffold; a link to an already-placed contig is skipped. Contigs not reached by any followed link
each become their own single-contig scaffold, emitted in ascending index order. Case is preserved;
no alphabet normalization is performed.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `contigs`/`links` non-null; empty contigs → empty result.
2. Index the forward links out of each contig (in-range, distinct endpoints), preserving input order.
3. For each contig not yet placed, in ascending index order, start a new scaffold with that contig.
4. Follow the link path: append the first unplaced successor preceded by a gap run of `gapSize`
   characters (`gapSize > 0`) or 100 characters (`gapSize ≤ 0`); mark it placed; repeat.
5. Emit the assembled scaffold string.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Gap length rule: `gapLength = gapSize > 0 ? gapSize : 100`. The constant 100 is the
  GenBank/EMBL/DDBJ standard length for a gap of unknown size [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Scaffold` | O(n + L + S) | O(n + L + S) | n = contigs, L = links (indexed once), S = total output length including gaps |

This is below O(n²); no property-based-test / benchmark requirement applies (Definition of Done
mandates those only for O(n²)+ algorithms). The dominant term is producing the output string `S`.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.Scaffold(contigs, links, gapCharacter)`: joins ordered contigs into scaffolds with gap-character runs.

### 5.2 Current Behavior

Links are indexed into a dictionary keyed by `contig1`, preserving input order so the first declared
forward link wins on ties. The path is followed greedily from each unplaced contig; the first
successor not yet placed is appended. This is a search/matching-adjacent unit only in the loose
sense of "find the next link"; it performs **no substring search** over sequence data (it indexes by
integer contig id), so the repository suffix tree is **not applicable** and is not used. Overlap
resolution for negative gap estimates is not performed; such estimates fall back to the unknown-gap
placeholder.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Contigs along a link path are concatenated, interspersed with a run of the gap character of length
  equal to the (positive) distance estimate [1].
- Non-positive / unknown gaps emit the GenBank/EMBL/DDBJ standard length of 100 characters [2].
- Each contig appears in exactly one scaffold; unlinked contigs become length-1 scaffolds [4].

**Intentionally simplified:**

- Negative-estimate overlap resolution: a negative estimate indicates the contigs overlap [1][3];
  **consequence:** instead of merging the contigs, a 100-character placeholder gap is emitted, so an
  overlapping pair is not collapsed here (use `MergeContigs` for known overlaps).
- Gap-distance estimation: the maximum-likelihood / EM estimator [1][3] is upstream; **consequence:**
  the caller supplies `gapSize`; this unit does not derive it from read alignments.

**Not implemented:**

- Contig orientation (reverse-complementing contigs on `-` strand links); **users should rely on:**
  pre-oriented contigs supplied in the correct strand.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Non-positive estimate → 100-N placeholder instead of overlap merge | Assumption | Overlapping pairs not collapsed; gap length is the AGP unknown-size standard, not the true (negative) distance | accepted | Constant source-backed [2]; overlap resolution out of scope (use `MergeContigs`) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty contigs | Empty result | No contigs → no scaffolds |
| No links | One single-contig scaffold per contig | Length-1 paths [4] |
| `gapSize = 0` | 100 gap characters | Zero-length gap invalid → unknown default [2] |
| `gapSize < 0` | 100 gap characters | Negative = overlap, unresolved → unknown default [1][2] |
| Out-of-range / self link | Ignored | Indices must reference distinct existing contigs |
| Link to placed contig | Skipped | A contig cannot appear twice (INV-04) |
| Null `contigs` / `links` | `ArgumentNullException` | Input validation |

### 6.2 Limitations

No contig orientation/strand handling; no overlap resolution for negative gaps; gap estimation is
the caller's responsibility. Cycles in links terminate naturally because a contig is never placed
twice. Output depends on link input order for tie-breaking.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var contigs = new[] { "ACGT", "TTGG", "CCAA" };
var links = new[] { (0, 1, 3), (1, 2, 2) };
IReadOnlyList<string> scaffolds = SequenceAssembler.Scaffold(contigs, links);
// scaffolds[0] == "ACGTNNNTTGGNNCCAA"  (4 + 3 + 4 + 2 + 4 = 17 chars)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_Scaffold_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_Scaffold_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ASSEMBLY-SCAFFOLD-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md)
- Related algorithms: [Contig_Merging](Contig_Merging.md)

## 8. References

1. Jackman SD, Vandervalk BP, Mohamadi H, Chu J, Yeo S, Hammond SA, Jahesh G, Khan H, Coombe L, Warren RL, Birol I. 2017. ABySS 2.0: resource-efficient assembly of large genomes using a Bloom filter. *Genome Research* 27:768–777. https://genome.cshlp.org/content/27/5/768
2. NCBI. AGP Specification v2.1 — Accessioned Golden Path file format. National Center for Biotechnology Information. https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
3. Sahlin K, Street N, Lundeberg J, Arvestad L. 2012. Improved gap size estimation for scaffolding algorithms. *Bioinformatics* 28(17):2215–2222. https://academic.oup.com/bioinformatics/article/28/17/2215/246308
4. Pop M, Kosack DS, Salzberg SL. 2004. Hierarchical Scaffolding With Bambus. *Genome Research* 14(1):149–159. https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics)
