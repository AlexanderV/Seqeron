# Scaffolding (Joining Ordered Contigs with N-Gaps)

| Field | Value |
|-------|-------|
| Algorithm Group | Extended Assembly |
| Test Unit ID | ASSEMBLY-SCAFFOLD-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Implemented (greedy end-to-end path merging; see §5.3) |
| Last Reviewed | 2026-11-25 |

## 1. Overview

Scaffolding orders and orients a set of contigs into longer sequences (scaffolds) using
paired-end / mate-pair link information, bridging the unknown stretches between contigs with runs
of the gap character `N` [1][4]. This unit consumes pre-computed links — `(contig1, contig2, gapSize)`
triples, or oriented, weighted `ScaffoldLink`s — selects a consistent set of them with the greedy
path-merging rule [5], and emits the scaffolds: oriented contigs along a path are concatenated,
interspersed with a run of `N` whose length is the distance estimate [1], and a negative estimate
whose overlap is confirmed by the sequences is collapsed [1][6]. It is deterministic and specification-driven; the upstream maximum-likelihood
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
one k-mer length [3]. "If such an overlap is indeed found in the contig overlap graph, the two
contigs are merged" [1]; ABySS `MergeContigs` merges when the last −d bases of the scaffold so far
and the first −d bases of the next contig agree (equal ignoring case, or one is `N`) [6]. When the
overlap is not found, no positive run of `N` represents it. The
governing file-format convention requires gap lengths to be positive — "Negative gaps and gap lines
with zero length are not valid" — and for "negative gaps, or gaps of unknown size, use … 100 as the
gap size, since 100 is the GenBank/EMBL/DDBJ standard for gaps of unknown size" [2].

### 2.3 Greedy path merging and orientation

The contig-mate-pair graph has two nodes per contig (its two ends) joined by a contig edge; a link
(mate-edge) joins one end of each of two contigs and carries a weight (number of mate pairs) and a
length (the gap; negative = overlap) [5 §2]. Huson, Reinert & Myers' Algorithm 5.1 [5 §5]: select
all contig edges; "for each mate-edge e in descending order of weight", if its end nodes lie on
different selected paths P1 ≠ P2, merge them, maintaining "the invariant that every node is
adjacent to at most two selected edges". Bambus likewise "joins together contigs with the most links
first" [4]. A selected path is a scaffold: an ordering *and orientation* of its contigs. A link
`A(o1) → B(o2), d` (ABySS notation `A+ B-,d`; AGP orientation `+`/`-`) joins the exit end of A
(3′ if `+`, 5′ if `-`) to the entry end of B (5′ if `+`, 3′ if `-`); read from the other strand it is
`B(¬o2) → A(¬o1), d`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A `gapSize > 0` link emits a contiguous run of exactly `gapSize` gap characters between the two contigs | Run length = distance estimate [1] |
| INV-02 | `gapSize = 0`, or `gapSize < 0` whose −gapSize overlap does not agree (or exceeds a contig / the bases since the last gap), emits exactly 100 gap characters; an agreeing overlap is collapsed | [2]; [1][6] |
| INV-03 | Scaffold length over a path = Σ(contig lengths) + Σ(emitted gap lengths) − Σ(collapsed overlaps) | Concatenation + gap runs [1], overlap merge [6] |
| INV-04 | Each contig appears in exactly one scaffold, in one orientation | Selected paths cover the nodes disjointly [5] |
| INV-05 | Each contig (or its reverse complement, for `-`) is preserved verbatim, in path order | "sequences … are concatenated" [1] |
| INV-06 | Every contig end carries at most one selected link; selected links form no cycle | Algorithm 5.1 invariant, P1 ≠ P2 [5] |
| INV-07 | A heavier link is considered before a lighter one; equal weights in input order | "descending order of weight" [5]; most links first [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `contigs` | `IReadOnlyList<string>` | required | Contigs indexed by position | non-null; indices referenced by links |
| `links` | `IReadOnlyList<(int contig1, int contig2, int gapSize)>` | required | Forward links (both `+`, weight 1) | non-null; out-of-range / self indices ignored |
| `links` (overload) | `IReadOnlyList<ScaffoldLink>` | required | `ScaffoldLink(Contig1, Contig2, GapSize, Contig1Reversed=false, Contig2Reversed=false, Weight=1)` | non-null; `Weight ≥ 1` else `ArgumentOutOfRangeException` |
| `gapCharacter` | `char` | `'N'` | Character used to fill gaps | any char |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IReadOnlyList<string>` | One string per scaffold, listed by the scaffold's lowest contig index |

### 3.3 Preconditions and Validation

Null `contigs` or `links` raise `ArgumentNullException`. An empty `contigs` list returns an empty
result. Link indices are 0-based positions into `contigs`; a link whose `contig1` or `contig2` is
out of range, or whose endpoints are equal, is ignored. Contigs without a selected link are
single-contig scaffolds. Case is preserved; reverse complements use the Core IUPAC complement.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `contigs`/`links` non-null and weights ≥ 1; empty contigs → empty result.
2. Sort links by descending weight (stable: input order breaks ties) [5].
3. For each link (in-range, distinct contigs): let `x` = exit end of contig1, `y` = entry end of
   contig2. Select it iff `x` and `y` carry no selected link and the contigs are on different paths
   (union–find); record the gap on the join. Otherwise it is a conflict (branch) or a cycle and is
   dropped.
4. For each path (by lowest contig index): read it from a terminal contig (free end first), each
   contig in the orientation the joins imply; keep the reading with more forward contigs (tie: the
   lowest-index contig forward).
5. Spell it: `gap > 0` → `gap` gap characters; `gap < 0` with an agreeing −gap overlap → collapse
   via `MergeContigs(tail, next, −gap, resolveOverlapConflicts: true)`; otherwise 100 gap characters.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Selection: descending weight, end nodes free, no cycle (Algorithm 5.1 restricted to end-to-end
  merges) [5].
- Orientation of the reported strand: majority of contigs forward (SSPACE reverses a scaffold that
  is wholly reverse-complemented [7]); tie → lowest-index contig forward. Forward-only links always
  yield forward scaffolds.
- Overlap: agreement = equal ignoring case, or either base `N` (ABySS `createConsensus`) [6]; the
  overlap may not exceed the next contig or the bases placed since the previous gap.
- Gap length: positive estimate verbatim [1]; otherwise 100 — the GenBank/EMBL/DDBJ standard length
  for a gap of unknown size, also mandated for negative gaps [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Scaffold` | O(n + L log L + S) | O(n + L + S) | n = contigs, L = links (sorted once), S = total output length including gaps |

This is below O(n²). The dominant term is producing the output string `S`.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.Scaffold(contigs, IReadOnlyList<(int,int,int)> links, gapCharacter)`: forward links, weight 1 (delegates).
- `SequenceAssembler.Scaffold(contigs, IReadOnlyList<ScaffoldLink> links, gapCharacter)`: oriented, weighted links.
- MCP `scaffold_contigs` (Alignment server) delegates; `ContigLinkInput` carries the optional orientation flags and weight.

### 5.2 Current Behavior

See §4. Reverse complement: Core `SequenceExtensions.TryGetReverseComplement`; overlap collapse:
`SequenceAssembler.MergeContigs`. No substring search over sequence data is performed (links are
integer contig ids), so the suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented:** greedy selection in descending weight with the ≤ 2-selected-edges-per-node and
P1 ≠ P2 rules [5]; orientation via bidirected end nodes [5][6]; gap runs of the estimated length [1];
overlap collapse when the −d bases agree [1][6]; 100-character unknown gap otherwise [2].

**Not implemented (documented):**

- Path interleaving ("zippering" P1 into P2 when a link hits a path interior) and the happy/unhappy
  weight test of Algorithm 5.1 [5]: both need each link's standard deviation σ (happy ⇔ |Δ| ≤ 3σ),
  which the input does not carry; such a link is treated as a conflict and not used.
- Edge bundling / transitive reduction [5 §2.2–2.3] and SSPACE's link-ratio repeat filter [7]: links
  are taken as already-bundled contig-level edges with a weight.
- ABySS's fallback gapped overlap alignment and IUPAC-subset consensus [6]: a non-agreeing overlap
  gets the 100-N unknown gap.
- Gap-distance estimation (ML / EM) [1][3] is upstream; the caller supplies `gapSize`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Unresolved negative / zero estimate → 100 N | Spec choice | AGP v2.1 mandates 100 for negative/unknown gaps [2]; ABySS and SSPACE emit a single `n` | accepted | AGP is the submission format |
| 2 | Orientation of the reported strand (majority forward, tie lowest-index forward) | Convention | Either strand is the same scaffold | accepted | SSPACE reverses wholly-reversed scaffolds [7] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty contigs | Empty result | No contigs → no scaffolds |
| No links | One single-contig scaffold per contig | Length-1 paths |
| Link into a lower-index contig, e.g. `(1,0,2)` over `AA,TT` | `TTNNAA` | Paths are read from their terminal, not from index 0 |
| Two successors at one end `(0,1)`, `(0,2)` | Heavier (else first) wins; the other contig is a singleton | ≤ 1 link per end [5] |
| Cycle `0→1→2→0` | Cycle-closing link dropped | P1 ≠ P2 [5] |
| `gapSize = 0` | 100 gap characters | Zero-length gap invalid [2] |
| `gapSize < 0`, overlap agrees | Collapsed (`ACGTAC`+`TACGGA`, −3 → `ACGTACGGA`) | [1][6] |
| `gapSize < 0`, overlap disagrees or too long | 100 gap characters | [2] |
| Out-of-range / self link | Ignored | Indices must reference distinct existing contigs |
| `Weight < 1` | `ArgumentOutOfRangeException` | Weight = number of mate pairs |
| Null `contigs` / `links` | `ArgumentNullException` | Input validation |

### 6.2 Limitations

See §5.3 "Not implemented". Gap estimation is the caller's responsibility.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var contigs = new[] { "ACGT", "TTGG", "CCAA" };
var links = new[] { (0, 1, 3), (1, 2, 2) };
IReadOnlyList<string> scaffolds = SequenceAssembler.Scaffold(contigs, links);
// scaffolds[0] == "ACGTNNNTTGGNNCCAA"  (4 + 3 + 4 + 2 + 4 = 17 chars)

var oriented = new[]
{
    new SequenceAssembler.ScaffoldLink(0, 1, 2),
    new SequenceAssembler.ScaffoldLink(0, 2, 3, Contig2Reversed: true, Weight: 4),
};
SequenceAssembler.Scaffold(new[] { "AAAA", "CCCC", "GGTT" }, oriented);
// ["AAAANNNAACC", "CCCC"]: the weight-4 link takes contig 0's 3' end; 0→1 is a conflict.
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_Scaffold_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAssembler_Scaffold_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ASSEMBLY-SCAFFOLD-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md)
- Related algorithms: [Contig_Merging](Contig_Merging.md)

## 8. References

1. Jackman SD, Vandervalk BP, Mohamadi H, Chu J, Yeo S, Hammond SA, Jahesh G, Khan H, Coombe L, Warren RL, Birol I. 2017. ABySS 2.0: resource-efficient assembly of large genomes using a Bloom filter. *Genome Research* 27:768–777. https://genome.cshlp.org/content/27/5/768
2. NCBI. AGP Specification v2.1 — Accessioned Golden Path file format. National Center for Biotechnology Information. https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
3. Sahlin K, Street N, Lundeberg J, Arvestad L. 2012. Improved gap size estimation for scaffolding algorithms. *Bioinformatics* 28(17):2215–2222. https://academic.oup.com/bioinformatics/article/28/17/2215/246308
4. Pop M, Kosack DS, Salzberg SL. 2004. Hierarchical Scaffolding With Bambus. *Genome Research* 14(1):149–159. https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics)
5. Huson DH, Reinert K, Myers EW. 2002. The greedy path-merging algorithm for contig scaffolding. *Journal of the ACM* 49(5):603–615. doi:10.1145/585265.585267 (PDF: https://publications.mpi-cbg.de/Huson_2002_6349.pdf)
6. ABySS source, `MergePaths/MergeContigs.cpp` (`mergeContigs`, `createConsensus`). https://github.com/bcgsc/abyss
7. Boetzer M, Henkel CV, Jansen HJ, Butler D, Pirovano W. 2011. Scaffolding pre-assembled contigs using SSPACE. *Bioinformatics* 27(4):578–579; SSPACE-Basic source `SSPACE_Basic.pl`, `bin/PairingAndScaffolding.pl` (https://github.com/nsoranzo/sspace_basic).
