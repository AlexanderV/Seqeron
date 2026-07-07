# Overlap-Layout-Consensus (OLC) Assembly

> **Baseline / reference method.** OLC is an important classical assembler paradigm; de Bruijn approaches are often more practical for many read regimes. See [Legacy / Baseline Methods](../CANONICAL_MAP.md).

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-OLC-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Overlap-Layout-Consensus (OLC) is a *de novo* genome-assembly paradigm that reconstructs
long contiguous sequences (contigs) from many shorter reads. It proceeds in three stages:
compute pairwise suffix-prefix **overlaps** between reads, arrange the reads into a
**layout** (a path through the resulting overlap graph), and read off a **consensus**
sequence for each layout [1][2]. OLC is a heuristic: exact layout is equivalent to
finding a Hamiltonian path through the overlap graph, which is NP-complete [1], so
practical implementations (including this one) use greedy overlap chaining. It is exact
for unambiguous (non-repeat) read tilings and degrades — splitting or mis-resolving — on
repeats longer than the read length [3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Shotgun sequencing yields many reads sampled from unknown positions of a genome. When a
suffix of one read matches a prefix of another, the two reads likely overlap in the
genome ("first law of assembly") [3]. Assembly stitches such overlapping reads back into
the original sequence.

### 2.2 Core Model

For reads `A` and `B`, the overlap `overlap(A → B)` is the length of the **longest**
suffix of `A` that equals a prefix of `B`, provided that length is ≥ a minimum threshold
`l` [2 (p.5, p.10)][3 (p.16)]. The **overlap graph** has one node per read and a directed
edge `A → B` weighted by `overlap(A → B)` whenever that overlap is ≥ `l` [1][2 (p.20)][3 (p.23)].

The **layout** is a path through the overlap graph visiting reads in genome order; in the
ideal noise-free case it is a Hamiltonian path, and traversing it while merging
overlapping reads reconstructs the genome [1]. Because finding a Hamiltonian path is
NP-complete [1], real assemblers (and this implementation) approximate the layout
greedily and via transitive-edge reduction, emitting one contig per non-branching stretch
[2 (p.21–25)]. The **consensus** assigns each column of the read pile-up its
majority-vote nucleotide [2 (p.28)].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Reads are error-free for the exact-overlap path (identity 1.0) | With sequencing error, true overlaps fall below the identity threshold or spurious branches appear; `minIdentity` < 1.0 admits mismatches but can also admit false overlaps [2 (p.11–15, p.26)] |
| ASM-02 | Repeats are shorter than the read length | Repeats ≥ read length collapse; the layout cannot determine which in-paths pair with which out-paths, so contigs split or mis-join [3 (p.58–62)] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The overlap graph contains no self-edge (`ReadIndex1 ≠ ReadIndex2`) | A read is not overlapped against itself [1][2 (p.5)] |
| INV-02 | Every reported overlap length is ≥ `minOverlap` and ≤ min(len(A), len(B)) | Threshold `l` and the suffix-prefix definition [2 (p.5, p.10)] |
| INV-03 | The reported overlap for an ordered pair is the single longest qualifying suffix-prefix match | "report only the longest suffix/prefix match" [2 (p.10)] |
| INV-04 | An assembled contig is a superstring of its constituent reads; its length ≤ Σ read lengths and ≥ the longest single read | Merging along a layout path concatenates non-overlapping suffixes [1][3 (p.26)] |
| INV-05 | When the overlap graph is edgeless, each read becomes its own contig | No edge ⇒ no chaining; reads pass through as singletons [1] |

### 2.5 Comparison with Related Methods

| Aspect | OLC (this algorithm) | de Bruijn graph (DBG) |
|--------|----------------------|------------------------|
| Graph node | one read | one (k−1)-mer |
| Layout target | Hamiltonian path (NP-complete) | Eulerian path (polynomial) [1] |
| Best fit | longer reads, fewer reads | very many short reads [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `reads` | `IReadOnlyList<string>` | required | Sequence reads (overlap-graph nodes) | null/empty → empty result |
| `parameters.MinOverlap` | `int` | 20 | Minimum overlap length `l` | ≥ 1; longer = stricter |
| `parameters.MinIdentity` | `double` | 0.9 | Minimum fractional identity over the overlap window | 0.0–1.0; 1.0 = exact match |
| `parameters.MinContigLength` | `int` | 100 | Contigs shorter than this are discarded | ≥ 0 |
| `minOverlap` (FindAllOverlaps) | `int` | 20 | Edge threshold `l` | ≥ 1 |
| `minIdentity` (FindAllOverlaps) | `double` | 0.9 | Edge identity threshold | 0.0–1.0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AssemblyResult.Contigs` | `IReadOnlyList<string>` | Assembled contigs (after the min-length filter) |
| `AssemblyResult.TotalReads` | `int` | Number of input reads |
| `AssemblyResult.N50` | `double` | N50 of the contig length distribution |
| `AssemblyResult.LongestContig` | `int` | Longest contig length |
| `AssemblyResult.TotalLength` | `int` | Sum of contig lengths |
| `FindAllOverlaps` → `Overlap` | record | Edge `(ReadIndex1 → ReadIndex2, OverlapLength, Position1, Position2)` |

### 3.3 Preconditions and Validation

Reads are 0-indexed. `AssembleOLC(null)` and `AssembleOLC([])` return an empty
`AssemblyResult` (no exception). Overlap matching is **case-insensitive** (identity is
computed on upper-cased characters). The alphabet is unconstrained (any `char` strings;
DNA/RNA/protein all work). `Position2` of an overlap is always 0 (prefix of the second
read). Identity is the fraction of matching positions over the overlap window; an overlap
is accepted only if identity ≥ `minIdentity`.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Overlap:** for every ordered pair (i, j), i ≠ j, find the longest suffix-of-read[i] /
   prefix-of-read[j] match of length ≥ `MinOverlap` and identity ≥ `MinIdentity`; record it
   as a weighted edge (`FindAllOverlaps` / `FindOverlap`).
2. **Layout:** sort edges by descending overlap length; assign each read its best (longest)
   successor; reads with no incoming best-edge are chain starts. Walk each chain, merging
   `successor[overlap:]` onto the growing contig (greedy chaining; `BuildContigsFromOverlaps`).
3. **Consensus / emit:** the merged superstring of each chain is the contig; unused reads
   become singleton contigs. Discard contigs shorter than `MinContigLength`, then compute
   statistics (N50, total length, longest).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Edge weight = overlap length (no scoring matrix; no penalties). The only thresholds are
  the caller-supplied `MinOverlap` (`l`) and `MinIdentity` [2 (p.5, p.10)].
- Layout tie-break: edges processed in descending overlap length; the first best successor
  recorded per read wins (greedy).
- No biological constants are involved; there are no fixed numeric tables in the OLC path.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindAllOverlaps` (all-pairs) | O(N²) = O(d²n²) | O(a) | d reads of length n, N = dn, a = #edges; matches Langmead OLC p.16 |
| `FindOverlap` (one pair) | O(n²) | O(1) | scans overlap lengths from longest down, O(n) identity each |
| `AssembleOLC` | O(N² + a log a) | O(N) | dominated by overlap detection; layout sorts a edges |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.AssembleOLC(reads, parameters)`: full OLC pipeline (overlap → greedy layout → chain merge → stats).
- `SequenceAssembler.FindAllOverlaps(reads, minOverlap, minIdentity)`: builds the overlap graph (edges).
- `SequenceAssembler.FindAllOverlaps(reads, minOverlap, minIdentity, CancellationToken, IProgress)`: cancellable/progress-reporting variant of the above (delegates to the same `FindOverlap`).
- `SequenceAssembler.FindOverlap(seq1, seq2, minOverlap, minIdentity)`: longest suffix-prefix overlap for one ordered pair (internal building block).

### 5.2 Current Behavior

- **Suffix-tree reuse evaluated, not used.** Overlap detection is a "find suffix-prefix
  overlaps" operation for which a generalized suffix tree gives O(N + a) [2 (p.6–10)].
  However, this implementation supports approximate overlaps via `MinIdentity` < 1.0
  (mismatches), which the exact-match suffix tree cannot model; Langmead notes that the
  exact suffix-tree method is faster but the all-pairs scan "is more flexible, allowing
  mismatches" [2 (p.16)]. The all-pairs O(N²) scan is therefore used. A future exact-match
  fast path could front the suffix tree for the `MinIdentity == 1.0` case.
- The layout is a **greedy best-successor chaining** rather than full transitive reduction
  + Hamiltonian search; this resolves unambiguous chains exactly but does not optimally
  resolve repeats (see 5.3).
- Overlap matching is case-insensitive; `Position2` is always 0.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Suffix-prefix overlap definition with minimum length `l`, reporting only the longest match per ordered pair [2 (p.5, p.10)].
- Overlap graph as a directed graph with one node per read and edges weighted by overlap length [1][2 (p.20)][3 (p.23)].
- Three-stage Overlap → Layout → Consensus structure; singleton reads with no edge pass through as their own contigs [1][2 (p.4, p.25)].

**Intentionally simplified:**

- Layout uses greedy best-successor chaining instead of full transitive-edge reduction + Hamiltonian-path search; **consequence:** repeat-containing inputs may be split or mis-resolved rather than optimally laid out (exact layout is NP-complete) [1][2 (p.21–25)].
- Consensus is realized as concatenation along the chain (the merged superstring), not a per-column majority vote over a multiple read pile-up; **consequence:** for exact-overlap chains the result is identical to the majority-vote consensus, but mismatch resolution within an accepted overlap is not performed [2 (p.28)].

**Not implemented:**

- Spurious-branch (sequencing-error) pruning of the overlap graph [2 (p.26)]; **users should rely on:** upstream read error-correction (`SequenceAssembler.ErrorCorrectReads`) before assembly, or an external assembler.
- Suffix-tree-accelerated overlap detection [2 (p.6–10)]; **users should rely on:** the O(N²) all-pairs path here (required for the approximate-identity feature).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Greedy layout vs Hamiltonian-path layout | Deviation | Repeat-containing genomes may not be optimally reconstructed | accepted | Exact layout NP-complete [1]; see ASM-02 |
| 2 | Empty read set → empty result | Assumption | Defines unspecified behavior | accepted | No source specifies; trivial identity case |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `reads` null or empty | empty `AssemblyResult` (TotalReads 0, TotalLength 0) | Contract; no source-specified behavior |
| No pair overlaps ≥ threshold | each read is a singleton contig | Edgeless overlap graph [1] (INV-05) |
| Overlap exactly = `minOverlap` | accepted; one shorter rejected | Threshold `l` [2 (p.5)] (INV-02) |
| Lowercase reads | same overlaps as uppercase | case-insensitive identity |
| Unambiguous chain | single merged contig = superstring of all reads | INV-04 [1][3 (p.26)] |
| Repeat ≥ read length | contigs may split / not collapse below longest read | ASM-02 [3 (p.58–62)] |

### 6.2 Limitations

Does not perform error-aware branch pruning, does not model paired-end/scaffolding within
`AssembleOLC`, and uses a greedy (suboptimal) layout — greedy maximal-overlap merging can
produce a longer-than-minimal superstring [3 (p.57)]. Repeats longer than the read length
are not resolvable from overlaps alone [3 (p.58–62)]. Overlap detection is O(N²) and not
suited to datasets with hundreds of millions of reads [2 (p.29)].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reads = new[] { "AAAAACCCCC", "CCCCCGGGGG", "GGGGGTTTTT" };
var result = SequenceAssembler.AssembleOLC(
    reads,
    new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 1.0, MinContigLength: 10));
// result.Contigs == [ "AAAAACCCCCGGGGGTTTTT" ]  (one contig, length 20)
```

**Overlap-graph walk-through:** the 6 distinct 6-mers of `GTACGTACGAT`
(`GTACGT, TACGTA, ACGTAC, CGTACG, GTACGA, TACGAT`) with `minOverlap = 4` yield 12
directed edges with weights 4 and 5 — e.g. `GTACGT → TACGTA` (5), `GTACGT → ACGTAC` (4),
`CGTACG → TACGAT` (4), `GTACGA → TACGAT` (5) — reproducing the edge weights on Langmead
SCS p.24–25 [3].

### 7.2 Performance Baseline

`AssembleOLC` is dominated by the O(N²) all-pairs overlap scan (§4.3). Measured on random
100-base reads tiling a 5,000-base genome (Release build, single run, seed 42):

| # reads | `AssembleOLC` time |
|---------|--------------------|
| 100 | 251 ms |
| 200 | 678 ms |
| 400 | 2,737 ms |

The ~4× time increase per read-count doubling confirms the quadratic scaling in read count
and the practical limit noted in §6.2 (not suited to datasets of hundreds of millions of reads).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_AssembleOLC_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAssembler_AssembleOLC_Tests.cs) — covers `INV-01`–`INV-05`, `ASM-02`.
- Evidence: [ASSEMBLY-OLC-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-OLC-001-Evidence.md)

## 8. References

1. Compeau PEC, Pevzner PA, Tesler G. 2011. How to apply de Bruijn graphs to genome assembly. *Nature Biotechnology* 29(11):987–991. https://doi.org/10.1038/nbt.2023
2. Langmead B. Overlap Layout Consensus assembly (lecture notes, Johns Hopkins University). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf
3. Langmead B. Assembly & Shortest Common Superstring (lecture notes, Johns Hopkins University). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/16_assembly_scs_v2.pdf
