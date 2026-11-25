# De Bruijn Graph (DBG) Assembly

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-DBG-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Reference |
| Last Reviewed | 2026-11-25 (audit round 1, WP4: A2-4, A2-5) |

## 1. Overview

De Bruijn graph (DBG) assembly is a *de novo* genome-assembly paradigm that builds
contiguous sequences (contigs) from short reads by decomposing them into fixed-length
k-mers. Every distinct (k-1)-mer becomes a node and every k-mer becomes a directed edge
from its prefix (k-1)-mer to its suffix (k-1)-mer [1][3][4]. Because real reads leave gaps
in k-mer coverage, overlap each other and contain repeats, a single Eulerian walk through
the whole graph is neither guaranteed to exist nor unique [1]; assemblers therefore output
**contigs — the strings spelled by the maximal non-branching paths (unitigs) of the
graph**, which "must be present in any assembly with a given k-mer composition" [4]. Velvet
merges exactly these chains ("whenever a node A has only one outgoing arc that points to
another node B that has only one ingoing arc, the two nodes are merged") [5].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Shotgun sequencing yields many short reads sampled from an unknown genome. Assembly
recovers the genome (or long fragments of it) from these reads. The DBG approach trades
the per-read overlap graph of OLC for a per-k-mer graph: nodes are (k-1)-mers and edges
are k-mers, so the number of nodes is bounded by the k-mer space rather than by the number
of reads [3].

### 2.2 Core Model

For a read `r` and k-mer length `k`, each length-k substring `r[i..i+k]` (`0 ≤ i ≤ |r|−k`)
is an **edge** whose endpoints are its left (k-1)-mer `r[i..i+k−1]` (prefix) and right
(k-1)-mer `r[i+1..i+k]` (suffix); the edge is directed prefix → suffix [1]. The graph
`G(V, E)` is a directed **multigraph**: a k-mer that occurs more than once contributes a
repeated edge [1].

A node is **1-in-1-out** when `indegree(v) = outdegree(v) = 1`. A path is
**non-branching** when all its internal nodes are 1-in-1-out; a **maximal non-branching
path** cannot be extended, i.e. its first and last nodes are not 1-in-1-out. A connected
component consisting only of 1-in-1-out nodes is an **isolated cycle**; the book's
`MaximalNonBranchingPaths` emits each isolated cycle as an additional path [4]:

```
MaximalNonBranchingPaths(Graph)
  Paths ← empty list
  for each node v in Graph
    if v is not a 1-in-1-out node and out(v) > 0
      for each outgoing edge (v, w) from v
        NonBranchingPath ← the path consisting of the single edge (v, w)
        while w is a 1-in-1-out node
          extend NonBranchingPath by the outgoing edge (w, u) from w;  w ← u
        add NonBranchingPath to Paths
  for each isolated cycle Cycle in Graph
    add Cycle to Paths
  return Paths
```

A path `p = (p₀, p₁, …, p_m)` is spelled as

```
contig = p₀ · (last char of p₁) · (last char of p₂) · … · (last char of p_m)
```

i.e. emit the first (k-1)-mer in full, then append the final character of each subsequent
node [1][4].

**Multiplicity.** `BuildDeBruijnGraph` builds the textbook `DeBruijn(Patterns)`
multigraph of a k-mer *composition*, in which a repeated k-mer is a repeated edge [1][4].
Read data cannot distinguish a genomic repeat from coverage, so the assembler uses the set
of **distinct** k-mers: "the fundamental data structure in the de Bruijn graph is based on
k-mers, not reads, thus high redundancy is naturally handled by the graph without affecting
the number of nodes" (Velvet [5]); unitigs are defined over the k-mer set (BCALM 2 /
Chikhi et al. 2016 [6]).

**Strands (default: both).** Each k-mer is entered with its reverse complement ("each k-mer is
recorded simultaneously to its reverse complement ... this ensures that overlaps between
reads from opposite strands are taken into account" [5]); the resulting graph is
reverse-complement symmetric, so each unitig's reverse complement is also a unitig and each
such pair is reported once, as the ordinally smaller of the two strings. Velvet always builds this
double-strand graph, so `AssembleDeBruijn(reads, parameters)` does too; `includeReverseComplement: false`
selects the single-strand textbook model (Rosalind BA3K/BA3M), which the oracle values of the BA3M tests use.

**Non-ACGT symbols (Velvet read rule).** Velvet's read loader `velvetifySequence` (`src/readSet.c`) keeps
A/C/G/T (either case) and turns every other symbol into `N`; `writeNucleotide` (`src/tightString.c`) stores `N`
as Adenine. `AssembleDeBruijn` therefore upper-cases each read and replaces every symbol outside A/C/G/T —
`N`, IUPAC codes, `U`, `-`, any non-ASCII character — by `A` before the k-mers are taken, so a k-mer that
contains such a symbol is the k-mer with an A in that position, the contigs are over {A, C, G, T}, and
`AssembledReads` is counted on the folded reads. The textbook primitives `BuildDeBruijnGraph` /
`FindMaximalNonBranchingPaths` take arbitrary symbols verbatim (Langmead / Rosalind text examples).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every edge runs from the (k-1)-prefix to the (k-1)-suffix of some input k-mer | construction rule [1] |
| INV-02 | `BuildDeBruijnGraph`: total #edges = total #k-mers chopped from all reads (multigraph) | one edge per k-mer occurrence [1] |
| INV-03 | A read of length < k contributes zero edges | the chop bound `i ∈ [0, |r|−k]` is empty [1] |
| INV-04 | A sequence with no repeated (k-1)-mer is a single simple path ⇒ one contig equal to the sequence (for any read coverage of it) | maximal non-branching path [4] + distinct k-mers [5] |
| INV-05 | Every distinct input k-mer lies in exactly one contig and no contig contains a k-mer absent from the input (unitigs partition the k-mer set) | every edge lies on exactly one maximal non-branching path [4][6] |
| INV-06 | `TotalLength` = Σ contig lengths and `LongestContig` = max contig length | statistics definition |
| INV-07 | Contigs do not depend on read order or on how many times a k-mer is covered | distinct k-mer set, ordinal ordering |

### 2.5 Comparison with Related Methods

| Aspect | De Bruijn (unitigs) | Overlap-Layout-Consensus (Hamiltonian) |
|--------|---------------------|-----------------------------------------|
| Graph node | (k-1)-mer | read |
| Graph edge | k-mer | pairwise suffix-prefix overlap |
| Contig definition | maximal non-branching path — linear time in \|E\| [4] | greedy layout path (Hamiltonian path NP-complete) [3] |
| Scales with | k-mer diversity | number of reads (all-pairs overlap) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `reads` | `IReadOnlyList<string>` | required | input reads | null/empty → empty result; null reads and reads of length < k skipped |
| `parameters` | `AssemblyParameters?` | `DefaultParameters` | null → k = 31, MinContigLength 100 | `new AssemblyParameters()` zero-initialises the struct — pass null or `DefaultParameters` |
| `parameters.KmerSize` (k) | `int` | 31 | k-mer length | ≥ 2, else `ArgumentOutOfRangeException` (checked before the reads) |
| `parameters.MinContigLength` | `int` | 100 | minimum emitted contig length | contigs shorter than this are discarded |
| `includeReverseComplement` | `bool` | **true** (both overloads of the 2-argument form; MCP default true) | add the reverse complement of every k-mer (Velvet); false = forward strand only (BA3M) | odd k recommended (no k-mer is its own reverse complement) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Contigs` | `IReadOnlyList<string>` | upper-case spelled maximal non-branching paths and isolated cycles, after length filtering |
| `TotalReads` | `int` | number of input reads |
| `AssembledReads` | `int` | reads (folded) with at least one k-mer — on either strand when double-stranded — in a reported contig (Velvet `usedReads`) |
| `N50` / `LongestContig` / `TotalLength` | numeric | `CalculateStats` over the emitted contigs (ASSEMBLY-STATS-001) |

`BuildDeBruijnGraph(reads, k)` returns the out-adjacency multimap
`Dictionary<string, List<string>>` of the k-mer composition (case preserved).
`FindMaximalNonBranchingPaths(graph)` returns the paths (node lists) of any out-adjacency
(multi)graph, e.g. Rosalind BA3M/BA3K.

### 3.3 Preconditions and Validation

Null or empty `reads` returns an empty `AssemblyResult`. k < 2 throws
`ArgumentOutOfRangeException` in both `AssembleDeBruijn` and `BuildDeBruijnGraph`.
`AssembleDeBruijn` folds each read (upper-case, non-ACGT -> A, Velvet) and counts k-mers through the canonical
`SequenceExtensions.CountKmersSpan` (upper-case keys); `BuildDeBruijnGraph` takes k-mers
verbatim (any alphabet). `FindMaximalNonBranchingPaths(null)` throws `ArgumentNullException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Distinct k-mers:** fold each read (Velvet: upper-case, non-ACGT -> A); for each read of length ≥ k, collect its
   k-mers (`CountKmersSpan`); by default add each k-mer's reverse complement.
2. **Build graph:** one edge per distinct k-mer from its (k-1)-prefix to its (k-1)-suffix;
   k-mers are inserted in ordinal order.
3. **Unitigs:** `FindMaximalNonBranchingPaths` (BA3M, above), including isolated cycles.
4. **Spell:** first node, then the last character of each subsequent node.
5. **Strand pairs (double-stranded, the default):** skip a path whose first k-mer is already covered by an
   emitted path or its reverse complement; emit the ordinally smaller of the spelled path and
   its reverse complement.
6. **Filter:** keep contigs of length ≥ `MinContigLength`; compute statistics.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Start nodes are visited in ordinal order and their out-edges in stored (ordinal) order; an
isolated cycle starts and ends at its ordinally smallest node. A 1-in-1-out node has a single
out-edge, so "visited" is tracked per 1-in-1-out node. Output order is therefore
deterministic and independent of read order.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| k-mer collection | O(N·k) | O(D·k) | N = total read length, D = distinct k-mers |
| graph build | O(D log D · k) | O(D·k) | ordinal sort for determinism |
| `FindMaximalNonBranchingPaths` | O(\|V\| log \|V\| + \|E\|) | O(\|V\| + \|E\|) | sorted node set |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.AssembleDeBruijn(reads, parameters)` and
  `AssembleDeBruijn(reads, parameters, includeReverseComplement)`.
- `SequenceAssembler.BuildDeBruijnGraph(reads, k)`: k-mer-composition multigraph.
- `SequenceAssembler.FindMaximalNonBranchingPaths(graph)`: BA3M.
- Canonical helpers used: `SequenceExtensions.CountKmersSpan` (k-mer counting),
  `SequenceExtensions.TryGetReverseComplement` (reverse complement).
- MCP: `assemble_de_bruijn` (Alignment server) delegates to the overload above.

### 5.2 Current Behavior

Not a substring-search unit, so the repository suffix tree is not applicable. A repeat
longer than k-2 creates branch nodes, so repeated regions become separate contigs; a
self-loop (homopolymer / tandem repeat of period 1) makes its node branching.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):** nodes = (k-1)-mers, edges = k-mers
[1]; contigs = maximal non-branching paths incl. isolated cycles [4]; distinct k-mers, reverse-complement k-mers (default) and the non-ACGT -> A read rule as in Velvet [5][7].

**Not implemented:** read error correction (see `ErrorCorrectReads`, ASSEMBLY-CORRECT-001),
tip / bubble removal (Velvet Tour Bus), coverage cut-offs, and repeat resolution with
paired-end or long reads; **users should rely on:** production assemblers (SPAdes, Velvet)
for real sequencing data.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty/null reads → empty result | Assumption | no exception on empty input | accepted | ASSUMPTION-2 |
| 2 | Case-insensitive k-mers, upper-case contigs | Convention | soft-masked bases join | accepted | canonical counter; Velvet stores 2-bit bases |
| 3 | Double-stranded by default (A2-4, 2026-11-25) | Default flip | reads of both strands join | resolved | Velvet always builds the double-strand graph [5]; `includeReverseComplement: false` is the BA3M single-strand model |
| 4 | Non-ACGT symbol -> A (A2-5, 2026-11-25) | Sourced rule | k-mers with N/IUPAC/U are the A-substituted k-mer; contigs never contain N | resolved | Velvet `velvetifySequence` + `writeNucleotide` [7] (an N is not a break: Velvet reads it as A) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / null reads | empty `AssemblyResult`; empty graph | ASSUMPTION-2 |
| Read length < k | contributes no k-mers | chop bound is empty (INV-03) [1] |
| `k < 2` | `ArgumentOutOfRangeException` | (k-1)-mer nodes must be non-empty [1] |
| Duplicate reads / overlapping reads | shared k-mers counted once | Velvet [5] |
| Isolated cycle | emitted as one contig starting at its smallest node | BA3M [4] |
| Disconnected graph (coverage gap) | separate contigs | [1][4] |

### 6.2 Limitations

Contigs stop at every branch node, so repeats ≥ k-1, sequencing errors (tips, bubbles) and
heterozygosity fragment the assembly; no error-graph cleaning is done.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var result = SequenceAssembler.AssembleDeBruijn(
    new[] { "AAACCCA" },
    new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1),
    includeReverseComplement: false);          // single-strand textbook model (BA3M)
// result.Contigs == ["AAA", "AACC", "CCA", "CCC"]
```

`AAACCCA` (Langmead's text example `AAABBBA` with B -> C, because the assembler reads every non-ACGT symbol as A), k = 3. Edges: `AA→AA, AA→AC, AC→CC, CC→CC, CC→CA`. AA (in 2, out 2) and CC
(in 2, out 2) are not 1-in-1-out, AC is. Paths: `AA→AA` (AAA), `AA→AC→CC` (AACC),
`CC→CA` (CCA), `CC→CC` (CCC). Overlapping reads `ATGGCGT`, `GGCGTGCA` (k = 4, single strand) give the single
contig `ATGGCGTGCA`; with the default double-stranded graph `GATTACAGGT` + `TAAGGACCTG` (k = 5) join into `GATTACAGGTCCTTA`, and a read with an `N` (`ACGTNGGTCA`, k = 4) is read as `ACGTAGGTCA`. Rosalind BA3K sample (k-mer multiset `ATG ATG TGT TGG CAT GGA GAT AGA`)
through `FindMaximalNonBranchingPaths(BuildDeBruijnGraph(...))` gives
`AGA ATG ATG CAT GAT TGGA TGT`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_AssembleDeBruijn_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAssembler_AssembleDeBruijn_Tests.cs)
- Evidence: [ASSEMBLY-DBG-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-DBG-001-Evidence.md)
- Related algorithms: [Overlap_Layout_Consensus](../Assembly/Overlap_Layout_Consensus.md)

## 8. References

1. Langmead B. De Bruijn Graph assembly (lecture notes, Johns Hopkins University). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf
2. Jones NC, Pevzner PA. 2004. *An Introduction to Bioinformatics Algorithms*. MIT Press. ISBN 0-262-10106-8 (Theorems 8.1, 8.2).
3. Compeau PEC, Pevzner PA, Tesler G. 2011. How to apply de Bruijn graphs to genome assembly. *Nature Biotechnology* 29(11):987-991. https://doi.org/10.1038/nbt.2023
4. Compeau P, Pevzner P. *Bioinformatics Algorithms: An Active Learning Approach*, ch.3; Rosalind BA3K (Generate Contigs from a Collection of Reads) https://rosalind.info/problems/ba3k/ and BA3M (Generate All Maximal Non-Branching Paths in a Graph) https://rosalind.info/problems/ba3m/
5. Zerbino DR, Birney E. 2008. Velvet: algorithms for de novo short read assembly using de Bruijn graphs. *Genome Research* 18:821-829. https://doi.org/10.1101/gr.074492.107
6. Chikhi R, Limasset A, Medvedev P. 2016. Compacting de Bruijn graphs from sequencing data quickly and in low memory. *Bioinformatics* 32(12):i201-i208 (unitigs = maximal non-branching paths).
7. Velvet 1.2.10 source (dzerbino/velvet): `src/readSet.c` `velvetifySequence` (non-ACGT -> N), `src/tightString.c` `writeNucleotide` (N -> Adenine), `src/kmer.c`. https://github.com/dzerbino/velvet
