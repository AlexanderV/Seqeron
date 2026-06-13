# De Bruijn Graph (DBG) Assembly

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-DBG-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Reference |
| Last Reviewed | 2026-06-13 |

## 1. Overview

De Bruijn graph (DBG) assembly is a *de novo* genome-assembly paradigm that reconstructs
contiguous sequences (contigs) from short reads by decomposing them into fixed-length
k-mers. Every distinct (k-1)-mer becomes a node and every input k-mer becomes a directed
edge from its prefix (k-1)-mer to its suffix (k-1)-mer, producing a directed multigraph
[1][3]. Assembly is then framed as an **Eulerian walk** — a walk that traverses every
edge exactly once — whose node sequence spells out the genome [1][2]. Unlike the
overlap/Hamiltonian-path formulation (NP-complete), an Eulerian walk can be found in time
linear in the number of edges, which is what makes the DBG formulation attractive [2][3].
The reconstruction is exact when the Eulerian walk is unique (no repeated (k-1)-mer); when
a (k-1)-mer repeats, multiple Eulerian walks exist and only one is the true genome, so the
emitted contig is not guaranteed to be the original sequence [1].

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

A node `v` is **balanced** if `indegree(v) = outdegree(v)` and **semi-balanced** if
`|indegree(v) − outdegree(v)| = 1` [1][2]. **Euler's theorem** (Jones & Pevzner, Theorem
8.1): a connected graph is Eulerian (has an Eulerian cycle) iff every node is balanced [2].
**Theorem 8.2**: a connected graph has an Eulerian path iff it has at most two
semi-balanced nodes and all others are balanced [2]. With perfect, even-coverage
sequencing of a linear genome the DBG construction satisfies Theorem 8.2: the genome's
first (k-1)-mer has one extra outgoing edge, the last has one extra incoming edge, and all
other nodes are balanced [1].

The genome is recovered from an Eulerian walk `p = (p₀, p₁, …, p_m)` by spelling:

```
genome = p₀ · (last char of p₁) · (last char of p₂) · … · (last char of p_m)
```

i.e. emit the first (k-1)-mer in full, then append the final character of each subsequent
node [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every edge runs from the (k-1)-prefix to the (k-1)-suffix of some input k-mer | construction rule [1] |
| INV-02 | Total #edges = total #k-mers chopped from all reads (multigraph) | one edge per k-mer occurrence [1] |
| INV-03 | A read of length < k contributes zero edges | the chop bound `i ∈ [0, |r|−k]` is empty [1] |
| INV-04 | If a component's Eulerian walk is unique, the spelled contig is the source genome | spelling rule + Theorem 8.2 [1][2] |
| INV-05 | Every input k-mer's string is a substring of some emitted contig | each k-mer is an edge traversed by the walk [1] |
| INV-06 | `TotalLength` = Σ contig lengths and `LongestContig` = max contig length | statistics definition |

### 2.5 Comparison with Related Methods

| Aspect | De Bruijn (Eulerian) | Overlap-Layout-Consensus (Hamiltonian) |
|--------|----------------------|-----------------------------------------|
| Graph node | (k-1)-mer | read |
| Graph edge | k-mer | pairwise suffix-prefix overlap |
| Layout problem | Eulerian path — linear time in \|E\| [2][3] | Hamiltonian path — NP-complete [3] |
| Scales with | k-mer diversity | number of reads (all-pairs overlap) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `reads` | `IReadOnlyList<string>` | required | input reads | null/empty → empty result; reads of length < k skipped |
| `parameters.KmerSize` (k) | `int` | 31 | k-mer length | must be ≥ 2 (`BuildDeBruijnGraph` throws `ArgumentOutOfRangeException` below 2) |
| `parameters.MinContigLength` | `int` | 100 | minimum emitted contig length | contigs shorter than this are discarded |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Contigs` | `IReadOnlyList<string>` | one spelled superstring per weakly-connected component, after length filtering |
| `TotalReads` | `int` | number of input reads |
| `N50` / `LongestContig` / `TotalLength` | numeric | assembly statistics over the emitted contigs |

`BuildDeBruijnGraph(reads, k)` returns the out-adjacency multimap
`Dictionary<string, List<string>>` from each (k-1)-mer to its successor (k-1)-mers.

### 3.3 Preconditions and Validation

Null or empty `reads` returns an empty `AssemblyResult`. `BuildDeBruijnGraph` returns an
empty graph for null/empty reads and throws `ArgumentOutOfRangeException` when `k < 2`.
Matching is case-sensitive (k-mers are taken verbatim); no reverse-complement strand is
considered. Reads shorter than k are skipped (INV-03).

## 4. Algorithm

### 4.1 High-Level Steps

1. **Chop:** for each read, emit every k-mer `r[i..i+k]`.
2. **Build graph:** add a directed edge from each k-mer's (k-1)-prefix to its (k-1)-suffix
   (multigraph; repeated k-mers → repeated edges).
3. **Find components:** split the graph into weakly-connected components.
4. **Eulerian walk per component:** start at the semi-balanced source (out−in = 1) if one
   exists, else any node; traverse every edge exactly once (Hierholzer).
5. **Spell:** emit the first node, then the last character of each subsequent node.
6. **Filter:** keep contigs of length ≥ `MinContigLength`; compute statistics.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The graph is stored as an out-adjacency multimap (`Dictionary<string, List<string>>`),
keyed by (k-1)-mer; edge multiplicity is preserved by list duplication. Hierholzer's
algorithm uses a per-node cursor into this list so each edge is consumed exactly once.
Component order and the start node within a balanced component are chosen as the
lexicographically smallest node, making the output deterministic and order-independent.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `BuildDeBruijnGraph` | O(N·k) | O(N·k) | N = total read length; k-mer substring extraction |
| Eulerian walk (Hierholzer) | O(\|E\|) | O(\|E\|) | \|E\| = number of k-mers [2] |
| `AssembleDeBruijn` (overall) | O(N·k) | O(N·k) | dominated by graph construction |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.AssembleDeBruijn(reads, parameters)`: builds the DBG and reconstructs
  one contig per connected component via an Eulerian walk.
- `SequenceAssembler.BuildDeBruijnGraph(reads, k)`: public; returns the (k-1)-mer
  out-adjacency multigraph.
- `SequenceAssembler.ReconstructContigs` / `SpellEulerianWalk` (private): component
  detection + Hierholzer Eulerian walk + spelling.

### 5.2 Current Behavior

This is not a substring-search / pattern-matching unit — it constructs and traverses a
graph over k-mers rather than locating occurrences of a query in a text — so the
repository suffix tree (`SuffixTree`) is **not applicable** and is not used.

Walk selection among multiple Eulerian walks is deterministic (lexicographically smallest
start node, edges consumed in insertion order). When the Eulerian walk is unique the
output equals the source genome; when it is not (a repeated (k-1)-mer), the emitted contig
is one valid Eulerian spelling, not necessarily the original sequence. A graph that is
disconnected (coverage gaps) yields one contig per component.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Nodes = (k-1)-mers, edges = k-mers, prefix→suffix direction, multigraph [1].
- Eulerian-walk reconstruction by spelling `p₀ + last-char(p_i)` [1].
- Hierholzer linear-time Eulerian traversal; semi-balanced start per Theorem 8.2 [1][2].

**Intentionally simplified:**

- Walk selection among multiple Eulerian walks: deterministic single choice;
  **consequence:** for inputs with unresolvable repeats the contig may differ from the true
  genome (any of several Eulerian walks is "correct" per the model) [1].
- No error correction, bubble/tip removal, or coverage-based edge weighting;
  **consequence:** real (errored, uneven-coverage) data produces non-Eulerian or fragmented
  graphs and therefore fragmented or mis-assembled output [1].

**Not implemented:**

- Read error correction, repeat resolution via paired-end/long reads, reverse-complement
  strand handling; **users should rely on:** production assemblers (SPAdes, Velvet) for
  real sequencing data.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Walk selection among multiple Eulerian walks | Assumption | repeat-containing inputs may not reconstruct the true genome | accepted | sources do not prescribe a walk; deterministic choice (ASSUMPTION-1) |
| 2 | Empty/null reads → empty result | Assumption | no exception on empty input | accepted | trivial identity (ASSUMPTION-2) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / null reads | empty `AssemblyResult`; empty graph | trivial identity (ASSUMPTION-2) |
| Read length < k | contributes no edges | chop bound is empty (INV-03) [1] |
| `k < 2` | `BuildDeBruijnGraph` throws `ArgumentOutOfRangeException` | (k-1)-mer nodes must be non-empty [1] |
| Repeated k-mer | multiedge in the graph | multigraph definition [1] |
| Disconnected graph (coverage gap) | one contig per component | components are individually Eulerian [1] |

### 6.2 Limitations

Reconstruction is exact only for perfect-coverage, error-free reads whose Eulerian walk is
unique. Repeats ≥ k-1 create multiple valid walks (only one is the true genome); coverage
gaps fragment the assembly; sequencing errors and uneven coverage make the graph
non-Eulerian [1]. This is a reference implementation for teaching the Eulerian-path
formulation, not a production assembler.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reads = new[] { "AAABBBA" };
var result = SequenceAssembler.AssembleDeBruijn(
    reads,
    new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1));
// result.Contigs[0] == "AAABBBA"
```

**Numerical / biological walk-through:**

`AAABBBA`, k = 3. 3-mers: `AAA, AAB, ABB, BBB, BBA`. Nodes (2-mers): `AA, AB, BB, BA`.
Edges: `AA→AA, AA→AB, AB→BB, BB→BB, BB→BA`. Degrees: AA (in 1, out 2 — source), BA (in 1,
out 0 — sink), AB and BB balanced ⇒ exactly two semi-balanced nodes ⇒ Eulerian. Walk
`AA → AA → AB → BB → BB → BA`; spelling `AA + A + B + B + B + A = AAABBBA` [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_AssembleDeBruijn_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_AssembleDeBruijn_Tests.cs) — covers `INV-01`…`INV-06`
- Evidence: [ASSEMBLY-DBG-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-DBG-001-Evidence.md)
- Related algorithms: [Overlap_Layout_Consensus](../Assembly/Overlap_Layout_Consensus.md)

## 8. References

1. Langmead B. De Bruijn Graph assembly (lecture notes, Johns Hopkins University). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf
2. Jones NC, Pevzner PA. 2004. *An Introduction to Bioinformatics Algorithms*. MIT Press. ISBN 0-262-10106-8 (Theorems 8.1, 8.2; §8.8-8.9). https://eclass.uoa.gr/modules/document/file.php/NURS565/BioinformaticsAlgsBook.pdf
3. Compeau PEC, Pevzner PA, Tesler G. 2011. How to apply de Bruijn graphs to genome assembly. *Nature Biotechnology* 29(11):987-991. https://doi.org/10.1038/nbt.2023
