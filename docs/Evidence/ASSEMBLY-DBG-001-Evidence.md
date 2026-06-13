# Evidence Artifact: ASSEMBLY-DBG-001

**Test Unit ID:** ASSEMBLY-DBG-001
**Algorithm:** De Bruijn graph genome assembly — graph construction (`BuildDeBruijnGraph`) and Eulerian-walk reconstruction (`AssembleDeBruijn`)
**Date Collected:** 2026-06-13

---

## Online Sources

### Langmead (JHU) — "De Bruijn Graph assembly" lecture notes

**URL:** https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf
**Accessed:** 2026-06-13
**Authority rank:** 3 (teaching notes by an established author; the notes themselves cite the Jones & Pevzner textbook for the Euler theorems)

**Retrieved by:** WebSearch query `Langmead de Bruijn graph assembly lecture notes k-1-mer node edge Eulerian path ZABCDABEFABY example`, then WebFetch of the PDF; binary PDF extracted locally page-by-page with `pypdf` (`PdfReader(...).pages[i].extract_text()`).

**Key Extracted Points:**

1. **k-mer / (k-1)-mer (p.4-5):** "'k-mer' is a substring of length k." "AAB is a k-mer (k = 3). AA is its left k-1-mer, and AB is its right k-1-mer." The "k-1-mer" refers to a substring of length k-1.
2. **Nodes and edges (p.6-7):** "Take each length-3 input string and split it into two overlapping substrings of length 2. … Let 2-mers be nodes in a new graph. Draw a directed edge from each left 2-mer to corresponding right 2-mer." "An edge corresponds to an overlap (of length k-2) between two k-1 mers. More precisely, it corresponds to a k-mer from the input." Each edge corresponds to exactly one length-k input string.
3. **Multigraph (p.8-9):** Adding "one more B to our input string: AAABBBBA" produces a **multiedge**. "De Bruijn graph is a directed multigraph G(V, E)". "Node's indegree = # incoming edges. Node's outdegree = # outgoing edges."
4. **Eulerian walk definitions (p.10):** "Eulerian walk visits each edge exactly once." "Node is balanced if indegree equals outdegree." "Node is semi-balanced if indegree differs from outdegree by 1." "A directed, connected graph is Eulerian if and only if it has at most 2 semi-balanced nodes and all other nodes are balanced." (Cites "Jones and Pevzner section 8.8".)
5. **Perfect-sequencing graph is Eulerian (p.15):** "With perfect sequencing, this procedure always yields an Eulerian graph." The left-end k-1-mer is semi-balanced with one more outgoing than incoming edge; the right-end k-1-mer is semi-balanced with one more incoming than outgoing; all other nodes are balanced (unless the genome is circular).
6. **Reconstruction by spelling the walk (p.18-19):** the Eulerian walk is a node sequence; the genome is recovered as `superstring = path[0] + ''.join(x[-1] for x in path[1:])` — i.e. emit the first (k-1)-mer in full, then append the **last character** of every subsequent node. Worked example: `DeBruijnGraph(["a_long_long_long_time"], 5)` → walk `['a_lo','_lon','long','ong_','ng_l','g_lo','_lon','long','ong_','ng_l','g_lo','_lon','long','ong_','ng_t','g_ti','_tim','time']`, spelling `a_long_long_long_time`.
7. **`to_every…` example (p.19, p.22):** `DeBruijnGraph(["to_every_thing_turn_turn_turn_there_is_a_season"], 4)` → superstring `to_every_thing_turn_turn_turn_there_is_a_season` (k = 4 reconstructs exactly); at k = 3 it does NOT ("to_every_turn_turn_thing_turn_there_is_a_season") "Due to repeats that are unresolvable at k = 3".
8. **Multiple Eulerian walks (p.21):** graph for `ZABCDABEFABY`, k = 3, has two alternative Eulerian walks; "graph can have multiple Eulerian walks, only one of which corresponds to original superstring"; AB is a repeat joining two edge-disjoint cycles.
9. **Construction implementation (p.16):** `chop(st,k)` yields `st[i:i+k]` for `i in range(0, len(st)-(k-1))`; for each k-mer `km1L, km1R = kmer[:-1], kmer[1:]`; nodes are created on demand and the edge `(nodeL → nodeR)` is appended to a multimap. (Confirms iteration bound and prefix/suffix split.)
10. **Eulerian-walk complexity (p.17):** "For Eulerian graph, Eulerian walk can be found in O(|E|) time. |E| is # edges" (Hierholzer-style recursive procedure).

### Jones & Pevzner (2004) — *An Introduction to Bioinformatics Algorithms*, MIT Press

**URL:** https://eclass.uoa.gr/modules/document/file.php/NURS565/BioinformaticsAlgsBook.pdf (full-text PDF; ISBN 0-262-10106-8)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed textbook; the primary cited by the Langmead notes for the Euler theorems)

**Retrieved by:** WebSearch query `Jones Pevzner "Introduction to Bioinformatics Algorithms" Eulerian path theorem balanced semi-balanced node de Bruijn section 8.8`, then WebFetch of the PDF; text extracted locally with `pypdf` and grepped for "Theorem", "balanced", "semibalanced".

**Key Extracted Points:**

1. **Euler's theorem — Eulerian cycle (Theorem 8.1):** "A connected graph is Eulerian if and only if each of its vertices is balanced." (A vertex is *balanced* when indegree = outdegree.) The construction proof is implementable "in time linear in the number of edges in the graph."
2. **Eulerian path (Theorem 8.2):** "A connected graph has an Eulerian path if and only if it contains at most two semibalanced vertices and all other vertices are balanced." A vertex v is *semibalanced* if `|indegree(v) − outdegree(v)| = 1`.
3. **Eulerian path → Eulerian cycle reduction (§8.9):** "If a graph has an Eulerian path starting at vertex s and ending at vertex t, then all its vertices are balanced, with the possible exception of s and t, which may be semibalanced." The Eulerian-path problem reduces to the Eulerian-cycle problem; fragment assembly (SBH) reduces to finding an Eulerian path.

### Compeau, Pevzner & Tesler (2011) — "How to apply de Bruijn graphs to genome assembly", Nature Biotechnology

**URL:** https://doi.org/10.1038/nbt.2023 (abstract page https://www.nature.com/articles/nbt.2023; scanned/open-access mirror https://www.cs.umb.edu/~rvetro/vetroBioComp/CancerGenome/Assembly/Compeau2011%20How%20to%20apply%20de%20Bruijn%20graphs%20to%20genome%20assembly.pdf)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)

**Retrieved by:** WebSearch query `Compeau Pevzner Tesler "How to apply de Bruijn graphs to genome assembly" Nature Biotechnology de Bruijn graph k-mer node edge Eulerian path`; the Nature page and two PDF mirrors were fetched. **Note:** both PDF mirrors are image/scanned (no extractable text layer; `pypdf` returned zero readable characters), so only the publication metadata and abstract obtained via WebSearch/Nature are cited here, not page-level quotes. The de Bruijn / Eulerian-path theory used for behavior is taken from the two text-extractable sources above (Langmead notes and Jones & Pevzner textbook).

**Key Extracted Points:**

1. **Publication identity:** Compeau PEC, Pevzner PA, Tesler G. 2011. "How to apply de Bruijn graphs to genome assembly." *Nature Biotechnology* 29(11):987-991. DOI 10.1038/nbt.2023. (Confirmed via the Nature article page and Semantic Scholar / ResearchGate listings returned by the search.)
2. **Thesis (abstract):** representing reads as a de Bruijn graph (k-mers as edges between (k-1)-mer nodes) "turns the formidable challenge of assembling a contiguous genome from billions of short sequencing reads into a tractable computational problem" — i.e. assembly becomes an Eulerian-path problem (solvable in polynomial time), versus the NP-complete Hamiltonian-path overlap formulation. The paper cites Pevzner, Tang & Waterman (2001) "An Eulerian path approach to DNA fragment assembly", *PNAS*.

---

## Documented Corner Cases and Failure Modes

### From Langmead DBG notes

1. **Multiple Eulerian walks at a repeat (p.21):** a (k-1)-mer that occurs more than once (a repeat ≥ k-1) becomes a shared node joining edge-disjoint cycles; several Eulerian walks exist, only one of which is the true genome. The reconstruction is therefore not guaranteed unique when a repeat is unresolvable at the chosen k.
2. **Repeat unresolvable at small k (p.22):** the same input reconstructs correctly at k = 4 but incorrectly at k = 3. Increasing k can resolve repeats shorter than k-1.
3. **Gaps in coverage → disconnected graph (p.24-25):** omitting a k-mer (`ong_t`) splits the graph; "Connected components are individually Eulerian, overall graph is not." Assembly then yields multiple contigs.
4. **Coverage differences → non-Eulerian (p.26):** an extra copy of a k-mer creates 4 semi-balanced nodes, so the graph is no longer Eulerian (Theorem 8.2 violated).
5. **Sequencing errors → non-Eulerian / disconnected (p.27):** an error turning `long_` into `lxng_` makes the largest component non-Eulerian.
6. **Eulerian framing is appealing but not practical at scale (p.28):** uneven coverage, errors and repeats make real graphs non-Eulerian; the De Bruijn Superwalk Problem is NP-hard (cites Medvedev et al. 2007; Kingsford et al. 2010).

### From Jones & Pevzner

1. **Existence preconditions (Theorems 8.1/8.2):** an Eulerian walk exists only if the graph is connected and has at most two semibalanced vertices, all others balanced. Disconnected or unbalanced graphs have no single Eulerian walk.

---

## Test Datasets

### Dataset: AAABBBA de Bruijn graph, k = 3 (Langmead DBG p.5-11)

**Source:** Langmead, "De Bruijn Graph assembly", p.5-11; node/edge set re-derived from the prefix/suffix split and confirmed computationally (`pypdf`-extracted procedure on p.16 reimplemented in Python).

The 3-mers of `AAABBBA` are `AAA, AAB, ABB, BBB, BBA`. Splitting each into left/right 2-mers gives nodes `{AA, AB, BB, BA}` and directed multigraph edges:

| Edge (k-mer) | Left node (prefix) | Right node (suffix) |
|--------------|--------------------|---------------------|
| AAA | AA | AA |
| AAB | AA | AB |
| ABB | AB | BB |
| BBB | BB | BB |
| BBA | BB | BA |

Degrees: AA in 1 / out 2 (semi-balanced, source), BA in 1 / out 0 (semi-balanced, sink), AB in 1 / out 1 (balanced), BB in 2 / out 2 (balanced). Exactly two semi-balanced nodes ⇒ Eulerian (Theorem 8.2). Eulerian walk `AA → AA → AB → BB → BB → BA` (Langmead p.11). Spelling the walk recovers `AAABBBA` (length 7).

### Dataset: a_long_long_long_time, k = 5 (Langmead DBG p.12-18)

**Source:** Langmead DBG p.18 — printed Eulerian-walk output of `DeBruijnGraph(["a_long_long_long_time"], 5)`.

| Input | k | Reconstructed superstring |
|-------|---|---------------------------|
| `a_long_long_long_time` | 5 | `a_long_long_long_time` (verbatim) |

The walk is 18 nodes long; spelling it (`path[0]` then last char of each subsequent node) reproduces the 21-character input. Verified computationally with a Hierholzer Eulerian-walk reimplementation (matches the published output exactly).

### Dataset: to_every_thing_turn_turn_turn_there_is_a_season, k = 4 vs k = 3 (Langmead DBG p.19, p.22)

**Source:** Langmead DBG p.19, p.22 — printed outputs.

| k | Reconstruction |
|---|----------------|
| 4 | `to_every_thing_turn_turn_turn_there_is_a_season` (correct) |
| 3 | `to_every_turn_turn_thing_turn_there_is_a_season` (mis-ordered by the `turn` repeat) |

Demonstrates that increasing k resolves a repeat unresolvable at smaller k. Verified computationally for k = 4.

### Dataset: ATGGCGTGCA, k = 4 (DNA exact-reconstruction smoke case)

**Source:** derived directly from the construction + spelling rule (Langmead p.6-7, p.18); a repeat-free DNA string whose unique Eulerian path reconstructs it.

| Input | k | Reconstruction |
|-------|---|----------------|
| `ATGGCGTGCA` | 4 | `ATGGCGTGCA` |

All 7 of its 4-mers are distinct and no 3-mer repeats, so the Eulerian path is unique. Verified computationally.

---

## Assumptions

1. **ASSUMPTION: single-read tie-break / walk selection among multiple Eulerian walks.** When the graph admits more than one Eulerian walk (repeat ≥ k-1, Langmead p.21), the sources do not prescribe which walk a deterministic implementation must emit — they only state the true genome is one of them. Tests therefore assert exact reconstruction ONLY on inputs whose Eulerian walk is unique (no (k-1)-mer repeats: `AAABBBA`, `a_long_long_long_time` at k=5, `to_every…` at k=4, `ATGGCGTGCA` at k=4). Non-unique inputs are asserted only on source-guaranteed invariants (e.g. each input k-mer is spelled by some contig; total reconstructed length over all contigs equals the sum over edges), not on a specific walk.
2. **ASSUMPTION: empty / null read set → empty result.** No source specifies behavior for an empty read collection; the repository returns an empty `AssemblyResult`. Treated as the trivial identity edge case (mirrors `AssembleOLC`).
3. **ASSUMPTION: reads shorter than k contribute no k-mers.** Per the `chop` bound `range(0, len(st)-(k-1))` (Langmead p.16), a read of length < k yields zero k-mers and is silently skipped; this follows directly from the construction rule and is not a separate modelling choice.

---

## Recommendations for Test Coverage

1. **MUST Test:** `BuildDeBruijnGraph(["AAABBBA"], 3)` yields nodes `{AA,AB,BB,BA}` and exactly the 5 directed edges `AA→AA, AA→AB, AB→BB, BB→BB, BB→BA` (multigraph; BB→BB present once) — Evidence: Langmead DBG p.5-11.
2. **MUST Test:** `AssembleDeBruijn` on the 3-mers of `AAABBBA` (k=3) reconstructs a single contig `AAABBBA` — Evidence: Langmead DBG p.11 Eulerian walk.
3. **MUST Test:** `AssembleDeBruijn` on the 5-mers of `a_long_long_long_time` (k=5) reconstructs `a_long_long_long_time` — Evidence: Langmead DBG p.18.
4. **MUST Test:** `AssembleDeBruijn` on the 4-mers of `to_every_thing_turn_turn_turn_there_is_a_season` (k=4) reconstructs the input exactly — Evidence: Langmead DBG p.19/p.22.
5. **MUST Test:** `BuildDeBruijnGraph` produces a multiedge for `AAABBBBA` (k=3): node BB has two outgoing edges to BB and BA, and an extra BB→BB self-loop appears once more than in `AAABBBA` — Evidence: Langmead DBG p.8.
6. **SHOULD Test:** a read shorter than k yields no k-mers / empty graph — Rationale: `chop` bound (Langmead p.16); ASSUMPTION-3.
7. **SHOULD Test:** empty and null read sets return an empty `AssemblyResult` — Rationale: ASSUMPTION-2 (mirrors OLC).
8. **SHOULD Test (invariant):** each distinct input k-mer corresponds to exactly one graph edge; |edges| over all nodes equals the total number of k-mers chopped from the reads — Evidence: Langmead p.7 ("each edge corresponds to a length-k input string").
9. **COULD Test (property, O(n·k)):** for a unique-walk DNA input, every input read is a substring of the single reconstructed contig and the contig length equals (#distinct (k-1)-nodes on the path) — Rationale: spelling rule (Langmead p.18).
10. **COULD Test:** the graph for an input with a repeated (k-1)-mer (`to_every…`, k=3) contains a branch node (out-degree ≥ 2), the structural cause of multiple Eulerian walks and unresolvable repeats — Evidence: Langmead DBG p.21-22. (Note: the *specific* mis-reconstruction printed by Langmead's `eulerianWalkOrCycle` depends on an arbitrary walk choice; under this library's deterministic Hierholzer walk the k=3 input may still be recovered, so the test asserts the documented branching structure, not a specific wrong string — see ASSUMPTION-1.)

---

## References

1. Langmead B. De Bruijn Graph assembly (lecture notes, Johns Hopkins University). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf
2. Jones NC, Pevzner PA. 2004. *An Introduction to Bioinformatics Algorithms*. MIT Press. ISBN 0-262-10106-8. (Theorems 8.1, 8.2; §8.8-8.9.) Full text: https://eclass.uoa.gr/modules/document/file.php/NURS565/BioinformaticsAlgsBook.pdf
3. Compeau PEC, Pevzner PA, Tesler G. 2011. How to apply de Bruijn graphs to genome assembly. *Nature Biotechnology* 29(11):987-991. https://doi.org/10.1038/nbt.2023
4. Pevzner PA, Tang H, Waterman MS. 2001. An Eulerian path approach to DNA fragment assembly. *PNAS* 98(17):9748-9753. https://doi.org/10.1073/pnas.171285098 (cited by [3]; identified via the search of [3]).

---

## Change History

- **2026-06-13**: Initial documentation.
