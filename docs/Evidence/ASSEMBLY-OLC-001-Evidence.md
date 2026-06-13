# Evidence Artifact: ASSEMBLY-OLC-001

**Test Unit ID:** ASSEMBLY-OLC-001
**Algorithm:** Overlap-Layout-Consensus (OLC) genome assembly — overlap detection (`FindAllOverlaps`) and OLC assembly (`AssembleOLC`)
**Date Collected:** 2026-06-13

---

## Online Sources

### Compeau, Pevzner & Tesler (2011) — "How to apply de Bruijn graphs to genome assembly", Nature Biotechnology

**URL:** https://www.cs.umb.edu/~rvetro/vetroBioComp/CancerGenome/Assembly/Compeau2011%20How%20to%20apply%20de%20Bruijn%20graphs%20to%20genome%20assembly.pdf (open-access mirror of Nat Biotechnol 29:987–991, DOI 10.1038/nbt.2023)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)

**Retrieved by:** WebSearch query `Compeau Pevzner Tesler "How to apply de Bruijn graphs to genome assembly" Nature Biotechnology 2011 overlap graph Hamiltonian path reads`, then WebFetch of the PDF mirror above.

**Key Extracted Points:**

1. **Overlap graph:** Each sequencing read becomes a node in a directed graph; an edge from read A to read B indicates "the suffix of A matches the prefix of B, with the overlap length exceeding a minimum threshold."
2. **Layout = Hamiltonian path:** "The genome assembly problem via overlap graphs reduces to finding a Hamiltonian path — a path visiting each node exactly once. If such a path exists, traversing it and merging overlapping reads in order reconstructs the original genomic sequence."
3. **Computational difficulty:** "finding a Hamiltonian path is NP-complete"; no known polynomial-time algorithm exists for the general case. This motivates heuristic layout and the shift to de Bruijn / Eulerian formulations.

### Langmead (JHU) — "Overlap Layout Consensus assembly" lecture notes

**URL:** https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf
**Accessed:** 2026-06-13
**Authority rank:** 3 (teaching notes by an established author; cite the primaries they reflect)

**Retrieved by:** WebSearch query `Ben Langmead overlap layout consensus shortest common superstring greedy assembly example reads GTACGTACGTA`, then WebFetch of the PDF; binary PDF extracted locally with `pypdf` (page-by-page `extract_text`).

**Key Extracted Points:**

1. **Three stages (p.4):** "**Overlap** — Build overlap graph. **Layout** — Bundle stretches of the overlap graph into contigs. **Consensus** — Pick most likely nucleotide sequence for each contig."
2. **Overlap definition (p.5):** for strings X, Y and minimum length `l`, find the length-`l` substring at the start of Y in X "going right-to-left", then "extend to left; … we confirm that a length-6 prefix of Y matches a suffix of X." I.e. the overlap is the longest suffix of X equal to a prefix of Y, length ≥ `l`.
3. **Longest match only (p.10):** "Assume for given string pair we report only the longest suffix/prefix match."
4. **Overlap-graph edge weight (p.20–25):** edges are labeled with the overlap length; the `to_every_thing_turn_turn_turn_there_is_a_season` example (`l = 4, k = 7`) shows edge weights 4, 5 and 6.
5. **Layout — transitive reduction (p.21–24):** "Remove transitively-inferrible edges, starting with edges that skip one node" (then two nodes); "Some edges can be inferred (transitively) from other edges."
6. **Layout — contig emission (p.25):** "Emit contigs corresponding to the non-branching stretches"; branching at an unresolvable repeat splits the layout into multiple contigs (`Contig 1` / `Contig 2`).
7. **Consensus (p.28):** "At each position, ask: what nucleotide (and/or gap) is here? … Take consensus, i.e. majority vote."
8. **Complexity (p.10, p.16):** suffix-tree overlap is `O(N + a)`; all-pairs dynamic-programming overlap is `O(d²n²) = O(N²)` where `d` = #reads of length `n`, `N = dn`, `a` = #overlapping pairs (worst case `a = O(d²)`).

### Langmead (JHU) — "Assembly & Shortest Common Superstring" lecture notes

**URL:** https://www.cs.jhu.edu/~langmea/resources/lecture_notes/16_assembly_scs_v2.pdf
**Accessed:** 2026-06-13
**Authority rank:** 3 (teaching notes by an established author)

**Retrieved by:** same WebSearch as above; WebFetch of the PDF; binary PDF extracted locally with `pypdf`.

**Key Extracted Points:**

1. **First law of assembly (p.16):** "If a suffix of read A is similar to a prefix of read B … then A and B might overlap in the genome."
2. **Overlap graph (p.23–25):** "Each node is a read. Draw edge A → B when suffix of A overlaps prefix of B." Worked example: "Nodes: all 6-mers from `GTACGTACGAT`. Edges: overlaps of length ≥ 4" — the displayed edge weights are 4 and 5 (verified by re-derivation, see Test Datasets).
3. **Greedy-SCS (p.45):** "Greedy-SCS: in each round, merge pair of strings with maximal overlap. Stop when there's 1 string left. `l` = minimum overlap." Worked trace (`l = 1`) on input `AAA AAB ABB BBB BBA` → `AAABBBA` (length 7).
4. **Greedy is not optimal (p.57):** a different merge order on the same input gives `AAABBABBB` (length 9); "Greedy answer isn't necessarily optimal."
5. **Repeats foil assembly (p.58–62):** assembling length-6 substrings of `a_long_long_long_time` is "Foiled by repeat!"; using length-8 substrings (one substring `g_long_l` spans all three `long`s) recovers the whole string. Repeats longer than the read length cannot be resolved.

---

## Documented Corner Cases and Failure Modes

### From Compeau, Pevzner & Tesler (2011)

1. **NP-completeness of exact layout:** finding a Hamiltonian path through the overlap graph is NP-complete, so an exact polynomial OLC layout does not exist in general; practical assemblers use heuristics.

### From Langmead OLC notes

1. **Unresolvable repeats split contigs:** a repeat longer than the read length creates branching in the overlap graph; the layout breaks into multiple contigs rather than one (p.25).
2. **Spurious subgraphs from sequencing error:** error-induced mismatches create dead-end branches that should be pruned (p.26).
3. **Report only the longest overlap per ordered pair (p.10).**

### From Langmead SCS notes

1. **Greedy layout is suboptimal (p.57):** greedy maximal-overlap merging can yield a longer-than-minimal superstring; it is a heuristic, not an exact SCS solver.
2. **Repeats below resolution length collapse (p.58–60):** k-mers/reads shorter than the repeat period cannot distinguish repeat copies.

---

## Test Datasets

### Dataset: GTACGTACGAT 6-mer overlap graph (Langmead SCS notes p.24–25)

**Source:** Langmead, "Assembly & Shortest Common Superstring", p.24–25; edge weights re-derived from the suffix-prefix definition (longest suffix of A equal to a prefix of B, length ≥ 4).

The 6 distinct 6-mers of `GTACGTACGAT` are `GTACGT, TACGTA, ACGTAC, CGTACG, GTACGA, TACGAT`. With minimum overlap 4, the directed overlap graph has the following 12 edges (A → B : overlap length):

| A | B | Overlap |
|---|---|---------|
| GTACGT | TACGTA | 5 |
| GTACGT | ACGTAC | 4 |
| TACGTA | ACGTAC | 5 |
| TACGTA | CGTACG | 4 |
| ACGTAC | GTACGT | 4 |
| ACGTAC | CGTACG | 5 |
| ACGTAC | GTACGA | 4 |
| CGTACG | GTACGT | 5 |
| CGTACG | TACGTA | 4 |
| CGTACG | GTACGA | 5 |
| CGTACG | TACGAT | 4 |
| GTACGA | TACGAT | 5 |

The source slide displays exactly these edge weights (4 and 5). The graph is symmetric in count: 12 directed edges total.

### Dataset: Greedy chain merge (consensus of an unambiguous tiling)

**Source:** Langmead OLC p.5 (longest suffix-prefix overlap) + consensus by concatenation along the chain.

| Reads (5-base overlaps) | Expected single contig |
|-------------------------|------------------------|
| `AAAAACCCCC`, `CCCCCGGGGG`, `GGGGGTTTTT` | `AAAAACCCCCGGGGGTTTTT` (length 20) |

Each adjacent pair shares a length-5 suffix-prefix overlap; merging along the chain (`A + B[overlap:]`) yields a 20-base superstring of all three reads.

### Dataset: Single suffix-prefix overlap (Langmead OLC p.5)

**Source:** Langmead OLC p.5.

| X | Y | l | Longest overlap (suffix X = prefix Y) |
|---|---|---|----------------------------------------|
| `CTCTAGGCC` | `TAGGCCCTC` | 3 | 6 (`TAGGCC`) |

---

## Assumptions

1. **ASSUMPTION: exact-match (identity = 1.0) overlap for canonical numeric cases.** The numeric datasets above assume error-free reads (identity 1.0). The sources also discuss approximate overlaps (mismatch/gap via dynamic programming, OLC p.11–15), but the published numeric examples (GTACGTACGAT graph, greedy traces) are stated for exact matches. The repository `minIdentity` parameter generalizes this; canonical tests use identity 1.0 where exact integer overlaps are published, and a separate case exercises the identity threshold.
2. **ASSUMPTION: empty-input → empty result.** No source explicitly specifies the behavior for an empty read set; the repository returns an empty `AssemblyResult`. This is the natural identity case and is treated as a trivially-correct edge, not a source-derived value.

---

## Recommendations for Test Coverage

1. **MUST Test:** `FindAllOverlaps` on the `GTACGTACGAT` 6-mers (minOverlap 4) returns exactly the 12 directed edges with the overlap lengths in the table — Evidence: Langmead SCS p.24–25 (re-derived).
2. **MUST Test:** `FindAllOverlaps` never emits self-overlaps (`ReadIndex1 != ReadIndex2`) and the reported `OverlapLength` is the longest suffix-prefix match ≥ minOverlap — Evidence: Langmead OLC p.5, p.10.
3. **MUST Test:** `AssembleOLC` on an unambiguous 5-overlap tiling (`AAAAACCCCC`,`CCCCCGGGGG`,`GGGGGTTTTT`) reconstructs the single contig `AAAAACCCCCGGGGGTTTTT` — Evidence: Langmead OLC p.5 + consensus by chain merge.
4. **MUST Test:** `AssembleOLC` on three non-overlapping reads returns 3 singleton contigs (no edges in the overlap graph) — Evidence: Compeau/Pevzner overlap-graph definition (no edge below threshold) + Langmead OLC layout.
5. **SHOULD Test:** identity threshold gates overlap acceptance (7/8 = 0.875 accepted at 0.85, rejected at 0.95) — Rationale: `minIdentity` is a correctness-affecting parameter of overlap detection (OLC p.11–15 generalization).
6. **SHOULD Test:** `MinOverlap` boundary — an overlap exactly equal to `minOverlap` is accepted, one below is rejected — Rationale: threshold semantics (OLC p.5 `l`).
7. **SHOULD Test:** invariant — every contig length ≤ sum of read lengths and ≥ longest single read; total result length is recomputed consistently — Rationale: a merged superstring cannot exceed the concatenation length and cannot be shorter than its longest constituent.
8. **COULD Test:** repeat limitation — assembling reads with an internal repeat does not falsely collapse into one contig shorter than the source genome — Rationale: documented repeat-foils-assembly behavior (SCS p.58–62).

---

## References

1. Compeau PEC, Pevzner PA, Tesler G. 2011. How to apply de Bruijn graphs to genome assembly. *Nature Biotechnology* 29(11):987–991. https://doi.org/10.1038/nbt.2023 (open-access mirror: https://www.cs.umb.edu/~rvetro/vetroBioComp/CancerGenome/Assembly/Compeau2011%20How%20to%20apply%20de%20Bruijn%20graphs%20to%20genome%20assembly.pdf)
2. Langmead B. Overlap Layout Consensus assembly (lecture notes, Johns Hopkins University). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf
3. Langmead B. Assembly & Shortest Common Superstring (lecture notes, Johns Hopkins University). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/16_assembly_scs_v2.pdf

---

## Change History

- **2026-06-13**: Initial documentation.
