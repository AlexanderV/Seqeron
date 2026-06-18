# Validation Report: ASSEMBLY-DBG-001 — De Bruijn graph assembly

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `SequenceAssembler.BuildDeBruijnGraph(reads, k)` (public (k-1)-mer
  out-adjacency multigraph) and `SequenceAssembler.AssembleDeBruijn(reads, parameters)`
  (Eulerian-walk reconstruction into contigs + stats); private helpers
  `ReconstructContigs` / `ChooseEulerStart` / `SpellEulerianWalk` (Hierholzer).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (test-coverage gap found and fixed in session)

## Stage A — Description

### Sources opened & what they confirm
- **Langmead B., "De Bruijn Graph assembly" (JHU lecture notes)** — fetched the PDF this
  session (`https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf`) and
  extracted text page-by-page with `pypdf`. Every Stage-A claim was read directly from the
  slides:
  - p.4-7: "'k-mer' is a substring of length k"; "I'll use 'k-1-mer' to refer to a substring
    of length k-1"; "Take each length-3 input string and split it into two overlapping
    substrings of length 2 … Draw a directed edge from each left 2-mer to corresponding right
    2-mer"; "Each edge … corresponds to a k-mer from the input."
  - p.8-9: AAABBBBA "we get a multiedge"; "De Bruijn graph is a directed multigraph G(V, E)";
    indegree/outdegree definitions.
  - p.10: "Node is balanced if indegree equals outdegree"; "semi-balanced if indegree differs
    from outdegree by 1"; "A directed, connected graph is Eulerian if and only if it has at
    most 2 semi-balanced nodes and all other nodes are balanced" (cites Jones & Pevzner §8.8);
    "Eulerian walk visits each edge exactly once."
  - p.11: AAABBBA k=3 — "AA → AA → AB → BB → BB → BA"; "AA and BA are semi-balanced, AB and BB
    are balanced ⇒ Yes [Eulerian]."
  - p.15: "With perfect sequencing, this procedure always yields an Eulerian graph" — left-end
    k-1-mer one extra outgoing, right-end one extra incoming, all others balanced.
  - p.16: `chop`: `for i in xrange(0, len(st)-(k-1)): yield st[i:i+k]`; `km1L,km1R =
    kmer[:-1], kmer[1:]`; edge `G.setdefault(nodeL,[]).append(nodeR)` — confirms iteration
    bound i ∈ [0, |r|−k] and prefix/suffix split.
  - p.18-19/22: spelling rule `superstring = path[0] + ''.join(map(lambda x: x[-1], path[1:]))`;
    `DeBruijnGraph(["a_long_long_long_time"],5)` → 18-node walk spelling
    `a_long_long_long_time`; `to_every_thing_turn_turn_turn_there_is_a_season` reconstructs
    exactly at k=4, mis-orders at k=3 ("Due to repeats that are unresolvable at k = 3").
  - p.21: ZABCDABEFABY k=3 has two Eulerian walks ("only one … corresponds to original
    superstring"); a repeated (k-1)-mer joins edge-disjoint cycles.
  - p.24: "Gaps in coverage can lead to disconnected graph … Connected components are
    individually Eulerian, overall graph is not."
- **Jones & Pevzner (2004), Theorems 8.1/8.2** — already extracted in the unit's Evidence;
  Theorem 8.2 (Eulerian path iff ≤2 semibalanced vertices, all others balanced) is exactly the
  condition Langmead restates and the code's `ChooseEulerStart` realizes.
- **Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987-991** — publication identity and
  thesis (k-mers as edges between (k-1)-mer nodes; assembly = Eulerian-path, polynomial vs the
  NP-complete Hamiltonian/OLC formulation) confirmed; used for framing only.

### Formula check
- Node = (k-1)-mer; edge = k-mer, directed prefix→suffix; multigraph — matches Langmead p.5-9.
- Chop bound i ∈ [0, |r|−k]; #k-mers per read = |r|−k+1 — matches p.16.
- balanced/semi-balanced and Euler condition — matches p.10 and J&P Thm 8.1/8.2 verbatim.
- Spelling `p₀ + last-char(pᵢ for i>0)` — matches p.18-19.

### Edge-case semantics
- Read length < k → 0 k-mers (chop bound empty): sourced (p.16). ✓
- Empty/null reads → empty result/graph: ASSUMPTION-2 (no source prescribes; trivial identity).
- k < 2 → reject (nodes would be empty): follows directly from the (k-1)-mer node definition. ✓
- Repeated k-mer → multiedge: sourced (p.8). ✓
- Disconnected graph → one Eulerian walk per component (one contig each): sourced (p.24). ✓
- Multiple Eulerian walks at a repeat ≥ k-1 → reconstruction not unique; only invariants
  asserted, not a specific walk: ASSUMPTION-1, consistent with p.21.

### Independent cross-check (numbers)
Re-derived an independent Hierholzer DBG assembler in Python (separate from the repo code) and
reproduced every published value:
- `AAABBBA` k=3 → `AAABBBA`; edges AA→AA, AA→AB, AB→BB, BB→BB, BB→BA (5 = 7−3+1 k-mers).
- `AAABBBBA` k=3 → 6 edges (8−3+1), BBB occurs twice ⇒ two BB→BB edges (multiedge).
- `a_long_long_long_time` k=5 → `a_long_long_long_time` (matches Langmead p.18 walk verbatim).
- `to_every_thing_turn_turn_turn_there_is_a_season` k=4 → input exactly (Langmead p.22).
- `ATGGCGTGCA` k=4 → `ATGGCGTGCA` (all 4-mers distinct, no 3-mer repeat ⇒ unique walk).

### Findings / divergences
None. The description (TestSpec, Evidence, algorithm doc) is faithful to the authoritative
sources. Note on `to_every` at k=3: Langmead's printed mis-ordered result comes from his
particular arbitrary walk; under a deterministic insertion-order Hierholzer walk (this library's
choice) the k=3 input is recovered. The spec correctly anticipates this (ASSUMPTION-1) and the
k=3 test (C1) asserts only the branching structure, not a string. Stage A = **PASS**.

## Stage B — Implementation

### Code path reviewed (file:line)
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs`
- `AssembleDeBruijn` L104-125 (null/empty guard, build → reconstruct → MinContigLength filter → stats)
- `BuildDeBruijnGraph` L366-402 (k<2 throw, null guard, chop bound `i ≤ read.Length−k`, prefix/suffix split, multimap append)
- `ReconstructContigs` L413-455 (degree tallies, weakly-connected components, per-component walk)
- `BuildUndirectedAdjacency` L457-474, `CollectComponent` L476-498
- `ChooseEulerStart` L506-524 (semi-balanced source out−in=1, else lexicographically-smallest node)
- `SpellEulerianWalk` L532-569 (iterative Hierholzer, edges consumed in insertion order, reverse, spell)

### Formula realised correctly? (evidence)
- Chop bound `for (int i = 0; i <= read.Length - k; i++)` = Langmead's `range(0, len−(k−1))`. ✓
- `prefix = read.Substring(i, k-1)`, `suffix = read.Substring(i+1, k-1)` = `kmer[:-1], kmer[1:]`. ✓
- Multimap append preserves edge multiplicity (multigraph). ✓
- `ChooseEulerStart` returns the unique out−in=1 node (Thm 8.2 source of an Eulerian path); in a
  valid Eulerian-path graph there is exactly one, so DFS-order iteration is safe. ✓
- `SpellEulerianWalk` builds `walk[0] + concat(last char of walk[i>0])`. ✓
- `adjacency` is a per-call mutable copy; edges are consumed by a per-node cursor; components
  are edge-disjoint, so sharing one `adjacency` across the per-component loop is correct. ✓

### Cross-verification table recomputed vs code
Ran the real `AssembleDeBruijn` (temporary console harness against the live assembly) and the
NUnit fixture; all match the independent reference:

| Input | k | Reference (Python) | Code output | Source |
|-------|---|--------------------|-------------|--------|
| AAABBBA | 3 | AAABBBA | AAABBBA | Langmead p.11 |
| a_long_long_long_time | 5 | a_long_long_long_time | a_long_long_long_time | Langmead p.18 |
| to_every…season | 4 | input verbatim | input verbatim | Langmead p.22 |
| ATGGCGTGCA | 4 | ATGGCGTGCA | ATGGCGTGCA | construction+spelling |
| {ATGGCGTGCA, GATTACAGGTC} | 4 | two contigs = each read | two contigs = each read | Langmead p.24 (disconnected) |
| same, MinContigLength=11 | 4 | {GATTACAGGTC} | {GATTACAGGTC} | contract §3.1 (filter) |
| AAABBBA edges | 3 | 5 edges | 5 edges | Langmead p.7 |
| AAABBBBA BB→BB | 3 | 2 (multiedge), 6 edges | 2, 6 | Langmead p.8 |

### Variant/delegate consistency
`AssembleDeBruijn` consumes `BuildDeBruijnGraph`; no `*Fast`/instance variants. Determinism test
confirms identical output across runs (lexicographically-smallest start + insertion-order edges).

### Test quality audit (HARD gate)
- **Sourced expectations, not code echoes:** every exact expected value (M1-M7 strings/edges,
  multiedge count, edge counts, disconnected contigs, filter) traces to Langmead (verified from
  the fetched PDF this session) or to the independent Hierholzer reference — not to the code's
  own output. A deliberately-wrong walk/spelling would fail M4-M7 and the new disconnected test.
- **No green-washing:** assertions use exact `Is.EqualTo` on strings, edge lists, node sets, and
  counts; no Greater/AtLeast/range substitutions where an exact value is known. S3 (INV-05) uses
  `Contains` because INV-05 is *defined* as a substring property, not because an exact value was
  weakened. No skips, no widened tolerances, no expected values bent to match output.
- **Cover all the logic — GAP FOUND AND FIXED:** the original fixture (17 tests) covered all M/S/C
  spec rows but did **not** exercise three documented Stage-A behaviours/branches:
  1. **Disconnected graph → one contig per weakly-connected component** (Langmead p.24; algorithm
     step 3-4; the `BuildUndirectedAdjacency`/`CollectComponent`/per-component loop). Untested.
  2. **`MinContigLength` filter** (contract §3.1; the `Where(c => c.Length >= MinContigLength)`
     branch). Untested (all tests used MinContigLength:1).
  3. **`BuildDeBruijnGraph(null, k)`** null guard (contract §3.3). Untested.
  Added three tests with externally-derived exact values:
  - `AssembleDeBruijn_DisconnectedGraph_OneContigPerComponent` — `{ATGGCGTGCA, GATTACAGGTC}` k=4
    → exactly those two contigs (each has all-distinct 3-mers and they share no 3-mer node, so
    the value was hand-derived from the unique-walk + components rule, then cross-checked).
  - `AssembleDeBruijn_MinContigLength_FiltersShortContigs` — same reads, MinContigLength=11 →
    only `GATTACAGGTC` (length 11) survives, the length-10 contig is dropped.
  - `BuildDeBruijnGraph_NullReads_ProducesEmptyGraph` — empty graph, no exception.
- **Honest green:** FULL unfiltered suite = **6497 passed, 0 failed**, 1 pre-existing skip
  (`MFE_Benchmark_AllScenarios`, unrelated). DBG fixture now 20 tests, all pass. `dotnet build`
  0 errors (the 4 warnings are pre-existing in `ApproximateMatcher_EditDistance_Tests.cs`,
  untouched by this session).

### Findings / defects
- **Test-coverage gap (Stage B defect, FIXED in session):** disconnected-component output,
  `MinContigLength` filtering, and `BuildDeBruijnGraph` null handling were untested branches with
  defined, sourced behaviour. Fixed by adding three exact-value tests (above). No implementation
  bug found — the code already produced the correct sourced values for all three.

## Verdict & follow-ups
- **Stage A: PASS.** Description faithful to Langmead + Jones & Pevzner + Compeau et al.
- **Stage B: PASS-WITH-NOTES.** Code realises the validated formula exactly for every cross-check;
  the only issue was a test-coverage gap (three documented branches), now closed with
  sourced-value tests.
- **End-state: CLEAN.** No implementation defect; the test gap was completely fixed; full suite
  green (Failed: 0), build warning-free on touched files.
- No follow-ups.
