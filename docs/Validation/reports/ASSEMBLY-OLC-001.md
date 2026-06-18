# Validation Report: ASSEMBLY-OLC-001 — Overlap-Layout-Consensus assembly

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `SequenceAssembler.AssembleOLC(reads, parameters)`,
  `SequenceAssembler.FindAllOverlaps(reads, minOverlap, minIdentity)` (+ cancellable overload),
  `SequenceAssembler.FindOverlap(seq1, seq2, minOverlap, minIdentity)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm
- **Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987–991, DOI 10.1038/nbt.2023** — overlap
  graph (one node per read; directed edge A→B when a suffix of A matches a prefix of B above a
  minimum threshold); layout = Hamiltonian path; finding a Hamiltonian path is NP-complete, hence
  heuristic layout. Confirmed against the open-access PDF cited in the Evidence file.
- **Langmead (JHU), "Overlap Layout Consensus assembly" lecture notes** — three stages
  Overlap/Layout/Consensus (p.4); suffix-prefix overlap with minimum length `l` (p.5); "report
  only the longest suffix/prefix match" (p.10); transitive reduction + non-branching contig
  emission (p.21–25); consensus = per-column majority vote (p.28); all-pairs O(N²) vs suffix-tree
  O(N+a) (p.10, p.16).
- **Langmead (JHU), "Assembly & Shortest Common Superstring" lecture notes** — first law of
  assembly (p.16); GTACGTACGAT 6-mer overlap graph, edges of length ≥ 4 with weights 4 and 5
  (p.23–25); greedy-SCS and its sub-optimality (p.45, p.57); repeats foil assembly (p.58–62).

The cited sources are authoritative and say what the Evidence/TestSpec claim. Formulae are
parameter-free (edge weight = overlap length; thresholds `MinOverlap` `l` and `MinIdentity` are
caller-supplied); there are no scoring tables or biological constants to source-check.

### Formula / definition check
- Overlap = longest suffix of A equal to a prefix of B, length ≥ `l`, report only the longest —
  matches Langmead OLC p.5/p.10 and SCS p.16/p.23.
- Coordinates: 0-indexed reads; `Position2` always 0 (prefix of seq2) — internally consistent and
  matches the suffix-prefix definition.
- Edge-case semantics: edgeless graph → singletons (Compeau overlap-graph definition); empty input
  → empty result (explicitly an ASSUMPTION, no source — acceptable trivial identity).

### Independent cross-check (numbers I recomputed this session)
Hand-derived (Python, from the suffix-prefix definition only — **not** from the repo code):

- **GTACGTACGAT 6-mers, minOverlap 4** → exactly the 12 directed edges with weights {4,5}; every
  edge and length matches the Evidence table and Langmead SCS p.24–25 (re-derived). 12 edges total.
- `FindOverlap("CTCTAGGCC","TAGGCCCTC")` → longest suffix-prefix = 6 (`TAGGCC`), pos1 = 3
  (Langmead OLC p.5).
- Chain `AAAAACCCCC`,`CCCCCGGGGG`,`GGGGGTTTTT` merged at 5-overlaps → `AAAAACCCCCGGGGGTTTTT`,
  length 20.
- Identity case: suffix `ACGTACGT` vs prefix `ACGTACCT` = 7/8 = 0.875 (1 mismatch), and the longer
  windows (len 9–12) all score < 0.85, so the length-8 overlap is the one reported.
- Boundary: `AAAACGTT`/`CGTTGGGG` longest suffix-prefix = 4 (`CGTT`); `ACGTACGT`/`CGTAAAAA` longest
  exact match = 3 (so no edge at minOverlap 4).

All cross-checks reproduce the Evidence values independently.

### Findings / divergences
Stage A is sound. The description is honest about the two intentional simplifications (greedy
best-successor chaining instead of transitive-reduction + Hamiltonian search; consensus realized
as chain concatenation rather than column majority vote), both correctly justified by NP-completeness
and exact-overlap equivalence respectively.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs`:
- `AssembleOLC` (L61–83) — overlap → `BuildContigsFromOverlaps` → MinContigLength filter (L77) → stats.
- `FindAllOverlaps` (L142–165) and cancellable overload (L177–212) — all ordered pairs i≠j, delegate
  to `FindOverlap`.
- `FindOverlap` (L221–242) — scans overlap lengths from `min(len1,len2)` down to `minOverlap`,
  returns the first window meeting `minIdentity` (i.e. the longest), pos1 = len1−L, pos2 = 0.
- `CalculateIdentity` (L247–260) — case-insensitive fraction of matching positions.
- `BuildContigsFromOverlaps` (L262–346) — greedy best-successor chaining; singleton fallback gated
  by MinContigLength.

### Formula realised correctly?
Yes. `FindOverlap` realises the longest-suffix-prefix definition exactly (descending scan returns
the longest qualifying window). `CalculateIdentity` is the matching-fraction over the window,
case-insensitive (matches doc §3.3). `FindAllOverlaps` excludes self-overlaps and reports one edge
per ordered pair. The greedy layout merges `successor[overlap:]`, which equals the superstring merge.

### Cross-verification table recomputed vs code (test run)
| Case | External value | Code (test) |
|------|----------------|-------------|
| GTACGTACGAT 6-mers @minOv4 | 12 edges, weights {4,5} (exact set) | exact match (M1 green) |
| FindOverlap CTCTAGGCC/TAGGCCCTC | len 6, pos1 3, pos2 0 | match (M3 green) |
| Chain merge | `AAAAACCCCCGGGGGTTTTT` (20) | match (M4 green) |
| 7/8 identity gate | accept @0.85, reject @0.95 | match (S1 green) |
| minOverlap boundary | accept @4, none @5 | match (S2 green) |
| below threshold (3<4) | no edges | match (S3 green) |
| MinContigLength discard | drop 5-base, keep 10-base | match (M5b green, added this session) |

### Variant / delegate consistency
Cancellable `FindAllOverlaps` (CancellationToken.None) yields the identical edge set (S5 green) and
throws `OperationCanceledException` on a pre-cancelled token (S6 green). Both delegate to the same
`FindOverlap`.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1–M5, S1–S3, M3 assert exact externally-derived values (12-edge set,
  exact contig strings, exact overlap length+positions, exact accept/reject). A deliberately-wrong
  implementation (e.g. shortest-overlap, off-by-one positions, wrong merge) would fail these. PASS.
- **No green-washing:** no weakened assertions on the numeric cases; exact equality used where exact
  values are known. The two property tests (S4 bounds, C1 ≥-longest) are inequalities *by design* —
  they assert the documented heuristic-limitation contract (INV-04 / ASM-02), which is the strongest
  honestly-assertable statement for a greedy layout per the sources (Langmead SCS p.57: greedy is not
  optimal). Acceptable.
- **Coverage:** all three public methods + the cancellable overload + all Stage-A branches are
  exercised (exact graph, no-self, below-threshold, identity gate, boundary, unambiguous chain,
  singletons, empty, null, case-insensitive, repeat, cancellation). **Gap found & fixed:** the
  `MinContigLength` discard branch (doc §3.2) had no test; added **M5b**
  (`AssembleOLC_MinContigLength_DiscardsShortContigs`) locking the contract (5-base read discarded,
  10-base kept, TotalLength 10, TotalReads 2).
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6494` (was 6493 + the new M5b);
  `dotnet build` 0 warnings / 0 errors.

Result: **test-quality gate PASS** after adding M5b.

### Findings / defects
No code defect. The implementation faithfully realises the validated description. The PASS-WITH-NOTES
on Stage B reflects only the documented (sourced) intentional simplifications, not a divergence:
greedy layout (not Hamiltonian-path optimal) and concatenation-consensus (equivalent to majority
vote for exact-overlap chains). One test-coverage gap (MinContigLength) was found and fixed in
session.

## Verdict & follow-ups
- **Stage A:** PASS — description and formulae independently confirmed against Compeau et al. (2011)
  and Langmead's OLC/SCS notes; the 12-edge GTACGTACGAT graph and all overlap/merge values re-derived
  by hand.
- **Stage B:** PASS-WITH-NOTES — code is correct for its stated contract; intentional greedy-layout /
  concatenation-consensus simplifications are documented and sourced.
- **End-state:** ✅ CLEAN — no defect; the single test-coverage gap (MinContigLength filter) was
  completely fixed in session with a sourced exact-value test; full suite green.
- **Test-quality gate:** PASS (after adding M5b).
