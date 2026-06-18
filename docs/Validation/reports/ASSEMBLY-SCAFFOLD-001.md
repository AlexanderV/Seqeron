# Validation Report: ASSEMBLY-SCAFFOLD-001 — Scaffolding (joining ordered contigs with N-gaps)

- **Validated:** 2026-06-15   **Area:** Assembly (Extended Assembly)
- **Canonical method(s):** `SequenceAssembler.Scaffold(IReadOnlyList<string> contigs, IReadOnlyList<(int,int,int)> links, char gapCharacter='N')`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of the repo)

1. **NCBI AGP Specification v2.1** — https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ (WebFetch, 2026-06-15).
   Confirmed verbatim: *"Gap lengths must be positive. Negative gaps and gap lines with zero length are not valid."* and
   *"For negative gaps, or gaps of unknown size, use U as the component_type and 100 as the gap size, since 100 is the
   GenBank/EMBL/DDBJ standard for gaps of unknown size."*
2. **Jackman et al. (2017), ABySS 2.0 scaffold paper** — https://github.com/sjackman/abyss-scaffold-paper/blob/master/scaffold.md
   (WebFetch, 2026-06-15). Confirmed verbatim: *"The sequences of the vertices in a path are concatenated, interspersed
   with gaps represented by a run of the character N, whose length corresponds to the estimate of the distance between
   those two contigs."* and the negative-estimate/overlap rule: *"It is possible that the distance estimate is negative,
   indicating that the two contigs should in fact overlap. If such an overlap is indeed found in the contig overlap graph,
   the two contigs are merged."* The paper text states **no** minimum gap size when contigs are placed close together.
3. **Sahlin et al. (2012)** — https://academic.oup.com/bioinformatics/article/28/17/2215/246308 — corroborates that a
   negative gap (overlap) is a common, expected input class for de Bruijn assemblers.
4. **Pop et al. (2004), Bambus** (via Wikipedia primary citation) — scaffold = a path of distinct contigs; greedy
   scaffolders join "most links first".

### Formula / model check
- Scaffold construction (concatenate path contigs, interspersing an N-run of length = distance estimate) matches Jackman
  et al. (2017) verbatim. INV-01, INV-03, INV-05 follow directly.
- Non-positive (zero/negative) estimate → 100 gap characters matches NCBI AGP v2.1 (unknown-gap = 100). INV-02 follows.
- Each contig in exactly one scaffold (path of distinct contigs) matches Pop et al. (2004). INV-04 follows.

### Edge-case semantics check
- Zero gap → invalid per AGP → unknown default (100 N): sourced and correct.
- Negative gap → overlap, unresolved here → unknown default (100 N): the *constant 100* is source-backed (AGP); the
  *scoping decision* to fall back rather than resolve the overlap is a documented assumption (Assumption Register #1).
  This is correctly disclosed in both the algorithm doc (§5.3 "Intentionally simplified") and the TestSpec.
- Out-of-range / self link ignored, null inputs throw: standard, consistent with the sibling `MergeContigs`.

### Independent cross-check (hand-computed against sources, not against code)
- M1: `["ACGT","TTGG","CCAA"]`, `[(0,1,3),(1,2,2)]` → `ACGT`+`NNN`+`TTGG`+`NN`+`CCAA` = `ACGTNNNTTGGNNCCAA`,
  length 4+3+4+2+4 = 17 (Jackman rule). Matches.
- M3: gap −5 → 100 N → `AAAA`+100N+`TTTT`, length 108 (AGP). Matches.
- M4: gap 0 → invalid → 100 N (AGP). Matches.

### Findings / divergences (Stage A)
- The only divergence from the cited theory is documented and scoped: Bambus joins "most links first", whereas this unit
  follows **input link order** (it consumes pre-ordered links and does not perform link-count ranking or overlap
  resolution / contig orientation). All three are explicitly listed as intentional simplifications. No biological or
  mathematical error in the description.

**Stage A verdict: PASS.**

## Stage B — Implementation

- **Code path reviewed:** `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:679-744` plus the
  `UnknownGapLength = 100` constant at line 651.
- **Formula realised correctly:** gap length rule `gapLength = gapSize > 0 ? gapSize : UnknownGapLength` (line 733) is
  exactly INV-01/INV-02. The path-follow loop appends `gapCharacter` × gapLength then `contigs[nextIndex]` (lines 733-736),
  realising the concatenate-with-N-run model verbatim. `used` HashSet enforces INV-04. Single-contig scaffolds for
  unreached contigs in ascending index order (outer `for`, line 709). Null guards (lines 684-685); empty → empty (687).
- **Cross-verification recomputed vs code:** traced M1, M3, M4, S1, and the two added branch cases by hand against the
  code — all match the externally sourced values above.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1/M2 lock exact strings/run-lengths from the Jackman rule; M3/M4 lock exactly 100 N from
  the AGP standard (a wrong `Math.Max(1, gapSize)` implementation — the pitfall called out in the TestSpec — would fail
  these). These tests would fail a deliberately-wrong implementation.
- **No green-washing:** all assertions are exact (`Is.EqualTo`), no widened tolerances, nothing skipped/ignored.
- **Coverage of logic / branches:** original suite (13 tests) covered M1–M7, S1–S5, C1. Two real logic branches were
  **uncovered**: (a) the inner "skip an already-placed successor and stop the path" branch when a contig has a forward
  link to a used contig mid-path, and (b) the multi-forward-link tie-break (first declared link wins). I added two tests:
  - `Scaffold_SuccessorAlreadyPlaced_SkippedKeepingContigsDistinct` — `["AA","CC","GG"]`, `[(0,1,2),(1,2,3),(2,1,4)]` →
    `"AANNCCNNNGG"` (1 scaffold, length 11), derived from INV-03/04/05.
  - `Scaffold_MultipleForwardLinks_FollowsFirstDeclared` — `["AA","CC","GG"]`, `[(0,1,2),(0,2,3)]` → `["AANNCC","GG"]`,
    locking the documented input-order tie-break and the unreached-contig length-1 scaffold.
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6531** (was 6529 + 2 new); `dotnet build` 0 errors
  (4 pre-existing warnings in unrelated files, none introduced here).

### Findings / defects (Stage B)
- No defect. The implementation faithfully realises the validated description; the only gap was missing branch coverage
  in the tests, now fixed with sourced/invariant-derived expectations.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS.** End-state: **CLEAN.**
- Test-quality gate: **PASS** (sourced exact assertions; all logic branches now exercised; full suite green).
- No outstanding defects. Documented simplifications (link-count ranking, overlap resolution, contig orientation) are
  out of scope and clearly disclosed; they are not correctness defects for the stated contract.
