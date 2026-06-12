# Validation Report: RESTR-DIGEST-001 — Restriction Digest Simulation

- **Validated:** 2026-06-12   **Area:** MolTools
- **Canonical method(s):** `RestrictionAnalyzer.Digest(seq, enzymeNames)`, `GetDigestSummary`, `CreateMap`, `AreCompatible`, `FindCompatibleEnzymes`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs`)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_Digest_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Restriction digest"** — confirms fragments are produced by cleavage at recognition
  sites, and the **fragment-sum invariant**, quoting Addgene-style verification verbatim:
  *"The sum of the individual fragments should equal the size of the original fragment, and each
  digest's fragments should also sum up to be the same size as each other."* (Invariant #2).
- **Wikipedia "Restriction map"** — confirms maps are built from single and double (multi-enzyme)
  digests and site positions; supports CreateMap's site-position model. (Does not itself state
  the unique-cutter or linear/circular fragment-count rules.)
- **Independent biology sources (search consensus, e.g. life.illinois.edu digest assignment,
  ScienceDirect "Restriction Digest" overview, Quora/Doubtnut/Brainly worked answers)** — confirm
  the canonical fragment-count rules:
  - **Linear molecule:** n cuts → **n + 1 fragments** ("ribbon cut at 3 points → 4 segments").
  - **Circular molecule:** n cuts → **n fragments** (no ends; need 2 cuts for 2 fragments).

### Formula / semantics check
- Fragment generation = partition the sequence at each forward-strand cut position; boundaries at
  the molecule ends (0 and L). Fragment k spans `[cut_{k-1}, cut_k)`; sizes sum to L. ✔
- Forward-strand-only counting for palindromic Type II enzymes (avoid double-counting the same
  physical site that matches on both strands) — standard convention. ✔
- Overhang compatibility (blunt↔blunt always; sticky requires same overhang **type** AND same
  overhang **sequence**; 5′ and 3′ overhangs never compatible even with identical strings) matches
  Wikipedia "Sticky and blunt ends". ✔

### Worked examples (hand-computed)
1. **Linear, 1 cut — EcoRI (G↓AATTC).** Sequence `AAAGAATTCAAA` (L=12). Site at index 3, cut at
   3+1 = 4. Fragments: `[0,4)`=`AAAG` (4), `[4,12)`=`AATTCAAA` (8). Count = 1 cut + 1 = **2**;
   sum = 4+8 = **12** = L. ✔
2. **Linear, 2 cuts — EcoRI.** Sequence `GAATTCAAAGAATTCAAA` (L=18). Sites at 0 and 9, cuts at 1
   and 10. Fragments: `[0,1)`=1, `[1,10)`=9, `[10,18)`=8. Count = 2+1 = **3**; sum = 1+9+8 = **18**
   = L. ✔
3. **Circular, 1 cut → 1 linear fragment = full length.** This topology is **not implemented** (see
   notes). The library models only linear molecules.

### Findings / divergences (Stage A)
- **NOTE (scope):** The TestSpec, invariants, and implementation cover **linear molecules only**.
  The circular-molecule rule (n cuts → n fragments; 1 cut on a plasmid → 1 linear fragment of full
  length) is **not modeled**. This is a documented scope limitation, *not* a correctness defect:
  every invariant in the spec (Invariant #1: "k cut positions → k+1 fragments") is explicitly the
  linear rule, and no test or API claims circular behavior. Stage A passes with this note.

## Stage B — Implementation

### Code path reviewed
`RestrictionAnalyzer.Digest` (RestrictionAnalyzer.cs:237–296):
- Collects forward-strand cut positions into a `SortedSet<int>` (dedup + sorted) — line 244–257.
- No cuts → returns the whole sequence as one fragment (LeftEnzyme/RightEnzyme null) — line 259–270
  (Invariant #7). ✔
- Builds `cuts = [0] + sortedCutPositions + [Length]` and emits fragment `[cuts[i], cuts[i+1])` —
  line 273–295. n cut positions ⇒ `cuts.Count = n+2` ⇒ **n+1 fragments** (linear rule). ✔
- Fragment sizes are contiguous half-open intervals spanning `[0, L]`, so they **sum to L** by
  construction. ✔
- First fragment `LeftEnzyme=null`, last `RightEnzyme=null` (Invariant #5); `FragmentNumber=i+1`
  sequential (Invariant #3); start positions strictly increasing (Invariant #4). ✔
- `length > 0` guard (line 283) means a cut at the very start (position 0) or end (position L)
  yields no zero-length fragment (Invariant #6). ✔

### Cross-verification (recomputed vs code, matches tests)
| Input | Enzyme(s) | Cuts | Fragments (len) | Sum | Test |
|-------|-----------|------|-----------------|-----|------|
| AAAGAATTCAAA (12) | EcoRI | {4} | AAAG(4), AATTCAAA(8) | 12 | M1/M2 ✔ |
| GAATTCAAAGAATTCAAA (18) | EcoRI | {1,10} | 1,9,8 | 18 | M4 ✔ |
| GAATTC…×3 (27) | EcoRI | {1,10,19} | 1,9,9,8 | 27 | props ✔ |
| AAAGAATTCAAAGGATCCAAA (21) | EcoRI+BamHI | {4,13} | 4,9,8 | 21 | M9 ✔ |
| AAATCGATCAAA (12) | TaqI+MboI | {4,5} | 4,1,7 | 12 | S3 (adjacent) ✔ |
| AAAAAAAAAAAA (12) | EcoRI | {} | whole(12) | 12 | M3 ✔ |

`GetDigestSummary`: sizes sorted descending, bounds, average, enzyme list — verified
(`GetDigestSummary_ReturnsCorrectSummaryWithInvariants`, sizes [9,8,1], avg 6.0). ✔
`CreateMap`: forward-strand-only TotalSites and UniqueCutters; NonCutters; search-all-when-empty —
verified. ✔
Compatibility (blunt/sticky/cross-type/symmetry/enzyme DB vs Wikipedia 20-case table) — verified. ✔

### Edge cases
- No cut site (linear) → 1 fragment = whole sequence ✔ (M3, S2).
- Adjacent sites → 1-bp middle fragment, all lengths > 0 ✔ (S3).
- Multi-enzyme digest → merged sorted cut set, correct enzyme attribution ✔ (M9).
- Null sequence / no enzymes → `ArgumentNullException` / `ArgumentException` ✔ (M10, M11).
- Cut exactly at start/end → `length > 0` guard prevents empty fragments. ✔

### Findings / defects (Stage B)
- No correctness defect: cut-position-based fragmentation, linear n+1 count, and the
  sum-equals-length invariant are all correctly realized.
- **NOTE (scope, same as Stage A):** No circular/plasmid topology. `Digest` always frames cuts
  with `[0 … L]`, i.e. linear-only. A circular digest would instead join the last and first
  segments around the origin (n cuts → n fragments). Not a bug against the stated spec, but a
  capability gap relative to the general biology; recorded as a known limitation, not fixed here
  (would require an API/contract change — a `topology`/`isCircular` parameter and new fragment-
  joining logic plus tests — which is out of scope for a validation pass and not requested by the
  TestSpec).

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (description correct for linear digests; circular topology out of
  scope and not claimed by the spec).
- **Stage B:** PASS-WITH-NOTES (implementation faithfully realizes the linear digest; same circular
  scope note).
- **End-state:** CLEAN — no defect found; no code change required. All linear-digest invariants
  (fragment count n+1, sum = L, boundary enzymes, ordering, positive length, no-cut whole-sequence)
  hold and are locked by exact-value tests.
- Tests: `--filter FullyQualifiedName~Digest` → 67 passed / 0 failed. Full suite: **4461 passed,
  0 failed** (baseline preserved; no code changed).

### Sources
- https://en.wikipedia.org/wiki/Restriction_digest
- https://en.wikipedia.org/wiki/Restriction_map
- https://en.wikipedia.org/wiki/Restriction_enzyme
- https://en.wikipedia.org/wiki/Sticky_and_blunt_ends
- https://www.addgene.org/protocols/restriction-digest/
- Linear/circular fragment-count consensus: life.illinois.edu digest assignment; ScienceDirect
  "Restriction Digest" overview; standard molecular-biology worked examples.
