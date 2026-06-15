# Validation Report: GENOMIC-MOTIFS-001 — Known Motif Search

- **Validated:** 2026-06-16   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindKnownMotifs(DnaSequence, IEnumerable<string>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Summary

GENOMIC-MOTIFS-001 ("Known Motif Search") is the classical **exact set-matching** problem:
given a subject sequence T and a set of fixed query motifs, report for each motif the 0-based
start positions of **all** its occurrences, including overlapping ones. It is *not* de novo
motif discovery / the (l,d) motif problem — those concern enriched/approximate motifs and are
covered by other units (MOTIF-DISCOVER-001, MOTIF-SHARED-001, etc.). This unit is exact
substring matching of a curated motif set against one reusable suffix-tree index.

## Stage A — Description

### Sources opened & what they confirm

1. **Rosalind "SUBS — Finding a Motif in DNA"** (https://rosalind.info/problems/subs/, retrieved
   this session via WebFetch). The problem is "find all starting positions where a substring
   appears within a larger DNA string", reporting **every** occurrence including overlapping ones.
   Sample: motif `ATAT` in `GATATATGCATATACTT` → output `2 4 10` (1-based). Positions 2 and 4 are
   overlapping. This directly confirms INV-01 (exact match), INV-02 (overlaps all reported), and
   the "set of all start positions" model. (Rosalind uses 1-based; the repo uses 0-based — see
   convention note below.)
2. **Wikipedia "EcoRI"** (https://en.wikipedia.org/wiki/EcoRI, retrieved this session). Quote:
   "The nucleic acid recognition sequence where the enzyme cuts is G↓AATTC". Confirms the
   biological motif `GAATTC` (length 6) used in dataset M2. (The Evidence file cited the
   "Restriction site" article; I cross-checked the dedicated EcoRI article and it agrees.)
3. **Gusfield exact-matching definition** (cited in Evidence via Tufts COMP 150GEN). The
   exact-matching problem is "find all occurrences of a pattern string P ... in a text string T";
   `P=aaa`, `T=aaaaa` → three overlapping occurrences at 0-based 0,1,2. Independently corroborated
   by the Rosalind overlapping example above.

### Formula / model check

The model is `{ i : T[i..i+n-1] == P }`, 0-based, overlaps included, per-motif independent sets,
absent motifs omitted, keys upper-cased. Every clause matches the cited sources. No log base,
normalization, or scoring is involved (exact match, not scored/approximate), so there is no
formula subtlety to mis-state.

### Definitions & conventions

- **Coordinate base:** 0-based start positions. The doc/spec state this explicitly and
  consistently. Rosalind is 1-based, but that is a presentation choice; the underlying set
  {1,3,9} (0-based) = {2,4,10} (1-based) is identical. Internal consistency holds.
- **Case handling:** `DnaSequence` upper-cases on construction (src line 30) and validates ACGT;
  `FindKnownMotifs` upper-cases each motif and keys the result by the upper-cased form
  (Biopython processes DNA upper-cased). Consistent and sourced.
- **Overlap convention:** all occurrences reported (not greedy non-overlapping). Sourced.

### Edge-case semantics

- Empty motif set → empty dictionary (no motif to search). Defined.
- Absent motif → omitted (empty occurrence set). Defined, INV-04.
- Empty/whitespace motif → skipped. This is an **assumption** (documented): no authoritative
  source defines an empty-pattern occurrence set as meaningful; a suffix-tree `FindAllOccurrences("")`
  would return every position, which is not a motif. The repo policy (skip) is reasonable and
  matches `FindMotif`. Accepted as an API-shape policy, not an algorithmic claim.
- `motifs` null → `ArgumentNullException`. Standard input validation.

### Independent cross-check (numbers)

Reference implementation: a naive overlapping exact-match scan
(`[i for i in range(len(T)-len(P)+1) if T[i:i+len(P)]==P]`), run this session:

| Dataset | Motif | Reference (0-based) | Spec expected |
|---------|-------|---------------------|----------------|
| `AAAAA` | `AAA` | `[0,1,2]` | `[0,1,2]` ✓ |
| `GAATTCAAAGAATTC` | `GAATTC` | `[0,9]` | `[0,9]` ✓ |
| `ACGTACGTAA` | `ACGT` | `[0,4]` | `[0,4]` ✓ |
| `ACGTACGTAA` | `AA` | `[8]` | `[8]` ✓ |
| `ACGTACGTAA` | `TTT` | `[]` | omitted ✓ |
| `AAAAA` | `A` | `[0,1,2,3,4]` | `[0,1,2,3,4]` ✓ |
| `ACAA` | `AA` | `[2]` | `[2]` ✓ |
| `GAATTC` | `GAATTC` | `[0]` | `[0]` ✓ |
| `ACGT` | `TTTT` | `[]` | omitted ✓ |
| `GATATATGCATATACTT` (Rosalind) | `ATAT` | `[1,3,9]` = 1-based `2 4 10` | matches Rosalind ✓ |

All values trace to externally-retrieved sources (Rosalind sample output; EcoRI site) and an
independent reference implementation, not to the repo's own code.

### Findings / divergences

None. Description is biologically and mathematically correct. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed

- `GenomicAnalyzer.FindKnownMotifs` — `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:190-228`.
- `SuffixTree.FindAllOccurrences` — `src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs:81-93`
  (and span overload :148). Collects all leaves under the matched node via `CollectLeaves`
  (DFS), so all start positions (overlapping included) are returned, unsorted.
- `DnaSequence` ctor — `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:22-31`
  (upper-case + ACGT validation); `SuffixTree` property :48 (built once, cached).

### Formula realised correctly?

Yes. For each non-empty motif: upper-case → `tree.FindAllOccurrences(normalized)` (all occurrences,
overlaps included) → if non-empty, `OrderBy(p => p)` (ascending sort for deterministic contract,
INV-03) → store under upper-cased key. Null guard throws `ArgumentNullException`. Empty/whitespace
skip via `string.IsNullOrWhiteSpace` (also covers a null element). Duplicate upper-cased keys
collapse (`result.ContainsKey` guard). Every Stage-A clause is realised.

### Cross-verification table recomputed vs code

The 12 unit tests assert the exact-value table above and all pass (see below). The values were
sourced independently first, then confirmed equal to actual code output — not the reverse.

### Variant/delegate consistency

Only one public entry point (`FindKnownMotifs`); no `*Fast`/instance variant. It delegates to the
shared `SuffixTree` used across the library. Consistent.

### Numerical robustness

Integer positions only; no floating point, overflow, or div-by-zero. Empty subject → empty suffix
tree → no occurrences → empty result. Non-ACGT motif characters: blocked upstream because the test
inputs are ACGT; a motif with a non-ACGT char would simply not match (returns no positions), per
documented limitation. Robust on the stated domain.

### Test quality audit (HARD gate)

File: `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindKnownMotifs_Tests.cs` (12 tests).

- **Sourced, not code-echoes:** every expected value matches the externally-sourced /
  independently-recomputed table above (Rosalind overlap definition, EcoRI site, naive reference).
  A deliberately-wrong implementation (e.g. greedy non-overlapping, or DFS-unsorted output) would
  FAIL M1/S4 (overlap) and M4 (ascending) respectively. Not code echoes.
- **No greenwashing:** all assertions use exact `Is.EqualTo(new[]{…})` or `Is.Empty` /
  `ContainsKey(...) Is.False`. No `Greater`/`AtLeast`/`Contains`-only/ranges. No widened
  tolerances, no skip/ignore. M4 pairs `Is.Ordered.Ascending` with an exact-value check.
- **Coverage of all branches:** overlapping (M1, S4), multi-motif per-motif sets (M3), absent
  omitted (M5 + TTT in M3), ascending sort (M4), empty motif set (S1), case normalization +
  upper-cased key (S2), empty/whitespace skip (S3), duplicate-key dedup (C1), null throws (Edge).
  This exercises every branch of the method body. INV-01..INV-05 each have a covering test.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6573`; FindKnownMotifs subset
  `Failed: 0, Passed: 12`; project builds with **0 errors** (4 pre-existing NUnit2007 warnings are
  in an unrelated file, ApproximateMatcher_EditDistance_Tests.cs — not touched this session).

**Test-quality gate: PASS.**

### Findings / defects

None. Code faithfully realises the validated description; tests are sourced, exact, and complete.
No code or test change was required this session.

## Verdict & follow-ups

- **Stage A: PASS** — description is correct against Rosalind SUBS, Wikipedia EcoRI, and Gusfield
  exact-matching definition; overlapping-occurrence and 0-based conventions are standard and
  internally consistent.
- **Stage B: PASS** — implementation computes the exact-matching set with overlaps via the suffix
  tree and sorts for determinism; all 12 tests pass; full suite green (6573/0); build clean.
- **End-state: CLEAN** — no defect found; algorithm fully functional. No follow-ups.
