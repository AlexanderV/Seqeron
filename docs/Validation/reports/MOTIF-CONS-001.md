# Validation Report: MOTIF-CONS-001 — Consensus Sequence from a Multiple Alignment

- **Validated:** 2026-06-15   **Area:** Matching
- **Canonical method(s):** `MotifFinder.CreateConsensusFromAlignment(IEnumerable<string>)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (this session)

1. **Rosalind — Consensus and Profile (CONS)**, https://rosalind.info/problems/cons/ (WebFetch 2026-06-15).
   - Consensus string = "A string of length n formed from our collection by taking the most common
     symbol at each position."
   - Profile matrix = "A 4 × n matrix P in which P₁,j represents the number of times that 'A' occurs
     in the jth position … P₂,j … C … and so on."
   - Sample input (7 strings, length 8): `ATCCAGCT`, `GGGCAACT`, `ATGGATCT`, `AAGCAACC`, `TTGGAACT`,
     `ATGCCATT`, `ATGGCACT`.
   - Sample output consensus: **`ATGCAACT`**; profile A=`5 1 0 0 5 5 0 0`, C=`0 0 1 4 2 0 6 1`,
     G=`1 1 6 3 0 1 0 0`, T=`1 5 0 0 0 1 1 6`.
   - Ties: "If several possible consensus strings exist, then you may return any one of them."
2. **Wikipedia — Consensus sequence**, https://en.wikipedia.org/wiki/Consensus_sequence (WebFetch 2026-06-15).
   - "The consensus sequence is the calculated sequence of most frequent residues … found at each
     position in a sequence alignment." Limitation noted: reduces variability to a single residue.
3. **Los Alamos HIV DB — Advanced Consensus Maker explanation**,
   https://hfv.lanl.gov/content/sequence/CONSENSUS/AdvConExplain.html (WebFetch 2026-06-15).
   - Consensus per column = most frequent character. Tie-breaking options documented: IUPAC ambiguity
     code, a *specified residue order*, or '?'. No single universal rule.
4. **Web search — alphabetical tie-break** (WebSearch 2026-06-15). The phrasing "In the event of a tie,
   the residue letter occurring earlier in the alphabet was chosen" appears in a USPTO patent describing
   a consensus-motif method; current Geneious Prime instead uses ambiguity codes. So *alphabetical*
   tie-break is a real, documented convention but not universal.

### Formula check
Profile column maximum = most-frequent residue (Rosalind point 2; Wikipedia definition). Matches the
spec/description exactly. No log/normalisation/units involved (integer column counts).

### Edge-case semantics check
- Equal-length precondition: explicit in Rosalind ("collection of equal-length DNA strings"). ✓
- Empty / single / identical: derivable from the definition (single → its own bases; identical → that
  sequence). ✓
- Tie: source-permitted to be any tied symbol; library fixes alphabetical (A<C<G<T) for **determinism**.
  This is an *assumption*, correctly registered in the TestSpec (Assumption #1) and Evidence. It is also
  the de-facto behaviour of Biopython `Bio.motifs` `.consensus` (scan alphabet A→C→G→T, keep first max).

### Independent cross-check (numbers — hand-computed this session from the 7 strings)
| col | bases | counts | max → base |
|----|-------|--------|-----------|
| 0 | A G A A T A A | A=5 G=1 T=1 | A |
| 1 | T G T A T T T | T=5 G=1 A=1 | T |
| 2 | C G G G G G G | G=6 C=1 | G |
| 3 | C C G C G C G | C=4 G=3 | C |
| 4 | A A A A A C C | A=5 C=2 | A |
| 5 | G A T A A A A | A=5 G=1 T=1 | A |
| 6 | C C C C C T C | C=6 T=1 | C |
| 7 | T T T C T T T | T=6 C=1 | T |

Consensus = **ATGCAACT**, profile rows identical to Rosalind's published profile. ✓

### Findings / divergences
- Description is biologically/mathematically correct. The only non-source-mandated element is the
  **alphabetical tie-break**, which the description honestly labels an assumption (Rosalind allows any
  tied symbol). Hence **PASS-WITH-NOTES** rather than PASS. The note is documented, defensible
  (matches Biopython `.consensus`), and does not affect the rank-5 Rosalind dataset (no decisive ties).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs:392-437` plus the alphabet table at
line 371 (`ConsensusAlphabet = {'A','C','G','T'}`).

### Formula realised correctly? (evidence)
- Per column, a 4-element count array indexed A,C,G,T is built (lines 409-422); invalid characters
  throw `ArgumentException` (415-419).
- Maximum found by scanning indices 1..3 with strict `>` against `bestIndex` initialised to 0
  (426-431); because the scan visits A→C→G→T and only a strictly-greater count displaces the
  incumbent, the **alphabetically-earliest** maximum wins on a tie. Matches the validated description
  and the Biopython convention.
- Null → `ArgumentNullException` (394); empty → `""` (397); unequal length → `ArgumentException`
  (400-402). Uppercasing via `ToUpperInvariant` (396) gives case-insensitivity.

### Cross-verification table recomputed vs code
Ran the suite: M1 returns `ATGCAACT` (matches the hand-computed/Rosalind value). Tie case `AT`,`GT`
→ `AT` (col0 A vs G tie resolves to A). All consistent with the per-column table above.

### Variant/delegate consistency
Single public method; no `*Fast`/delegate variants. The sibling `GenerateConsensus` (IUPAC-degenerate)
is a different unit and is explicitly distinguished in the XML docs.

### Test quality audit (gate)
File `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_CreateConsensusFromAlignment_Tests.cs` (10 tests):
- **Sourced, exact values, not code-echoes:** M1 `ATGCAACT` traces to the Rosalind sample (re-fetched +
  hand-verified this session). M3 (`AT`,`GT`→`AT`) is **discriminating**: a reversed tie-break would
  return `GT`, so the test would fail against a wrong implementation — not a tautology.
- **No green-washing:** every assertion is exact `Is.EqualTo` or a specific exception type
  (`ArgumentNullException`/`ArgumentException`); no Greater/AtLeast/Contains/ranges/skips.
- **Coverage:** the single public method's every branch is exercised — happy path (M1,M2,M4,C2),
  case-insensitivity (M5), tie-break (M3), null (S1), empty (S2), unequal length (S3), invalid char
  (C1). All Stage-A edge cases covered. C2 (`AA`,`AA`,`CC`→`AA`) tests a strict majority distinct from
  a tie.
- **Honest green:** full unfiltered suite **Passed: 6570, Failed: 0** (10 s); `dotnet build` 0 errors.
  Pre-existing NUnit2007 warnings live only in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`,
  not in this unit's files; this session changed no code.

**Test-quality gate: PASS.**

### Findings / defects
None. Code faithfully realises the validated description; tests lock the sourced values and cover all
branches.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** (alphabetical tie-break is an honestly-documented assumption beyond the
  Rosalind "any tied symbol" allowance; consistent with Biopython `.consensus`).
- **Stage B: PASS.**
- **End-state: CLEAN** — no defect found; no code/test changes required.
