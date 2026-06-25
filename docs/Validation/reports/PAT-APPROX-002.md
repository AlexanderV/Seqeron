# Validation Report: PAT-APPROX-002 — Approximate Matching (Edit / Levenshtein Distance)

- **Validated:** 2026-06-24   **Area:** Pattern Matching
- **Canonical method(s):** `ApproximateMatcher.EditDistance(string, string)`, `ApproximateMatcher.FindWithEdits(string, string, int)`, `FindWithEdits(DnaSequence, string, int)` (wrapper)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Distance model (and contrast with PAT-APPROX-001)

PAT-APPROX-002 is **edit distance / Levenshtein** with **unit costs** for insertion,
deletion, and substitution (match = 0). Transpositions are NOT a single operation
(that would be Damerau–Levenshtein, which is not claimed). This contrasts with
PAT-APPROX-001 (**Hamming distance** = substitutions only, equal-length strings,
no indels). `EditDistance` is case-sensitive (characters are distinct symbols per the
standard definition); the `FindWithEdits` search uppercases text+pattern by design for
DNA matching. Convention: returns a non-negative integer count of edits; `FindWithEdits`
positions are 0-based start offsets in the sequence.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Levenshtein distance"**: the recurrence below with base cases
  `lev(i,0)=i`, `lev(0,j)=j`; ins/del/sub each cost 1, match 0; worked examples
  `kitten→sitting = 3` and `flaw→lawn = 2`; transpositions not free in standard Levenshtein.
- **Rosetta Code — "Levenshtein distance"**: canonical vectors (`rosettacode→raisethysword=8`,
  `saturday→sunday=3`, `stop→tops=2`, `sleep→fleeting=5`). Reproduced independently below
  (a from-scratch DP) rather than trusting the citation label.
- **Navarro (2001), "A Guided Tour to Approximate String Matching"**: edit-distance foundation,
  consistent with unit-cost Levenshtein.

### Formula check
```
lev(i,0) = i ;  lev(0,j) = j
lev(i,j) = min( lev(i-1,j)+1,                          // deletion
                lev(i,j-1)+1,                          // insertion
                lev(i-1,j-1) + (a[i-1]==b[j-1]?0:1) )  // substitution / match
```
Substitution cost = 1 (Levenshtein), not 2 (LCS-style). Confirmed verbatim against Wikipedia.

### Edge-case semantics check
empty vs empty = 0 ✔; empty vs length-n = n ✔; identical = 0 ✔; transposition NOT free
(`stop→tops = 2`) ✔; case-sensitive (`A` vs `a` = 1) ✔. Null → ArgumentNullException;
negative maxEdits → ArgumentOutOfRangeException.

### Independent cross-check (numbers)
A from-scratch full-matrix Wagner–Fischer reference (unit costs) reproduced **every** spec
vector exactly:

| s1 | s2 | reference | spec |
|----|----|-----------|------|
| kitten | sitting | 3 | 3 |
| rosettacode | raisethysword | 8 | 8 |
| saturday | sunday | 3 | 3 |
| flaw | lawn | 2 | 2 |
| stop | tops | 2 | 2 |
| sleep | fleeting | 5 | 5 |
| "" | abc | 3 | 3 |
| abc | "" | 3 | 3 |
| "" | "" | 0 | 0 |
| A | a | 1 | 1 |
| Kitten | kitten | 1 | 1 |
| ABC | abc | 3 | 3 |

Invariants (symmetry, identity, empty-string = length, triangle inequality,
`|len a − len b| ≤ d ≤ max(len)`) are genuine metric properties and hold.

### Findings / divergences
None.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs`
- `EditDistance` — lines 195–236
- `FindWithEdits(string,…)` — lines 118–155; `FindWithEdits(DnaSequence,…)` wrapper — lines 160–164

### Formula realised correctly? (evidence)
`EditDistance` is a two-row (rolling) Wagner–Fischer DP:
- base row `prev[j]=j` (212–213) and `curr[0]=i` (217) → `lev(0,j)=j`, `lev(i,0)=i`;
- `cost = s1[i-1]==s2[j-1] ? 0 : 1` (221) → match 0 / mismatch 1;
- `min(prev[j]+1 /*del*/, curr[j-1]+1 /*ins*/, prev[j-1]+cost /*sub*/)` (222–228) →
  exactly the validated recurrence;
- empty fast paths `if(m==0) return n; if(n==0) return m;` (204–205) match base cases;
- `EditDistance` does NOT uppercase → correctly case-sensitive per the standard definition.
- Null → ArgumentNullException (197–198); two-row buffer cannot overflow; `int` cannot
  overflow (distance ≤ max length).

### Cross-verification table recomputed vs code
The 24-test suite asserts the exact sourced vectors and all five invariants via
`Is.EqualTo(<constant>)`; all match the Stage-A reference numbers.

### Variant/delegate consistency
`FindWithEdits(DnaSequence,…)` forwards `sequence.Sequence` to the string overload
(verified by C01). `FindWithEdits` uses `EditDistance` as its variable-length window
scorer (139) with a `maxEdits ≥ 0` guard (124–125); empty inputs yield empty (121–122).

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcher_EditDistance_Tests.cs` — 24 tests
(M01–M13 minus none, S01–S10, C01; M01–M13 map to the spec). Builds with **0 warnings**;
all assertions check exact sourced values, deterministic, cover every Stage-A edge case.
Test run: **Failed: 0, Passed: 24**.

### Findings / defects
None.

## Verdict & follow-ups
Stage A PASS, Stage B PASS. State: CLEAN. No code changes. No follow-ups.
