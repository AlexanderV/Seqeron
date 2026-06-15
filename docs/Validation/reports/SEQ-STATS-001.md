# Validation Report: SEQ-STATS-001 — Sequence Composition Statistics

- **Validated:** 2026-06-15   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateNucleotideComposition(string)`; delegates `SummarizeNucleotideSequence(string?)`, `CalculateAminoAcidComposition(string)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (with extracted numbers)

1. **Biopython `Bio.SeqUtils` source** — https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py (fetched 2026-06-15).
   - `gc_fraction`: numerator counts C, G, S (case-insensitive); default `ambiguous="remove"` denominator counts only C, G, S, A, T, W, U. **Returns 0 for an empty sequence.**
   - `GC_skew`: formula `(g - c) / (g + c)`, counts G and C upper+lowercase; **returns 0.0 when g+c = 0** (catches `ZeroDivisionError`); default `window=100`.
2. **Wikipedia "GC skew"** — https://en.wikipedia.org/wiki/GC_skew (fetched 2026-06-15).
   - Modern GC skew = **(G − C)/(G + C)**; AT skew = **(A − T)/(A + T)**.
   - Lobry's original 1996 convention was (C − G)/(C + G); modern implementations flip it.
   - Under the modern flipped definition: **positive GC skew = G-rich, negative = C-rich.**
3. **Lobry (1996)** Mol Biol Evol 13(5):660–665 (DOI 10.1093/oxfordjournals.molbev.a025626) — primary source for strand compositional asymmetry; cited via Wikipedia.

### Formula check
- GC content = (G+C)/(A+T+G+C+U), float in [0,1] — matches Biopython `gc_fraction` over the standard alphabet (numerator G+C, denominator the unambiguous bases). ✓
- GC skew = (G−C)/(G+C), 0 when G+C=0 — matches Biopython/Wikipedia exactly, including the modern sign convention. ✓
- AT skew = (A−T)/(A+T), 0 when A+T=0 — matches Wikipedia. ✓
- Sign interpretation (positive GC skew = G-rich) matches Wikipedia. ✓

### Edge-case semantics
- Empty/null → all-zero composition, GC content 0: sourced to Biopython "returns zero for an empty sequence". ✓
- Zero-denominator skew → 0: sourced to Biopython `ZeroDivisionError` handling. ✓
- Case-insensitive: Biopython counts lowercase; impl upper-cases first. ✓

### Independent cross-check (hand computation, all from the sourced formulas)
| Input | A T G C U | GC content | AT content | GC skew | AT skew |
|-------|-----------|-----------|-----------|---------|---------|
| ATGC | 1 1 1 1 0 | 2/4 = 0.5 | 2/4 = 0.5 | 0/2 = 0 | 0/2 = 0 |
| GGGC | 0 0 3 1 0 | 4/4 = 1.0 | 0 | 2/4 = 0.5 | 0 (A+T=0) |
| GCCC | 0 0 1 3 0 | 4/4 = 1.0 | 0 | −2/4 = −0.5 | 0 (A+T=0) |
| AAAT | 3 1 0 0 0 | 0 | 1.0 | 0 (G+C=0) | 2/4 = 0.5 |
| AAUUGGCC | 2 0 2 2 2 | 4/8 = 0.5 | (2+0+2)/8 = 0.5 | 0/4 = 0 | (2−0)/2 = 1.0 |

### Findings / divergences (Stage A → PASS-WITH-NOTES)
- **N1 (documented divergence):** Degenerate IUPAC codes (S, W, R, Y, …) are counted as `Other` and excluded from GC/AT totals, whereas Biopython's `gc_fraction` partially counts S toward GC and W toward the length. For the unit's declared scope ({A,T,G,C,U}) the two agree exactly. Explicitly documented in the algorithm doc §5.3 and the assumption register. Acceptable.
- **N2 (documented convention):** `AtContent` uses (A+T+U)/total — it counts U with A/T for RNA. No external source defines an "AT content" for RNA, so this is a defensible, documented convention (algorithm doc §5.2 / §7). `AtSkew` correctly uses the DNA-specific (A−T)/(A+T) without U, matching Lobry/Wikipedia. No correctness issue.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs`
- `CalculateNucleotideComposition` lines 49–94.
- `CalculateAminoAcidComposition` lines 99–139.
- `SummarizeNucleotideSequence` lines 990–1019.

### Formula realised correctly?
- Null/empty → all-zero record (lines 51–54). ✓
- Single-pass switch counting A/T/G/C/U/N, everything else → Other (lines 56–70). ✓
- `total = a+t+g+c+u`; `gcContent = gc/total` guarded by `total>0` (lines 72–76). ✓
- `gcSkew = (g-c)/(g+c)` guarded by `(g+c)>0`, else 0 (line 78). ✓ Matches Biopython/Wikipedia incl. sign.
- `atSkew = (a-t)/(a+t)` guarded by `(a+t)>0`, else 0 (line 79). ✓
- `Length = sequence.Length` (line 82) — includes N/Other, so INV-03 partition holds. ✓
- `SummarizeNucleotideSequence` delegates GcContent/counts to the canonical method (line 996, 1015). ✓
- `CalculateAminoAcidComposition` counts letters case-insensitively; Length = sum of counts (lines 106–115). ✓

### Cross-verification table recomputed vs code
All five sequences above were checked against the code logic and the canonical test assertions; every value matches the externally-sourced hand computation (GC content 0.5/1.0, GC skew 0.5/−0.5/0, AT skew 0.5/1.0/0).

### Variant/delegate consistency
`SummarizeNucleotideSequence` returns the same GcContent and per-base counts as the canonical method (verified by code path and test C1). `CalculateAminoAcidComposition` residue counts verified exact (test C2). ✓

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** Canonical file `SequenceStatistics_CalculateNucleotideComposition_Tests.cs` asserts exact values traced to Biopython/Wikipedia (GC content 0.5/1.0, GC skew ±0.5 exact, AT skew 0.5 exact, empty→0, zero-denominator→0). These would fail a deliberately-wrong implementation. ✓
- **No green-washing:** No weakened assertions in the canonical file — all use `Is.EqualTo(...)` within 1e-10, not Greater/Less/ranges. ✓
- **Coverage:** All Stage-A branches exercised — counts partition (M1, S5), GC content (M2/M3/S4), AT content incl. U (added M4b), GC skew ±/zero (M5/M6/S2), AT skew/zero (M7/S3), AT skew RNA U-exclusion (added M7b), empty (M8), null (M9), case-insensitivity (S1), N/Other (S4), delegates (C1/C2/C2b). ✓
- **Legacy weak tests:** `SequenceStatisticsTests.cs` retains pre-existing permissive SEQ-STATS-001 assertions (`GcSkew Is.GreaterThan(0)`, Summary `GreaterThan(0)`). These are true-but-weak and explicitly superseded by the exact canonical file (TestSpec §5.2/§5.3); they are not green-wash substitutes since the exact assertions exist independently. Left in place (they also support neighbouring concerns) — no defect.

### Tests added this session (lock the documented RNA conventions)
- **M4b** `CalculateNucleotideComposition_RnaSequence_AtContentIncludesUracil` — `AAUUGGCC` AtContent = 0.5 exact (locks the (A+T+U)/total branch, sourced from the Evidence worked example).
- **M7b** `CalculateNucleotideComposition_RnaSequence_AtSkewExcludesUracil` — `AAUUGGCC` AtSkew = 1.0 exact (locks the DNA-specific (A−T)/(A+T) skew formula, sourced from Wikipedia + Evidence table).

### Build & test result
- `dotnet build` — 0 errors (4 pre-existing warnings, none in SEQ-STATS-001 files).
- Full unfiltered `dotnet test` — **Failed: 0, Passed: 6512, Skipped: 0** (confirmed on two consecutive runs). One transient single-test flake observed on an intermediate run (not in any SEQ-STATS-001 test; isolated class run was 19/19, trx confirmed no SEQ-STATS failures); both clean full runs ended Failed: 0.

### Findings / defects
None. No code change required; two exact-value tests added to lock the documented RNA conventions.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (two documented, sourced-or-defensible divergences: IUPAC handling N1, RNA AtContent convention N2 — both correct within scope).
- **Stage B:** PASS — implementation faithfully realises the validated formulas; tests assert exact externally-sourced values and cover all branches.
- **Test-quality gate:** PASS.
- **End-state:** CLEAN — fully functional; no defect found; coverage strengthened.
- No defect logged in FINDINGS_REGISTER (N1/N2 are documented intentional conventions, not defects).
