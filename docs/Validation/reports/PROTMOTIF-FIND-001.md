# Validation Report: PROTMOTIF-FIND-001 — Protein Motif Search

- **Validated:** 2026-06-12   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.FindMotifByPattern(sequence, regexPattern, motifName, patternId)`, `ProteinMotifFinder.FindCommonMotifs(proteinSequence)`, `ProteinMotifFinder.CommonMotifs` (dictionary)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

Source file: `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs`
Test file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_MotifSearch_Tests.cs` (34 tests)

---

## Stage A — Description

### Sources opened & what they confirm

1. **PROSITE User Manual** (https://prosite.expasy.org/prosuser.html) — confirmed the full pattern syntax:
   - `-` separates pattern elements; standard IUPAC one-letter codes.
   - `[ABC]` = any one of A/B/C (alternatives).
   - `{ABC}` = any amino acid EXCEPT A/B/C (exclusion).
   - `x` = any amino acid; `x(n)` = `x` repeated n times; `x(m,n)` = between m and n.
   - `<` = N-terminus anchor, `>` = C-terminus anchor; a period ends the pattern.
   - IUPAC codes are uppercase; matching is treated case-insensitively in practice.
2. **Wikipedia "Sequence motif"** (https://en.wikipedia.org/wiki/Sequence_motif) — confirmed the N-glycosylation motif: "*Asn, followed by anything but Pro, followed by either Ser or Thr, followed by anything but Pro*" written `N{P}[ST]{P}`; confirmed `x`, `{...}`, `[XY]`, and `-` semantics.
3. **PROSITE PS00001** (https://prosite.expasy.org/PS00001) — official pattern `N-{P}-[ST]-{P}.` (ASN_GLYCOSYLATION). Matches implementation `N[^P][ST][^P]`.
4. **PROSITE PS00007** (https://prosite.expasy.org/PS00007) — official pattern `[RK]-x(2)-[DE]-x(3)-Y` (fixed repeats). Matches implementation regex `[RK].{2}[DE].{3}Y` (the documented PS00007 fix is in place).
5. **PROSITE PS00018** (https://prosite.expasy.org/PS00018) — official pattern `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]`. Matches implementation regex exactly, including `[^W]` at position 2 and trailing `[LIVMFYW]` (the documented PS00018 fix is in place).

### Semantics covered by this unit

This unit covers **regex-based motif search** (`FindMotifByPattern`) and a **curated `CommonMotifs` dictionary** of pre-converted PROSITE/literature patterns scanned by `FindCommonMotifs`. The PROSITE-string→regex parser (`ConvertPrositeToRegex`) is the subject of the sibling unit **PROTMOTIF-PROSITE-001** and is out of scope here.

- **Position convention:** 0-based; `Start`/`End` are inclusive indices into the (uppercased) sequence. `Sequence == input[Start..End+1]`.
- **All occurrences incl. overlapping:** discovered via lookahead wrapper `(?=(pattern))`, consistent with ScanProsite (De Castro et al. 2006).
- **Case:** case-insensitive (`RegexOptions.IgnoreCase` + uppercasing).
- **Edge cases:** empty/null sequence → empty; empty/invalid regex → empty (no throw); motif = whole sequence → single match.

### Independent cross-check (hand computation)

- **RGD** in `MRGDKLARGDPMRGD` (regex `RGD`): occurrences at 0-based indices 1, 7, 12 → 3 matches, each `Sequence="RGD"`. Matches spec M1/M2.
- **N-glycosylation** in `AAAANFTAAAA` (regex `N[^P][ST][^P]`): N@4, F(≠P), T(∈ST), A(≠P) → match `NFTA` at [4,7]; exactly 1. Matches spec M3.
- **Proline exclusion**: `AAAANPSAAAAANPTAAA` — N followed by P fails `[^P]` → 0 matches. Matches spec M4.
- **Overlapping**: `P..P` on `PPPPP` → matches at 0 and 1 (2 matches), confirming lookahead overlap discovery.
- **Score (IC)**: RGD = 3 fixed positions × log₂(20) ≈ 12.97 bits; PKC `[ST].[RK]` = 2×log₂(10) ≈ 6.64 bits. Both consistent with Schneider & Stephens (1990) information content.

### Findings / divergences

None. All cited PROSITE patterns (PS00001, PS00005, PS00006, PS00004, PS00007, PS00008, PS00009, PS00016, PS00017, PS00018, PS00028, PS00029) and the 5 non-PROSITE literature patterns (NLS1, NES1, SIM1, WW1, SH3_1) match the authoritative definitions byte-for-byte.

---

## Stage B — Implementation

### Code path reviewed

- `FindMotifByPattern` (lines 161–198): uppercases input, wraps the user regex in `(?=(...))` for overlapping discovery, returns `MotifMatch(Start=match.Index, End=Index+len-1, Sequence=captured.Value, …)`. Try/catch around `new Regex(...)` yields empty on malformed pattern.
- `FindCommonMotifs` (lines 140–154): null/empty guard → empty; iterates `CommonMotifs.Values`, delegating to `FindMotifByPattern`.
- `CommonMotifs` dictionary (lines 69–131): each entry carries both the PROSITE-style string and the equivalent regex, with inline source URLs.
- Scoring: `CalculateMotifScore` (IC = Σ log₂(20/allowed)) and `CalculateEValue` (E = (N−L+1)·2^(−IC)), with `ParseRegexAllowedCounts` reading character classes and `{n}` quantifiers.

### Formula realised correctly?

Yes. Lookahead overlapping discovery, 0-based inclusive positions, case-insensitive matching, and IC/E-value scoring all match the validated description. The substring invariant `Sequence == upper[Start..End+1]` holds by construction (`captured.Value` is the captured group, `End=Index+len−1`).

### Cross-verification table recomputed vs code (via tests)

| Case | Input | Expected | Result |
|------|-------|----------|--------|
| M1/M2 RGD | `MRGDKLARGDPMRGD` | 3 at 1,7,12 | PASS |
| M3 N-glyc | `AAAANFTAAAA` | 1 × `NFTA` [4,7] | PASS |
| M4 Pro-excl | `AAAANPSAAAAANPTAAA` | 0 | PASS |
| M5 PKC | `AAAAASARKAAA` | 1 × `SAR` [5,7] | PASS |
| M6 CK2 | `AAAASAAEASDEDAAA` | 2 (`SAAE`[4,7],`SDED`[9,12]) | PASS |
| M7 P-loop | `AAAAAGXXXXGKSAAAA` | 1 × `GXXXXGKS` [5,12] | PASS |
| S4 overlap | `PPPPP` / `P..P` | 2 at 0,1 | PASS |
| M17 IC | RGD | ≈12.97 bits | PASS |

### Variant/delegate consistency

`FindMotifByProsite` delegates through `ConvertPrositeToRegex` → `FindMotifByPattern` (PROSITE-syntax conversion validated under PROTMOTIF-PROSITE-001). `FindCommonMotifs` is a thin loop over `FindMotifByPattern` — consistent.

### Test quality audit

34 canonical tests assert exact counts, positions, and sequences (not "no throw"/tautologies), and lock the official PROSITE/literature pattern strings (M12–M14, S5). Edge cases (empty/null/invalid-regex/empty-pattern/overlapping/case) all covered with exact expectations.

### Findings / defects

- **Cosmetic (fixed):** test `FindMotifByPattern_EValue_UsesProperProbability` cited "Altschul et al. 1990" in its assertion message, but the Evidence doc explicitly removed the BLAST citation in favor of Schneider & Stephens (1990) direct combinatorial probability. Updated the comment to match the authoritative source. No behavioral change; formula and assertion were already correct.

---

## Verdict & follow-ups

- **Stage A: PASS** — every pattern and semantic convention independently confirmed against PROSITE (ExPASy) and Wikipedia.
- **Stage B: PASS** — code faithfully realises the validated description; all worked examples recomputed; tests are strong.
- **End state: CLEAN.** One cosmetic test-comment fix applied. Build green; `ProteinMotifFinder_MotifSearch_Tests` 34/34 pass; full suite 4486/4486 pass.

**Files changed:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_MotifSearch_Tests.cs` (one assertion-message correction).
