# Validation Report: PROTMOTIF-FIND-001 — Protein Motif Search

- **Validated:** 2026-06-24   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.FindMotifByPattern(sequence, regexPattern, motifName, patternId)`, `ProteinMotifFinder.FindCommonMotifs(proteinSequence)`, `ProteinMotifFinder.CommonMotifs` (dictionary)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End state:** CLEAN

Source file: `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs`
Canonical test file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_MotifSearch_Tests.cs` (34 tests)

---

## Stage A — Description

### Motif model

The unit covers **regex-based protein motif search** plus a curated `CommonMotifs` dictionary of
pre-converted PROSITE / literature patterns:

- Motifs are specified as **.NET regexes** that are the regex translation of PROSITE PA-line
  patterns (the PROSITE→regex string converter `ConvertPrositeToRegex` is the sibling unit
  PROTMOTIF-PROSITE-001 and is out of scope here).
- **Alphabet:** the 20 standard one-letter amino-acid codes; non-standard residues simply do not
  match the constrained classes.
- **Position convention:** 0-based; `Start`/`End` are inclusive indices into the uppercased
  sequence, so `Sequence == input[Start..End+1]`.
- **All occurrences incl. overlapping** are discovered via the lookahead wrapper `(?=(pattern))`,
  consistent with ScanProsite (De Castro et al. 2006).
- **Case:** case-insensitive (`RegexOptions.IgnoreCase` + uppercasing) per PROSITE convention.

### Sources opened & what they confirm (this session, fetched live)

Every PROSITE PA-line was fetched verbatim from prosite.expasy.org and compared to both the
`CommonMotifs.Pattern` string and its `RegexPattern`:

| Accession | Official PA line (fetched) | Impl. regex | Match |
|-----------|----------------------------|-------------|-------|
| PS00001 | `N-{P}-[ST]-{P}.` | `N[^P][ST][^P]` | ✓ |
| PS00004 | `[RK](2)-x-[ST]` | `[RK]{2}.[ST]` | ✓ |
| PS00005 | `[ST]-x-[RK]` | `[ST].[RK]` | ✓ |
| PS00006 | `[ST]-x(2)-[DE]` | `[ST].{2}[DE]` | ✓ |
| PS00007 | `[RK]-x(2)-[DE]-x(3)-Y` | `[RK].{2}[DE].{3}Y` | ✓ |
| PS00008 | `G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}` | `G[^EDRKHPFYW].{2}[STAGCN][^P]` | ✓ |
| PS00009 | `x-G-[RK]-[RK]` | `.G[RK][RK]` | ✓ |
| PS00016 | `R-G-D` | `RGD` | ✓ |
| PS00017 | `[AG]-x(4)-G-K-[ST]` | `[AG].{4}GK[ST]` | ✓ |
| PS00018 | `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]` | exact (see code l.109-110) | ✓ |
| PS00028 | `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` | `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H` | ✓ |
| PS00029 | `L-x(6)-L-x(6)-L-x(6)-L` | `L.{6}L.{6}L.{6}L` | ✓ |

The 5 non-PROSITE literature patterns (NLS1, NES1, SIM1, WW1, SH3_1) carry inline citations
(Chelsky 1989 / Dingwall & Laskey 1991; la Cour 2004; Hecker 2006; Chen & Sudol 1995; Mayer 2001)
matching the Evidence doc; all 12 PROSITE patterns above are the priority external cross-check.

### Edge-case semantics

Empty/null sequence → empty; empty/invalid/pathological regex → empty (no throw, bounded 2 s
match timeout); zero-width captures are skipped so the `0 ≤ Start ≤ End ≤ n−1` invariant holds;
overlapping occurrences are all returned. All sourced (PROSITE docs / robustness).

### Independent cross-check (hand computation, this session)

`MRGDKLARGDPMRGD` indexed: M0 R1 G2 D3 K4 L5 A6 R7 G8 D9 P10 M11 R12 G13 D14
→ RGD at **0-based 1, 7, 12** (`Sequence="RGD"` each). The code and canonical tests use these
exact values.

Other hand checks confirmed against code/tests: N-glyc `NFTA`@[4,7] in `AAAANFTAAAA`;
0 matches for proline exclusion `AAAANPSAAAAANPTAAA`; overlapping `P..P` on `PPPPP` → matches at
0 and 1; IC scores RGD = 3·log₂20 ≈ 12.97 bits, PKC `[ST].[RK]` = 2·log₂10 ≈ 6.64 bits
(Schneider & Stephens 1990).

### Findings / divergences

- **TestSpec doc typo (FIXED):** `tests/TestSpecs/PROTMOTIF-FIND-001.md` rows M1/M2 listed the
  middle RGD at position **6** (“positions 1, 6, 12” / “Start=6,End=8”). The correct 0-based index
  for the second RGD in `MRGDKLARGDPMRGD` is **7** (R is at index 7, not 6). The code, the
  canonical tests, and the prior validation report all use 7. Corrected the TestSpec to
  “1, 7, 12” / “Start=7,End=9”. Documentation-only change; no behavioural impact. This is the sole
  reason Stage A is PASS-WITH-NOTES rather than PASS.

---

## Stage B — Implementation

### Code path reviewed

`FindMotifByPattern` (l.178-240): lookahead-wrapped regex with `IgnoreCase` + 2 s match timeout;
swallows compile errors and `RegexMatchTimeoutException` → empty; skips zero-width captures;
emits `(Start=match.Index, End=Index+len-1, Sequence=captured.Value, …)`. `FindCommonMotifs`
(l.150-164) uppercases once and scans every `CommonMotifs` entry. Scoring `CalculateMotifScore`
(l.1208) and `CalculateEValue` (l.1228) implement IC = Σ log₂(20/allowed) and E = (N−L+1)·2^−IC.

### Formula realised correctly?

Yes. Position convention (0-based inclusive, substring invariant), overlapping discovery,
case-insensitivity, IC/E-value formulas, and all 12 PROSITE regex strings match the validated
description and the live-fetched PROSITE definitions.

### Cross-verification vs code

The 34 canonical tests encode the hand-checked values (RGD@1,7,12; NFTA@[4,7]; proline exclusion
→ 0; CK2 two sites; P-loop@[5,12]; overlapping P..P → 2; IC = 12.97 / 6.64 bits; E-value
= 7·2^−IC) and all pass. PS00007/PS00018 fixes (exact repeats, `[^W]` at pos 2, trailing
`[LIVMFYW]`) are present and asserted byte-for-byte.

### Variant/delegate consistency

`FindCommonMotifs` delegates to `FindMotifByPattern`; `FindMotifByProsite` delegates via
`ConvertPrositeToRegex` (separate unit). `FindDomains` reuses `FindMotifByPattern`. Consistent.

### Test quality audit

Strong: exact counts, exact Start/End/Sequence, exact pattern+regex strings vs official PROSITE,
exact IC/E-value with `.Within`. Edge cases (empty/null/invalid/empty-pattern/overlapping) covered
with deterministic asserts. No tautologies.

### Findings / defects

None in code.

---

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** — description and all 12 PROSITE patterns independently re-confirmed
  verbatim against prosite.expasy.org; one TestSpec position typo (6→7) found and fixed.
- **Stage B: PASS** — code faithfully realises the validated model; 34/34 canonical tests pass.
- **End state: CLEAN.** Code unchanged (correct already); only the TestSpec documentation was
  corrected.
