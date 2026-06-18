# Validation Report: SEQ-REPLICATION-001 — Replication Origin Prediction (cumulative GC-skew minimum)

- **Validated:** 2026-06-16   **Area:** Composition
- **Canonical method(s):** `GcSkewCalculator.PredictReplicationOrigin(DnaSequence)`, `GcSkewCalculator.PredictReplicationOrigin(string)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session (retrieved, not trusted by label)

| Source | URL | What it confirms |
|--------|-----|------------------|
| Rosalind, Minimum Skew Problem (BA1F) | https://rosalind.info/problems/ba1f/ | skew = (#G − #C); "All integer(s) i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to \|Genome\|)"; sample input = the 100-nt string used in the tests (verbatim); sample output = `53 97`. |
| Wikipedia, GC skew | https://en.wikipedia.org/wiki/GC_skew | "the maximum value of the cumulative skew corresponds to the terminal, and the minimum value corresponds to the origin of replication." Formula GC skew = (G − C)/(G + C). Cites Lobry 1996 and Grigoriev 1998. |

(Grigoriev 1998 / Lobry 1996 are the cited primaries behind the min=origin/max=terminus convention; the convention itself was confirmed verbatim from Wikipedia's quote of them, and the numeric worked example from Rosalind BA1F.)

### Formula check
- Per-nucleotide cumulative skew: Skew_0 = 0; G → +1, C → −1, A/T → 0; prefix indices i ∈ [0, n]. Matches Rosalind BA1F definition exactly.
- Origin = global minimum of the cumulative diagram; terminus = global maximum. Matches Wikipedia's quoted sentence (citing Grigoriev/Lobry).

### Edge-case semantics
- **Ties:** BA1F returns ALL minimizers (`53 97`); the API returns a single position, with a documented deterministic tie-break of "first (smallest) minimizing index". Defensible and documented (TestSpec §7, Evidence corner-case 1).
- **Flat diagram (no G/C asymmetry):** amplitude 0, origin = terminus = 0; documented as "not resolvable" (Evidence corner-case 3). Sourced behaviour, not implementation-defined.
- **IsSignificant:** ASSUMPTION 1 — no authoritative numeric cutoff exists; redefined as `max > min` (non-zero amplitude). The previous invented `0.01 × count` threshold was removed. This is the weakest non-invented predicate and is honestly flagged as an assumption (INV-5). Accepted as PASS (not a divergence from any external source, because no source defines significance).

### Independent cross-check (numbers re-derived this session, Python)
| Sequence | Min value @ first pos | Max value @ first pos | Source it matches |
|----------|----------------------|----------------------|-------------------|
| BA1F (100 nt) | −4 @ **53** (also 97) | +2 @ 16 | Rosalind sample output `53 97`; min value −4 |
| `CCGGGG` | −2 @ 2 | +2 @ 6 | BA1F definition (re-derived) |
| `GGGCCC` | 0 @ 0 | +3 @ 3 | max = terminus |
| `CCGGCC` | −2 @ 2 (also 6) | 0 @ 0 | tie-break first index |
| `GAATTG` | 0 @ 0 | +2 @ 6 | A/T ignored |
| `AAAATTTT`/`AATT` | 0 @ 0 | 0 @ 0 | flat diagram |
| `G` | 0 @ 0 | +1 @ 1 | single base |

The published Rosalind output `53 97` (min −4) was reproduced exactly from the verbatim sample string — confirming both the formula and the prefix-indexing convention against an external worked example.

### Invariants
INV-1..INV-4 and INV-6 are genuine mathematical properties of the cumulative-skew diagram (confirmed by the re-derivation). INV-3 (OriginSkew ≤ 0 ≤ TerminusSkew) holds because Skew_0 = 0 is always in range. INV-5 is the flagged assumption above.

**Stage A findings:** none. Description is biologically and mathematically correct and matches the retrieved sources verbatim.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs:236-286` (`PredictReplicationOrigin` overloads + `PredictReplicationOriginCore`).

### Formula realised correctly?
Yes. `PredictReplicationOriginCore` accumulates G:+1 / C:−1 (constants at :215-216), leaves A/T/other unchanged, starts cumulative at 0 with `minPos = maxPos = 0` (Skew_0), and uses `prefixIndex = i + 1`. Min/max are updated with strict `<` / `>`, which preserves the first-occurrence (smallest-index) tie-break. This is the exact per-nucleotide cumulative-skew model validated in Stage A — not a windowed/approximate model. (Note: a separate legacy windowed `PredictReplicationOrigin(windowSize)` overload and its invented `0.01×count` threshold were removed in the authoring session; confirmed absent.)

Overloads: `(DnaSequence)` null-guards then delegates to core; `(string)` returns the zero prediction for null/empty and otherwise uppercases (case-insensitive) before delegating. Consistent.

### Cross-verification table recomputed vs code
The 15 tests in the canonical file assert exactly the externally-sourced numbers in the table above; the full suite passes (below). M1 / string-M1 assert `53` and `−4` from the published Rosalind output — values a windowed or off-by-one implementation could not produce, so these are genuine sourced checks, not code echoes.

### Variant/delegate consistency
`(string)` overload reproduces the BA1F result identically to `(DnaSequence)` (tested, S3 + string-M1). Confirmed.

### Numerical robustness
Integer cumulative counter (no float drift); widened to double only when stored in the record. No div-by-zero (no division). O(n), single pass.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1/string-M1 lock `53`/`−4` from Rosalind's published output; M2–M6, S-, C-cases lock hand-re-derived diagram values. ✅
- **No green-washing:** exact `Is.EqualTo(...)` with tight `Within(1e-10)` tolerances; INV-3 uses `LessThanOrEqualTo/GreaterThanOrEqualTo` legitimately (it asserts a sign invariant, not a known exact value); no skipped/ignored/weakened tests. ✅
- **Coverage:** every public overload exercised; Stage-A branches covered — per-nt increments (M2), max=terminus (M3), first-index tie-break (M4), A/T ignored (M5), sign invariant (M6), flat diagram + IsSignificant=false (S1), IsSignificant=true (S2), bounds (S4), case-insensitive (S3), null DnaSequence throws (C1), null/empty string → zero prediction (C2), single base (C3). ✅
- **Honest green:** full unfiltered suite **Passed: 6607, Failed: 0, Skipped: 0**; `dotnet build` 0 errors. ✅

Minor coverage observation (not a defect): the exact mirror of `GGGCCC` — an all-negative diagram where the terminus stays pinned at position 0 — is not asserted as its own case, but the symmetric origin-side behaviour is covered by `GAATTG`/`GGGCCC`, and the implementation is symmetric. The spec §5.4 reports "16 tests" while the file contains 15; this is a documentation count slip, not a missing required case (every M/S/C row maps to a present test).

**Stage B findings:** none requiring a code or test change. No defect.

## Verdict & follow-ups
- **Stage A: PASS.** Description matches Rosalind BA1F (formula, indexing, sample `53 97`/−4) and Wikipedia's min=origin/max=terminus quote (citing Grigoriev 1998 / Lobry 1996), all retrieved this session.
- **Stage B: PASS.** Code faithfully realises the validated per-nucleotide cumulative-skew model; tests assert externally-sourced values and cover all branches; full suite green.
- **End-state: ✅ CLEAN** — no defect found; algorithm fully functional.
- Non-blocking note: TestSpec §5.4 says "16 tests" vs 15 present (count slip only).
