# Test Specification: GENOMIC-TANDEM-001

**Test Unit ID:** GENOMIC-TANDEM-001
**Area:** Analysis
**Algorithm:** Exact tandem repeats (`GenomicAnalyzer.FindTandemRepeats`; MCP `find_tandem_repeats`)
**Status:** ☑ Complete (re-validated and re-implemented 2026-11-18, B09 F35–F38)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-18

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Benson, G. (1999) *Tandem Repeats Finder*, *Nucleic Acids Research* 27(2):573–580 | 1 | https://doi.org/10.1093/nar/27.2.573 | 2026-06-14 (definition re-confirmed by search snippets 2026-11-18; publisher page returns 403) |
| 2 | Kolpakov & Kucherov (1999) *Finding maximal repetitions in a word in linear time*, FOCS, 596–604 | 1 | https://dl.acm.org/doi/10.5555/795665.796470 | 2026-11-18 (search snippets) |
| 3 | Bannai, I, Inenaga, Nakashima, Takeda, Tsuruta (2017) *The "Runs" Theorem*, SIAM J. Comput. 46(5) | 1 | https://arxiv.org/abs/1406.0263 | 2026-11-18 (full text opened: Definition 1, Lemmas 3/6/7/8, §3, §4) |
| 4 | Thiel et al. (2003) MISA, `misa.pl` ("reject false type motifs") | 2 | https://raw.githubusercontent.com/cfljam/SSR_marker_design/master/misa.pl | 2026-11-18 (opened; line 110) |
| 5 | Wikipedia "Tandem repeat" (cited primaries) | 4 | https://en.wikipedia.org/wiki/Tandem_repeat | 2026-06-14 |

### 1.2 Key Evidence Points

1. A tandem repeat is "two or more contiguous … copies of a pattern of nucleotides" (k ≥ 2; copies directly adjacent) — Benson 1999 [1].
2. Worked example: `ATTCG ATTCG ATTCG` — unit `ATTCG` (period 5) repeated 3 times, contiguous — Wikipedia [5].
3. A **run** (maximal repetition) is an interval whose smallest period `p` satisfies length ≥ 2p and which cannot be extended left or right with the same period; its exponent is length/p — Kolpakov & Kucherov [2], Bannai et al. Definition 1 [3].
4. A microsatellite/tandem unit is **primitive**: MISA drops "false type motifs [e.g. (TT)6 or (ACAC)5]" [4]. A period-`p` interval of length ≥ 2p whose unit is primitive has smallest period `p` (Fine–Wilf), so "run with primitive unit" and "run with smallest period p" coincide [3].
5. A string of length `n` has fewer than `n` runs (runs theorem) and all of them are computable from the longest Lyndon words of the two alphabet orders with LCE queries [3].
6. Detection of tandem repeats in strings can use suffix trees/arrays [5].

### 1.3 Documented Corner Cases

- Fewer than two contiguous copies is not a tandem repeat (k ≥ 2 floor) — Benson 1999 [1].
- One report per maximal run and unit length: rotations (`TA` in `ATATATA`) are not repeated; non-primitive multiples (`ATAT` in `ATATATAT`) are not reported; runs of different (non-multiple) periods may overlap and are each reported — [3][4].
- A trailing partial copy is not counted (`ATATATA` → `AT`×3) — run exponent ⌊length/p⌋ convention of MISA/pytrf (REP-STR-001).

### 1.4 Known Failure Modes / Pitfalls

1. Benson's TRF detects *approximate* copies; this detector is exact only. Approximate repeats: `RepeatFinder.FindApproximateTandemRepeats` (TRF model, REP-APPROX-001) / MCP `find_approximate_tandem_repeats`.
2. The pre-2026-11 implementation reported non-primitive units (`ATATATAT`, minUnit 1, minReps 2 → `AT`×4@0 **and** `ATAT`×2@0) — fixed (F35).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindTandemRepeats(sequence, minUnitLength = 2, minRepetitions = 2)` | GenomicAnalyzer | **Canonical** | Facade over `RepeatFinder.FindMicrosatellites` / run enumeration; no upper unit length. |
| `FindTandemRepeats(sequence, minUnitLength, minRepetitions, maxUnitLength)` | GenomicAnalyzer | **Canonical** | Restricted unit-length range (additive overload). |
| `FindTandemRepeats(string, min, minReps[, maxUnitLength])` | GenomicAnalyzer | **Canonical** (A1-21) | N/IUPAC-tolerant string overloads: canonical `RepeatFinder.FindMicrosatellites(string)` contract (units over A/C/G/T only, case-insensitive; `null` → `ArgumentNullException`, empty → none) |
| `RepeatFinder.FindMicrosatellitesWithEngine(..., Auto/Scan/Runs)` | RepeatFinder (internal) | Engine hook | Run enumeration (`MaximalRepetitions.FindRuns`) vs the per-unit-length scan, forced one by one in tests; since 2026-11-19 (R15, F78) the engine choice (> 256 unit lengths -> runs) lives in `RepeatFinder.FindMicrosatellites`, `GenomicAnalyzer.FindTandemRepeats` is a pure facade. |
| `find_tandem_repeats` | MCP Analysis | Wrapper | Delegates; optional `maxUnitLength`. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported repeat has Repetitions ≥ `minRepetitions` ≥ 2 and unit length in `[minUnitLength, maxUnitLength]` | Yes | Benson 1999 [1] |
| INV-2 | `TotalLength = Unit.Length × Repetitions` | Yes | period × copy number = repeat length [1] |
| INV-3 | `Position + Unit.Length × Repetitions ≤ sequence.Length` | Yes | copies contiguous and inside the sequence [1][5] |
| INV-4 | Unit is primitive; the run is maximal | Yes | [2][3][4] |
| INV-5 | Order = unit length ascending, then start ascending; both engines identical | Yes | contract |
| INV-6 | String input (A1-21): a window containing any non-A/C/G/T symbol never forms a repeat (unit over A/C/G/T only, MISA `[acgt]`); an A/C/G/T string = the `DnaSequence` overload; both engines identical on N/IUPAC strings | Yes | `generate_tandem_oracle.py --ambiguous` (892 rows), differential run engine vs per-unit scan |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Exact trinucleotide unit | `ATGATGATG`, minUnit 3, minReps 3 | `ATG`, pos 0, 3 copies | [1][5] |
| M2 | Worked example | `ATTCGATTCGATTCG`, minUnit 5, minReps 2 | `ATTCG`, pos 0, 3 copies, total length 15 | [5] |
| M3 | minRepetitions floor | sequence with only 2 copies, minReps 3 | not reported | k ≥ 2 / threshold [1] |
| M4 | No repeat / empty | `ACGT` / `""` | empty result | definition |
| M5 | Primitive unit | `ATATATAT`, minUnit 1, minReps 2 | `AT`×4@0 only (not `ATAT`×2); `AAAAAA` → `A`×6 only | [3][4], oracle |
| M6 | Python oracle | 1156 rows | exactly the oracle hit strings (unit:pos:copies, order) | oracle (3 independent Python methods) |
| M7 | Run enumeration = scan | oracle rows + 1500 random + 62 000 exhaustive + 60 long strings | identical sets/order | Kolpakov–Kucherov runs, brute force |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Position is 0-based | tandem starting at index 0 | Position == 0 | contract |
| S2 | TotalLength invariant | any reported repeat | Unit.Length × Repetitions | INV-2 |
| S3 | Left-end phase | `GCAGCAGCAGCAGG` minUnit 3 | `GCA`×4@0; `TCAGCAGCAGCAGT` → `CAG`×4@1 | MISA/pytrf convention |
| S4 | Partial copy | `ATATATA` | `AT`×3@0 | exponent ⌊len/p⌋ |
| S5 | Order / overlap | `AAAGAAAGAAAG` minUnit 1 | `A`×3@0/4/8 then `AAAG`×3@0 | contract |
| S6 | maxUnitLength | oracle seed `ACGACGACGACGTGTTAGCTTAGCTTAGC` | exact lists for max 1/3/5; `minUnit > maxUnit` throws (eager) | oracle |
| S7 | Parameter validation | minUnit 0, minReps 1, null | `ArgumentOutOfRangeException` / `ArgumentNullException`, eager | contract |
| S8 | MCP | `find_tandem_repeats` incl. `maxUnitLength` | sourced values | wrapper |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Telomere motif | `TTAGGG` ×4 | unit `TTAGGG`, 4 copies | vertebrate telomere repeat |
| C2 | 1 Mbp | bounded (2–60, ≥ 5 copies), unbounded, homopolymer, period-8 | planted repeat found, ≤ 6 s (Debug) | complexity |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Tests.cs` — canonical sourced fixture (27 tests: worked examples, thresholds, invariants); expectations unchanged by the contract fix (every expected hit is a primitive maximal run).
- `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Oracle_Tests.cs` — new (B09 F35–F38): Python oracle (M6), brute force, both engines (M7), contract S3–S7, complexity C2.
- `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_ExactTandems_Oracle_Tests.cs` — new (B09 A2-34): MUMmer `exact-tandems` oracle (3630 rows, per-start smallest-unit mapping) and the canonical-rotation option (least rotation == brute force on every oracle hit; `OriginalUnit` / `FullSequence` unchanged; string overload with N).
- `Fuzzing/GenomicTandemFuzzTests.cs`, `Fuzzing/RepeatsFuzzTests.cs` (row 14) — the fuzz oracle was the old greedy scan (non-primitive units) → replaced by the maximal-run definition; seeded-tandem test now seeds primitive roots.
- `tests/Seqeron/Seqeron.Mcp.Analysis.Tests/FindTandemRepeatsTests.cs` — `AAAAA` no longer reports the `AA`×2 tiling; `maxUnitLength` binding added.

### 5.2 Coverage Classification

All M/S/C items ✅ covered (see §4); no ❌/⚠ remain.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | By default the reported unit is the motif phase at the run's left end (MISA/pytrf convention, REP-STR-001); `canonicalRotation: true` (A2-34) reports its least rotation (Booth 1980). | §1.3, algorithm doc §5.4 |

---

## 7. Open Questions / Decisions

1. **DECISION — canonical engine.** The former duplicate-Registry resolution ("same method as REP-TANDEM-001, nothing to do") hid a real defect: the method was an independent O(n²·m) scan reporting non-primitive units (B04 request, DUP_MAP §17). It now delegates to the canonical `RepeatFinder.FindMicrosatellites` (maximal runs, primitive units, ACGT) and uses the run enumeration of the runs theorem for wide unit-length ranges. REP-TANDEM-001 keeps the summary helpers; this unit owns the `FindTandemRepeats` contract.
2. **DECISION — approximate repeats are not wrapped.** `RepeatFinder.FindApproximateTandemRepeats` (TRF) and MCP `find_approximate_tandem_repeats` are the canonical approximate detector; a pass-through in `GenomicAnalyzer` would duplicate a public API with a different result type, and neither the Benson paper nor the unit contract calls for it.
