# Validation Report: GENOMIC-TANDEM-001 — Tandem Repeat Detection (GenomicAnalyzer)

- **Validated:** 2026-06-16; re-validated and re-implemented 2026-11-18 (B09 review campaign, F35–F38)
- **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int minUnitLength = 2, int minRepetitions = 2)` and the additive `(…, int maxUnitLength)` overload → `IEnumerable<TandemRepeat>`; MCP `find_tandem_repeats`
- **Stage A verdict:** PASS-WITH-NOTES (the 2026-06 description accepted every period interpretation of a region as a separate tandem repeat, including non-primitive multiples; the sourced definition — Kolpakov–Kucherov runs / MISA primitive motifs — reports one maximal run per unit length with a primitive unit)
- **Stage B verdict:** FAIL → fixed (non-primitive units reported; O(n²·m) substring scan; no unit-length bound) → PASS
- **End-state:** CLEAN (fixed)
- **Relationship:** the 2026-06 session classified this unit a "duplicate of REP-TANDEM-001, resolved by consolidation". The shared method was nevertheless a stand-alone implementation that disagreed with the canonical `RepeatFinder.FindMicrosatellites` (B04 request, DUP_MAP §17); REP-TANDEM-001 keeps the summary helpers, this unit owns the `FindTandemRepeats` contract.

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES. **Stage B:** FAIL -> PASS. **State:** FIXED.
- **Fixes (B09.md, F-numbers):** F35 facade over `RepeatFinder.FindMicrosatellites` (non-primitive units no longer reported), F36 `maxUnitLength` overload + run-enumeration engine (Bannai et al. Lyndon runs, `MaximalRepetitions`), F37 documentation, F38 independent oracle, F50 N / IUPAC strings, F78 (R15, 'B09 for B04') the engine choice moved into the canonical `RepeatFinder.FindMicrosatellites`; `GenomicAnalyzer.FindTandemRepeats` is a pure facade.
- **Audit round 2 (A2-34 / A2-35, F98 / F99):** MUMmer `exact-tandems` oracle (3 630 rows) and the `canonicalRotation` option (least rotation, Booth 1980); string overloads of `FindTandemRepeats` with a rotation option.
- **Oracle numbers:** three-method Python oracle 1 156 rows + 62 000 exhaustive strings; 1 Mbp unbounded ~ 3 s.
- **Open items:** registered (B11, not finished): `GenomeAnnotator.cs:979` private tandem scan -> `GenomicAnalyzer.FindTandemRepeats(seq, 1, minCopies, 60)` (QUEUE R16).
- The unit's TestSpec (`tests/TestSpecs/GENOMIC-TANDEM-001.md`), Evidence (`docs/Evidence/GENOMIC-TANDEM-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

## Stage A — Description

### Sources opened (2026-11-18)

1. **Bannai, I, Inenaga, Nakashima, Takeda, Tsuruta, "The 'Runs' Theorem"** — https://arxiv.org/abs/1406.0263, abstract + full PDF. Definition 1 (run: smallest period p, length ≥ 2p, not extendable left/right), Lemma 3/6/7/8, Theorem 9 (ρ(n) < n), §4.1–4.2 (stack computation of the longest Lyndon words; candidate extension by LCE on w and its reverse; "$ ≺0 a (and thus a ≺1 $)").
2. **MISA `misa.pl`** (mirror https://raw.githubusercontent.com/cfljam/SSR_marker_design/master/misa.pl) line 110: "reject false type motifs [e.g. (TT)6 or (ACAC)5]".
3. **Kolpakov & Kucherov, FOCS 1999** — full paper opened 2026-11-19 (`monge.univ-eiffel.fr/~koutcher/PAPERS/KolpakovKucherovFOCS99.pdf`): period = smallest period, repetition = exponent >= 2, maximal repetition = extension by one letter in either direction gives a bigger period; linear number of runs. (Earlier: search snippets only.)
4. **Benson 1999** — the paper's pages stay unreachable (PMC reCAPTCHA, OUP 403, Europe PMC 403); the TRF author's site `tandem.bu.edu/trf` was opened 2026-11-19 and states the definition ("two or more adjacent, approximate copies of a pattern of nucleotides"; pattern size 1-2000); see GENOMIC-TANDEM-001-Evidence.md "Primary sources opened". Earlier: search snippets ("two or more contiguous, approximate copies of a pattern of nucleotides"). Wikipedia "Tandem repeat" (2026-06 fetch): worked example `ATTCG ATTCG ATTCG`.
5. EMBOSS `etandem` 6.6.0 run for information only (its scoring model is not an exact-run enumerator; see the Evidence file).

### Formula / model check

A run `(a, e, p)` is a maximal interval with period `p` (`S[x] = S[x+p]` for `a ≤ x < e−p`), length ≥ 2p; the reported repeat is the run's primitive unit `S[a..a+p)` with `⌊(e−a)/p⌋` complete copies (≥ `minRepetitions`), once per run, at its left end. Fine–Wilf: a period-p interval of length ≥ 2p whose unit is primitive has smallest period p, so "run with primitive unit" = "run with smallest period p" (Bannai Definition 1). This is also the per-motif-size convention of MISA and pytrf (REP-STR-001). The 2026-06 claim "reporting each unit-length interpretation is sound" is refuted for non-primitive multiples: `ATATATAT` has the single run `AT`×4.

### Edge-case semantics

Empty / shorter than `minUnitLength × minRepetitions` → empty; `minUnitLength ≥ 1`, `minRepetitions ≥ 2`, `maxUnitLength ≥ minUnitLength`, non-null sequence — eager `ArgumentOutOfRangeException` / `ArgumentNullException`; partial last copy not counted; rotations not repeated; runs of different non-multiple periods may overlap and are each reported; order = unit length, then start.

### Independent cross-check (numbers)

Three independent Python implementations of the definition (brute-force Kolpakov–Kucherov runs, boolean-array runs, regex), asserted equal on 1156 `(sequence, minUnit, minReps, maxUnit)` rows (`scripts/analysis/generate_tandem_oracle.py`, `TestData/Analysis/tandem_oracle.tsv`). Worked values: `ATATATAT` (1,2) → `AT:0:4`; `AAAAA` (1,2) → `A:0:5`; `ATATATA` (2,2) → `AT:0:3`; `GCAGCAGCAGCAGG` (3,2) → `GCA:0:4`; `AAAGAAAGAAAG` (1,2) → `A:0:3;A:4:3;A:8:3;AAAG:0:3`; `ACGACGACGACGTGTTAGCTTAGCTTAGC` (1,2) → `T:14:2;T:19:2;T:24:2;GT:11:2;ACG:0:4;TTAGC:14:3`.

**Stage A: PASS-WITH-NOTES.**

## Stage B — Implementation

### Code path reviewed

Before: `GenomicAnalyzer.FindTandemRepeatsCore` — for each unit length, for each start, `Substring` compare, skip to the end of the block (O(n²·m)).
After: `GenomicAnalyzer.FindTandemRepeats` (validation) → `RepeatFinder.FindMicrosatellites(seq, min, max, reps)` (canonical per-unit-length maximal-run scan, REP-STR-001) when ≤ 256 unit lengths are admissible, else `MaximalRepetitions.FindRuns` (Bannai et al. Lyndon-array runs algorithm, O(n log² n)). MCP `find_tandem_repeats` delegates (optional `maxUnitLength`).

### Cross-verification recomputed against the code

* Old code vs oracle: 276 / 1156 rows differed, **all** by extra non-primitive units (487 extra hits; on 20 000 further random strings 26 747 extras; 0 missing, 0 non-maximal) — the DUP_MAP/B04 remark "not left-maximal" does not hold for the primitive hits, the defect is primitivity only.
* New code (both engines) vs oracle: 1156 / 1156 rows identical (also vs `RepeatFinder.FindMicrosatellites` directly); vs in-test brute force: 1500 random strings (bounded and unbounded), 62 000 exhaustive strings (`{A,C}` ≤ 14, `{A,C,G}` ≤ 9) for the run enumeration, 60 long low-entropy strings (2.5 kb) engine-vs-engine.
* Complexity (1 Mbp random, Debug build): units 2–60 / ≥ 5 copies ≈ 1 s (scan, ≈ 17 ms per unit length); unbounded unit length ≈ 3 s (run enumeration; the scan would need up to 500 000 passes); worst case 1 Mbp period-8 string ≈ 6 s. Old code: cubic.

### Variant / delegate consistency

Scan engine = run-enumeration engine = oracle (set and order); the MCP wrapper maps the same values; the 3-argument overload = the 4-argument overload with `int.MaxValue`.

### Test quality audit

Fuzz oracle (`GenomicTandemFuzzTests.Oracle`) replayed the old greedy scan → replaced by the run definition (brute force); the seeded-tandem fuzz test now seeds primitive roots; the MCP `AAAAA` test no longer expects the `AA`×2 tiling. All new expectations are sourced oracle values; no assertion was loosened.

### Findings / defects

- **F35** non-primitive units reported; method a separate copy of the canonical scan (B04 / DUP_MAP §17) → facade over `RepeatFinder.FindMicrosatellites`.
- **F36** no `maxUnitLength`, O(n²·m) (cubic) complexity, unbounded default impractical → additive overload + MCP parameter; run-enumeration engine (Bannai et al. 2017) for wide ranges.
- **F37** docs: algorithm docs said approximate tandem repeats need external tools (TRF); they are `RepeatFinder.FindApproximateTandemRepeats` / MCP `find_approximate_tandem_repeats` (B04 request); MCP description/docs of `find_tandem_repeats` (AA tiling example).
- **F38** independent three-method oracle + engine differential tests.

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES. Stage B: PASS after fixes (FAIL before).** End-state **CLEAN**.
- No `GenomicAnalyzer` pass-through for approximate repeats: the canonical `RepeatFinder.FindApproximateTandemRepeats` (TRF 4.10.0 model) and MCP tool already exist with a different result type; the unit contract and Benson (1999) do not call for a second wrapper.
- Registered for B04 (finished batch; behaviour/performance change of a public method, so not done here): `RepeatFinder.FindMicrosatellites` with a very wide `maxUnitLength` could use `MaximalRepetitions.FindRuns` instead of one scan per unit length (QUEUE.md "Requests to finished batches").

- **Audit r1 (A1-21):** `FindTandemRepeats(string, …)` overloads for N/IUPAC input over the canonical `RepeatFinder.FindMicrosatellites(string)`; 892-row ambiguous oracle; MCP accepts IUPAC (B09.md F50).
- **Audit r3 (A3-8):** exact ports of EMBOSS `etandem` / `equicktandem` (`GenomicAnalyzer.FindTandemRepeatsEtandem` / `FindTandemRepeatsEquicktandem`, MCP `find_tandem_repeats_emboss`); 3000 / 3000 oracle rows equal (B09.md F136).
