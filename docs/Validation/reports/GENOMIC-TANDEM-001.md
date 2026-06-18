# Validation Report: GENOMIC-TANDEM-001 — Tandem Repeat Detection (GenomicAnalyzer)

- **Validated:** 2026-06-16   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int minUnitLength = 2, int minRepetitions = 2)` → `IEnumerable<TandemRepeat>` (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:100`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **Relationship:** Documented duplicate of REP-TANDEM-001 (same method, same class). Resolved by consolidation; this session re-validated the shared implementation and its canonical fixture against externally-retrieved sources and closed one test-coverage gap.

## Duplicate determination

The TestSpec (§7) and algorithm doc declare GENOMIC-TANDEM-001 a duplicate Registry entry for the
identical `GenomicAnalyzer.FindTandemRepeats` already shipped under REP-TANDEM-001, with the shared
canonical fixture `GenomicAnalyzer_TandemRepeat_Tests.cs`. I confirmed this is a genuine duplicate:
there is exactly one implementation of `FindTandemRepeats` and one canonical test fixture. Per the
duplicate-elimination rule, no second test file was created. However, "duplicate" does not exempt the
algorithm/tests from validation, so I independently validated the shared method and fixture against
sources retrieved THIS session (not the repo's own artifacts).

## Stage A — Description

### Sources opened this session (URLs + extracted numbers)

1. **Wikipedia, "Tandem repeat"** — https://en.wikipedia.org/wiki/Tandem_repeat (accessed 2026-06-16).
   - Definition (verbatim): "tandem repeats occur in DNA when a pattern of one or more nucleotides is
     repeated and the repetitions are directly adjacent to each other."
   - Worked example (verbatim): "ATTCG ATTCG ATTCG … the sequence ATTCG is repeated three times." →
     unit `ATTCG`, period 5, copy number 3, total length 15, 0-based start 0.
   - Classification: dinucleotide (2 nt), trinucleotide (3 nt); microsatellites/STRs (short units),
     minisatellites (10–60 nt), macrosatellites (~1000 nt).
   - Detection: "Tandem repeats in strings … can be efficiently detected using suffix trees or suffix
     arrays."
2. **Wikipedia, "Microsatellite"** — https://en.wikipedia.org/wiki/Microsatellite (accessed 2026-06-16).
   - Repeat-unit range "generally ten nucleotides or less"; "typically repeated 5–50 times".
   - Vertebrate telomeric repeat: "the hexanucleotide repeat motif TTAGGG in vertebrates."
   - Forensic STRs "are all tetra- or penta-nucleotide repeats."
   - Trinucleotide repeat disorders incl. Huntington's disease.
3. **Benson, G. (1999), Tandem Repeats Finder, NAR 27(2):573–580** —
   https://academic.oup.com/nar/article/27/2/573/1061099 (accessed 2026-06-16). DOI 10.1093/nar/27.2.573.
   - Formal definition (verbatim): "A tandem repeat in DNA is two or more contiguous, approximate
     copies of a pattern of nucleotides." → minimum copy number k ≥ 2.
   - Period size = "the most common matching distance between corresponding characters in the
     alignment" (= unit length for an exact repeat); copy number = number of pattern copies; consensus
     size may differ from period size.

### Formula / model check

The algorithm doc §2.2 states: a tandem repeat exists when `S[p .. p+k|U|) = U^k` with `k ≥ 2`, where
period = `|U|` and copy number = `k`. This matches Benson's definition (restricted to *exact* copies)
and the Wikipedia worked example. The exact-only restriction (vs Benson's approximate copies) is an
explicitly documented simplification (algorithm doc §5.3); over exact tandems the two definitions
coincide, so the simplification is correct (a subset, not a wrong rule). INV-1/2/3 are genuine
mathematical properties (k ≥ 2; TotalLength = period × copies; block within bounds).

### Edge-case semantics

Empty → no hits; no-tandem → empty; whole-sequence tandem → one spanning result; same region under
multiple periods → reported once per qualifying period interpretation. Each interpretation is itself a
valid tandem repeat under the formal definition, so reporting all of them is sound (and documented).

### Independent cross-check (hand computation, traced to sources)

- `ATTCGATTCGATTCG` (Wikipedia worked example): unit `ATTCG`, period 5, copies 3, total 15, start 0.
  Hand trace of the implementation matches the source exactly (confirmed by execution, below).
- `AAAA`: periods 1 (A×4) and 2 (AA×2) both satisfy the definition — both genuine tandem repeats.

**Stage A verdict: PASS.** Description, formulae, conventions (0-based, k ≥ 2, period = unit length),
and edge cases are correct and source-backed. One minor doc nit (the TestSpec/Evidence quote Benson's
definition as "two or more contiguous … copies" eliding "approximate"; the elision is harmless because
the implementation is deliberately exact-only and that simplification is separately documented).

## Stage B — Implementation

### Code path reviewed

`GenomicAnalyzer.FindTandemRepeats` (`GenomicAnalyzer.cs:100–126`) and `struct TandemRepeat`
(`:540–554`). Brute-force scan: for each unit length `minUnitLength .. seq.Length/minRepetitions`,
for each start, extract candidate unit, count consecutive exact copies, yield if count ≥
`minRepetitions`, then `start = pos - unitLen` to skip past the block within the current unit-length
pass. `TotalLength = Unit.Length × Repetitions`; `FullSequence = Unit repeated Repetitions times`.

### Formula realised correctly?

Yes. Direct execution (temporary observation harness, since removed) against externally-derived
expectations:

| Input | minUnit/minReps | Source expectation | Code output | Match |
|-------|-----------------|--------------------|-------------|-------|
| `ATTCGATTCGATTCG` | 5 / 2 | ATTCG, pos 0, 3 copies, len 15 (Wikipedia) | ATTCG, 0, 3, 15 | ✅ |
| `ATTCGATTCGATTCG` | 2 / 2 | only ATTCG period-5 (no smaller exact sub-period) | ATTCG, 0, 3, 15 | ✅ |
| `AAAA` | 1 / 2 | A×4 and AA×2 (both valid periods) | A/0/4 + AA/0/2 | ✅ |
| `ATATAT` | 1 / 2 | only AT×3 (no contiguous A or T runs) | AT, 0, 3 | ✅ |
| `ATGATGCCATGATGATG` | 3 / 2 | ATG×2 @0 and ATG×3 @8 (skip logic) | as expected | ✅ |
| `` (empty) | 2 / 2 | empty | empty | ✅ |

### Variant / delegate consistency

No `*Fast`/delegate variants for this method; `TotalLength`/`FullSequence` derived properties verified
consistent with `Unit`/`Repetitions`.

### Numerical robustness

Integer indexing only; loops bounded by `seq.Length`; no overflow/precision concerns on stated ranges.
Thresholds are not guarded (documented accepted assumption); the public tests only exercise sensible
(≥1) thresholds, which is the supported contract.

### Test quality audit (against sources, not code)

Canonical fixture `GenomicAnalyzer_TandemRepeat_Tests.cs`. All correctness assertions use **exact**
sourced values (no Greater/AtLeast/Contains weakening, no widened tolerances, no skips). Public surface
fully exercised: `Unit`, `Position`, `Repetitions`, `TotalLength`, `FullSequence`; both `minUnitLength`
and `minRepetitions` filters; empty/no-repeat boundaries; 0-based position; long repeat; whole-sequence;
adjacent distinct runs; case normalization; and three invariant property tests. The lone range-style
assertion is a performance bound (`< 30 s`), not a correctness assertion — legitimate.

**Gap found and fixed (Stage-B defect under the test-quality gate):** the TestSpec §5.2 claimed M2 —
the literal Wikipedia worked example `ATTCGATTCGATTCG` — was "covered", but the fixture only used
*substitute* motifs (telomere/tetra/penta). The canonical sourced example itself was never asserted.
I added `FindTandemRepeats_WikipediaWorkedExample_ATTCG`, locking the externally-retrieved values
(unit `ATTCG`, pos 0, 3 copies, total 15) directly from the Wikipedia article fetched this session.

### Findings / defects

- **F1 (fixed this session):** missing test for the verbatim Wikipedia worked example asserted as
  covered by the spec. Closed by adding the sourced exact-value test. No production-code defect found.

## Verdict & follow-ups

- **Stage A: PASS.** **Stage B: PASS-WITH-NOTES** (one missing sourced test, now added; minor doc nit
  on the elided "approximate" wording — harmless given the documented exact-only simplification).
- **Test-quality gate: PASS** after fix. Full unfiltered suite: **Failed: 0, Passed: 6619** (1
  benchmark skipped by design). Changed test file builds warning-free.
- **End-state: CLEAN** (duplicate confirmed; shared algorithm independently validated; the one
  coverage gap completely fixed in-session).
