# Validation Report: REP-STR-001 — Microsatellite / Short Tandem Repeat (STR) Detection

- **Validated:** 2026-06-25 (fresh re-validation; supersedes 2026-06-24)   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindMicrosatellites(DnaSequence|string, int minUnitLength=1, int maxUnitLength=6, int minRepeats=3)` (+ `CancellationToken`/`IProgress<double>` overloads); `GetTandemRepeatSummary` integration; **`RepeatFinder.FindApproximateTandemRepeats(DnaSequence|string, int minPeriod=1, int maxPeriod=6, int minScore=50)`** (opt-in approximate / imperfect / interrupted detector, TRF model — added 2026-06-24); **`RepeatFinder.ComputeBernoulliStatistics(string, int period, double expectedMatchProbability=0.80)`** (opt-in TRF Bernoulli statistical-significance measures, Benson 1999 — added 2026-06-24).
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`
- **Test file(s):** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_Microsatellite_Tests.cs` (34 tests), `RepeatFinder_ApproximateTandemRepeats_Tests.cs` (13 tests), `RepeatFinderTests.cs` (summary), `Properties/RepeatFinderProperties.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (code correct; sole note = TRF k-tuple **seeding** / genome-scale-performance residual, out of scope)
- **State:** ✅ CLEAN

> **2026-06-25 — fresh independent re-validation (REP-STR-001).** This unit was reset to ⬜ pending
> after the limitation-elimination campaign added the TRF Bernoulli significance scoring
> (`ComputeBernoulliStatistics`) and the opt-in approximate detector (`FindApproximateTandemRepeats`).
> Re-validated this session against authoritative EXTERNAL sources retrieved 2026-06-25 (Wikipedia
> Microsatellite; TRF definitions page tandem.bu.edu/trf/trf.definitions.html for the Bernoulli PM/PI
> model and adjacent-copy statistics scope; Richard, Kerrest & Dujon 2008 / Benson 1999 for STR and TRF
> definitions). All cross-check numbers were re-derived by hand AND reproduced by running the actual
> code in an independent harness (not by trusting the repo tests). Full suite green
> (Seqeron.Genomics.Tests = 18780 passed, 0 failed). No defect found → **CLEAN**.
>
> **2026-06-24 — limitation fix (REP-STR-001).** The "perfect-repeat-only" scope note has been
> superseded: an **opt-in approximate (imperfect / interrupted) tandem-repeat detector** was added per
> the Tandem Repeats Finder model (Benson 1999) — see §"Approximate detection" below. The default
> perfect-repeat path is unchanged.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Microsatellite** (fetched 2026-06-24) — confirmed: repeat unit "typically ten nucleotides or less" (library uses 1–6 bp), "repeated 5–50 times"; verbatim examples **"TATATATATA is a dinucleotide microsatellite"** and **"GTCGTCGTCGTCGTC is a trinucleotide microsatellite"**; history: first microsatellite "a polymorphic **GGAT** repeat in the human **myoglobin** gene" by **Weller, Jeffreys** et al. (1984). All match the TestSpec evidence table verbatim.
- **MISA — MIcroSAtellite identification tool** (literature search, MISA-web PMC5870701) — MISA default per-motif thresholds: 10 mono / 5 di / 4 tri / 3 tetra-penta-hexa, **configurable**. Establishes that a per-length minimum-count threshold is a *tool parameter convention*, not a fixed biological constant — so the library's single configurable `minRepeats` (default 3, enforced ≥2) is a legitimate, standard choice.
- **Tandem Repeats Finder (Benson 1999) / microsatellite literature** — perfect vs imperfect/compound distinction confirmed; this detector targets **perfect** (uninterrupted, exact) tandem repeats, in-scope per the TestSpec (all M/S cases use perfect tracts).

### Definitions & conventions verified
- **Motif length range:** 1–6 bp (mono–hexa) — matches defaults and `ClassifyRepeatType`.
- **Minimum copy threshold:** single configurable `minRepeats` (default 3), enforced **≥2** (a repeat requires ≥2 occurrences). Sourced as a tool parameter; no per-length minimum mandated. No divergence.
- **Perfect vs imperfect:** perfect-only. In-scope.
- **Coordinates:** `Position` **0-based**; span half-open `[Position, Position+TotalLength)`; `TotalLength = RepeatUnit.Length × RepeatCount`; copy number = `RepeatCount`. Standard, internally consistent.
- **Canonicalisation / overlap:** smallest non-redundant unit preferred (`AAAA`→`A×4`, not `AA×2`); a result whose span is fully contained in an already-reported span is suppressed; reading-frame rotations that are not mutually contained (e.g. `AC×5` and `CA×5`) are both reported. Documented and tested.

### Independent cross-check (hand computation)
- `CACACACA` → unit `CA`, 4 copies, span [0,8) (end index 7). ✔
- `AAAA` → `A×4` mono (canonicalisation, `IsRedundantUnit("AA")`=true). ✔
- `TATATATATA` → `TA×5`, pos 0, len 10 (Wikipedia di example). ✔
- `GTCGTCGTCGTCGTC` → `GTC×5`, pos 0, len 15 (Wikipedia tri example). ✔
- Designed `(GATA)n`: `AAGATAGATAGATAGATAAA` → `AGAT×4` pos 1 (span [1,17)) + `GATA×4` pos 2 (span [2,18)); neither contains the other → both reported. ✔ Matches the half-open containment policy.
- Single occurrence is NOT a repeat: `minRepeats≥2` enforced; `ATAT` (AT×2) with minRepeats=3 → empty. ✔

### Independent cross-check reproduced THIS session (2026-06-25, live code in a standalone harness)
Microsatellite detector (`FindMicrosatellites`, exact values from running the code):
- `CACACACA` (1–6,3) → `CA×4 @0 len8 Dinucleotide` (smallest unit, not CACA×2). ✔
- `AAAA` (1–6,3) → `A×4 @0 len4 Mononucleotide` (AA rejected by `IsRedundantUnit`). ✔
- `TATATATATA` (2,2,3) → `TA×5 @0 len10` (Wikipedia di example). ✔
- `GTCGTCGTCGTCGTC` (3,3,3) → `GTC×5 @0 len15` (Wikipedia tri example). ✔
- `AAGATAGATAGATAGATAAA` (4,4,3) → `AGAT×4 @1` + `GATA×4 @2` (both reported; neither span contains the other). ✔
- `ATAT` (2,2,3) → empty (AT×2 < minRepeats 3, sub-threshold). ✔
- `CAGCAGCAGTAGCAGCAG` (1–6,3) → `CAG×3 @0` only (perfect detector fragments at the interruption). ✔

Bernoulli statistics (`ComputeBernoulliStatistics`, adjacent-copy, hand-derived AND reproduced):
- `CACACACACA` p2 → trials=8, matches=8, PM=1.000000, PI=0, E[matches]=8, meets=True. ✔
- `CAGCAGCAGTAGCAGCAG` p3 → trials=15, matches=13, mismatches=2, PM=0.866667 (=13/15), PI=0, E[matches]=13, meets=True. ✔
- `CACACATACACA` p2 → trials=10, matches=8, mismatches=2, PM=0.800000 (=8/10, exactly the default threshold, inclusive ⇒ meets=True), E[matches]=8. ✔
- `ACACTGTG` p4 → trials=4, matches=0, PM=0, meets=False (below default PM 0.80). ✔

**Stage A finding:** Description faithful to authoritative sources. No defect.

## Stage B — Implementation

### Code path reviewed
`RepeatFinder.cs:25–260` — `FindMicrosatellites` overloads → `FindMicrosatellitesCore` (174–219) / `FindMicrosatellitesCancellable` (101–172); `IsRedundantUnit` (221–246); `ClassifyRepeatType` (248–260); `MicrosatelliteResult` record (630–641).

### Realises the validated description? (evidence)
- Outer loop `unitLen = min..max`; inner greedy exact extension (`seq.Substring(j,unitLen) == unit`) → perfect tandem counting. ✔
- Validation: `minUnitLength≥1`, `maxUnitLength≥minUnitLength`, `minRepeats≥2`, on **all four** public overloads (string + DnaSequence, plain + cancellable) — M16 parity holds. ✔
- 0-based `Position`, `TotalLength = repeats*unitLen`, `RepeatCount`, `FullSequence = unit×count`. ✔
- `IsRedundantUnit` rejects composite units (`AA`, `ATAT`) by testing every divisor sub-length → smallest-unit canonicalisation. ✔
- Containment suppression via `reported` half-open span set (`r.Start<=i && r.End>=end`). ✔
- `ClassifyRepeatType`: 1→Mono … 6→Hexa, else Complex. ✔

### Cross-verification table recomputed vs code (tests executed)
| Case | Expected (source) | Code | Match |
|------|-------------------|------|-------|
| `ACGTAAAAAACGT` (1–6,3) | A×6 pos 4 | A×6 pos 4 | ✔ |
| `TATATATATA` (2,3) | TA×5 pos 0 | TA×5 pos 0 | ✔ |
| `GTCGTCGTCGTCGTC` (3,3) | GTC×5 pos 0 | GTC×5 pos 0 | ✔ |
| `ATGCAGCAGCAGCAGCAGTGA` (3,3) | GCA×5 pos2 + CAG×5 pos3 | same | ✔ |
| `AAGATAGATAGATAGATAAA` (4,4,3) | AGAT×4 pos1 + GATA×4 pos2 | same | ✔ |
| `AAGGATGGATGGATGGATAA` (4,4,3) | GGAT×4 pos 2 | GGAT×4 pos 2 | ✔ |
| `AAAGAATTCGAATTCGAATTCAAA` (6,3) | GAATTC×3 | GAATTC×3 | ✔ |
| `ATATATAT` (2–4,2) | AT×4 (not ATAT×2) | AT×4 | ✔ |
| `ATAT` (2,3) | empty | empty | ✔ |
| `ATATAT` (2,3) | AT×3 | AT×3 | ✔ |
| `""` (1–6,3) | empty | empty | ✔ |

### Variant/delegate consistency
String and DnaSequence overloads produce identical results (`StringOverload_ProducesSameResults`); both uppercase input; all four overloads validate parameters identically. Cancellable variant mirrors core logic (adds progress/cancellation only). ✔

### Test quality audit
44 tests in scope (34 microsatellite + summary) all assert exact sourced values (unit, count, position, type, total length), not "no-throw" tautologies. Stage-A edge cases covered: empty, too-short, exactly/below threshold, entire-sequence-repeat, canonicalisation, rotations, lowercase, N-handling (DnaSequence rejects / string accepts), null/invalid params on both surfaces, cancellation smoke. Deterministic.

### Findings / notes
- **Perfect-repeat-only (historical scope note — RESOLVED 2026-06-24):** interrupted/imperfect or compound microsatellites were previously reported only as separate perfect tracts. The opt-in `FindApproximateTandemRepeats` now reports them as one approximate locus (see below). The default perfect detector is unchanged.

### Approximate detection (added 2026-06-24 — TRF model, Benson 1999)
- **Source (retrieved this session):** Benson G (1999) "Tandem repeats finder", Nucleic Acids Res 27(2):573–580, https://doi.org/10.1093/nar/27.2.573; TRF README/definitions (Benson-Genomics-Lab/TRF, tandem.bu.edu). Captured verbatim: approximate-repeat definition; reported statistics (period size, copy number, consensus size, percent matches/indels between adjacent copies overall, alignment score); recommended scoring (match +2, mismatch −7, indel −7) and report threshold Minscore = 50.
- **Method:** for each period, align the observed window against a whole number of tandem copies of the majority-rule consensus (reusing `SequenceAligner.GlobalAlign` with a flat-indel TRF scoring matrix); read the alignment columns for the statistics; report best-score-first when score ≥ minScore.
- **Hand-verified worked examples (all pass):** `CACACACACA` → CA, 5 copies, 100% matches, 0% indels, score 20; `CAGCAGCAGTAGCAGCAG` (one substitution) → CAG, 6 copies, 94.4̄% matches, score 27 (perfect detector fragments to CAG×3); `CACACATACACA` → 91.6̄% matches, score 15; `(CAG)×10`-with-one-deletion → 9.6̄ copies, 96.6̄% matches, 3.3̄% indels, score 51.
- **Residual (honest):** candidate discovery is a deterministic exhaustive (start, period) scan (not whole-genome scale). The per-repeat statistics follow Benson (1999) exactly. Recorded in `LIMITATIONS.md` and the algorithm doc §5.3.

### Bernoulli statistical-significance scoring (added 2026-06-24 — TRF model, Benson 1999)
- **Source (retrieved this session):** Benson (1999) NAR 27(2):573–580 (https://doi.org/10.1093/nar/27.2.573); TRF detailed-description / definitions pages (tandem.bu.edu/trf/trf.desc.html, trf.definitions.html) and Benson-Genomics-Lab/TRF README. Captured **verbatim**: "We model alignment of two tandem copies of a pattern of length n by a sequence of n independent Bernoulli trials"; "P(Heads), which we also call PM or matching probability, represents the average percent identity between the copies"; "PI or indel probability … the average percentage of insertions and deletions between the copies"; statistics are "between adjacent copies … not between the sequence and the consensus pattern"; defaults "PM = .80 and PI = .10". The sum-of-heads R(d,k,pM) and random-walk W(d,pI) variables (k-tuple seeding) were also captured verbatim and identified as the non-reproducible (simulation-table) residual.
- **Method added (opt-in):** `RepeatFinder.ComputeBernoulliStatistics(string repeatTract, int period, double expectedMatchProbability = 0.80)` → `TandemRepeatBernoulliStatistics` (PM, PI, BernoulliTrials, Matches/Mismatches/Indels, ExpectedMatches = PM·d, MeetsExpectedMatchProbability). Estimates PM/PI **between adjacent copies** (segment into period copies, align each adjacent pair with the recommended +2/−7/−7 scoring) — distinct from the consensus-based percent matches of `FindApproximateTandemRepeats`. `FindMicrosatellites` and `FindApproximateTandemRepeats` are unchanged.
- **Hand-verified worked examples (all pass — 10 new tests):** perfect `CACACACACA` → PM = 8/8 = 1.0, PI = 0, E[matches] = 8; `CAGCAGCAGTAGCAGCAG` (one substitution) → PM = 13/15 ≈ 0.8667 (adjacent-copy, not the 17/18 consensus value), PI = 0, E[matches] = 13; `CACACATACACA` → PM = 8/10 = 0.80 (exactly the default threshold, inclusive), E[matches] = 8; `ACACTGTG` → PM = 0, below threshold; indel `(CAG)×10`-with-deletion → PI > 0 with match/mismatch/indel partitioning the trials; exposed defaults PM = .80 / PI = .10.
- **Residual (honest, narrowed):** the per-repeat **Bernoulli statistical measures are now computed**; what remains unbundled is only TRF's **k-tuple seeding** (the R(d,k,pM) 95% sum-of-heads percentile cut-off and the W(d,pI) random-walk band), whose percentile values come from TRF's non-redistributable simulation tables — i.e. a **whole-genome-scale performance** index, not a per-repeat-correctness gap. The deterministic exhaustive scan already finds every candidate a seed would; seeding is a performance heuristic, not a capability.

## Verdict & follow-ups
- **Stage A: PASS** — description matches Wikipedia (Microsatellite: motif "generally ten nucleotides or less", repeated 5–50×; verbatim TATATATATA di / GTCGTCGTCGTCGTC tri examples; GGAT myoglobin / Weller & Jeffreys 1984), MISA convention, Richard et al. 2008, and the TRF Bernoulli model (TRF definitions page: PM=.80, PI=.10, statistics "between adjacent copies … not between the sequence and the consensus pattern", match +2 / mismatch & indel −7 stringent, Minscore 50). All worked examples hand-verified.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the validated perfect-STR definition, the opt-in TRF approximate detector, and the opt-in TRF Bernoulli statistical measures; coordinates (0-based, half-open), thresholds (`minRepeats≥2`), copy counts, canonicalisation/overlap, TRF statistics, and PM/PI/expected-matches all correct and independently reproduced this session. All public methods/overloads exercised (string + DnaSequence, plain + cancellable `FindMicrosatellites`; string + DnaSequence `FindApproximateTandemRepeats`; `ComputeBernoulliStatistics`). Tests assert exact sourced/hand-computed values (no code-echo, no no-throw tautologies) and cover every Stage-A edge case (empty, single-char, all-N, sub-threshold, interrupted/imperfect, motif-longer-than-seq, each motif length 1–6). Sole note: TRF k-tuple **seeding** (genome-scale-performance) residual — out of scope per the campaign boundary; the per-repeat Bernoulli significance is fully implemented and validated.
- **State (2026-06-25): ✅ CLEAN.** Fresh independent re-validation of the campaign-added methods (`ComputeBernoulliStatistics`, `FindApproximateTandemRepeats`) plus the canonical `FindMicrosatellites`. All cross-check numbers re-derived from primary sources and reproduced by running the live code. Full suite green: Seqeron.Genomics.Tests = 18780 passed, 0 failed; 0 warnings (no code changed this session). Root registry flipped `☐ → ☑`.
- **No defects logged.**
