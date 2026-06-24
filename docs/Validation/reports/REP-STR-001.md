# Validation Report: REP-STR-001 — Microsatellite / Short Tandem Repeat (STR) Detection

- **Validated:** 2026-06-24   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindMicrosatellites(DnaSequence|string, int minUnitLength=1, int maxUnitLength=6, int minRepeats=3)` (+ `CancellationToken`/`IProgress<double>` overloads); `GetTandemRepeatSummary` integration; **`RepeatFinder.FindApproximateTandemRepeats(DnaSequence|string, int minPeriod=1, int maxPeriod=6, int minScore=50)`** (opt-in approximate / imperfect / interrupted detector, TRF model — added 2026-06-24).
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`
- **Test file(s):** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_Microsatellite_Tests.cs` (34 tests), `RepeatFinder_ApproximateTandemRepeats_Tests.cs` (13 tests), `RepeatFinderTests.cs` (summary), `Properties/RepeatFinderProperties.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (TRF probabilistic seeding / significance scoring residual; documented)

> **2026-06-24 — limitation fix (REP-STR-001).** The "perfect-repeat-only" scope note has been
> superseded: an **opt-in approximate (imperfect / interrupted) tandem-repeat detector** was added per
> the Tandem Repeats Finder model (Benson 1999) — see §"Approximate detection" below. The default
> perfect-repeat path is unchanged. Per the validation campaign protocol, the unit Status in the root
> registry is **reset to `☐`** pending independent re-validation of the new method.

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
- **Residual (honest):** TRF's probabilistic k-tuple distance-list **seeding** and **sum-of-Bernoulli statistical-significance** scoring are not bundled; candidate discovery is a deterministic exhaustive (start, period) scan (not whole-genome scale). The per-repeat statistics follow Benson (1999) exactly. Recorded in `LIMITATIONS.md` and the algorithm doc §5.3.

## Verdict & follow-ups
- **Stage A: PASS** — description matches Wikipedia (Microsatellite), MISA convention, and TRF/literature; all worked examples hand-verified.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the validated perfect-STR definition plus the opt-in TRF approximate detector; coordinates, thresholds, copy counts, canonicalisation/overlap, and TRF statistics all correct. Sole note: TRF probabilistic seeding / significance-scoring residual (out of scope).
- **State (2026-06-24):** limitation fix landed — opt-in `FindApproximateTandemRepeats` added (TRF model). Status **reset to `☐`** in the root registry pending independent re-validation of the new method. Build green (0 warnings); 13 new approximate-detector tests + 44 perfect-STR/summary tests pass; full suite green.
- **No defects logged.**
