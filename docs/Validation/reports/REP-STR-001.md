# Validation Report: REP-STR-001 — Microsatellite / Short Tandem Repeat (STR) Detection

- **Validated:** 2026-06-12   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindMicrosatellites(DnaSequence|string, int minUnitLength=1, int maxUnitLength=6, int minRepeats=3)` (+ `CancellationToken`/`IProgress<double>` overloads), with `GetTandemRepeatSummary` integration.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`
- **Test file(s):** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_Microsatellite_Tests.cs` (34 tests), `RepeatFinderTests.cs` (summary), `Properties/RepeatFinderProperties.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (perfect-repeat-only scope; documented below)

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Microsatellite** (https://en.wikipedia.org/wiki/Microsatellite) — confirmed verbatim:
  - Repeat unit length spans single nucleotides up to ~6 bp (some authors extend the upper bound), classified as **mono- (1 bp), di- (2 bp), tri- (3 bp), tetra- (4 bp), penta- (5 bp), hexa-nucleotide (6 bp)**.
  - "typically repeated **5–50 times**."
  - Exact examples: **"the sequence TATATATATA is a dinucleotide microsatellite, and GTCGTCGTCGTCGTC is a trinucleotide microsatellite."**
  - Forensic markers "are all **tetra- or penta-nucleotide repeats**."
  - History: "The first microsatellite was characterised in 1984 ... by **Weller, Jeffreys** and colleagues as a polymorphic **GGAT** repeat in the human **myoglobin** gene."
- **Wikipedia: Trinucleotide repeat disorder** (https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder) — Huntington CAG: **normal 6–35, pathogenic 36–250**; trinucleotide repeats are "a subset of a larger class of unstable microsatellite repeats." Confirms the spec's medical-relevance framing.
- **Richard GF et al. (2008) MMBR** (cited) — comprehensive review establishing perfect/imperfect/compound repeat categories; consistent with the spec's perfect-repeat definition used here.

### Definitions & conventions verified
- **Motif length range:** 1–6 bp (mono–hexa). Matches implementation defaults (`minUnitLength=1`, `maxUnitLength=6`).
- **Minimum copy threshold:** The library uses a single configurable `minRepeats` (default **3**), enforced as **≥ 2** (a "repeat" requires ≥2 occurrences). This is a *tool parameter*, not a biological constant; Wikipedia's "5–50" is a typical-occurrence statement, not a calling threshold, so no per-motif-length minimum is mandated by the sources. No divergence.
- **Perfect vs imperfect:** Detector finds **perfect** (uninterrupted, exact) tandem repeats only. In-scope per the test spec (all M/S cases use perfect tracts).
- **Coordinates:** `Position` is **0-based**; span is `[Position, Position+TotalLength)` (half-open); `TotalLength = RepeatUnit.Length × RepeatCount`; copy number = `RepeatCount`. Standard and internally consistent.
- **Canonicalization (overlap/nesting/reading-frame):** smallest non-redundant unit preferred (`AAAA` → mono `A×4`, not `AA×2`); a result fully contained in an already-reported span is suppressed; different reading-frame rotations of one tract that are NOT mutually contained (e.g. `AC×5` and `CA×5`) are both reported. This policy is explicitly documented and tested by the spec.

### Independent cross-check (hand computation)
- `CACACACA` (8 bp) → unit `CA`, 4 copies, span 0–7 (end index 7). **Matches** the prompt's worked example.
- `AAAA` → `A×4` (mono), not `AA×2`. Canonicalization confirmed.
- `TATATATATA` → `TA×5`, pos 0, len 10 (Wikipedia di example).
- `GTCGTCGTCGTCGTC` → `GTC×5`, pos 0, len 15 (Wikipedia tri example).
- `ATGCAGCAGCAGCAGCAGTGA` → `GCA×5` pos 2 + `CAG×5` pos 3 (CAG expansion).

All hand computations match both the spec's expected values and a standalone Python re-implementation of the algorithm.

**Stage A finding:** Description is faithful to authoritative sources. No defect.

## Stage B — Implementation

### Code path reviewed
`RepeatFinder.cs:25–246` — `FindMicrosatellites` overloads → `FindMicrosatellitesCore` / `FindMicrosatellitesCancellable`; `IsRedundantUnit` (221–246); `ClassifyRepeatType` (248–260); `MicrosatelliteResult` record (560–571).

### Realises the validated description? (evidence)
- Motif-length loop `unitLen = minUnitLength..maxUnitLength` with greedy exact extension (`seq.Substring(j,unitLen) == unit`) → perfect tandem repeat counting. ✔
- Threshold guards: `minUnitLength≥1`, `maxUnitLength≥minUnitLength`, `minRepeats≥2`. ✔ (string overloads validate identically — M16/parity confirmed.)
- 0-based `Position`, `TotalLength = repeats*unitLen`, `RepeatCount = repeats`, `FullSequence = unit×count`. ✔
- `IsRedundantUnit` skips composite units (`AA`, `ATAT`) → smallest-unit canonicalization. ✔
- Containment suppression via `reported` set → nested/overlapping policy as specified. ✔
- `ClassifyRepeatType`: 1→Mono … 6→Hexa, else Complex. ✔

### Cross-verification table recomputed vs code (tests executed)
| Case | Expected (source) | Code result | Match |
|------|-------------------|-------------|-------|
| `ACGTAAAAAACGT` (1–6,3) | A×6 pos 4 | A×6 pos 4 | ✔ |
| `TATATATATA` (2,3) | TA×5 pos 0 | TA×5 pos 0 | ✔ |
| `GTCGTCGTCGTCGTC` (3,3) | GTC×5 pos 0 | GTC×5 pos 0 | ✔ |
| `ATGCAGCAGCAGCAGCAGTGA` (3,3) | GCA×5 pos2 + CAG×5 pos3 | same | ✔ |
| `AAGGATGGATGGATGGATAA` (4,3) | GGAT×4 pos 2 | GGAT×4 pos 2 | ✔ |
| `AAAGAATTCGAATTCGAATTCAAA` (6,3) | GAATTC×3 | GAATTC×3 | ✔ |
| `ATATATAT` (2–4,2) | AT×4 (not ATAT×2) | AT×4 | ✔ |
| `ATAT` (2,3) | empty | empty | ✔ |
| `ATATAT` (2,3) | AT×3 | AT×3 | ✔ |
| `""` (1–6,3) | empty | empty | ✔ |

### Variant/delegate consistency
String and DnaSequence overloads produce identical results (test `StringOverload_ProducesSameResults`); both uppercase input; both validate parameters identically. Cancellable variant mirrors core logic. ✔

### Test quality audit
Tests assert exact sourced values (unit, count, position, type, total length), cover all Stage-A edge cases (empty, too-short, exactly/below threshold, entire-sequence-repeat, motif at boundaries via positioned cases, canonicalization, lowercase, N-handling, null/invalid params, cancellation smoke). 63 microsatellite-related tests pass.

### Findings / notes
- **Perfect-repeat-only (scope note, not a defect):** the detector does not find imperfect/interrupted or compound microsatellites (e.g. `(CA)₅TG(CA)₅` reported as two perfect tracts, not one imperfect locus). This is consistent with the test spec, which defines and tests only perfect repeats. Imperfect-repeat detection would be a feature extension, not a correctness bug.

## Verdict & follow-ups
- **Stage A: PASS** — description matches Wikipedia (Microsatellite; Trinucleotide repeat disorder) and Richard et al. (2008); all worked examples hand-verified.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the validated perfect-STR definition; coordinates, thresholds, copy counts, and canonicalization all correct. Sole note: perfect-repeat-only scope (in-spec).
- **State: CLEAN** — no defect found; no code change required. Build green, 63 microsatellite tests pass.
- **No defects logged.**
