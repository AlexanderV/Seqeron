# Validation Report: CHROM-ALPHASAT-001 — Alpha-Satellite Monomer Detection

- **Validated:** 2026-06-25   **Area:** Chromosome
- **Canonical method(s):** `ChromosomeAnalyzer.DetectAlphaSatellite(string)`, `ChromosomeAnalyzer.FindCenpBBoxes(string)`; sourced constants `AlphaSatelliteMonomerLength` (171), `CenpBBoxConsensus` (`YTTCGTTGGAARCGGGA`)
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs:630` (`DetectAlphaSatellite`), `:692` (`FindCenpBBoxes`), `:32` / `:41` (constants).
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_AlphaSatellite_Tests.cs`

> Scope note: suprachromosomal-family / chromosome-specific HOR family assignment is a documented
> **data-blocked boundary** — it requires curated reference HOR libraries (T2T/CHM13) not present in the
> repo. `DetectHigherOrderRepeat` (period/copy-number structure, no family naming) is out of this unit's
> scope (CHROM-CENT-001 / HOR). This unit covers monomer-periodicity + AT-richness + CENP-B box detection.

## Stage A — Description

### Sources opened this session
- **Masumoto H, Masukata H, Muro Y, Nozaki N, Okazaki T (1989), J Cell Biol 109(4):1963–1973** — origin of the 17-bp CENP-B box.
- **PMC4843215** (Kugou et al., New World monkey CENP-B box) — quotes the consensus verbatim:
  *"The CENP-B box was described as a 17 bp motif (YTTCGTTGGAARCGGGA), based on sequence information derived from humans … in which the underlined nucleotides form the core recognition sequence."* Core = `TTCG…CGGG`.
- **PMC6121732** (Hartley & O'Neill 2019, *Alpha satellite DNA biology*) — *"Alpha satellite DNA is composed of fundamental 171bp monomeric repeat units."* and *"individual monomers within a HOR unit have 50–70% identity."*
- WebSearch corroboration: Willard 1985 / Waye & Willard 1987 for the 171-bp AT-rich monomer; CENP-B core "NTTCGNNNNANNCGGGN" (9 underlined nt sufficient).

### Facts confirmed
- **Monomer length = 171 bp** (Willard 1985; Waye & Willard 1987; PMC6121732). ✓ matches `AlphaSatelliteMonomerLength = 171`.
- **CENP-B box = 17-bp consensus `YTTCGTTGGAARCGGGA`**, Y = C/T, R = A/G (Masumoto 1989; PMC4843215). ✓ matches `CenpBBoxConsensus`.
- **AT-richness**: the alphoid monomer is widely described as AT-rich; threshold AT > 0.5 is a reasonable above-balance gate.

### Divergence noted (does NOT affect verdict)
- PMC6121732 renders the CENP-B box as `5'-T/CTCGTTGGAAA/GCGGGA-3'`, which expands to `YTCGTTGGAARCGGGA` = **16 bp** — it drops one `T` after the leading Y. This is a typo in that review; the **canonical Masumoto 1989 / PMC4843215 form is the 17-bp `YTTCGTTGGAARCGGGA`**, which is what the code uses. The code's doc-comment cites both but implements the correct 17-bp form. No action needed.

### Independent cross-check (hand-/Python-computed, no dependence on the C# code)
- Consensus length = 17. ✓
- All four IUPAC corner instances match at index 0: `CTTCGTTGGAAACGGGA` (Y=C,R=A), `TTTCGTTGGAAACGGGA` (Y=T,R=A), `CTTCGTTGGAAGCGGGA` (Y=C,R=G), `TTTCGTTGGAAGCGGGA` (Y=T,R=G). ✓
- `Y` rejects `A` (leading base A → no match). ✓
- Fixed consensus position 2 `T`→`A` → no match. ✓
- Box embedded after a 50-bp flank → reported at 0-based offset 50. ✓
- 171-bp tandem array (60 A + 40 T + 36 C + 35 G), 20 copies: periodicity@171 = **1.0**; AT content = **100/171 = 0.5847953…**. ✓
- Box-carrying monomer (77 A + 17-bp box + 77 T), 10 copies: **10** CENP-B boxes. ✓

## Stage B — Implementation

### Code path reviewed
- `DetectAlphaSatellite` (`:630`): too-short guard (< 171+5+1); AT content over ACGT bases only; periodicity = best self-similarity over period window [166,176]; CENP-B count; `IsAlphaSatellite = periodicity ≥ 0.50 AND AT > 0.50`.
- `FindCenpBBoxes` (`:692`) and `CountCenpBBoxes` (`:855`) → `MatchesIupac` (`:872`): Y→C/T, R→A/G, all other positions exact; 0-based ascending positions.

### Realises the description?
Yes. Periodicity is exactly "fraction of bases identical to the base `period` upstream"; AT content excludes non-ACGT from the denominator; CENP-B matching implements the Y/R degeneracy and exact match at fixed positions. All sourced constants present and correct.

### Cross-verification recomputed vs code
Ran the fixture (filtered) and the full suite: every fixture expected value (period 171, periodicity 1.0, AT 100/171, box count, IUPAC corners, negatives, offsets, empty/null/short) matches the independent hand computation above — values trace to the literature consensus, not to code echoes.

### Variant/delegate consistency
`CountCenpBBoxes` and `FindCenpBBoxes` share `MatchesIupac`, so count == positions.Count (verified by the box-count test = 10 == FindCenpBBoxes hits). Mixed-case input matches uppercase (test).

### Test-quality audit (HARD gate)
Covers every public surface and every Stage-A edge case: 171-bp tandem detected; CENP-B boxes located + position; IUPAC Y/R both resolutions + wrong-base + fixed-position-violation; non-satellite negatives (random, AT-rich-non-repetitive, GC-rich 16-bp period mismatch); empty / null / too-short; mixed case; sourced constants (171, 17-bp consensus). Expected numbers are literature- or hand-derived, not code echoes.

**Gap found & fixed (Stage-B test gap, not a code defect):** the prompt's edge-case list includes "non-ACGT" but the fixture had no test exercising non-ACGT input. Added `DetectAlphaSatellite_NonAcgtBases_ExcludedFromAtContentDenominator`: monomer 60 A + 40 T + 36 C + 30 G + 5 N → AT content = **100/166** (N excluded from the ACGT denominator, hand-derived), periodicity 1.0, still called alpha-satellite. Locks the defined non-ACGT semantics.

### Findings / defects
None in the code. The CENP-B consensus, monomer length, AT-richness gate, periodicity, and IUPAC degeneracy are all correct and faithfully implemented.

## Verdict & follow-ups
- **Stage A:** ✅ PASS — 171-bp monomer, 17-bp `YTTCGTTGGAARCGGGA` consensus, and AT-richness confirmed against Masumoto 1989 / PMC4843215 / PMC6121732. PMC6121732's 16-bp rendering is a typo; code uses the correct 17-bp form.
- **Stage B:** ✅ PASS — code realises the description; all cross-check values match the independent computation; one test gap (non-ACGT) closed.
- **State:** ✅ CLEAN. Family/HOR-library assignment remains a documented data-blocked boundary outside this unit's scope.
- Full unfiltered `dotnet test Seqeron.sln -c Debug`: Failed 0 (Seqeron.Genomics.Tests 18779 passed).
