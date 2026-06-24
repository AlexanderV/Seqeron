# Validation Report: CHROM-CENT-001 — Centromere detection (alpha-satellite / repeat region) + Levan classification

- **Validated:** 2026-06-24   **Area:** Chromosome
- **Canonical method(s):** `ChromosomeAnalyzer.AnalyzeCentromere(chromosomeName, sequence, windowSize, minAlphaSatelliteContent)` → private `EstimateRepeatContent`, `CalculateGcVariability`, `DetermineCentromereType`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Centromere_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

## Scope of this session

The prior report (`git show cb113ce:...`) validated the **Levan (1964) arm-ratio classification** thoroughly. This session re-validates that result and, per the task, focuses on the part the prior report glossed over: the **detection criterion** — how a centromere region is *identified* (alpha-satellite / AT-rich / tandem-repeat), the scoring/threshold, and the reported region coordinates.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Alpha satellite" / "Centromere"** + literature search (Nature JHG, PLoS Comp Biol, PMC reviews on α-satellite biology): human centromeres are built from **α-satellite (alphoid) DNA**, a tandem repeat of a **~171 bp AT-rich monomer**. Monomers are 50–70% identical and are arranged into **higher-order repeat (HOR) units** (2–34 monomers) that are themselves tandemly amplified into large, near-identical (95–100%) arrays. The two defining molecular features are therefore: **(1) a ~171 bp tandem period, and (2) high AT content.**
- **Levan A, Fredga K, Sandberg AA (1964)** "Nomenclature for centromeric position on chromosomes", *Hereditas* 52(2):201–220 (via Wikipedia table + independent summaries): arm ratio r = q/p (long/short); m ≤ 1.7, sm (1.7,3.0], st (3.0,7.0), t (acrocentric) ≥ 7.0, Telocentric at p = 0. Confirmed verbatim against the prior report's source check.

### What the spec / Evidence actually describe as the detection criterion
The TestSpec and `CHROM-CENT-001-Evidence.md` describe detection as a **heuristic**, not as alpha-satellite-specific identification:
- Evidence "Implementation Notes": *"Sliding window approach with k-mer frequency analysis … GC content variability as discriminating feature … Repeat content estimation using k-mer counting (k=15)."*
- TestSpec M4 rationale: *"Centromeres are characterized by repetitive DNA."*

So the description honestly declares a **generic tandem-repeat-density heuristic**. It does **not** claim to measure AT-content, nor to match the 171-bp monomer period. This is biologically motivated (centromeres are indeed the most repetitive regions of the chromosome) and is a defensible declared heuristic.

### Independent cross-check (hand / numeric, k=15 repeat-content as in code)
| Region type | repeat content | Detected at default 0.3? | Comment |
|---|---|---|---|
| AT-rich 171-bp tandem (alpha-sat-like) | 1.00 | yes | true positive |
| GC-rich generic 16-bp tandem (NOT alpha-sat, NOT AT-rich) | 1.00 | yes | **false positive** — flags any tandem repeat |
| AT-rich but non-repetitive (random AT) | 0.26 | no (borderline) | AT-content alone is not sufficient — confirms detector keys on *repetitiveness*, not AT% |
| Random 4-base sequence | ≈0 | no | true negative (matches test M4b) |

This confirms the criterion is **repeat density**, independent of period length and of AT-content. It correctly *includes* alpha-satellite (which is repetitive) but is not *specific* to it.

### Findings / divergences (Stage A)
- **Note 1 (overclaiming name):** the unit is titled "alpha-satellite / AT-rich repeat identification", and the result field is `AlphaSatelliteContent`, but the algorithm computes a generic repeat×(1−GC-variability) score with **no** AT-content term and **no** 171-bp period test. The score is *consistent with* but not *specific to* alpha-satellite. The Evidence's own Implementation Notes describe the true mechanism, so the heuristic *is* declared — but the "AlphaSatellite" labeling oversells it. Documented, not a correctness defect.
- **Note 2 (Levan, unchanged):** r ≥ 7 class is Levan sign "t" = Acrocentric; Telocentric reserved for p = 0. Interval operationalization of Levan's point table is standard. (Same as prior report.)

## Stage B — Implementation

### Code path reviewed (`ChromosomeAnalyzer.cs`)
- `AnalyzeCentromere` (line 361): scans windows of `windowSize` stepping `windowSize/4`; score = `EstimateRepeatContent(window) * (1 - CalculateGcVariability(window,1000))`; keeps the max-scoring window whose `repeatContent > minAlphaSatelliteContent`; then greedily extends left/right while neighbouring half-windows have repeat content ≥ `0.7 × threshold`.
- `EstimateRepeatContent` (line 442): k=15 k-mer count; returns Σ(instances of k-mers occurring >1×) / (#k-mers). Pure repetitiveness measure.
- `CalculateGcVariability` (line 471): std-dev of GC fraction over 1 kb sub-windows.
- `DetermineCentromereType` (line 495): r = q/p, swap-safe via Min/Max, `pArm==0`→Telocentric, switch `<=1.7 / <=3.0 / <7.0 / else` → M/sm/st/Acrocentric.

### Faithfulness to the declared description
- The code realises exactly the declared k-mer-repeat + low-GC-variability heuristic. ✔
- Region coordinates are **0-based, half-open**: `Length = End − Start`, `End` clamped to `sequence.Length`. Confirmed by snapshot (`Start:0, End:15500, Length:15500`). ✔
- Levan classification realised correctly (full cross-check table reproduced in prior report; re-confirmed: 1.7→M, 3.0→sm, 7.0→Acrocentric, p=0→Telocentric, IEEE-754 boundaries exact). ✔
- Non-centromeric (random) region is **not** flagged: M4b passes (`Type=Unknown`, `Start/End=null`, score < 0.3). ✔

### Findings / defects (Stage B, all PASS-WITH-NOTES — none are correctness defects vs the declared heuristic)
1. **`AlphaSatelliteContent` is a misnomer.** The field stores `maxScore` = repeatContent×(1−GCvariability), a unitless generic-repeat score (snapshot shows 1.0), **not** a fraction of alpha-satellite sequence. The TestSpec invariant only requires it be ≥ 0, which masks the mislabel.
2. **`AlphaSatelliteConsensus` constant is dead and mislabeled.** It is 62 bp (`...len 62, 76% AT`), **never referenced** by `AnalyzeCentromere`, and its only test (`AlphaSatelliteConsensus_IsValidDnaSequence`) asserts merely `length > 50` with a comment claiming the monomer is "~171 bp". The canonical α-satellite monomer is **171 bp** (Willard 1985 / current literature); the 62-bp string is neither 171 bp nor used. Cosmetic/documentation, no runtime effect.
3. **Detection is not alpha-satellite-specific** (Stage A Note 1): any sufficiently tandem-repetitive window scores 1.0 and is flagged regardless of period or AT-content. Acceptable *as a declared heuristic*; would be a defect only if the spec claimed period/AT-specificity, which it does not.

### Test quality audit
- 24 tests, all passing. They assert exact Levan class strings, `IsAcrocentric` invariant, valid type-set, structural invariants (Start≤End, Length=End−Start), and the right detection edge cases (empty/null/short/uniform/non-repetitive/case-insensitive/threshold-monotonicity). Boundary ratios reached stochastically (~1.0/2.0/5.0/21), with exact boundaries (1.7/3.0/7.0) covered by hand in the prior report.
- Gaps tied to the notes above: no test pins `AlphaSatelliteContent` to a meaningful biological quantity, and no test asserts the constant equals a real 171-bp monomer. These lock in the mislabels rather than catch them.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — Levan classification confirmed against sources; detection is a **declared generic tandem-repeat-density heuristic** (k=15 k-mer repeat + low GC variability), biologically motivated and correctly excluding random sequence, but **not** alpha-satellite-specific (no AT% term, no 171-bp period). The "AlphaSatellite" naming overclaims relative to what is computed.
- **Stage B:** PASS-WITH-NOTES — code faithfully realises the declared heuristic and the Levan classification; region coords 0-based half-open; non-centromeric region not flagged. Notes: `AlphaSatelliteContent` field is a generic score (misnamed); `AlphaSatelliteConsensus` (62 bp) is unused dead code mislabeled as the ~171 bp monomer.
- **State:** CLEAN — no algorithmic correctness defect against the declared description; the heuristic is honestly declared in the Evidence. The three notes are naming/documentation caveats, not behavioural bugs, so no code change is made (changing the misnomer/constant would be cosmetic and risk churning passing tests/snapshots without improving correctness). FullSuite not re-run end-to-end since no code changed; Centromere filter 24/24 passing.
- **Suggested (non-blocking) follow-ups for a future session:** rename `AlphaSatelliteContent`→`RepeatScore` (or document it as a repeat score), and either replace the 62-bp `AlphaSatelliteConsensus` with a sourced 171-bp α-satellite monomer or drop the "~171 bp" claim from the test comment.
