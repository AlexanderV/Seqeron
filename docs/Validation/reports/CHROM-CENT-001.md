# Validation Report: CHROM-CENT-001 — Centromere position classification (Levan 1964 centromere index / arm ratio)

- **Validated:** 2026-06-25   **Area:** Chromosome
- **Canonical method(s) (this unit's OWN surface):**
  - `ChromosomeAnalyzer.CalculateArmRatio(centromerePosition, chromosomeLength)` — arm-length ratio (p/q)
  - `ChromosomeAnalyzer.ClassifyChromosomeByArmRatio(armRatio)` — public Levan category classifier
  - `ChromosomeAnalyzer.AnalyzeCentromere(chrName, seq, windowSize, minAlphaSatelliteContent)` → private
    `DetermineCentromereType` (the Levan classifier applied to a detected region), `EstimateRepeatContent`,
    `CalculateGcVariability`
  - record `CentromereResult`; constant `AlphaSatelliteConsensus`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`
- **Test files:** `ChromosomeAnalyzer_Centromere_Tests.cs`, `ChromosomeAnalyzer_MutationKillers_Tests.cs`,
  `ChromosomeAnalyzerTests.cs` (Levan classification rows)
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ❌→✅ FAIL-THEN-FIXED (one real defect in `ClassifyChromosomeByArmRatio`, fully fixed this session)
- **State:** ✅ CLEAN

> **Scope note.** This re-validation targets CHROM-CENT-001's OWN canonical surface — centromere
> position / centromere index / **Levan classification**. The alpha-satellite 171-bp / CENP-B detection
> (`DetectAlphaSatellite`, `FindCenpBBoxes`) and HOR detection (`DetectHigherOrderRepeat`) added during
> the limitation-elimination campaign are **separate, already-CLEAN units** (CHROM-ALPHASAT-001 /
> CHROM-HOR-001) and are referenced but not re-litigated here. Suprachromosomal-family / specific
> α-satellite-family assignment is a documented data-blocked boundary (needs curated reference HOR
> libraries the library does not embed).

## Stage A — Description

### Sources opened this session & what they confirm
- **Levan A, Fredga K, Sandberg AA (1964)** "Nomenclature for centromeric position on chromosomes",
  *Hereditas* 52(2):201–220 (Wiley `10.1111/j.1601-5223.1964.tb01953.x`). Located via SCIRP / Wiley /
  Google Scholar. The classification is defined by the **arm ratio r = L/S** (long arm / short arm) and
  equivalently by the **centromeric index ci = 100·p/(p+q)** (short-arm fraction).
- **Standard cytogenetics operationalisation of Levan's point table** (Tandfonline Caryologia
  `10.1080/00087114.2015.1032614`; Vedantu CBSE; Wikipedia *Centromere*), retrieved this session, agree
  on the numeric boundaries:

  | Category | Arm ratio r = L/S | Centromeric index ci = 100·p/(p+q) |
  |---|---|---|
  | **Metacentric (m)** | 1.0 ≤ r ≤ 1.7 | 37.5 – 50.0 |
  | **Submetacentric (sm)** | 1.7 < r ≤ 3.0 | 25.0 – 37.5 |
  | **Subtelocentric (st)** | 3.0 < r < 7.0 | 12.5 – 25.0 |
  | **Acrocentric (a/t)** | r ≥ 7.0 | < 12.5 |
  | **Telocentric (T)** | r = ∞ (one arm absent, p = 0) | 0 |

### Formula / threshold check
- **Centromere index = short-arm fraction × 100** and **arm ratio = long/short** — confirmed verbatim
  against Levan 1964 (CI = 100·p/(p+q); r = q/p). The four defining numeric cut-points **1.7, 3.0, 7.0**
  (and ci 37.5/25.0/12.5) are the canonical Levan boundaries.
- **Boundary convention.** The shared endpoints (1.7, 3.0, 7.0) are conventionally assigned to the more
  symmetric (lower) category at 1.7/3.0 and to acrocentric at 7.0 (Wikipedia "1.0–1.7", "≥7"). This is a
  measure-zero tie-break in real karyotyping; the repo's `DetermineCentromereType` uses exactly that
  reading (`≤1.7→m, ≤3.0→sm, <7.0→st, else a`). Documented, acceptable.

### Independent cross-check (hand computation, exact numbers)
With r normalised to long/short (≥1):

| Input (p/q or q/p) | r = long/short | Levan category (hand) |
|---|---|---|
| p/q 1.0 | 1.00 | Metacentric |
| p/q 0.7 → q/p 1.43 | 1.43 | Metacentric (< 1.7) |
| q/p 1.70 (boundary) | 1.70 | Metacentric |
| q/p 1.71 | 1.71 | Submetacentric |
| p/q 0.5 → q/p 2.0 | 2.00 | Submetacentric |
| q/p 3.0 (boundary) | 3.00 | Submetacentric |
| q/p 3.01 / p/q 0.3 (q/p 3.33) | 3.01 / 3.33 | Subtelocentric |
| p/q 0.2 → q/p 5.0 | 5.00 | Subtelocentric |
| q/p 6.99 | 6.99 | Subtelocentric |
| q/p 7.0 (boundary) | 7.00 | Acrocentric |
| p/q 0.1 → q/p 10.0 | 10.0 | Acrocentric |
| one arm absent (p = 0) | ∞ | Telocentric |

End-to-end (centromere position c of length 100 → p=c, q=100−c, r=max/min):
c=50→r=1.00 m; c=40→r=1.50 m; c=33→r=2.03 sm; c=20→r=4.00 st; c=10→r=9.00 a. (All re-derived in Python
this session, not read from code.)

### Stage A verdict
**PASS.** Centromere index formula, arm-ratio definition, the four boundaries (1.7/3.0/7.0; ci
37.5/25.0/12.5), the five categories, and the telocentric single-arm semantics are all confirmed against
Levan 1964 and standard cytogenetics references retrieved this session.

## Stage B — Implementation

### Code paths reviewed (`ChromosomeAnalyzer.cs`)
- `DetermineCentromereType` (`:568`) — r = qArm/pArm with qArm=Max, pArm=Min (so r≥1); `pArm==0`→Telocentric;
  switch `≤1.7 m / ≤3.0 sm / <7.0 st / else a`. **Correct vs Levan** (verified IEEE-754 boundaries
  1.7→m, 3.0→sm, 7.0→a). Reached only via `AnalyzeCentromere`.
- `CalculateArmRatio` (`:1299`) — returns p/q where p = centromerePosition, q = length−centromerePosition;
  guards `centromerePosition≤0 || length≤0 || q≤0` → 0. Correct (a ratio, not a classifier).
- `ClassifyChromosomeByArmRatio` (`:1313`) — **DEFECT FOUND.**

### Defect (Stage B FAIL → FIXED)
The public `ClassifyChromosomeByArmRatio` originally classified a **p/q** ratio with ad-hoc cuts:
`[0.9,1.1]→m, [0.5,0.9)→sm, [0.2,0.5)→a, <0.2→T`, plus a mirrored `(1.1,2.0]→sm, (2.0,5.0]→a, >5.0→T`.
Compared to Levan 1964 this is wrong across most of the range:

| p/q | q/p (= r) | Levan (correct) | Original code | 
|---|---|---|---|
| 0.7 | 1.43 | **Metacentric** | Submetacentric ✗ |
| 0.49 | 2.04 | **Submetacentric** | Acrocentric ✗ |
| 0.30 | 3.33 | **Subtelocentric** | Acrocentric ✗ |
| 0.20 | 5.00 | **Subtelocentric** | Acrocentric ✗ |
| 0.10 | 10.0 | **Acrocentric** | Telocentric ✗ |

Three independent errors: (1) the **Subtelocentric category never appears**; (2) the thresholds do not
correspond to Levan's 1.7/3.0/7.0; (3) a finite arm ratio is labelled "Telocentric" (telocentric is the
single-arm p=0 case only). This also made the two public Levan classifiers in the class
(`DetermineCentromereType` vs `ClassifyChromosomeByArmRatio`) **mutually inconsistent**.

The two existing tests (`ClassifyChromosomeByArmRatio_BoundaryTable`,
`ClassifyChromosomeByArmRatio_ClassifiesCorrectly`) **green-washed the defect** — their expected strings
were echoes of the wrong cuts (e.g. they asserted p/q 0.3 → "Acrocentric").

**Fix.** Rewrote `ClassifyChromosomeByArmRatio` to (a) treat `armRatio ≤ 0` as Telocentric (degenerate
single-arm), (b) normalise any p/q or q/p input to `r = max(armRatio, 1/armRatio) ≥ 1`, and (c) apply the
Levan switch `≤1.7 m / ≤3.0 sm / <7.0 st / else a` — now identical in convention to
`DetermineCentromereType`. Added the Levan 1964 citation to the XML doc.

### Tests rewritten / added (Levan-sourced, hand-computed — no code echoes)
- `ChromosomeAnalyzer_MutationKillers_Tests.ClassifyChromosomeByArmRatio_BoundaryTable` — replaced the
  echoing table with exact boundaries r=1.0/1.7/1.71/3.0/3.01/5.0/6.99/7.0/7.01/21.0 (q/p form) **and**
  the reciprocal p/q forms (0.7/0.5/0.3/0.2/0.1), each expected value hand-derived from Levan.
- New `ClassifyChromosomeByArmRatio_DegenerateSingleArm_IsTelocentric` (ratio 0.0 / −1.0 → Telocentric).
- New `ArmRatioPipeline_ClassifiesPerLevan` — end-to-end `CalculateArmRatio→ClassifyChromosomeByArmRatio`
  for centromere 50/40/33/20/10 of length 100 → m/m/sm/st/a (hand-derived r = 1.00/1.50/2.03/4.00/9.00).
- `ChromosomeAnalyzerTests.ClassifyChromosomeByArmRatio_ClassifiesCorrectly` — corrected the echoing rows
  (0.7→Metacentric, 0.3→Subtelocentric, 0.1→Acrocentric, 1.5→Metacentric, 3.0→Submetacentric,
  10.0→Acrocentric).

### Other surfaces (verified, unchanged)
- `AnalyzeCentromere` + `DetermineCentromereType`: region detection is a declared generic
  tandem-repeat-density heuristic (k=15 k-mer repeat × low GC variability); region coords 0-based
  half-open (`Length = End − Start`); the Levan classification it applies is correct; non-repetitive
  random input is not flagged (`Type=Unknown`). 24-test Centromere fixture all green. (Detection scope and
  the `AlphaSatelliteContent` generic-score naming were validated in prior rounds; the alpha-satellite-/
  HOR-specific detectors are the separate CHROM-ALPHASAT-001 / CHROM-HOR-001 units.)

### Stage B verdict
**FAIL → FIXED.** One real correctness defect (`ClassifyChromosomeByArmRatio` diverged from Levan 1964 and
omitted Subtelocentric) found and fully fixed; the green-washing tests that locked the wrong values were
rewritten to Levan-sourced hand-computed expectations; the two public Levan classifiers are now consistent.

## Verdict & follow-ups
- **Stage A:** ✅ PASS — Levan 1964 centromere index / arm-ratio thresholds (1.7/3.0/7.0; ci 37.5/25.0/12.5)
  and the five categories confirmed against external sources retrieved this session.
- **Stage B:** ❌→✅ FAIL-THEN-FIXED — `ClassifyChromosomeByArmRatio` corrected to Levan; tests de-green-washed.
- **State:** ✅ **CLEAN** — defect fully fixed, tests trace to Levan/hand-computation, full unfiltered
  `dotnet test Seqeron.sln -c Debug` **Failed: 0** (Genomics 18817 passed), 0 warnings on changed files.
- **Documented boundary (acceptable):** suprachromosomal-family / specific α-satellite-family assignment is
  data-blocked (needs curated reference HOR libraries not embedded in the library).

---

## Historical note

Earlier rounds (PASS-WITH-NOTES) examined the **detection criterion** of `AnalyzeCentromere` (generic
tandem-repeat heuristic vs alpha-satellite specificity) and drove the additive, opt-in
`DetectAlphaSatellite` / `FindCenpBBoxes` / `DetectHigherOrderRepeat` capabilities (171-bp monomer,
CENP-B box `YTTCGTTGGAARCGGGA`, HOR period/copy-number) — now the separate CLEAN units
CHROM-ALPHASAT-001 / CHROM-HOR-001. Those additions did not touch the Levan classification surface that is
the subject of this report; this session re-validated that surface independently and fixed the
`ClassifyChromosomeByArmRatio` defect that the prior rounds did not examine.
