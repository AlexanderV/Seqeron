# Validation Report: ONCO-CHIP-001 — Clonal Hematopoiesis (CHIP) / variant-origin calling

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CallVariantOrigin` (Bolton 2020 strict matched-WBC rule — re-validation focus), with `IdentifyCHIPVariants`, `FilterCHIP`, `IsCanonicalChipGene` (gene+VAF heuristic, no-WBC fallback)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN (no defect found; no code changed)

This re-validation targets commit `eb5e0f7f` (limitations campaign), which added `CallVariantOrigin(WbcObservation)` implementing the Bolton 2020 matched-WBC origin rule. The pre-existing gene+VAF heuristic (`IdentifyCHIPVariants` / `FilterCHIP`) was validated in the prior report (2026-06-16, PASS-WITH-NOTES / CLEAN) and is unchanged; it is retained as the by-design no-WBC fallback.

## Stage A — Description

### Sources opened this session (external, not memory)

| # | Source | What it confirms (verbatim) | URL |
|---|--------|------------------------------|-----|
| 1 | Bolton et al. (2020) *Nat Genet* 52(11):1219–1226 (Methods) | CH-mutation read/VAF floor: *"We required a variant allele fraction of at least 2% and at least 10 supporting reads."* Blood-to-tumour fold rule: *"Variant calls that were present in the blood with a VAF of at least twice that in the tumor or 1.5 times the VAF if the tumor biopsy site was a lymph node were considered somatic."* Rationale: *"This ratio was chosen based on minimizing sensitivity and specificity of CH calls through simulations of leukocyte contamination in the tumor…"* | pmc.ncbi.nlm.nih.gov/articles/PMC7891089 |
| 2 | Steensma et al. (2015) *Blood* 126(1):9–16 | CHIP VAF gate: *"the mutant allele fraction must be ≥2% in the peripheral blood"*; CHIP = *"detectable somatic clonal mutations in genes recurrently mutated in hematologic malignancies…but who lack a known hematologic malignancy."* | pmc.ncbi.nlm.nih.gov/articles/PMC4624443 |
| 3 | Razavi 2019 / Arango-Argoty 2025 (per prior report + TestSpec) | matched cfDNA–WBC sequencing is the definitive tumour-vs-CH origin test; VAF–origin relationship alone "remains unclear" (motivates strict matched-WBC over gene+VAF) | (corroborated; carried from Evidence doc) |

### Formula / definition check

The implemented `CallVariantOrigin` rule (INV-7) — a variant is **Chip** iff a matched-WBC observation at the same locus has **WBC VAF ≥ τ_w (0.02) AND WBC reads ≥ ρ (10) AND WBC VAF ≥ φ·tumour VAF (φ = 2.0; 1.5 for lymph node)**, else **Tumor** — matches Bolton 2020 Methods exactly (#1). Each threshold maps to a verbatim source clause:

- `chipMinWbcVaf = 0.02`, inclusive ≥ ⟵ "variant allele fraction of at least 2%"
- `minWbcAltReads = 10`, inclusive ≥ ⟵ "at least 10 supporting reads"
- `wbcVafFold = 2.0`, inclusive ≥ ⟵ "at least twice that in the tumor"
- `LymphNodeWbcVafFold = 1.5` ⟵ "1.5 times the VAF if the tumor biopsy site was a lymph node"

"At least" → inclusive ≥ on all three thresholds: correct. Locus absent from WBC ⇒ no CH evidence ⇒ Tumor (INV-8): consistent with the matched-WBC design (no blood support ⇒ somatic/tumour).

### Edge-case semantics

All sourced and defined: boundary at exactly 2%/10 reads/2× (CHIP, inclusive); WBC VAF < 2% (Tumor); reads < 10 (Tumor); fold < φ (Tumor); lymph-node 1.5× changes the call; absent locus (Tumor); one call per input in input order (INV-9). No "implementation-defined" gaps.

### Independent cross-check (hand-computed vs the four prompt cases)

| Case | tumour VAF | WBC VAF | reads | φ | meetsVaf | meetsReads | meetsFold | Expected | Source |
|------|-----------|---------|-------|---|----------|-----------|-----------|----------|--------|
| WBC 3% / 12 reads / tumour 1% (3×) | 0.01 | 0.03 | 12 | 2.0 | 0.03≥0.02 ✓ | 12≥10 ✓ | 0.03 ≥ 0.02 ✓ | **CHIP** | Bolton 2020 |
| WBC 1% (below 2% gate) | — | 0.01 | — | 2.0 | 0.01≥0.02 ✗ | — | — | **not CHIP / Tumor** | Bolton/Steensma 2% gate |
| Borderline exactly 2× (OC4) | 0.01 | 0.02 | 10 | 2.0 | 0.02≥0.02 ✓ | 10≥10 ✓ | 0.02 ≥ 0.02 ✓ | **CHIP** (inclusive) | Bolton 2020 |
| Lymph-node 1.5× (OC7) | 0.25 | 0.40 | 40 | 2.0 / 1.5 | ✓ | ✓ | 0.40≥0.50 ✗ / 0.40≥0.375 ✓ | **Tumor (2.0) / CHIP (1.5)** | Bolton 2020 lymph-node |

All four match the Bolton rule.

### Findings / divergences (Stage A)

None. The rule, all thresholds, inclusivity, and the lymph-node exception are reproduced verbatim from Bolton 2020. The gene+VAF heuristic remains a labelled, sourced (Steensma 2015) fallback for the no-WBC case, by design.

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:9877–9944` (`CallVariantOrigin`); constants `ChipVafThreshold=0.02` (9599), `ChipMinWbcSupportingReads=10` (9783), `DefaultWbcVafFold=2.0` (9791), `LymphNodeWbcVafFold=1.5` (9798); `WbcObservation`/`VariantOrigin`/`VariantOriginCall` (9805–9850); `LocusKey` (9770–9775).

- **Formula realised correctly:** lines 9931–9937 compute `meetsVaf = wbc.Vaf >= chipMinWbcVaf`, `meetsReads = wbc.AltReads >= minWbcAltReads`, `meetsFold = wbc.Vaf >= wbcVafFold * variant.Vaf`, and set `Chip` iff all three; else default `Tumor`. This is the validated Bolton rule with inclusive ≥ on every threshold. Absent locus (no dictionary hit, 9924) leaves origin = `Tumor` and reports WbcVaf/WbcAltReads = 0 (INV-8).

- **Locus matching & robustness:** WBC observations indexed by exact `(chrom, 1-based pos, ref, alt)` tuple; on duplicate loci the highest-VAF observation is kept (9911) — conservative, never misses a confident CH call. Plain `double >=` comparisons, no precision-sensitive arithmetic, no div-by-zero (fold multiplies, never divides). One call per input in input order (9917–9941, INV-9).

- **Validation guards:** `variants`/`whiteBloodCellObservations` null → `ArgumentNullException`; `wbcVafFold < 1`, `chipMinWbcVaf ∉ (0,1]`, `minWbcAltReads < 1` → `ArgumentOutOfRangeException` (9884–9903). The `!(fold >= 1.0)` / `!(vaf > 0.0)` NaN-safe negated forms are correct.

- **Fallback consistency:** `CallVariantOrigin` deliberately does NOT apply the gene+VAF heuristic — a CH-driver-gene variant genuinely absent from matched WBC is called `Tumor`, resolving the prior report's "heuristic over-removal" note. `FilterCHIP`/`IdentifyCHIPVariants` unchanged and remain the no-WBC fallback.

### Cross-verification table recomputed vs code (tests OC1–OC9)

| Test | Inputs | Code result | Expected (Bolton) | Match |
|------|--------|-------------|-------------------|-------|
| OC1 | tumour 0.10, WBC 0.30, 40 | Chip | Chip | ✓ |
| OC2 | DNMT3A 0.40, no WBC | Tumor | Tumor | ✓ |
| OC3 | tumour 0.10, WBC 0.30, 9 reads | Tumor | Tumor (reads<10) | ✓ |
| OC4 | tumour 0.01, WBC 0.02, 10 | Chip | Chip (all boundary) | ✓ |
| OC5 | tumour 0.005, WBC 0.015, 50 | Tumor | Tumor (VAF<2%) | ✓ |
| OC6 | tumour 0.30, WBC 0.40, 40 | Tumor | Tumor (1.33×<2×) | ✓ |
| OC7 | tumour 0.25, WBC 0.40, 40 | Tumor (2.0) / Chip (1.5) | same | ✓ |
| OC8 | mixed CHIP+Tumor | order preserved, 1/input | same (INV-9) | ✓ |
| OC9 | null/empty/out-of-range | throws / empty | same (V7) | ✓ |

### Variant/delegate consistency

`CallVariantOrigin`, `FilterCHIP`, and `IdentifyCHIPVariants` share `LocusKey` and the same `ChipVafThreshold = 0.02`. Default parameters bind to the sourced public constants. No divergent constants.

### Test-quality audit (HARD gate)

| Check | Result |
|-------|--------|
| Sourced expectations, not code echoes | PASS — OC1–OC9 assert Bolton thresholds (2%/10/2×/1.5×) with rationale comments, not current output |
| No green-washing (exact values, no widened tolerances) | PASS — VAF asserted `Within(1e-10)`; boundary OC4 locks inclusive ≥; OC3 (9 reads) and OC6 (1.33×) would catch `>` regressions |
| All public methods/overloads covered | PASS — strict caller + lymph-node fold overload (OC7); all 3 guard params (OC9) |
| All Stage-A branches / edge / error cases | PASS — VAF floor, read floor, fold rule, lymph-node, absent-locus, order, invalid args all covered |
| Honest green (filtered suite Failed: 0) | PASS — 45/45 CHIP tests pass; project builds 0 warnings / 0 errors |

### Findings / defects (Stage B)

None. The implementation faithfully realises the Bolton 2020 rule; tests are evidence-based and cover every Stage-A branch.

## Verdict & follow-ups

- **Stage A:** PASS — Bolton 2020 read/VAF floors, 2×/1.5× fold rule, and inclusivity reproduced verbatim from the primary source; Steensma 2% gate confirmed. Four hand cross-checks (incl. exactly-2× boundary and lymph-node 1.5×) all match.
- **Stage B:** PASS — code at OncologyAnalyzer.cs:9877–9944 computes the validated rule exactly; OC1–OC9 recomputed against code all match; guards and order/multiplicity correct.
- **End-state:** CLEAN — no defect found; no code changed. The no-WBC gene+VAF heuristic fallback is by design.
- **Build/test env:** system dotnet 10.0.301 (net10.0); `--filter FullyQualifiedName~FilterCHIP` → 45 passed / 0 failed.

### Logged defects
- None.
