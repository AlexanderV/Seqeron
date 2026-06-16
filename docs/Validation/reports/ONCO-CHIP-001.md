# Validation Report: ONCO-CHIP-001 — Clonal Hematopoiesis (CHIP) Filtering for cfDNA Liquid Biopsy

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.IdentifyCHIPVariants`, `OncologyAnalyzer.FilterCHIP`, `OncologyAnalyzer.IsCanonicalChipGene`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** CLEAN (defects found were fully fixed in-session)

## Stage A — Description

### Sources opened this session (independent of the repo)

| # | Source | What it confirms | URL |
|---|--------|------------------|-----|
| 1 | Steensma et al. (2015) *Blood* 126(1):9–16 | "the mutant allele fraction must be ≥2% in the peripheral blood"; CHIP = somatic mutation in genes recurrently mutated in hematologic malignancies, no malignancy diagnosis | ashpublications.org/blood/article/126/1/9; pubmed 25931582 |
| 2 | Genovese et al. (2014) *NEJM* 371(26):2477–2487 (PMC4290021) | "Four genes (DNMT3A, TET2, ASXL1, and PPM1D) had disproportionately high numbers of somatic mutations"; DNMT3A 190, ASXL1 35, TET2 31; JAK2 V617F (24), SF3B1 K700E (9); mutant allele "present in less than 50% of the sequencing reads" | pmc.ncbi.nlm.nih.gov/articles/PMC4290021 |
| 3 | Razavi et al. (2019) *Nat Med* 25:1928–1937 (pubmed 31768066) | matched cfDNA + WBC sequencing (508 genes, >60,000×); 81.6% controls / 53.2% cancer-patient cfDNA mutations consistent with CH; matched-WBC presence is the origin test | nature.com/articles/s41591-019-0652-7 |
| 4 | Arango-Argoty et al. (2025) *NPJ Precis Oncol* 9:147 (PMC12092662) | "the exact relationship between VAF and variant origin remains unclear"; gene/VAF is a candidate flag; matched-WBC is definitive; top-3 DNMT3A/TET2/ASXL1 | pmc.ncbi.nlm.nih.gov/articles/PMC12092662 |
| 5 | Wikipedia *Clonal hematopoiesis* (cited primaries) | ≥2% VAF threshold; driver genes DNMT3A, TET2, ASXL1, JAK2, SF3B1, SRSF2, TP53, PPM1D (+ IDH1/2, FLT3, RUNX1) | en.wikipedia.org/wiki/Clonal_hematopoiesis |
| 6 | Sun, *Ann Transl Med* (Razavi commentary) | Origin call rests on **WBC presence/absence**, not gene identity; matched WBC "should also be analyzed" | atm.amegroups.org/article/view/33616/html |

### Formula / definition check

- **VAF threshold τ = 0.02, inclusive (≥):** confirmed verbatim by Steensma 2015 (#1, #5). INV-2 holds.
- **CHIP candidate flag `gene ∈ G ∧ VAF ≥ τ`:** matches the CHIP definition (#1). INV-1 holds.
- **Default panel {DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D}:** every gene is a documented canonical CH driver (#2, #5). The set is a labelled, overridable default (Assumption 1) and is a strict subset of the canonical driver list — defensible.
- **Matched-WBC subtraction (rule a):** a cfDNA variant present in matched WBC is CH-derived and removed regardless of gene (#3, #6). INV-4 holds.

### Edge-case semantics

VAF = 0.02 boundary (CHIP), non-CHIP gene at high VAF (not CHIP), empty-WBC, null/empty gene — all defined and sourced. The WBC alt-read presence cutoff (≥1, Wan 2020) is an assay-specific, labelled, configurable assumption — acceptable.

### Independent cross-check (numbers)

Hand-computed against sourced rules:

| Gene | VAF | In WBC? | IdentifyCHIP (source) | FilterCHIP kept? (documented contract) |
|------|-----|---------|----------------------|----------------------------------------|
| DNMT3A | 0.05 | no | CHIP (#1, ≥0.02) | removed (rule b) |
| DNMT3A | 0.01 | no | not CHIP (<0.02) | kept |
| EGFR | 0.30 | no | not CHIP (not a driver) | kept |
| EGFR | 0.30 | yes | n/a | removed (rule a, #3) |
| TET2 | 0.02 | — | CHIP (boundary, #1) | removed (rule b) |

All match the implementation.

### Findings / divergences (Stage A)

1. **DEFECT (doc, FIXED) — internal inconsistency in Evidence worked-example.** Evidence `Dataset: Worked classification cases` row 5 (TP53 0.40, not in WBC) listed FilterCHIP = **"kept (not in WBC ⇒ tumor)"**, which directly contradicts row 1 (DNMT3A 0.05, not in WBC → "removed") and the algorithm's own documented contract (`Clonal_Hematopoiesis_Filtering.md` §2.2: remove when `inWBC OR meetsChipHeuristic`). The same scenario (CH gene, VAF ≥ τ, absent from WBC) cannot be both removed and kept. Corrected row 5 to "removed (gene+VAF heuristic, rule b)" and added an origin caveat note explaining the strict-matched-WBC interpretation and the conservative heuristic.

2. **NOTE — heuristic over-removes vs strict matched-WBC.** Per Razavi 2019 / Arango-Argoty 2025 (#4, #6), gene identity is **not** a definitive origin test; the VAF–origin relationship is "unclear", and a CH-gene variant absent from WBC is not provably CH. `FilterCHIP` rule (b) deliberately removes such variants as a **conservative, labelled heuristic fallback** — this is honestly documented as a heuristic (not a definitive call), and callers can disable it via a custom/empty `chipGenes` panel so only matched-WBC subtraction applies. Acceptable as documented ⇒ PASS-WITH-NOTES rather than FAIL.

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:5834–6014`.
- **Formula realised correctly:** `IdentifyCHIPVariants` (5920) flags `IsCanonicalChipGene(gene) && Vaf >= minVaf` — exact, inclusive `>=` (INV-1/2). `FilterCHIP` (5966) builds a `HashSet` of WBC loci with `AltReads >= minWbcAltReads` and retains a variant iff `!inMatchedWbc && !meetsChipHeuristic` — matches §2.2 (INV-3/4/5). `IsCanonicalChipGene` (5886) is case-insensitive (`OrdinalIgnoreCase`) and treats null/empty gene as non-CHIP (INV-6).
- **Validation:** null `variants`/`whiteBloodCellVariants` → `ArgumentNullException`; `minVaf ∉ (0,1]` and `minWbcAltReads < 1` → `ArgumentOutOfRangeException` (all present in code).
- **Variant/delegate consistency:** `IdentifyCHIPVariants` and `FilterCHIP` both route gene membership through `IsCanonicalChipGene`; default `minWbcAltReads = DefaultMrdMinSupportingReads = 1` (consistent with the MRD convention, Wan 2020).
- **Numerical robustness:** locus key is exact tuple; VAF comparison is a plain double `>=` with no precision-sensitive arithmetic.

### Cross-verification recomputed vs code

Ran the full unfiltered suite: all CHIP cases pass; hand-traced table above matches code output exactly (incl. TET2 @0.02 boundary and WBC alt-read @cutoff inclusivity).

### Test-quality audit (HARD gate)

| Check | Result |
|-------|--------|
| Sourced expectations, not code echoes | PASS — assertions trace to Steensma/Genovese/Razavi values, not to current output |
| No green-washing (exact values, no widened tolerances/skips) | PASS — exact counts/identities; boundary `>=` tests would catch a `>` regression |
| Cover all public methods/overloads | PASS — both canonical methods + `IsCanonicalChipGene`; default & caller panel; custom minVaf |
| All Stage-A branches / edge / error cases | **Strengthened this session** |
| Honest green (FULL suite Failed: 0, warning-free changed files) | PASS — 6677 passed / 0 failed; CHIP files build with no warnings |

Defects fixed in the tests (Stage B):
- **M5** strengthened from a count-only assertion (`Has.Count.EqualTo(8)`) to `Is.EquivalentTo(genes)` — the count-only form would pass even if a gene were dropped and another duplicated.
- **Added `FilterCHIP_MinVafOutOfRange_Throws`** — the `minVaf ∉ (0,1]` error path on `FilterCHIP` was previously untested (only `IdentifyCHIPVariants` covered it).
- **Added `FilterCHIP_MinWbcAltReadsBelowOne_Throws`** — the documented `minWbcAltReads < 1 → ArgumentOutOfRangeException` contract (§3.3) had **no** test.
- **Added `FilterCHIP_WbcAltReadsExactlyAtCutoff_Removed`** — locks the inclusive (`>=`) WBC alt-read cutoff boundary; previously only the 0-read absence case (S2) was tested.

### Findings / defects (Stage B)

No code defect. The implementation faithfully realises the documented (and now self-consistent) description. Three test-quality gaps (one weak assertion, two missing documented error/boundary paths) were fixed in-session.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — biology/threshold/genes/matched-WBC rule all confirmed against primary sources; one documentation inconsistency (Evidence row 5) fixed; heuristic-over-removal noted as a documented, defensible design choice.
- **Stage B:** PASS-WITH-NOTES — code correct; test suite strengthened (4 changes) and full suite green.
- **End-state:** CLEAN — all defects fully fixed; `dotnet build` 0 errors, full suite Failed: 0 (6677 passed).
- **Test-quality gate:** PASS (after in-session fixes).

### Logged defects
- **FND:** Evidence ONCO-CHIP-001 worked-example row 5 contradicted the documented FilterCHIP contract and row 1 (doc inconsistency) — fixed.
