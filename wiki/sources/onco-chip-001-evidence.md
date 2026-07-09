---
type: source
title: "Evidence: ONCO-CHIP-001 (clonal hematopoiesis / CHIP filtering for cfDNA liquid biopsy)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-CHIP-001-Evidence.md
sources:
  - docs/Evidence/ONCO-CHIP-001-Evidence.md
source_commit: 90f75a142c015ef57f04ebf747b01f8b855634db
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-CHIP-001

The validation-evidence artifact for test unit **ONCO-CHIP-001** — **Clonal Hematopoiesis
(CHIP) filtering for cfDNA liquid biopsy**. The **sixth ingested unit of the Oncology
family** and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern. The distinct method is synthesized in
[[clonal-hematopoiesis-cfdna-filtering]]; [[test-unit-registry]] tracks the unit. It removes
the dominant liquid-biopsy false-positive class *before* the clinical-interpretation ONCO units
run, a QC sibling of [[sequencing-artifact-detection]].

## What this file records

- **Online sources (five primaries + one repo-convention citation; mutually consistent, no contradictions):**
  - **Steensma et al. (2015)** *Blood* 126(1):9–16, PMC4624443 (rank 1 — coined + formally defined
    CHIP): **VAF ≥ 2%** in peripheral blood, **somatic mutation in a driver gene recurrently mutated
    in hematologic malignancies**, and **absence of a diagnosed malignancy / MDS**. Driver genes
    (Fig 2A): DNMT3A, TET2, ASXL1, TP53, SF3B1, JAK2, PPM1D, BCORL1, GNAS, and others.
  - **Genovese et al. (2014)** *NEJM* 371(26):2477, PMC4290021 (rank 1 — 12,380-person whole-exome
    study): recurrent CH driver genes **DNMT3A (most, 190 mutations), ASXL1 (35), TET2 (31), PPM1D**,
    plus JAK2 V617F and SF3B1 K700E; somatic CH variants are **sub-clonal (VAF < 0.5)**.
  - **Razavi et al. (2019)** *Nat. Med.* 25:1928 (rank 1 — paywalled, facts via index + corroborating
    sources): **CH is the dominant cfDNA confounder — 81.6% of cfDNA mutations in controls and 53.2%
    in cancer patients** have CH features; the **matched-WBC design** (high-intensity sequencing of
    cfDNA AND matched white-blood-cell DNA) — a cfDNA variant also present in matched WBC is WBC/CH-
    derived, not tumor. This is the **definitive origin test**.
  - **Arango-Argoty et al. (2025)** *NPJ Precis Oncol* 9:147, PMC12092662 (rank 3 — operationalizes +
    cites Razavi): matched-WBC subtraction is the **gold-standard** origin call; the **top-3 CH genes
    DNMT3A/TET2/ASXL1** are removed; **VAF caveat** — "the exact relationship between VAF and variant
    origin remains unclear" ⇒ gene+VAF flags only *candidate* CHIP.
  - **Bolton et al. (2020)** *Nat. Genet.* 52(11):1219, PMC7891089 (rank 1 — 24,439-patient paired
    blood–tumour MSK-IMPACT study): **strict origin rule** — a WBC/CH call requires **WBC VAF ≥ 2%
    AND ≥ 10 supporting reads AND WBC VAF ≥ φ × tumour VAF** (φ = 2.0; **1.5 for a lymph-node biopsy
    site**); the fold ratio was chosen via leukocyte-contamination simulations (a sourced parameter,
    not invented).
  - **Wan et al. (2020)** *Sci. Transl. Med.* 12(548):eaaz8084 (repo-convention citation for the
    `IsVariantDetected`-style ≥1-alt-read matched-WBC presence test).

- **Three methods validated:**
  - **`IdentifyCHIPVariants`** — flag a variant as CHIP iff it is in a (default or caller-supplied)
    CHIP gene **and** VAF ≥ 0.02 (inclusive). Case-insensitive HGNC gene match.
  - **`FilterCHIP`** — remove cfDNA variants that are CH-derived. Rule (a) **matched-WBC subtraction**
    (present in matched WBC ⇒ removed, even a non-CHIP-gene variant); rule (b) **labelled gene+VAF
    heuristic fallback** (a CH-driver-gene variant at VAF ≥ τ removed even with no WBC evidence).
    Deliberately **conservative / over-removing**; callers pass an empty/custom `chipGenes` panel for
    the strict matched-WBC-only behaviour.
  - **`CallVariantOrigin`** — the strict Bolton-2020 origin call (Chip vs Tumor) from tumour VAF, WBC
    VAF, WBC read count, and the fold φ.

- **Documented corner cases / failure modes:** sub-2% driver mutation is **not** CHIP; driver mutation
  with a diagnosed malignancy is not CHIP (out of scope — assays non-diagnostic plasma); **without
  matched-WBC subtraction CH variants are mis-called as tumor** (the 81.6%/53.2% dominant FP source);
  a cfDNA variant **absent** from matched WBC is retained as a candidate tumour variant even in a CHIP
  gene.

- **Datasets (deterministic, hand-derived):**
  - **Canonical CHIP genes** {DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D}, VAF threshold
    0.02, top-3 DNMT3A/TET2/ASXL1.
  - **Worked classification cases** (gene / VAF / in-WBC → IdentifyCHIP / FilterCHIP): DNMT3A 0.05 →
    CHIP / removed; DNMT3A 0.01 → not CHIP / kept; EGFR 0.30 → not CHIP / kept; EGFR 0.30 in-WBC →
    removed (WBC-matched); TP53 0.40 not-in-WBC → CHIP candidate / removed (heuristic rule b).
  - **Strict matched-WBC origin (Bolton 2020, φ default 2.0 / 1.5 lymph node), 7 rows:** e.g. tumour
    0.10 / WBC 0.30 / 40 reads → **Chip**; tumour 0.40 / WBC absent → **Tumor**; WBC 9 reads → **Tumor**
    (9 < 10); tumour 0.01 / WBC 0.02 / 10 reads → **Chip** (all boundaries inclusive); WBC 0.015 →
    **Tumor** (< 2% floor); tumour 0.30 / WBC 0.40 → **Tumor** (0.40 < 2×0.30); tumour 0.25 / WBC 0.40
    at φ=1.5 → **Chip** (0.40 ≥ 0.375) but **Tumor** under default 2.0.

- **Coverage recommendations:** MUST — VAF≥0.02 driver → CHIP and sub-2% → not; non-CHIP-gene never
  flagged; caller-supplied panel honoured; `FilterCHIP` removes matched-WBC variants and retains
  absent-from-WBC ones; VAF exactly 0.02 boundary is CHIP; `CallVariantOrigin` Chip/Tumor per the
  ≥2%/≥10-reads/≥φ× rules incl. lymph-node 1.5×. SHOULD — case-insensitive gene match; null/empty
  handling. COULD — kept-variant order preservation.

## Deviations and assumptions

- **ASSUMPTION — canonical default gene set.** With no caller panel the default is the source-cited
  {DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D} (Steensma 2015 Fig 2A; Genovese 2014) — a
  *labeled* canonical set, not invented; overridable. Gene match case-insensitive.
- **ASSUMPTION — matched-WBC presence test.** "Present in matched WBC" = `IsVariantDetected`-style
  ≥1-alt-read evidence at the same locus (repo MRD convention, Wan 2020; matched-WBC subtraction,
  Razavi 2019); the universal alt-read cutoff is assay-specific and configurable.
- **In-file correction (2026-06-16):** the earlier "kept" outcome for the TP53-0.40-not-in-WBC row
  contradicted both the DNMT3A row and the documented `FilterCHIP` contract (rule b removes a CH-gene
  variant at VAF ≥ τ regardless of WBC evidence) and was corrected to "removed". The conservative
  heuristic over-removes vs the strict matched-WBC definition — a documented, intentional trade-off.

No source contradictions among the five primaries: the VAF-2% + driver-gene definition (Steensma /
Genovese), the matched-WBC origin gold standard (Razavi / Arango-Argoty), and the strict fold-ratio
origin rule (Bolton) address disjoint stages and reinforce one another.
