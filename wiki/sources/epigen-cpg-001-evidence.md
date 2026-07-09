---
type: source
title: "Evidence: EPIGEN-CPG-001 (CpG site & CpG-island detection)"
tags: [validation, epigenetics]
doc_path: docs/Evidence/EPIGEN-CPG-001-Evidence.md
sources:
  - docs/Evidence/EPIGEN-CPG-001-Evidence.md
source_commit: c066419cf49d3dfb54307f5f20853e66d89c0dbd
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: EPIGEN-CPG-001

The validation-evidence artifact for test unit **EPIGEN-CPG-001** — **CpG site detection**, the
canonical **CpG observed/expected ratio**, and **CpG-island detection** by the Gardiner-Garden &
Frommer length/GC%/O-E criteria. This is the **fourth ingested unit of the Epigenetics family** and
one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. The algorithm is synthesized in its own concept, [[cpg-island-detection]];
[[test-unit-registry]] tracks the unit. Sibling of [[bisulfite-methylation-calling]] — that caller
reuses this unit's `FindCpGSites` to enumerate the reference CpGs it scores.

## What this file records

- **Online sources (all four mutually consistent, no contradictions):**
  - **Gardiner-Garden M & Frommer M (1987)** "CpG islands in vertebrate genomes", *J Mol Biol*
    196(2):261–282 (authority rank 1, primary paper) — the canonical CpG-island criteria (**≥200 bp,
    GC% >50%, CpG O/E >0.6**) and the O/E formula `O/E = CpG_count / ((C_count × G_count) / L)` used by
    the UCSC Genome Browser and most tools; boundaries found by a **sliding-window** scan.
  - **Takai D & Jones PA (2002)** *PNAS* 99(6):3740–5 (rank 1) — **stricter** criteria (**≥500 bp, GC%
    >55%, O/E >0.65**) that reduce Alu-repeat false positives; validated on human chr 21/22. Their
    explicit "length ≥ 200 bp, ObsCpG/ExpCpG ≥ 0.6, %GC ≥ 50%" citation confirms the **≥ (inclusive)**
    threshold operators.
  - **Saxonov S, Berg P & Brutlag DL (2006)** *PNAS* 103(5):1412–1417 (rank 1) — the **alternative
    expected** formula `Expected = ((C + G)/2)² / L`; genome-wide analysis split promoters into two
    CpG-content classes.
  - **Wikipedia — "CpG site"** (rank 4, cites the three primaries) — CpG = 5'—C—phosphate—G—3'
    (cytosine then guanine, **not** GpC); islands typically 300–3,000 bp in mammals.

- **Documented corner cases / failure modes:** no CpG in sequence → O/E = 0; **zero C or G → expected
  = 0 → O/E returns 0 (division-by-zero guard)**; sequence shorter than the ≥200 bp minimum cannot be
  an island; **GpC ≠ CpG** (only 5'→3' C-then-G qualifies); mixed-case input must be handled
  (uppercase-normalize); adjacent CpGs `CGCG` = two distinct sites at 0 and 2 (each a separate 2-nt
  window, not overlapping); a length-1 input yields 0 sites.

- **Datasets (documented oracles, hand-derived from the CpG definition):**
  - `CGCGCGCGCGCGCGCGCGCG` (20 bp) — CpG count 10 at even positions, C=G=10, Expected =
    (10·10)/20 = 5.0, **O/E = 10/5.0 = 2.0**, GC 100%.
  - `AATTAATTAATTAATTAATT` (20 bp) — CpG count 0, **O/E = 0.0**.
  - `ACGTCGACG` (9 bp) — CpG at 1/4/7, count 3, C=G=3, Expected = 1.0, **O/E = 3.0**.
  - `ACGT` (4 bp) — CpG at 1, count 1, C=G=1, Expected = 0.25, **O/E = 4.0**.
  - **Island**: 400 bp of `CGCG` repeats → **is a CpG island** (length ≥200, GC 100% >50%, O/E 2.0
    >0.6).

- **Test-coverage traceability:** 24 tests in `EpigeneticsAnalyzer_CpGDetection_Tests.cs` — CpG =
  C-then-G in 5'→3' (M1/M2/M3/M6/M7/M18), Gardiner-Garden O/E (M9/M11/M12), island ≥200 bp/GC≥50%/O-E≥0.6
  (M15/M16/M17/S4), null/empty (M4a/M4b/M13/C2), case-insensitivity (M5/S2), edge cases single-char/no-G/minimal-CG
  (M8/M14/S1/C1), GpC≠CpG (M18), AT-only→O/E=0 (M10).

## Deviations and assumptions

- **ASSUMPTIONS: none** — the artifact records "None". All behaviour is formally defined by
  Gardiner-Garden & Frommer (1987) and confirmed by Takai & Jones (2002) and Wikipedia's cited
  primaries; the **≥ (inclusive)** threshold operators are explicitly confirmed by Takai & Jones's
  citation of the 1987 criteria.
- **Implementation note (from the algorithm doc, not the Evidence file):** the repository realises the
  **default** Gardiner-Garden thresholds (≥200 bp / GC ≥0.5 / O-E ≥0.6) as method defaults; the
  stricter Takai & Jones and alternative Saxonov formulas are **not** offered as named presets — a
  caller supplies custom `minLength`/`minGc`/`minCpGRatio` arguments. Island tuples use a 0-based
  inclusive `Start` / exclusive `End`.

No source contradictions.
