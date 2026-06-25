# Seqeron.Genomics — Current Limitations

**Library:** Seqeron.Genomics (mission-critical)   **Last reviewed:** 2026-06-25

This file lists **only what the library currently does NOT do** — genuine open limitations. Anything already
implemented (including opt-in alternatives, convention-parity modes, and resolved defects) is deliberately
**excluded**; those are recorded in the per-unit reports under `docs/Validation/reports/{UNIT}.md`, not here.

Every limitation below is `BY-DESIGN` and the unit is `✅ CLEAN` for its stated contract — these are honest scope
boundaries, never defects. They fall into three kinds:

- **Irreducible** — no algorithm or model can close it (physics / information theory).
- **Data-blocked** — needs a trained model / matrix / database that is gated, non-redistributable, or has never been measured / published. (Where the algorithm ships, the gated data is **caller-supplied** via a loader — the limitation is that nothing usable ships out-of-the-box.)
- **Scope** — a deliberate out-of-scope boundary; use the named reference tool, or supply the input yourself.

---

## 1. Irreducible (cannot be closed by any implementation)

| Unit | Not done | Why irreducible |
|------|----------|-----------------|
| PARSE-FASTQ-001 | Auto-disambiguation of Phred+33 vs Phred+64 when a read carries only overlap-range (ASCII 64–73) qualities. | The encodings overlap; the input is genuinely ambiguous (information-theoretic). Callers must decode with an explicit offset. |
| RNA-STRUCT-001 | Recovering tertiary-stabilised knots (e.g. BWYV / PDB 437D) as the MFE structure. | Not representable by *any* nearest-neighbour thermodynamic model — an energy-model floor, not an algorithm gap. |
| RNA-STRUCT-001 | Detecting pseudoknot classes outside the csr-PK grammar (kissing hairpins / loop–loop, triple-crossing / chained, non-canonical bulged helices). | The loop–loop interaction energy has never been experimentally measured (Sperschneider 2011 → heuristic estimate only); a faithful detector would need an unsourced energy constant. |

## 2. Data-blocked (needs a trained model / matrix / database that is gated, non-redistributable, or never measured)

| Unit | Not done | Blocking data |
|------|----------|---------------|
| ONCO-MHC-001 | A **bundled** trained HLA coefficient matrix (the predictor matrix is caller-supplied), and the pan-allele NetMHCpan/MHCflurry neural model. | No redistributable, cross-verifiable trained HLA matrix exists in the open: BIMAS files are served only by a now-defunct CGI, the Parker 1994 table is paywalled, IEDB SMM matrices are non-commercial / no-redistribution. The pan-allele neural model is a separate trained network. |
| ONCO-IMMUNE-001 | The **CIBERSORT-LM22-identical** signature matrix specifically (the 547-gene × 22-cell-type LM22 is caller-supplied; exact-CIBERSORT parity is not claimed). | LM22 is distributed by Stanford under a non-commercial **no-redistribution** licence and gated behind registration. A permissive signature matrix now ships out of the box: the **ABIS-Seq** matrix (Monaco et al., 2019, *Cell Reports*, CC BY 4.0; 1296 genes × 17 immune cell types) is bundled via `LoadBundledAbisSignatureMatrix()`, so immune deconvolution (the ν-SVR engine) works out-of-the-box — only LM22 itself remains caller-supplied. |
| META-BIN-001 | The **per-lineage-specific** CheckM marker refinement + the **reference genome tree** for lineage placement (the gated `checkm_data`). | Domain-level completeness/contamination now works out of the box: the Pfam (CC0) subsets of the GTDB bac120 (Bacteria, 6) and ar122 (Archaea, 35) universal single-copy marker sets ship bundled, feeding the CheckM formula. Only the per-lineage refinement + reference tree (`checkm_data`, a large gated trained DB) and the TIGRFAM-defined markers (CC BY-SA 4.0, not redistributable here) remain caller-supplied. |
| CHROM-CENT-001 | Suprachromosomal-family / specific α-satellite family (J1/J2/W/…) assignment for a detected HOR. | Naming the family needs curated chromosome-specific reference HOR (consensus) libraries — external curated data not embedded. |
| DISORDER-REGION-001 | A calibrated per-residue / per-region disorder **confidence** value. | No disorder predictor publishes a calibrated confidence standard (the region boundaries themselves follow the validated TOP-IDP threshold). |

## 3. Scope boundaries (deliberate — use the named reference tool, or supply the input)

| Unit | Not done | Use instead |
|------|----------|-------------|
| PROBE-DESIGN-001 | The quantitative **MGB (minor-groove binder) ΔTm**, and dual-quencher labelling. | The MGB ΔTm model (Kutyavin 2000 / MGB-Eclipse) is empirical/proprietary with no published closed form — use a chemistry-specific tool; dual-quencher has no Tm impact. |
| MIRNA-TARGET-001 | The context++ feature `PCT` (multi-species conservation); `TA_3UTR` / `SPS` / `Len_ORF` / `ORF8m` are caller-supplied. The score is a partial context++. | The TargetScan reference pipeline (a multi-species alignment for PCT, a transcriptome for TA). |
| MIRNA-PRECURSOR-001 | A **trained** natural-vs-background precursor classifier (read-stacking probabilistic model). | miRDeep2 / a trained miRNA caller. |
| PARSE-EMBL-001 | Fetching the **sequence** of a remote entry referenced in a location (the reference is parsed; only local parts are extracted). | Retrieve the remote accession from the source database. |

---

## How to read this

- This file lists capability/scope limitations only; pure **performance** optimisations (e.g. a genome-scale seeded index where an exhaustive scan already returns the same result) are not limitations and are not listed.
- For **research / pipeline** use, every row is a normal scope boundary of a from-first-principles library — the algorithm is validated correct for its stated contract.
- For **clinical / decision-grade** use, the §2 oncology rows mark layers that require an external validated predictor and clinical sign-off — the library computes the rule, not the trained model behind it.
- **Irreducible** rows can never be closed; **data-blocked** rows reopen only when the gated / unmeasured data becomes available (several already accept it via a caller-supplied loader); **scope** rows point to the reference tool, or accept the input from the caller.
- Each row traces to its per-unit validation report under `docs/Validation/reports/` (all `✅ CLEAN`).
