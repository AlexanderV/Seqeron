# Algorithms Roadmap — Prioritization of Proposed Units

**Date:** 2026-07-08
**Scope:** the 114 proposed (`☐ Not Started`) algorithm units added to
[ALGORITHMS_CHECKLIST_V2.md](ALGORITHMS_CHECKLIST_V2.md) (§58–78). This file sorts them by
**value ÷ effort**, applies a **mission-fit gate**, and records what is **dropped** and why.
It does **not** touch the 255 validated units.

---

## Mission-fit gate (from the README)

A unit is kept only if it clears every gate below; otherwise it is dropped as scope-creep.

1. **From-scratch, dependency-free** — implementable in pure .NET 10, no external binaries, no shipped
   trained neural nets, no mandatory reference database to function.
2. **Deterministic & validatable** — produces a number/structure that can be checked against primary
   literature or a reference tool (the library's core promise: *computed, never guessed*, with provenance).
3. **An algorithm, not a product subsystem** — not a whole read-aligner / GATK-class caller / assembly
   pipeline that would need its own I/O stack and data model.
4. **Cohesive** — extends one of the existing domains (alignment, annotation, variants, phylogenetics,
   pop/comparative genomics, metagenomics, transcriptomics, epigenetics, RNA structure, oncology,
   molecular design) rather than opening an unrelated field.

## Scoring rubric

- **Value (L/M/H)** — scientific impact × gap it fills × reuse by other modules.
- **Effort (L/M/H)** — implementation complexity **plus validation difficulty** (needing ground-truth
  data or a hard-to-reproduce reference raises effort).
- **Tier** = value ÷ effort. Tier 1 = do first; Tier 3 = high value but heavy/uncertain, do only with a
  concrete driver.

---

## Dropped (5 units) — do not implement

| Unit | Reason it fails the gate |
|------|--------------------------|
| **MAP-SEEDEXTEND-001** | A full short-read mapper (BWA/Bowtie class) — a product, not an analysis algorithm. Fast substring search is already served by the Ukkonen suffix tree. Gate 3. |
| **VARIANT-HAPLO-001** | GATK HaplotypeCaller-class local reassembly + PairHMM — a whole variant-calling subsystem needing a read/BAM data model the library isn't built around, and hard to validate deterministically. Gate 2 & 3. |
| **HIC-LOOP-001** | HiCCUPS loop calling — heavy, highly specialized, low incremental value over TAD/compartment, and validation needs large real contact matrices. Gate 2. |
| **HIC-SCAFFOLD-001** | Hi-C scaffolding (3D-DNA/LACHESIS class) — an assembly *pipeline* over contact data, overlapping the assembly domain and needing its own infra. Gate 3. |
| **IDX-WAVELET-001** | A rank/select data-structure detail *inside* the FM-index, not a user-facing biological algorithm unit. Fold into IDX-FM-001 if that is built. Gate 3. |

> Net: **114 → 109 kept.** Registry, section counts, and Quick Reference in the checklist are updated to match.

---

## Tier 1 — Start here (high value, low effort, easy to validate)

Small, self-contained, mostly closed-form or textbook DP; validate against a formula or reference tool.
Several are **shared infrastructure** many other modules will call.

| Unit | V | E | Why it's a quick win |
|------|---|---|----------------------|
| STAT-MWU-001 | H | L | Non-parametric two-group test reused across DE, abundance, QC. |
| STAT-KS-001 | H | L | Distribution comparison for fragment/quality/enrichment. |
| STAT-QVALUE-001 | H | L | More power than BH for genome-wide testing; tiny. |
| STAT-PERM-001 | H | L | General permutation/bootstrap backstop for any statistic. |
| STAT-BENJYEK-001 | M | L | Dependency-robust FDR + Holm; correlated-test correctness. |
| CLIN-ROC-001 | H | L | ROC/AUC + PR — evaluation metric every classifier reuses. |
| VARIANT-NORM-001 | H | L | Indel left-align/normalize — correctness prerequisite for merge/annotate. |
| PHYLO-DISTCORR-001 | H | L | TN93/F84/LogDet extend the existing distance engine; closed-form. |
| POP-NEUTRAL-001 | H | L | Fu & Li D / Fay & Wu H — formulas over the SFS. |
| POP-MK-001 | M | L | McDonald-Kreitman / HKA — 2×2 + neutrality index. |
| TRANS-NORM-001 | H | L | TMM / median-of-ratios — fixes a real normalization gap; small. |
| RNA-MEA-001 | H | L | MEA/centroid reuses the already-shipped McCaskill probabilities. |
| META-UNIFRAC-001 | H | L | Phylogenetic beta diversity; small over existing tree + tables. |
| META-MASH-001 | H | L | Reuses the existing MinHash sketch for genome distance/ANI. |
| ALIGN-MYERS-001 | H | L | Bit-vector edit distance — big speed-up, trivially validated. |
| ALIGN-HIRSCHBERG-001 | M | L | Linear-space global alignment; standard divide-and-conquer. |
| ALIGN-BANDED-001 | M | L | Banded DP for near-identical sequences. |
| ALIGN-SUBOPT-001 | M | L | Waterman-Eggert top-k local alignments. |
| PHYLO-PARSIMONY-001 | H | M | Fitch small-parsimony; foundational for tree search. |
| PHYLO-CONSENSUS-001 | M | L | Strict/majority-rule consensus over bootstrap trees. |
| PHYLO-ROOT-001 | M | L | Midpoint / outgroup rooting. |
| PROT-GOR-001 | M | L | GOR IV secondary structure — windowed propensities. |
| PROT-RSA-001 | M | L | Solvent accessibility — windowed prediction. |
| PROT-BEPITOPE-001 | M | L | B-cell epitope propensity scales. |
| ONCO-KATAEGIS-001 | H | L | Inter-mutation-distance segmentation; small, real oncology value. |
| ONCO-WGD-001 | H | L | Whole-genome-doubling rule from existing ASCAT ploidy. |
| ONCO-SIGID-001 | H | M | Indel (ID83) catalog reuses the existing NMF/NNLS engine. |
| ONCO-SIGDBS-001 | H | M | Doublet (DBS78) catalog reuses the same engine. |
| ONCO-MSIPRO-001 | M | M | Single-sample MSI from repeat-length distribution shift. |
| ASSEMBLY-UNITIG-001 | M | L | Non-branching-path compaction of the existing graph. |
| ASSEMBLY-SIMPLIFY-001 | H | M | Tip removal + bubble popping — core assembly-graph cleaning. |
| ASSEMBLY-KMERHIST-001 | H | M | Genome-size/heterozygosity from the k-mer spectrum. |
| ASSEMBLY-MISASSEMBLY-001 | H | M | NGA50 + misassembly break vs reference; validates against QUAST. |
| CRISPR-CFD-001 | H | L | CFD off-target matrix score — small upgrade to the guide designer. |
| CRISPR-EFF-001 | H | M | Doench Rule-Set-2 on-target efficiency ranking. |

## Tier 2 — High value, moderate effort (do next)

Worth building; needs more code or validation data than Tier 1. Grouped by domain.

- **Oncology signatures & HRD:** ONCO-SIGCN-001, ONCO-SIGSV-001, ONCO-SIGFIT-001, ONCO-HRDETECT-001, ONCO-CHORD-001.
- **Oncology instability & drivers:** ONCO-CHROMOTHRIPSIS-001, ONCO-TIMING-001, ONCO-DNDS-001, ONCO-GISTIC-001.
- **Oncology subclonal:** ONCO-PYCLONE-001, ONCO-SUBCLONALCN-001.
- **Oncology immuno:** ONCO-MHCII-001, ONCO-NEOFITNESS-001, ONCO-IMMUNOGEN-001 (implement the DAI/agretopicity/foreignness features; the DNN references are inspiration only, not a shipped model).
- **Oncology liquid biopsy:** ONCO-LOHHLA-001, ONCO-ICHOR-001, ONCO-FRAGMENT-001, ONCO-PHASED-001.
- **Variants & phasing:** VARIANT-GLK-001, VARIANT-DENOVO-001, PHASE-READ-001, PHASE-STAT-001.
- **Population genetics:** POP-KINSHIP-001, POP-PCA-001, POP-IBD-001, POP-NE-001.
- **Transcriptome:** TRANS-NBDE-001, TRANS-GSVA-001, TRANS-ESTIMATE-001, TRANS-BATCH-001, TRANS-WGCNA-001.
- **Clinical biostatistics (optional support module for oncology):** CLIN-KM-001, CLIN-COX-001, CLIN-CINDEX-001, CLIN-SUBTYPE-001.
- **Phylogenetics:** PHYLO-ME-001, PHYLO-ANCESTRAL-001, PHYLO-MODELSEL-001, PHYLO-ALRT-001, PHYLO-REARRANGE-001, PHYLO-LSFIT-001.
- **Assembly:** ASSEMBLY-STRINGGRAPH-001, ASSEMBLY-MULTIK-001, ASSEMBLY-MINIMIZER-001, ASSEMBLY-MINHASH-001, ASSEMBLY-POA-001, ASSEMBLY-GAPFILL-001, ASSEMBLY-BLOOMDBG-001.
- **3D genome (kept core):** HIC-NORM-001, HIC-TAD-001, HIC-COMPARTMENT-001.
- **Epigenomic signal:** EPIGEN-PEAK-001, EPIGEN-FOOTPRINT-001, EPIGEN-NUCLEO-001.
- **RNA structure:** RNA-SAMPLE-001, RNA-INTERACT-001.
- **Gene/ncRNA finding:** ANNOT-IMM-001, ANNOT-TRNA-001, ANNOT-CRISPRARRAY-001.
- **Metagenomics:** META-ASV-001, META-COABUND-001.
- **Text indexing (infra):** IDX-SUFFIXARRAY-001, IDX-LCP-001, IDX-BWT-001.
- **Molecular design:** CRISPR-PRIME-001, CLONE-ASSEMBLY-001.

## Tier 3 — High value but heavy / hard to validate (defer; build only with a concrete driver)

| Unit | Why it waits |
|------|--------------|
| PHYLO-ML-001 + PHYLO-SUBMODEL-001 | Full likelihood engine (matrix exponentiation + branch-length optimization). High value but a large, careful build; SUBMODEL is its prerequisite. |
| ASSEMBLY-REPEATGRAPH-001 | Flye-class repeat/A-Bruijn graph for long reads — heavy, needs realistic long-read test data. |
| ONCO-ECDNA-001 | AmpliconArchitect-class breakpoint-graph reconstruction — high scientific interest but complex and hard to validate deterministically. |
| ANNOT-INFERNAL-001 | Covariance-model CYK engine — powerful but a substantial SCFG build. |
| PROT-DCA-001 | Direct coupling analysis — MSA-scale maximum-entropy optimization; borderline for a genomics library. |
| RNA-SANKOFF-001 | Simultaneous fold-and-align — O(n⁴) even banded; niche. |
| IDX-FM-001 | FM-index backward search — useful infra but overlaps the existing suffix tree; build only if a mapping/search driver appears. |

---

## Summary

| Tier | Count | Meaning |
|------|-------|---------|
| Dropped | 5 | Out of mission — do not implement. |
| Tier 1 | 35 | Quick wins — high value, low effort, easy validation. |
| Tier 2 | 60 | High value, moderate effort — the main body of work. |
| Tier 3 | 9 | High value but heavy/uncertain — defer until a concrete need. |

The highest-leverage near-term batch is **Tier 1's shared infrastructure** (STAT-*, CLIN-ROC,
VARIANT-NORM, TRANS-NORM) plus the **small oncology wins** (kataegis, WGD, ID83/DBS78 signatures that
reuse the existing NMF/NNLS engine) — cheap to build, broadly reused, and each validatable against a
published formula or reference tool.
