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
| **HIC-LOOP-001** | HiCCUPS loop calling — heavy, highly specialized, low incremental value over TAD/compartment, and validation needs large real contact matrices. Gate 2. |
| **HIC-SCAFFOLD-001** | Hi-C scaffolding (3D-DNA/LACHESIS class) — an assembly *pipeline* over contact data, overlapping the assembly domain and needing its own infra. Gate 3. |
| **IDX-WAVELET-001** | A rank/select data-structure detail *inside* the FM-index, not a user-facing biological algorithm unit. Fold into IDX-FM-001 if that is built. Gate 3. |

> **2026-11-10 (owner):** MAP-SEEDEXTEND-001 and VARIANT-HAPLO-001 were removed from this table and reinstated as Tier 3 units of the neoantigen mRNA vaccine pipeline (section below). The original reasons (a product rather than an algorithm; a whole variant-calling subsystem needing a read/BAM data model) still describe the cost; the owner accepted it.
>
> Net at the time of dropping: **114 → 109 kept.** Registry, section counts, and Quick Reference in the checklist are updated to match.

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

## Neoantigen mRNA vaccine pipeline (owner decision 2026-11-10)

Source of the list: [docs/NEOANTIGEN_VACCINE_GAP_ANALYSIS.md](docs/NEOANTIGEN_VACCINE_GAP_ANALYSIS.md) (86 numbered items: 26 present, 12 partial, 48 missing). The owner decided on 2026-11-10 that
the implementation of all missing and partial items belongs in the plan. The review campaign (B01–B26) is extended to them as
new batches B27–B33 (`docs/Validation/review-2026-09/README.md` → "Batches B27–B33"). Tiers below are **proposed** (value ÷ effort for this pipeline) and
override the tier of an already listed unit for scheduling only; each unit still passes the Gates and the Stage A/B protocol.
"Gap items" are the numbers in the gap document. The references in the checklist rows are starting points written without
fetching the papers: Stage A re-sources and verifies every one of them (a row saying "Stage A to source" has no reference yet).
**Nothing here is implemented yet; no batch is started.** Clinical use stays out of
scope: all results are research-grade.

| Unit | Tier | Gap items | Depends on | Review batch | Status |
|------|------|-----------|------------|--------------|--------|
| IO-BAM-001 | 2 | 4,5 | — | B27 | new |
| IO-FAIDX-001 | 1 | 8 | — | B27 | new |
| QC-ADAPTER-001 | 1 | 3 | — | B27 | new |
| READ-DEDUP-001 | 2 | 9 | IO-BAM-001 | B27 | new |
| READ-BQSR-001 | 3 | 9 | IO-BAM-001, IO-FAIDX-001, DATA-RESOURCES-001 | B27 | new |
| IDX-FM-001 | 2 | 15 | IDX-SUFFIXARRAY-001, IDX-BWT-001 | B27 | existing, promoted |
| MAP-SEEDEXTEND-001 | 3 | 10 | IDX-FM-001, ALIGN-LOCAL-001, IO-BAM-001 | B27 | reinstated |
| MAP-RNA-001 | 3 | 11 | MAP-SEEDEXTEND-001, SPLICE-PREDICT-001 | B27 | new |
| MAP-QUANT-001 | 2 | 12 | MAP-RNA-001 | B27 | new |
| VARIANT-SOMATIC-BAM-001 | 3 | 17 | MAP-SEEDEXTEND-001, READ-DEDUP-001, IO-BAM-001 | B28 | new |
| VARIANT-HAPLO-001 | 3 | 17 | MAP-SEEDEXTEND-001, READ-DEDUP-001, IO-BAM-001 | B28 | reinstated |
| VARIANT-GLK-001 | 2 | 18 | IO-BAM-001 | B28 | existing, promoted |
| VARIANT-GERMLINE-001 | 3 | 18 | VARIANT-GLK-001, VARIANT-NORM-001 | B28 | new |
| VARIANT-NORM-001 | 2 | 18,19 | — | B28 | existing, promoted |
| PHASE-READ-001 | 2 | 20,37 | IO-BAM-001 | B28 | existing, promoted |
| ONCO-ARTIFACT-BAM-001 | 2 | 19 | IO-BAM-001, DATA-RESOURCES-001 | B28 | new |
| FUSION-READS-001 | 2 | 24 | MAP-RNA-001 | B28 | new |
| QC-IDENTITY-001 | 2 | 75 | IO-BAM-001, PARSE-VCF-001 | B28 | new |
| HLA-DB-001 | 1 | 31 | DATA-RESOURCES-001 | B29 | new |
| HLA-TYPE-I-001 | 3 | 27 | MAP-SEEDEXTEND-001, HLA-DB-001 | B29 | new |
| HLA-TYPE-II-001 | 3 | 28 | MAP-SEEDEXTEND-001, HLA-DB-001 | B29 | new |
| ONCO-LOHHLA-001 | 2 | 30 | HLA-TYPE-I-001, IO-BAM-001 | B29 | existing, promoted |
| MHC-WEIGHTS-001 | 1 | 44 | — | B29 | new |
| DATA-RESOURCES-001 | 2 | 83,84,85 | — | B29 | new |
| NEO-PROTEIN-001 | 2 | 33,34,41 | VARIANT-ANNOT-001, PARSE-GFF-001 | B30 | new |
| NEO-FRAMESHIFT-001 | 2 | 35 | NEO-PROTEIN-001 | B30 | new |
| NEO-INFRAME-001 | 2 | 36 | NEO-PROTEIN-001 | B30 | new |
| NEO-PHASED-001 | 2 | 37 | NEO-PROTEIN-001, PHASE-READ-001 | B30 | new |
| NEO-FUSION-001 | 2 | 38 | NEO-PROTEIN-001, ONCO-FUSION-001 | B30 | new |
| NEO-SPLICE-001 | 3 | 39 | MAP-RNA-001, SPLICE-PREDICT-001 | B30 | new |
| NEO-NONCANON-001 | 3 | 40 | MAP-RNA-001 | B30 | new |
| ONCO-MHCII-001 | 2 | 43,45 | MHC-WEIGHTS-001 (class II pack) | B31 | existing, promoted |
| NEO-PROCESS-001 | 3 | 48 | — | B31 | new |
| NEO-STABILITY-001 | 3 | 49 | — | B31 | new |
| ONCO-IMMUNOGEN-001 | 2 | 50 | ONCO-NEO-001 | B31 | existing, promoted |
| ONCO-NEOFITNESS-001 | 2 | 50,53 | ONCO-IMMUNOGEN-001 | B31 | existing, promoted |
| NEO-TCR-001 | 3 | 51 | — | B31 | new |
| NEO-SELF-001 | 1 | 52 | DATA-RESOURCES-001, protein-alphabet multi-string index (PLAN.md section 8) | B31 | new |
| NEO-EXPR-001 | 2 | 54 | MAP-QUANT-001 or supplied TPM | B31 | new |
| NEO-RANK-001 | 2 | 55,56 | ONCO-NEO-001 and every NEO-* above | B31 | new |
| NEO-EXTERNAL-001 | 2 | 45,47 | ONCO-NEO-001, ONCO-MHCII-001, MHC-WEIGHTS-001 | B31 | new |
| MRNA-EPITOPE-ORDER-001 | 2 | 57 | NEO-RANK-001, NEO-PROCESS-001 | B32 | new |
| MRNA-LINKER-001 | 1 | 58 | — | B32 | new |
| MRNA-SIGNAL-001 | 1 | 59,60 | PROTMOTIF-SP-001 | B32 | new |
| MRNA-CODON-STRUCT-001 | 3 | 64,65 | CODON-OPT-001, RNA-MFE-001 | B32 | new |
| MRNA-UTR-001 | 2 | 66 | — | B32 | new |
| MRNA-TAIL-CAP-001 | 1 | 67 | — | B32 | new |
| MRNA-MODNUC-001 | 3 | 68 | RNA-ENERGY-001 | B32 | new |
| MRNA-IVT-001 | 1 | 69 | RESTR-FIND-001 | B32 | new |
| MRNA-IMMUNOSTIM-001 | 2 | 71 | — | B32 | new |
| MRNA-STABILITY-001 | 3 | 72 | CODON-OPT-001, RNA-PARTITION-001 | B32 | new |
| MRNA-ASSEMBLE-001 | 2 | 74 + integration | every MRNA-* unit | B32 | new |
| TCR-ASSEMBLE-001 | 3 | 79 | MAP-SEEDEXTEND-001 | B33 | new |
| REPORT-VACCINE-001 | 2 | 76 | NEO-RANK-001, MRNA-ASSEMBLE-001 | B33 | new |
| BENCH-NEO-001 | 2 | 78 | NEO-RANK-001 | B33 | new |
| BENCH-MAP-001 | 2 | — (validation against a reference) | MAP-SEEDEXTEND-001, VARIANT-GERMLINE-001 | B33 | new |
| BENCH-SOMATIC-001 | 2 | — (validation against a reference) | VARIANT-SOMATIC-BAM-001, VARIANT-HAPLO-001 | B33 | new |
| BENCH-HLA-001 | 2 | — (validation against a reference) | HLA-TYPE-I-001, HLA-TYPE-II-001 | B33 | new |
| INFRA-STREAM-001 | 2 | 86 | IO-BAM-001 | B33 | new |
| MCP-ONCO-001 | 1 | — (new server Seqeron.Mcp.Oncology) | — | B34 | new (MCP) |
| MCP-GAP-001 | 1 | — (existing pipeline-relevant capabilities without an MCP tool) | — | B34 | new (MCP) |
| MCP-WORKSPACE-001 | 2 | — (data by reference, long jobs, provenance) | — | B34 | new (MCP) |
| MCP-VACCINE-FLOW-001 | 2 | — (end-to-end route for an LLM (skill + MCP prompts)) | MCP-ONCO-001, MCP-WORKSPACE-001, B27–B33 tools | B34 | new (MCP) |

63 units: 55 new in the checklist (of which 2 reinstated from Dropped) and 8 already listed as proposed and promoted
(IDX-FM-001 from Tier 3 to 2, the others stay). Items of the gap document already ✅ need no unit.

**Validation against references (added 2026-11-10 at the owner's request).** A port is not assumed equal to its reference; parity is
measured per stage. Three units join `BENCH-NEO-001` (TESLA) in B33: BENCH-MAP-001 (mapping and germline variants on the Genome in a
Bottle truth sets), BENCH-SOMATIC-001 (the SEQC2 HCC1395/HCC1395BL call set and the DREAM synthetic tumours) and BENCH-HLA-001 (samples
with independently typed HLA alleles). The data-set names and papers are from memory and are confirmed at Stage A. Rules, also in the
review README (rule 7 of B27–B33):
- **Pre-registered thresholds.** The numeric acceptance threshold of each stage (for example a lower bound for F1 against the reference
  tool, a concordance rate for HLA types, a tolerance on MFE against ViennaRNA) is proposed before any measurement and approved by the
  owner; it is fixed before the run and checked on data that was not used while writing the code. The numbers are the owner's to set.
- **Independent measurer.** Whoever measures is not whoever implemented.
- **Speed is a separate measurement.** Genome-scale mapping is compared with BWA-MEM2 on time and memory on a larger machine than the
  development one, as a hand-off package under the existing rule.
- **External tools stay a first-class path.** BAM and VCF produced by external tools (BWA, GATK and others) are accepted as input, like
  the NetMHCpan adapter, so the library's own mapper and caller are a validated alternative and never the only route.
- **Where parity is not reached** (closed models), the measured difference is written into `LIMITATIONS.md` as a number.

**NEO-EXTERNAL-001 (added 2026-11-10 at the owner's request).** NetMHCpan-4.1 (class I) and NetMHCIIpan-4.0 (class II) are the vendor's
(DTU Health Tech) trained models and cannot be embedded; the library only talks to a copy the **caller** has installed under the
vendor's licence. The unit has three parts: (1) an input writer (peptides and alleles in the vendor format) and an output parser for a
**local executable**, with the version read from the output and a caller attestation that the licence permits the intended use (the
library does not try to judge licence terms); (2) a funnel module: a permissive pre-filter on the open models (MHCflurry port
affinity), a threshold calibrated on a set where both answers are known, a comparison on the **percentile rank per allele** and not on
IC50 across models, and a union of several models rather than one; (3) no web-service path by default, because the peptides come from a
patient's tumour and must not leave the machine. Scoring everything with the vendor tool is also cheap (an estimate: about 70 000
predictions for 300 missense mutations over 6 alleles), so the funnel is a cross-check and a licence economy, not a speed need. Its MCP
tools ship with the unit like every other unit's (B34 rule).

### MCP exposure of the pipeline (owner decision 2026-11-10, second step)

The owner asked that an LLM can run the whole pipeline through MCP. Today it cannot even run the existing oncology code: the
library has 483 tools on 11 servers (`docs/mcp/MCP_STATUS.md`), none of them oncological (`OncologyAnalyzer`, `ImmuneAnalyzer` and
`MhcflurryAffinityPredictor` are reachable only from C#, as the `seqeron-oncology` skill says). Four cross-cutting units close that,
in the review batch **B34**; in addition **every unit of B27–B33 ships its own MCP tools as part of its Definition of Done**
(gold-standard binding, Schema+Binding test, `{tool}.md` and `{tool}.mcp.json`, exactly as in `MCP_STATUS.md`), so no separate wrapper
unit per algorithm exists.

- **MCP-ONCO-001.** A new server `Seqeron.Mcp.Oncology` for the code that already exists. The source has about 230 `public static`
  members in 15 files (constants and overloads included, so the tool count will be lower; Stage A decides what is a tool). Needs no
  new algorithm, so it does not wait for B27–B33. The guarded units (ONCO-MHC-001, ONCO-IMMUNE-001) must reach the model as a
  structured limitation error, never as a wrong number; the MHCflurry weight pack is passed as a workspace handle.
- **MCP-GAP-001.** An audit and fill-in: for every item of the gap document marked present or partial, a tool exists or the reason is
  recorded. Found so far by name search only (to be confirmed at Stage A): no tool for the RNA partition function and base-pair
  probabilities, for the SAM parser, or for `ReportGenerator`.
- **MCP-WORKSPACE-001.** The model must not carry sequencing data. Tools take and return **handles** to files under a configured
  root (path, size, content hash), long steps run as jobs (start, poll, cancel, progress), every call appends a provenance record
  (tool, version, hash of inputs, parameters), a handle cannot leave the root, and the servers stay local over stdio so patient data
  does not leave the machine. Other units' tools depend on it.
- **MCP-VACCINE-FLOW-001.** The route for the LLM: an agent skill `seqeron-vaccine` plus MCP prompts that chain the tools, with
  mandatory human confirmation points (HLA type, the final neoantigen list, the final construct), the stop rule at guarded units, and
  the research-grade, not-for-clinical-use statement on every result. Also updates the `seqeron-oncology` and `seqeron-discovery`
  skills, which now say there is no oncology MCP server.

MCP does not remove the missing algorithms: until B27–B30 land, an LLM can drive only the middle of the chain through MCP.

---

## Summary

| Tier | Count | Meaning |
|------|-------|---------|
| Dropped | 3 | Out of mission — do not implement. (2 reinstated 2026-11-10, see the vaccine section.) |
| Tier 1 | 35 | Quick wins — high value, low effort, easy validation. |
| Tier 2 | 60 | High value, moderate effort — the main body of work. |
| Tier 3 | 9 | High value but heavy/uncertain — defer until a concrete need. |
| Vaccine pipeline | 63 | Owner decision 2026-11-10; tiers inside the section; 8 of these units are also counted in Tier 1–3 above. |

The highest-leverage near-term batch is **Tier 1's shared infrastructure** (STAT-*, CLIN-ROC,
VARIANT-NORM, TRANS-NORM) plus the **small oncology wins** (kataegis, WGD, ID83/DBS78 signatures that
reuse the existing NMF/NNLS engine) — cheap to build, broadly reused, and each validatable against a
published formula or reference tool.
