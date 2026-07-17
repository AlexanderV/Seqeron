---
type: source
title: "Validation verdict registry (CLEAN per-unit reports)"
tags: [validation, governance, registry]
sources:
  - docs/Validation/reports/
created: 2026-07-11
updated: 2026-07-15
---

# Validation verdict registry

One row per **CLEAN** per-unit validation report (no defect, no page needed). Reports that found a defect get their own `wiki/sources/<unit>-report.md` and a correction to the concept instead. See [[validation-ledger]] / [[validation-protocol]] for the governance model and [[validation-and-testing]] for the campaign.

| Unit | Concept | State | Stage A/B | Validated | Tests |
|------|---------|-------|-----------|-----------|-------|
| EPIGEN-BISULF-001 | [[bisulfite-methylation-calling]] | CLEAN | PASS / PASS | 2026-06-15 | 6539 |
| EPIGEN-CHROM-001 | [[chromatin-state-prediction]] | CLEAN | PASS / PASS | 2026-06-15 |  |
| EPIGEN-CPG-001 | [[cpg-island-detection]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| EPIGEN-DMR-001 | [[differentially-methylated-regions]] | CLEAN | PASS / PASS | 2026-06-15 | 6482 |
| EPIGEN-METHYL-001 | [[methylation-context-classification]] | CLEAN | PASS / PASS | 2026-06-15 | 6478 |
| GENOMIC-COMMON-001 | [[longest-common-substring]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 | 6573 |
| GENOMIC-MOTIFS-001 | [[known-motif-search]] | CLEAN | PASS / PASS | 2026-06-16 | 6573 |
| GENOMIC-ORF-001 | [[open-reading-frame-detection]] | CLEAN | PASS / PASS | 2026-06-16 | 6619 |
| GENOMIC-SIMILARITY-001 | [[kmer-jaccard-similarity]] | CLEAN | PASS / PASS | 2026-06-16 |  |
| GENOMIC-TANDEM-001 | [[repetitive-element-detection]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6619 |
| IMMUNE-NUSVR-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 | 50 |
| KMER-ASYNC-001 | [[asynchronous-kmer-counting]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6607 |
| KMER-BOTH-001 | [[both-strand-kmer-counting]] | CLEAN | PASS / PASS | 2026-06-16 | 6607 |
| KMER-COUNT-001 | [[k-mer-counting]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| KMER-DIST-001 | [[k-mer-euclidean-distance]] | CLEAN | PASS / PASS | 2026-06-15 | 6570 |
| KMER-FIND-001 | [[k-mer-search]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| KMER-FREQ-001 | [[k-mer-frequency-analysis]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| KMER-GENERATE-001 | [[k-mer-generation]] | CLEAN | PASS / PASS | 2026-06-16 | 6607 |
| KMER-POSITIONS-001 | [[k-mer-positions]] | CLEAN | PASS / PASS | 2026-06-16 | 6607 |
| KMER-STATS-001 | [[k-mer-statistics]] | CLEAN | PASS / PASS | 2026-06-16 | 6607 |
| KMER-UNIQUE-001 | [[unique-and-mincount-kmers]] | CLEAN | PASS / PASS | 2026-06-16 | 6607 |
| META-ALPHA-001 | [[alpha-diversity]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| META-BETA-001 | [[beta-diversity]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| META-BIN-001 | [[metagenomic-binning]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-24 | 4484 |
| META-CHECKM-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 | 100 |
| META-CLASS-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| META-FUNC-001 | [[functional-prediction]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| META-PATHWAY-001 | [[pathway-enrichment-ora]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| META-PROF-001 | [[taxonomic-profile]] | CLEAN | PASS / PASS | 2026-06-24 | 50 |
| META-TAXA-001 | [[significant-taxa-detection]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| META-TETRA-001 | [[metagenomic-binning]] | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| MHC-MATRIX-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 | 18756 |
| MHC-NN-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| MIRNA-CLASSIFY-001 | [[pre-mirna-hairpin-detection]] | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| MIRNA-CLEAVAGE-001 | [[pre-mirna-hairpin-detection]] | CLEAN | PASS / 🟡 PASS-WITH-NOTES | 2026-06-25 |  |
| MIRNA-CONTEXT-001 | [[mirna-target-site-prediction]] | CLEAN | PASS / ✅ PASS | 2026-06-25 | 816 |
| MIRNA-PCT-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| MIRNA-PRECURSOR-001 | [[pre-mirna-hairpin-detection]] | CLEAN | PASS / ✅ PASS | 2026-06-24 |  |
| MIRNA-SEED-001 | [[seed-sequence-analysis]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-24 |  |
| MOTIF-CONS-001 | [[consensus-from-alignment]] | CLEAN | PASS / PASS | 2026-06-15 | 6570 |
| MOTIF-DISCOVER-001 | [[overrepresented-kmer-discovery]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 |  |
| MOTIF-GENERATE-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6606 |
| MOTIF-REGULATORY-001 | [[regulatory-element-detection]] | CLEAN | PASS / PASS | 2026-06-16 |  |
| MOTIF-SHARED-001 | [[shared-motifs]] | CLEAN | PASS / PASS | 2026-06-16 | 6606 |
| ONCO-ACTION-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6694 |
| ONCO-ANNOT-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6628 |
| ONCO-ARTIFACT-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 |  |
| ONCO-ASCAT-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| ONCO-CCF-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 180 |
| ONCO-CHIP-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| ONCO-CLONAL-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| ONCO-CNA-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6649 |
| ONCO-CNA-002 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| ONCO-CNA-003 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| ONCO-CTDNA-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6667 |
| ONCO-DRIVER-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6624 |
| ONCO-EXPR-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6694 |
| ONCO-FUSION-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6644 |
| ONCO-FUSION-002 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6648 |
| ONCO-FUSION-003 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6649 |
| ONCO-HETERO-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6681 |
| ONCO-HLA-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6691 |
| ONCO-HRD-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| ONCO-IMMUNE-001 | ? | CLEAN | PASS / PASS | 2026-06-25 (fresh re-validation; supersedes 2026-06-24) |  |
| ONCO-LOH-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 |  |
| ONCO-MHC-001 | ? | CLEAN | PASS / PASS | 2026-06-26 (fresh re-validation after the 8–11 → 8–14 length-window fix, commit 66c24491; supersedes 2026-06-25) |  |
| ONCO-MRD-001 | ? | CLEAN | PASS / PASS | 2026-06-25 (fresh re-validation; supersedes 2026-06-24) |  |
| ONCO-MSI-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| ONCO-NEO-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6661 |
| ONCO-PHYLO-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6677 |
| ONCO-PLOIDY-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| ONCO-PURITY-001 | ? | CLEAN | PASS / PASS | 2026-06-24 | 30 |
| ONCO-SIG-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6637 |
| ONCO-SIG-002 | ? | CLEAN | PASS / PASS | 2026-06-24 | 999237 |
| ONCO-SIG-003 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| ONCO-SIG-004 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6644 |
| ONCO-SOMATIC-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| ONCO-SV-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6694 |
| ONCO-TMB-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| ONCO-VAF-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6621 |
| PANGEN-CLUSTER-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 875 |
| PANGEN-HEAP-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| PANGEN-MARKER-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6548 |
| PARSE-BED-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| PARSE-EMBL-001 | [[insdc-feature-location]] | CLEAN | PASS / ✅ PASS | 2026-06-26 (re-validated fresh; unit reset to ⬜ pending after the limitation-fix `d40318fb` ADDED remote-aware location-sequence assembly) |  |
| PARSE-FASTA-001 | ? | CLEAN | PASS / PASS  _(see "Re-validation 2026-06-25" below; supersedes the earlier PASS-WITH-NOTES once the DNA-only scope note was addressed by the opt-in alphabets)_ | 2026-06-25 (re-validated; original 2026-06-24) |  |
| PARSE-FASTQ-001 | ? | CLEAN | PASS / PASS | 2026-06-24; | 18880 |
| PARSE-GFF-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| PARSE-VCF-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| PAT-APPROX-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| PAT-APPROX-002 | ? | CLEAN | PASS / PASS | 2026-06-24 | 24 |
| PAT-APPROX-003 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 | 6536 |
| PAT-EXACT-001 | [[known-motif-search]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| PAT-IUPAC-001 | [[regulatory-element-detection]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| PAT-PWM-001 | [[regulatory-element-detection]] | CLEAN | PASS / PASS | 2026-06-24 | 11 |
| PHYLO-BOOT-001 | [[phylogenetic-bootstrap-support]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| PHYLO-COMP-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-24 |  |
| PHYLO-DIST-001 | ? | CLEAN | PASS / PASS | 2026-06-24 | 13674 |
| PHYLO-NEWICK-001 | [[mirna-target-site-prediction]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| PHYLO-STATS-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6561 |
| PHYLO-TREE-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-24 (independent re-confirmation) |  |
| POP-ANCESTRY-001 | ? | CLEAN | PASS / PASS | 2026-06-15 |  |
| POP-DIV-001 | ? | CLEAN | PASS / PASS | 2026-06-24 | 500000 |
| POP-FREQ-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| POP-FST-001 | ? | CLEAN | PASS / PASS | 2026-06-24 | 16 |
| POP-HW-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| POP-LD-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| POP-ROH-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6552 |
| POP-SELECT-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 | 6548 |
| PRIMER-DESIGN-001 | [[primer-design]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| PRIMER-DIMER-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| PRIMER-HAIRPIN-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| PRIMER-NNTM-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| PRIMER-STRUCT-001 | [[primer-structure-qc-screens]] | CLEAN | PASS / PASS | 2026-06-24 | 86 |
| PRIMER-TM-001 | [[primer-dimer-thermodynamics-tm]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| PROBE-DESIGN-001 | [[hybridization-probe-design]] / [[taqman-probe-design-rules]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| PROBE-EVALUE-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| PROBE-LNATM-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| PROBE-VALID-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| PROTMOTIF-CC-001 | [[coiled-coil-prediction]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6582 |
| PROTMOTIF-COMMON-001 | [[common-protein-motifs]] | CLEAN | PASS / PASS | 2026-06-16 |  |
| PROTMOTIF-DOMAIN-001 | [[protein-domain-and-signal-peptide-prediction]] | CLEAN | PASS / PASS (✅) — one test-data defect found and fixed; no code change | 2026-06-25 (fresh independent re-validation; see top section) | 18213 |
| PROTMOTIF-FIND-001 | [[protein-motif-pattern-search]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| PROTMOTIF-HMM-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| PROTMOTIF-LC-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6583 |
| PROTMOTIF-PATTERN-001 | [[protein-motif-pattern-search]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 |  |
| PROTMOTIF-PROSITE-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| PROTMOTIF-SP-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6579 |
| QUALITY-PHRED-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6510 |
| QUALITY-STATS-001 | [[fastq-quality-statistics]] | CLEAN | PASS / PASS | 2026-06-15 | 6510 |
| REP-APPROX-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 | 10 |
| REP-DIRECT-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| REP-INV-001 | ? | CLEAN | PASS / PASS | 2026-06-24 | 700 |
| REP-PALIN-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| REP-STR-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-25 (fresh re-validation; supersedes 2026-06-24) |  |
| REP-TANDEM-001 | [[repetitive-element-detection]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| RESTR-DIGEST-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-24 |  |
| RESTR-FILTER-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| RESTR-FIND-001 | [[restriction-enzyme-filtering]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| RNA-ACCESS-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| RNA-DOTBRACKET-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| RNA-ENERGY-001 | ? | CLEAN | PASS / PASS | 2026-06-24 | 4486 |
| RNA-HAIRPIN-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 |  |
| RNA-INVERT-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6593 |
| RNA-MFE-001 | [[pre-mirna-hairpin-detection]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6591 |
| RNA-PAIR-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6586 |
| RNA-PSEUDOKNOT-001 | [[rna-pseudoknot-detection]] | CLEAN | PASS / PASS | 2026-06-16 |  |
| RNA-STEMLOOP-001 | [[rna-stem-loop-enumeration]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| RNA-STRUCT-001 | [[pre-mirna-hairpin-detection]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| SEQ-ATSKEW-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6607 |
| SEQ-CODON-FREQ-001 | [[codon-usage-comparison]] | CLEAN | PASS / PASS | 2026-06-16 | 6618 |
| SEQ-COMP-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| SEQ-COMPLEX-001 | [[linguistic-complexity]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| SEQ-COMPLEX-DUST-001 | [[dust-low-complexity-score]] | CLEAN | PASS / PASS | 2026-06-16 | 6598 |
| SEQ-COMPLEX-KMER-001 | [[k-mer-statistics]] | CLEAN | PASS / PASS | 2026-06-16 |  |
| SEQ-COMPLEX-WINDOW-001 | [[windowed-sequence-complexity-profile]] | CLEAN | PASS / PASS | 2026-06-16 |  |
| SEQ-DINUC-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6523 |
| SEQ-ENTROPY-001 | [[shannon-entropy]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| SEQ-ENTROPY-PROFILE-001 | [[entropy-profile]] | CLEAN | PASS / PASS | 2026-06-16 | 6617 |
| SEQ-GC-001 | [[base-composition]] | CLEAN | PASS / ✅ PASS | 2026-06-24 |  |
| SEQ-GC-ANALYSIS-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| SEQ-GC-PROFILE-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-16 | 6617 |
| SEQ-GCSKEW-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| SEQ-HYDRO-001 | [[hydrophobicity-gravy-and-profile]] | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| SEQ-MW-001 | [[molecular-weight]] | CLEAN | PASS / PASS | 2026-06-15 | 6516 |
| SEQ-PI-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6517 |
| SEQ-REPLICATION-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6607 |
| SEQ-REVCOMP-001 | ? | CLEAN | PASS / PASS | 2026-06-24 |  |
| SEQ-RNACOMP-001 | ? | CLEAN | PASS / PASS | 2026-06-16 | 6573 |
| SEQ-SECSTRUCT-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6523 |
| SEQ-STATS-001 | ? | CLEAN | PASS / PASS | 2026-06-24 | 31 |
| SEQ-SUMMARY-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| SEQ-THERMO-001 | ? | CLEAN | PASS / PASS | 2026-06-15 |  |
| SEQ-TM-001 | ? | CLEAN | PASS / PASS | 2026-06-16 |  |
| SEQ-VALID-001 | [[sequence-validation]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| SPLICE-ACCEPTOR-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 (fresh re-validation) | 18 |
| SPLICE-DONOR-001 | ? | CLEAN | PASS / PASS | 2026-06-25 (fresh re-validation; supersedes 2026-06-24) |  |
| SPLICE-MAXENT3-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| SPLICE-MAXENT5-001 | ? | CLEAN | PASS / ✅ PASS | 2026-06-25 |  |
| SPLICE-PREDICT-001 | [[gene-structure-prediction-intron-exon]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| SV-BREAKPOINT-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| TRANS-CODON-001 | [[genetic-code-translation]] | CLEAN | PASS / PASS | 2026-06-24 |  |
| TRANS-DIFF-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 | 6501 |
| TRANS-EXPR-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 |  |
| TRANS-PROT-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-24 | 18211 |
| TRANS-SIXFRAME-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 | 6529 |
| VARIANT-ANNOT-001 | ? | CLEAN | PASS / PASS | 2026-06-15 |  |
| VARIANT-CALL-001 | ? | CLEAN | PASS / PASS | 2026-06-24 (re-validated after commit 6e900e92) | 34 |
| VARIANT-INDEL-001 | ? | CLEAN | PASS / PASS | 2026-06-15 | 6482 |
| VARIANT-SNP-001 | ? | CLEAN | PASS / PASS-WITH-NOTES | 2026-06-15 | 6482 |
