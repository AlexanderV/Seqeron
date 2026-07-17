---
type: index
title: "Ingestion backlog — pending algorithm docs (per-domain)"
tags: [meta, coverage]
created: 2026-07-09
updated: 2026-07-17
---

# Ingestion backlog — pending algorithm docs

The per-domain pending tables (18 algorithm docs across 5 domains) split out of **[[backlog]]** (which keeps the *Covered via concept* table, queued source batches, and notes). No concept page synthesizes these algorithm docs yet; each resolves when a concept lists it in `sources:`. See [[backlog]] for the full reconciliation model. (The K-mer domain is now fully covered — `K-mer_Search.md` → [[k-mer-search]], resolved 2026-07-13. The Metagenomics domain is now fully covered — `PanGenome_Core_Accessory.md` → [[pan-genome-core-accessory-partition]], resolved 2026-07-13. The MolTools domain is now fully covered — `Restriction_Site_Detection.md` → [[restriction-site-detection]], resolved 2026-07-14. The Oncology domain is now fully covered — `Variant_Allele_Frequency.md` → [[variant-allele-frequency-and-binomial-ci]], resolved 2026-07-15. The PanGenome domain is now fully covered — `Phylogenetic_Marker_Selection.md` → [[phylogenetic-marker-selection]], resolved 2026-07-15. The Pattern_Matching domain is now fully covered — `Suffix_Tree.md` → [[suffix-tree]], resolved 2026-07-15, ingesting the last pending Pattern_Matching doc. The Phylogenetics domain is now fully covered — `Tree_Statistics.md` → [[tree-statistics]], resolved 2026-07-15, ingesting the last pending Phylogenetics doc. The Population_Genetics domain is now fully covered — `Runs_Of_Homozygosity.md` → [[runs-of-homozygosity-inbreeding]], resolved 2026-07-16, ingesting the last pending Population_Genetics doc. The ProteinMotif domain is now fully covered — `Transmembrane_Helix_Prediction.md` → [[transmembrane-helix-prediction]], resolved 2026-07-16, ingesting the last pending ProteinMotif doc. The ProteinPred domain is now fully covered — `MoRF_Prediction.md` → [[morf-prediction-dip-in-disorder]], resolved 2026-07-16, ingesting the last pending ProteinPred doc. The Quality domain is now fully covered — `Quality_Statistics.md` → [[fastq-quality-statistics]], resolved 2026-07-16 (REUSE: enriched the existing Evidence-derived concept with the `QualityScoreAnalyzer` entry points, the multi-read `PerPositionMeanQuality` delegate variant, and complexity), ingesting the last pending Quality doc. The Repeat_Analysis domain is now fully covered — `Tandem_Repeat_Detection.md` → [[repetitive-element-detection]], resolved 2026-07-16 (REUSE: enriched the existing repeats/tandem family anchor with the `RepeatFinder.GetTandemRepeatSummary` → `FindMicrosatellites(1,6,minRepeats)` aggregation helper and its `TandemRepeatSummary` totals plus the 1–6 bp-only / per-class-stops-at-tetranucleotide scope caveats), ingesting the last pending Repeat_Analysis doc. The RnaStructure domain is now fully covered — `Turner_McCaskill_Partition_Function.md` → [[turner-mccaskill-partition-function]], resolved 2026-07-17 (NEW concept: the full Turner-2004 nearest-neighbour McCaskill partition function is a genuinely-distinct engine from the base-pair-counting [[rna-partition-function-mccaskill]] — different entry points `CalculateUnpairedProbabilities`/`CalculateRegionUnpairedProbability`, distinct outputs (`p_unpaired`, ensemble ΔG, region accessibility), and the TargetScan SA feature; it fills the "full Turner-parameter partition function" the simplified page listed as not implemented), ingesting the last pending RnaStructure doc. The Sequence_Composition domain is now fully covered — `Shannon_Entropy.md` → [[shannon-entropy]], resolved 2026-07-17 (NEW concept: the base per-symbol Shannon entropy `H = −Σ p·log₂ p` — a genuinely-distinct scalar member of the `SEQ-COMPLEX-*` complexity/entropy family with two entry points, `SequenceComplexity.CalculateShannonEntropy` (canonical DNA, only A/T/G/C) and `SequenceStatistics.CalculateShannonEntropy` (general-alphabet), kept distinct from the k-mer k-entropy of [[k-mer-statistics]] and the windowed profile), ingesting the last pending Sequence_Composition doc. The Splicing domain is now fully covered — `Gene_Structure_Prediction.md` → [[gene-structure-prediction-intron-exon]], resolved 2026-07-17 (REUSE: the composite intron/exon predictor spec reconciled onto the existing splicing-family composite page — enriched with the `PredictGeneStructure`/`PredictIntrons` signatures + defaults, the `minScore*0.8` site threshold + non-canonical opt-in, the `[acceptor-50, acceptor-18]` branch-point window at min score `0.4`, the `(donor+acceptor+branch)/3` vs `/2` scoring, greedy non-overlap selection, the `GeneStructure`/`Exon`/`Intron` output records, invariants INV-01…INV-04, and the SplicedSequence-vs-exon-record coordinate caveat), ingesting the last pending Splicing doc and closing the domain.)

### Statistics (5)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Statistics/Isoelectric_Point.md` | `isoelectric-point` |
| `docs/algorithms/Statistics/Melting_Temperature.md` | `melting-temperature` |
| `docs/algorithms/Statistics/Molecular_Weight_Calculation.md` | `molecular-weight-calculation` |
| `docs/algorithms/Statistics/Secondary_Structure_Prediction.md` | `secondary-structure-prediction` |
| `docs/algorithms/Statistics/Sequence_Summary.md` | `sequence-summary` |

### StructuralVar (3)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/StructuralVar/Breakpoint_Detection.md` | `breakpoint-detection` |
| `docs/algorithms/StructuralVar/Copy_Number_Variation.md` | `copy-number-variation` |
| `docs/algorithms/StructuralVar/SV_Detection.md` | `sv-detection` |

### Transcriptome (3)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Transcriptome/Alternative_Splicing.md` | `alternative-splicing` |
| `docs/algorithms/Transcriptome/Differential_Expression.md` | `differential-expression` |
| `docs/algorithms/Transcriptome/Expression_Quantification.md` | `expression-quantification` |

### Translation (3)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Translation/Codon_Translation.md` | `codon-translation` |
| `docs/algorithms/Translation/Protein_Translation.md` | `protein-translation` |
| `docs/algorithms/Translation/Six_Frame_Translation.md` | `six-frame-translation` |

### Variants (4)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Variants/Indel_Detection.md` | `indel-detection` |
| `docs/algorithms/Variants/SNP_Detection.md` | `snp-detection` |
| `docs/algorithms/Variants/Variant_Annotation.md` | `variant-annotation` |
| `docs/algorithms/Variants/Variant_Detection.md` | `variant-detection` |
