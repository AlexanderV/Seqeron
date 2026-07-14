---
type: index
title: "Ingestion backlog — docs/algorithms reconciliation + queued sources"
tags: [meta, coverage]
created: 2026-07-09
updated: 2026-07-14
---

# Ingestion backlog

Coverage reconciliation for `docs/algorithms/**` (kept **in scope** by the
[coverage exclude policy](SCHEMA.md#coverage-exclude-policy)) plus source batches
queued for ingestion. Generated during the 2026-07-09 lint pass; regenerate when
concept pages are added or algorithm docs change.

The **pending** rows are a real coverage gap that folds into the main per-algorithm
ingest campaign (the same campaign advancing the `docs/Evidence/**` files) — not a
separate effort. A pending algorithm doc is resolved when a concept page lists it in
`sources:`; at that point it moves to the covered table.

Status at generation: **124** algorithm docs covered-via-concept, **120** pending across 18 domains
(K-mer, Metagenomics and MolTools domains now fully covered;
Oncology/Focal_Amplification_Detection → [[focal-amplification-detection]] resolved 2026-07-14;
Oncology/Driver_Mutation_Detection → [[driver-gene-classification-20-20-rule]] resolved 2026-07-14;
Oncology/CtDNA_Analysis → [[ctdna-detection-and-tumor-fraction]] resolved 2026-07-14;
Oncology/Copy_Number_Alteration_Classification → [[copy-number-alteration-classification]] resolved 2026-07-14;
Oncology/Complex_Rearrangement_Classification → [[chromothripsis-inference]] resolved 2026-07-14;
Oncology/Clonal_Hematopoiesis_Filtering → [[clonal-hematopoiesis-cfdna-filtering]] resolved 2026-07-14;
Oncology/Clinical_Actionability_Assessment → [[clinical-actionability-oncokb-levels]] resolved 2026-07-14;
Oncology/Cancer_Variant_Annotation → [[cancer-variant-tier-classification-amp-asco-cap]] resolved 2026-07-14;
Oncology/Cancer_Cell_Fraction_Estimation → [[cancer-cell-fraction-clonal-clustering]] resolved 2026-07-14;
Oncology/Clonal_Subclonal_Classification → [[clonal-subclonal-classification-ccf-posterior]] resolved 2026-07-14; MolTools/Restriction_Site_Detection →
[[restriction-site-detection]] resolved 2026-07-14, closing the MolTools domain; K-mer_Search and PanGenome_Core_Accessory
resolved 2026-07-13; DNA_Dimer_Tm, DNA_Hairpin_Folding_Tm, DNA_Hairpin_Special_Loop_Bonus,
LNA_Adjusted_Nearest_Neighbor_Tm and NearestNeighbor_Salt_Corrected_Tm →
[[primer-dimer-thermodynamics-tm]] resolved 2026-07-13;
Hybridization_Probe_Design → [[hybridization-probe-design]] resolved 2026-07-13;
MolTools/Melting_Temperature → [[melting-temperature]] resolved 2026-07-13;
Primer3_Penalty_Objective → [[primer3-weighted-penalty-objective]] resolved 2026-07-13;
Primer_Design → [[primer-design]] resolved 2026-07-13;
Primer_Structure_Analysis → [[primer-structure-qc-screens]] resolved 2026-07-13;
Probe_Validation → [[probe-offtarget-specificity-scan]] resolved 2026-07-13).

## Covered via concept (done)

Each algorithm doc below is already synthesized by a concept page that lists it in
`sources:` (added at commit `9ce49ba`; staleness-clean, so no re-sync needed).

| Algorithm doc | Concept page |
| --- | --- |
| `docs/algorithms/Alignment/Alignment_Statistics.md` | [[alignment-statistics]] |
| `docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md` | [[global-alignment-needleman-wunsch]] |
| `docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md` | [[local-alignment-smith-waterman]] |
| `docs/algorithms/Alignment/Multiple_Sequence_Alignment.md` | [[multiple-sequence-alignment]] |
| `docs/algorithms/Alignment/Semi_Global_Alignment.md` | [[semi-global-alignment-fitting]] |
| `docs/algorithms/Analysis/Open_Reading_Frame_Detection.md` | [[open-reading-frame-detection]] |
| `docs/algorithms/Analysis/Sequence_Similarity.md` | [[kmer-jaccard-similarity]] |
| `docs/algorithms/Annotation/GFF3_IO.md` | [[gff3-io]] |
| `docs/algorithms/Annotation/Gene_Prediction.md` | [[prokaryotic-gene-prediction-rbs]] |
| `docs/algorithms/Annotation/ORF_Detection.md` | [[open-reading-frame-detection]] |
| `docs/algorithms/Annotation/Promoter_Detection.md` | [[promoter-detection]] |
| `docs/algorithms/Annotation/Relative_Synonymous_Codon_Usage.md` | [[relative-synonymous-codon-usage]] |
| `docs/algorithms/Annotation/Repetitive_Element_Detection.md` | [[repetitive-element-detection]] |
| `docs/algorithms/Assembly/Assembly_Statistics.md` | [[assembly-statistics]] |
| `docs/algorithms/Assembly/Consensus_Computation.md` | [[consensus-sequence]] |
| `docs/algorithms/Assembly/Coverage_Calculation.md` | [[coverage-depth-calculation]] |
| `docs/algorithms/Assembly/De_Bruijn_Graph_Assembly.md` | [[de-bruijn-graph-assembly]] |
| `docs/algorithms/Assembly/Error_Correction.md` | [[kmer-spectrum-error-correction]] |
| `docs/algorithms/Assembly/Overlap_Layout_Consensus.md` | [[overlap-layout-consensus-assembly]] |
| `docs/algorithms/Assembly/Quality_Trimming.md` | [[quality-trimming-running-sum]] |
| `docs/algorithms/Complexity/DUST_Score.md` | [[dust-low-complexity-score]] |
| `docs/algorithms/Complexity/K-mer_Entropy.md` | [[k-mer-statistics]] |
| `docs/algorithms/Complexity/Lempel_Ziv_Complexity.md` | [[sequence-complexity-compression-lempel-ziv]] |
| `docs/algorithms/Complexity/Windowed_Complexity.md` | [[windowed-sequence-complexity-profile]] |
| `docs/algorithms/Chromosome_Analysis/Aneuploidy_Detection.md` | [[aneuploidy-detection]] |
| `docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md` | [[centromere-analysis]] |
| `docs/algorithms/Chromosome_Analysis/Higher_Order_Repeat_Detection.md` | [[centromere-analysis]] |
| `docs/algorithms/Chromosome_Analysis/Karyotype_Analysis.md` | [[karyotype-analysis]] |
| `docs/algorithms/Chromosome_Analysis/Synteny_Analysis.md` | [[synteny-and-rearrangement-detection]] |
| `docs/algorithms/K-mer/Asynchronous_K-mer_Counting.md` | [[asynchronous-kmer-counting]] |
| `docs/algorithms/K-mer/Both_Strand_Kmer_Counting.md` | [[both-strand-kmer-counting]] |
| `docs/algorithms/K-mer/K-mer_Counting.md` | [[k-mer-counting]] |
| `docs/algorithms/K-mer/K-mer_Euclidean_Distance.md` | [[k-mer-euclidean-distance]] |
| `docs/algorithms/K-mer/K-mer_Frequency_Analysis.md` | [[k-mer-frequency-analysis]] |
| `docs/algorithms/K-mer/K-mer_Generation.md` | [[k-mer-generation]] |
| `docs/algorithms/K-mer/K-mer_Positions.md` | [[k-mer-positions]] |
| `docs/algorithms/K-mer/K-mer_Search.md` | [[k-mer-search]] |
| `docs/algorithms/K-mer/K-mer_Statistics.md` | [[k-mer-statistics]] |
| `docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md` | [[unique-and-mincount-kmers]] |
| `docs/algorithms/Chromosome_Analysis/Telomere_Analysis.md` | [[telomere-analysis]] |
| `docs/algorithms/Codon/Codon_Usage_Statistics.md` | [[codon-usage-statistics]] |
| `docs/algorithms/Codon/Effective_Number_of_Codons.md` | [[effective-number-of-codons]] |
| `docs/algorithms/Codon/Relative_Synonymous_Codon_Usage.md` | [[relative-synonymous-codon-usage]] |
| `docs/algorithms/Codon_Optimization/CAI_Calculation.md` | [[codon-adaptation-index]] |
| `docs/algorithms/Codon_Optimization/Codon_Usage_Analysis.md` | [[codon-usage-comparison]] |
| `docs/algorithms/Codon_Optimization/Rare_Codon_Detection.md` | [[rare-codon-analysis]] |
| `docs/algorithms/Codon_Optimization/Sequence_Optimization.md` | [[codon-optimization]] |
| `docs/algorithms/Comparative_Genomics/Average_Nucleotide_Identity.md` | [[average-nucleotide-identity]] |
| `docs/algorithms/Comparative_Genomics/Conserved_Gene_Clusters.md` | [[conserved-gene-clusters-common-intervals]] |
| `docs/algorithms/Comparative_Genomics/Dot_Plot_Generation.md` | [[dot-plot-word-match]] |
| `docs/algorithms/Comparative_Genomics/Genome_Comparison.md` | [[genome-comparison-core-dispensable]] |
| `docs/algorithms/Comparative_Genomics/Genome_Rearrangement_Detection.md` | [[genome-rearrangement-breakpoint-distance]] |
| `docs/algorithms/Comparative_Genomics/Ortholog_Identification.md` | [[ortholog-detection-reciprocal-best-hits]] |
| `docs/algorithms/Comparative_Genomics/Reciprocal_Best_Hits.md` | [[ortholog-detection-reciprocal-best-hits]] |
| `docs/algorithms/Comparative_Genomics/Reversal_Distance.md` | [[genome-rearrangement-breakpoint-distance]] |
| `docs/algorithms/Comparative_Genomics/Synteny_Block_Detection.md` | [[synteny-and-rearrangement-detection]] |
| `docs/algorithms/Epigenetics/Bisulfite_Sequencing_Analysis.md` | [[bisulfite-methylation-calling]] |
| `docs/algorithms/Epigenetics/Chromatin_State_Prediction.md` | [[chromatin-state-prediction]] |
| `docs/algorithms/Epigenetics/CpG_Site_Detection.md` | [[cpg-island-detection]] |
| `docs/algorithms/Epigenetics/Differentially_Methylated_Regions.md` | [[differentially-methylated-regions]] |
| `docs/algorithms/Epigenetics/Epigenetic_Age_Estimation.md` | [[epigenetic-age-horvath-clock]] |
| `docs/algorithms/Epigenetics/Methylation_Analysis.md` | [[methylation-context-classification]] |
| `docs/algorithms/Extended_Annotation/Coding_Potential_Calculation.md` | [[coding-potential-hexamer-score]] |
| `docs/algorithms/Extended_Assembly/Contig_Merging.md` | [[contig-merge-overlap-collapse]] |
| `docs/algorithms/Extended_Assembly/Scaffolding.md` | [[scaffolding]] |
| `docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md` | [[nucleotide-composition-skew]] |
| `docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md` | [[windowed-gc-profile-and-variance]] |
| `docs/algorithms/FileIO/BED_Parsing.md` | [[bed-format-parsing]] |
| `docs/algorithms/FileIO/EMBL_Parsing.md` | [[insdc-feature-location]] |
| `docs/algorithms/FileIO/FASTA_Parsing.md` | [[fasta-parsing]] |
| `docs/algorithms/FileIO/FASTQ_Parsing.md` | [[fastq-parsing]] |
| `docs/algorithms/FileIO/GFF_Parsing.md` | [[gff-parsing]] |
| `docs/algorithms/FileIO/GenBank_Parsing.md` | [[insdc-feature-location]] |
| `docs/algorithms/FileIO/VCF_Parsing.md` | [[vcf-parsing]] |
| `docs/algorithms/Genomic_Analysis/Tandem_Repeat_Detection.md` | [[repetitive-element-detection]] |
| `docs/algorithms/Metagenomics/Antibiotic_Resistance_Detection.md` | [[antibiotic-resistance-gene-detection]] |
| `docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md` | [[pan-genome-core-accessory-partition]] |
| `docs/algorithms/Metagenomics/Alpha_Diversity.md` | [[alpha-diversity]] |
| `docs/algorithms/Metagenomics/Beta_Diversity.md` | [[beta-diversity]] |
| `docs/algorithms/Metagenomics/Functional_Prediction.md` | [[functional-prediction]] |
| `docs/algorithms/Metagenomics/Genome_Binning.md` | [[metagenomic-binning]] |
| `docs/algorithms/Metagenomics/Pathway_Enrichment_ORA.md` | [[pathway-enrichment-ora]] |
| `docs/algorithms/Metagenomics/Significant_Taxa_Detection.md` | [[significant-taxa-detection]] |
| `docs/algorithms/Metagenomics/Taxonomic_Classification.md` | [[taxonomic-classification]] |
| `docs/algorithms/Metagenomics/Taxonomic_Profile.md` | [[taxonomic-profile]] |
| `docs/algorithms/MiRNA/MiRNA_Target_Pairing.md` | [[rna-base-pairing]] |
| `docs/algorithms/MiRNA/Pre_miRNA_Detection.md` | [[pre-mirna-hairpin-detection]] |
| `docs/algorithms/MiRNA/Seed_Sequence_Analysis.md` | [[seed-sequence-analysis]] |
| `docs/algorithms/MiRNA/Target_Site_Prediction.md` | [[mirna-target-site-prediction]] |
| `docs/algorithms/Motif_Analysis/Known_Motif_Search.md` | [[known-motif-search]] |
| `docs/algorithms/Motif_Discovery/Overrepresented_Kmer_Discovery.md` | [[overrepresented-kmer-discovery]] |
| `docs/algorithms/Motif_Discovery/Regulatory_Elements.md` | [[regulatory-element-detection]] |
| `docs/algorithms/Motif_Discovery/Shared_Motifs.md` | [[shared-motifs]] |
| `docs/algorithms/Pattern_Matching/Consensus_From_Alignment.md` | [[consensus-from-alignment]] |
| `docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Consensus.md` | [[iupac-degenerate-consensus]] |
| `docs/algorithms/ProteinPred/Low_Complexity_Region_Detection.md` | [[protein-low-complexity-seg]] |
| `docs/algorithms/Repeat_Analysis/Repeat_Detection.md` | [[longest-repeated-substring]] |
| `docs/algorithms/Sequence_Comparison/Common_Region_Detection.md` | [[longest-common-substring]] |
| `docs/algorithms/MolTools/Guide_RNA_Design.md` | [[crispr-guide-rna-design]] |
| `docs/algorithms/MolTools/Off_Target_Analysis.md` | [[crispr-guide-rna-design]] |
| `docs/algorithms/MolTools/PAM_Site_Detection.md` | [[crispr-guide-rna-design]] |
| `docs/algorithms/MolTools/Hybridization_Probe_Design.md` | [[hybridization-probe-design]] |
| `docs/algorithms/MolTools/DNA_Dimer_Tm.md` | [[primer-dimer-thermodynamics-tm]] |
| `docs/algorithms/MolTools/DNA_Hairpin_Folding_Tm.md` | [[primer-dimer-thermodynamics-tm]] |
| `docs/algorithms/MolTools/DNA_Hairpin_Special_Loop_Bonus.md` | [[primer-dimer-thermodynamics-tm]] |
| `docs/algorithms/MolTools/LNA_Adjusted_Nearest_Neighbor_Tm.md` | [[primer-dimer-thermodynamics-tm]] |
| `docs/algorithms/MolTools/NearestNeighbor_Salt_Corrected_Tm.md` | [[primer-dimer-thermodynamics-tm]] |
| `docs/algorithms/MolTools/Melting_Temperature.md` | [[melting-temperature]] |
| `docs/algorithms/MolTools/Primer3_Penalty_Objective.md` | [[primer3-weighted-penalty-objective]] |
| `docs/algorithms/MolTools/Primer_Design.md` | [[primer-design]] |
| `docs/algorithms/MolTools/Primer_Structure_Analysis.md` | [[primer-structure-qc-screens]] |
| `docs/algorithms/MolTools/Probe_Validation.md` | [[probe-offtarget-specificity-scan]] |
| `docs/algorithms/MolTools/Restriction_Digest_Simulation.md` | [[restriction-digest-simulation]] |
| `docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md` | [[restriction-enzyme-filtering]] |
| `docs/algorithms/MolTools/Restriction_Site_Detection.md` | [[restriction-site-detection]] |
| `docs/algorithms/Oncology/Allele_Specific_Copy_Number_Derivation.md` | [[allele-specific-copy-number-ascat]] |
| `docs/algorithms/Oncology/Cancer_Cell_Fraction_Estimation.md` | [[cancer-cell-fraction-clonal-clustering]] |
| `docs/algorithms/Oncology/Cancer_Variant_Annotation.md` | [[cancer-variant-tier-classification-amp-asco-cap]] |
| `docs/algorithms/Oncology/Clinical_Actionability_Assessment.md` | [[clinical-actionability-oncokb-levels]] |
| `docs/algorithms/Oncology/Clonal_Hematopoiesis_Filtering.md` | [[clonal-hematopoiesis-cfdna-filtering]] |
| `docs/algorithms/Oncology/Clonal_Subclonal_Classification.md` | [[clonal-subclonal-classification-ccf-posterior]] |
| `docs/algorithms/Oncology/Complex_Rearrangement_Classification.md` | [[chromothripsis-inference]] |
| `docs/algorithms/Oncology/Copy_Number_Alteration_Classification.md` | [[copy-number-alteration-classification]] |
| `docs/algorithms/Oncology/CtDNA_Analysis.md` | [[ctdna-detection-and-tumor-fraction]] |
| `docs/algorithms/Oncology/Driver_Mutation_Detection.md` | [[driver-gene-classification-20-20-rule]] |
| `docs/algorithms/Oncology/Focal_Amplification_Detection.md` | [[focal-amplification-detection]] |

## Pending (fold into the ingest campaign)

The per-domain pending tables (121 algorithm docs across 18 domains, no concept page yet) live in **[[backlog-pending]]** to keep this hub under the page-size cap. A pending row is resolved when a concept page lists the algorithm doc in `sources:`, at which point it moves to the *Covered via concept* table above.

## Queued source batches (approved 2026-07-09)

Approved for ingestion in the 2026-07-09 lint triage; pending `/wiki:ingest`.

### Testing methodology checklists (10) — `docs/checklists/`

- `docs/checklists/01_PROPERTY_BASED_TESTING.md`
- `docs/checklists/02_METAMORPHIC_TESTING.md`
- `docs/checklists/03_FUZZING.md`
- `docs/checklists/04_MUTATION_TESTING.md`
- `docs/checklists/05_SNAPSHOT_TESTING.md`
- `docs/checklists/06_ALGEBRAIC_TESTING.md`
- `docs/checklists/07_ARCHITECTURE_TESTING.md`
- `docs/checklists/08_DIFFERENTIAL_TESTING.md`
- `docs/checklists/09_COMBINATORIAL_TESTING.md`
- `docs/checklists/10_CHARACTERIZATION_TESTING.md`

### Validation governance ledgers (4) — `docs/Validation/`

- `docs/Validation/FINDINGS_REGISTER.md`
- `docs/Validation/LIMITATIONS.md`
- `docs/Validation/VALIDATION_LEDGER.md`
- `docs/Validation/VALIDATION_PROTOCOL.md`

### MCP top-level docs (3) — `docs/mcp/`

- `docs/mcp/MCP_STATUS.md`
- `docs/mcp/README.md`
- `docs/mcp/traceability.md`

## Notes

- `docs/algorithms/README.md` and `docs/algorithms/CANONICAL_MAP.md` are index/map
  docs, not algorithm units. `CANONICAL_MAP.md` is ingested as the source page
  [[canonical-algorithm-map]] (canonical-identity map: alias→canonical IDs, folder
  buckets, legacy baselines) — the identity counterpart to this coverage ledger.
  `README.md` remains index-only; ingest if a navigational need arises.
- The `docs/Evidence/**` campaign (175 of 213 remaining) is the primary driver: each
  Evidence ingest typically creates or extends the concept that also covers the
  matching algorithm doc, clearing a pending row here as a side effect.
