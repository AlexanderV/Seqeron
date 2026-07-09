---
type: index
title: "Ingestion backlog — docs/algorithms reconciliation + queued sources"
tags: [meta, coverage]
created: 2026-07-09
updated: 2026-07-09
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

Status at generation: **56** algorithm docs covered-via-concept, **189** pending across 30 domains.

## Covered via concept (done)

Each algorithm doc below is already synthesized by a concept page that lists it in
`sources:` (added at commit `9ce49ba`; staleness-clean, so no re-sync needed).

| Algorithm doc | Concept page |
| --- | --- |
| `docs/algorithms/Alignment/Alignment_Statistics.md` | [[alignment-statistics]] |
| `docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md` | [[global-alignment-needleman-wunsch]] |
| `docs/algorithms/Alignment/Multiple_Sequence_Alignment.md` | [[multiple-sequence-alignment]] |
| `docs/algorithms/Alignment/Semi_Global_Alignment.md` | [[semi-global-alignment-fitting]] |
| `docs/algorithms/Analysis/Open_Reading_Frame_Detection.md` | [[open-reading-frame-detection]] |
| `docs/algorithms/Analysis/Sequence_Similarity.md` | [[kmer-jaccard-similarity]] |
| `docs/algorithms/Annotation/Relative_Synonymous_Codon_Usage.md` | [[relative-synonymous-codon-usage]] |
| `docs/algorithms/Annotation/Repetitive_Element_Detection.md` | [[repetitive-element-detection]] |
| `docs/algorithms/Assembly/Assembly_Statistics.md` | [[assembly-statistics]] |
| `docs/algorithms/Assembly/Consensus_Computation.md` | [[consensus-sequence]] |
| `docs/algorithms/Assembly/Coverage_Calculation.md` | [[coverage-depth-calculation]] |
| `docs/algorithms/Assembly/De_Bruijn_Graph_Assembly.md` | [[de-bruijn-graph-assembly]] |
| `docs/algorithms/Assembly/Error_Correction.md` | [[kmer-spectrum-error-correction]] |
| `docs/algorithms/Assembly/Overlap_Layout_Consensus.md` | [[overlap-layout-consensus-assembly]] |
| `docs/algorithms/Assembly/Quality_Trimming.md` | [[quality-trimming-running-sum]] |
| `docs/algorithms/Chromosome_Analysis/Aneuploidy_Detection.md` | [[aneuploidy-detection]] |
| `docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md` | [[centromere-analysis]] |
| `docs/algorithms/Chromosome_Analysis/Karyotype_Analysis.md` | [[karyotype-analysis]] |
| `docs/algorithms/Chromosome_Analysis/Synteny_Analysis.md` | [[synteny-and-rearrangement-detection]] |
| `docs/algorithms/K-mer/Asynchronous_K-mer_Counting.md` | [[asynchronous-kmer-counting]] |
| `docs/algorithms/K-mer/Both_Strand_Kmer_Counting.md` | [[both-strand-kmer-counting]] |
| `docs/algorithms/K-mer/K-mer_Euclidean_Distance.md` | [[k-mer-euclidean-distance]] |
| `docs/algorithms/K-mer/K-mer_Generation.md` | [[k-mer-generation]] |
| `docs/algorithms/K-mer/K-mer_Positions.md` | [[k-mer-positions]] |
| `docs/algorithms/K-mer/K-mer_Statistics.md` | [[k-mer-statistics]] |
| `docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md` | [[unique-and-mincount-kmers]] |
| `docs/algorithms/Chromosome_Analysis/Telomere_Analysis.md` | [[telomere-analysis]] |
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
| `docs/algorithms/Genomic_Analysis/Tandem_Repeat_Detection.md` | [[repetitive-element-detection]] |
| `docs/algorithms/Motif_Analysis/Known_Motif_Search.md` | [[known-motif-search]] |
| `docs/algorithms/ProteinPred/Low_Complexity_Region_Detection.md` | [[protein-low-complexity-seg]] |
| `docs/algorithms/Repeat_Analysis/Repeat_Detection.md` | [[longest-repeated-substring]] |
| `docs/algorithms/Sequence_Comparison/Common_Region_Detection.md` | [[longest-common-substring]] |

## Pending (fold into the ingest campaign)

No concept page synthesizes these yet. `Expected slug` is the anticipated concept
page name once ingested (subject to change if a shared anchor concept fits better).

### Alignment (1)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md` | `local-alignment-smith-waterman` |

### Annotation (4)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Annotation/GFF3_IO.md` | `gff3-io` |
| `docs/algorithms/Annotation/Gene_Prediction.md` | `gene-prediction` |
| `docs/algorithms/Annotation/ORF_Detection.md` | `orf-detection` |
| `docs/algorithms/Annotation/Promoter_Detection.md` | `promoter-detection` |

### Chromosome_Analysis (1)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Chromosome_Analysis/Higher_Order_Repeat_Detection.md` | `higher-order-repeat-detection` |

### Codon (1)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Codon/Codon_Usage_Statistics.md` | `codon-usage-statistics` |

### Complexity (4)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Complexity/DUST_Score.md` | `dust-score` |
| `docs/algorithms/Complexity/K-mer_Entropy.md` | `k-mer-entropy` |
| `docs/algorithms/Complexity/Lempel_Ziv_Complexity.md` | `lempel-ziv-complexity` |
| `docs/algorithms/Complexity/Windowed_Complexity.md` | `windowed-complexity` |

### Extended_GC_Skew_Analysis (2)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md` | `at-skew` |
| `docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md` | `comprehensive-gc-analysis` |

### FileIO (7)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/FileIO/BED_Parsing.md` | `bed-parsing` |
| `docs/algorithms/FileIO/EMBL_Parsing.md` | `embl-parsing` |
| `docs/algorithms/FileIO/FASTA_Parsing.md` | `fasta-parsing` |
| `docs/algorithms/FileIO/FASTQ_Parsing.md` | `fastq-parsing` |
| `docs/algorithms/FileIO/GFF_Parsing.md` | `gff-parsing` |
| `docs/algorithms/FileIO/GenBank_Parsing.md` | `genbank-parsing` |
| `docs/algorithms/FileIO/VCF_Parsing.md` | `vcf-parsing` |

### K-mer (3)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/K-mer/K-mer_Counting.md` | `k-mer-counting` |
| `docs/algorithms/K-mer/K-mer_Frequency_Analysis.md` | `k-mer-frequency-analysis` |
| `docs/algorithms/K-mer/K-mer_Search.md` | `k-mer-search` |

### Metagenomics (10)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Metagenomics/Alpha_Diversity.md` | `alpha-diversity` |
| `docs/algorithms/Metagenomics/Antibiotic_Resistance_Detection.md` | `antibiotic-resistance-detection` |
| `docs/algorithms/Metagenomics/Beta_Diversity.md` | `beta-diversity` |
| `docs/algorithms/Metagenomics/Functional_Prediction.md` | `functional-prediction` |
| `docs/algorithms/Metagenomics/Genome_Binning.md` | `genome-binning` |
| `docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md` | `pangenome-core-accessory` |
| `docs/algorithms/Metagenomics/Pathway_Enrichment_ORA.md` | `pathway-enrichment-ora` |
| `docs/algorithms/Metagenomics/Significant_Taxa_Detection.md` | `significant-taxa-detection` |
| `docs/algorithms/Metagenomics/Taxonomic_Classification.md` | `taxonomic-classification` |
| `docs/algorithms/Metagenomics/Taxonomic_Profile.md` | `taxonomic-profile` |

### MiRNA (4)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/MiRNA/MiRNA_Target_Pairing.md` | `mirna-target-pairing` |
| `docs/algorithms/MiRNA/Pre_miRNA_Detection.md` | `pre-mirna-detection` |
| `docs/algorithms/MiRNA/Seed_Sequence_Analysis.md` | `seed-sequence-analysis` |
| `docs/algorithms/MiRNA/Target_Site_Prediction.md` | `target-site-prediction` |

### MolTools (17)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/MolTools/DNA_Dimer_Tm.md` | `dna-dimer-tm` |
| `docs/algorithms/MolTools/DNA_Hairpin_Folding_Tm.md` | `dna-hairpin-folding-tm` |
| `docs/algorithms/MolTools/DNA_Hairpin_Special_Loop_Bonus.md` | `dna-hairpin-special-loop-bonus` |
| `docs/algorithms/MolTools/Guide_RNA_Design.md` | `guide-rna-design` |
| `docs/algorithms/MolTools/Hybridization_Probe_Design.md` | `hybridization-probe-design` |
| `docs/algorithms/MolTools/LNA_Adjusted_Nearest_Neighbor_Tm.md` | `lna-adjusted-nearest-neighbor-tm` |
| `docs/algorithms/MolTools/Melting_Temperature.md` | `melting-temperature` |
| `docs/algorithms/MolTools/NearestNeighbor_Salt_Corrected_Tm.md` | `nearestneighbor-salt-corrected-tm` |
| `docs/algorithms/MolTools/Off_Target_Analysis.md` | `off-target-analysis` |
| `docs/algorithms/MolTools/PAM_Site_Detection.md` | `pam-site-detection` |
| `docs/algorithms/MolTools/Primer3_Penalty_Objective.md` | `primer3-penalty-objective` |
| `docs/algorithms/MolTools/Primer_Design.md` | `primer-design` |
| `docs/algorithms/MolTools/Primer_Structure_Analysis.md` | `primer-structure-analysis` |
| `docs/algorithms/MolTools/Probe_Validation.md` | `probe-validation` |
| `docs/algorithms/MolTools/Restriction_Digest_Simulation.md` | `restriction-digest-simulation` |
| `docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md` | `restriction-enzyme-filtering` |
| `docs/algorithms/MolTools/Restriction_Site_Detection.md` | `restriction-site-detection` |

### Motif_Discovery (3)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Motif_Discovery/Overrepresented_Kmer_Discovery.md` | `overrepresented-kmer-discovery` |
| `docs/algorithms/Motif_Discovery/Regulatory_Elements.md` | `regulatory-elements` |
| `docs/algorithms/Motif_Discovery/Shared_Motifs.md` | `shared-motifs` |

### Oncology (37)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Oncology/Allele_Specific_Copy_Number_Derivation.md` | `allele-specific-copy-number-derivation` |
| `docs/algorithms/Oncology/Cancer_Cell_Fraction_Estimation.md` | `cancer-cell-fraction-estimation` |
| `docs/algorithms/Oncology/Cancer_Variant_Annotation.md` | `cancer-variant-annotation` |
| `docs/algorithms/Oncology/Clinical_Actionability_Assessment.md` | `clinical-actionability-assessment` |
| `docs/algorithms/Oncology/Clonal_Hematopoiesis_Filtering.md` | `clonal-hematopoiesis-filtering` |
| `docs/algorithms/Oncology/Clonal_Subclonal_Classification.md` | `clonal-subclonal-classification` |
| `docs/algorithms/Oncology/Complex_Rearrangement_Classification.md` | `complex-rearrangement-classification` |
| `docs/algorithms/Oncology/Copy_Number_Alteration_Classification.md` | `copy-number-alteration-classification` |
| `docs/algorithms/Oncology/CtDNA_Analysis.md` | `ctdna-analysis` |
| `docs/algorithms/Oncology/Driver_Mutation_Detection.md` | `driver-mutation-detection` |
| `docs/algorithms/Oncology/Focal_Amplification_Detection.md` | `focal-amplification-detection` |
| `docs/algorithms/Oncology/Fusion_Breakpoint_Analysis.md` | `fusion-breakpoint-analysis` |
| `docs/algorithms/Oncology/Fusion_Gene_Detection.md` | `fusion-gene-detection` |
| `docs/algorithms/Oncology/HLA_Nomenclature_And_Allele_Specific_LOH.md` | `hla-nomenclature-and-allele-specific-loh` |
| `docs/algorithms/Oncology/HRD_Score.md` | `hrd-score` |
| `docs/algorithms/Oncology/Homozygous_Deletion_Detection.md` | `homozygous-deletion-detection` |
| `docs/algorithms/Oncology/Immune_Infiltration_Estimation.md` | `immune-infiltration-estimation` |
| `docs/algorithms/Oncology/Known_Fusion_Database_Lookup.md` | `known-fusion-database-lookup` |
| `docs/algorithms/Oncology/Loss_Of_Heterozygosity.md` | `loss-of-heterozygosity` |
| `docs/algorithms/Oncology/MHC_Peptide_Binding_Classification.md` | `mhc-peptide-binding-classification` |
| `docs/algorithms/Oncology/MRD_Detection.md` | `mrd-detection` |
| `docs/algorithms/Oncology/Microsatellite_Instability_Detection.md` | `microsatellite-instability-detection` |
| `docs/algorithms/Oncology/Mutational_Process_Classification.md` | `mutational-process-classification` |
| `docs/algorithms/Oncology/Mutational_Signature_Exposure_Bootstrap.md` | `mutational-signature-exposure-bootstrap` |
| `docs/algorithms/Oncology/Mutational_Signature_Extraction_NMF.md` | `mutational-signature-extraction-nmf` |
| `docs/algorithms/Oncology/Mutational_Signature_Fitting.md` | `mutational-signature-fitting` |
| `docs/algorithms/Oncology/Neoantigen_Peptide_Generation.md` | `neoantigen-peptide-generation` |
| `docs/algorithms/Oncology/SBS96_Trinucleotide_Context_Catalog.md` | `sbs96-trinucleotide-context-catalog` |
| `docs/algorithms/Oncology/Sequencing_Artifact_Detection.md` | `sequencing-artifact-detection` |
| `docs/algorithms/Oncology/Somatic_Mutation_Calling.md` | `somatic-mutation-calling` |
| `docs/algorithms/Oncology/Tumor_Gene_Expression_Outlier.md` | `tumor-gene-expression-outlier` |
| `docs/algorithms/Oncology/Tumor_Heterogeneity_Analysis.md` | `tumor-heterogeneity-analysis` |
| `docs/algorithms/Oncology/Tumor_Mutational_Burden.md` | `tumor-mutational-burden` |
| `docs/algorithms/Oncology/Tumor_Phylogeny_Reconstruction.md` | `tumor-phylogeny-reconstruction` |
| `docs/algorithms/Oncology/Tumor_Ploidy_Estimation.md` | `tumor-ploidy-estimation` |
| `docs/algorithms/Oncology/Tumor_Purity_Estimation.md` | `tumor-purity-estimation` |
| `docs/algorithms/Oncology/Variant_Allele_Frequency.md` | `variant-allele-frequency` |

### PanGenome (3)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/PanGenome/Gene_Clustering.md` | `gene-clustering` |
| `docs/algorithms/PanGenome/Pan_Genome_Growth_Model.md` | `pan-genome-growth-model` |
| `docs/algorithms/PanGenome/Phylogenetic_Marker_Selection.md` | `phylogenetic-marker-selection` |

### Pattern_Matching (9)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md` | `approximate-matching-hamming` |
| `docs/algorithms/Pattern_Matching/Consensus_From_Alignment.md` | `consensus-from-alignment` |
| `docs/algorithms/Pattern_Matching/Edit_Distance.md` | `edit-distance` |
| `docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md` | `exact-pattern-search` |
| `docs/algorithms/Pattern_Matching/Frequent_Words_With_Mismatches.md` | `frequent-words-with-mismatches` |
| `docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Consensus.md` | `iupac-degenerate-consensus` |
| `docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md` | `iupac-degenerate-matching` |
| `docs/algorithms/Pattern_Matching/Position_Weight_Matrix.md` | `position-weight-matrix` |
| `docs/algorithms/Pattern_Matching/Suffix_Tree.md` | `suffix-tree` |

### Phylogenetics (6)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Phylogenetics/Bootstrap_Analysis.md` | `bootstrap-analysis` |
| `docs/algorithms/Phylogenetics/Distance_Matrix.md` | `distance-matrix` |
| `docs/algorithms/Phylogenetics/Newick_Format.md` | `newick-format` |
| `docs/algorithms/Phylogenetics/Tree_Comparison.md` | `tree-comparison` |
| `docs/algorithms/Phylogenetics/Tree_Construction.md` | `tree-construction` |
| `docs/algorithms/Phylogenetics/Tree_Statistics.md` | `tree-statistics` |

### Population_Genetics (8)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Population_Genetics/Allele_Frequency.md` | `allele-frequency` |
| `docs/algorithms/Population_Genetics/Ancestry_Estimation.md` | `ancestry-estimation` |
| `docs/algorithms/Population_Genetics/Diversity_Statistics.md` | `diversity-statistics` |
| `docs/algorithms/Population_Genetics/F_Statistics.md` | `f-statistics` |
| `docs/algorithms/Population_Genetics/Hardy_Weinberg_Test.md` | `hardy-weinberg-test` |
| `docs/algorithms/Population_Genetics/Integrated_Haplotype_Score.md` | `integrated-haplotype-score` |
| `docs/algorithms/Population_Genetics/Linkage_Disequilibrium.md` | `linkage-disequilibrium` |
| `docs/algorithms/Population_Genetics/Runs_Of_Homozygosity.md` | `runs-of-homozygosity` |

### ProteinMotif (10)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md` | `coiled-coil-prediction` |
| `docs/algorithms/ProteinMotif/Common_Motif_Finding.md` | `common-motif-finding` |
| `docs/algorithms/ProteinMotif/Domain_Prediction.md` | `domain-prediction` |
| `docs/algorithms/ProteinMotif/Low_Complexity_Region_Detection.md` | `low-complexity-region-detection` |
| `docs/algorithms/ProteinMotif/Motif_Search.md` | `motif-search` |
| `docs/algorithms/ProteinMotif/PROSITE_Pattern_Matching.md` | `prosite-pattern-matching` |
| `docs/algorithms/ProteinMotif/Pattern_Matching_Methods.md` | `pattern-matching-methods` |
| `docs/algorithms/ProteinMotif/Profile_HMM_Domain_Detection.md` | `profile-hmm-domain-detection` |
| `docs/algorithms/ProteinMotif/Signal_Peptide_Prediction.md` | `signal-peptide-prediction` |
| `docs/algorithms/ProteinMotif/Transmembrane_Helix_Prediction.md` | `transmembrane-helix-prediction` |

### ProteinPred (4)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/ProteinPred/Disorder_Prediction.md` | `disorder-prediction` |
| `docs/algorithms/ProteinPred/Disorder_Propensity.md` | `disorder-propensity` |
| `docs/algorithms/ProteinPred/Disordered_Region_Detection.md` | `disordered-region-detection` |
| `docs/algorithms/ProteinPred/MoRF_Prediction.md` | `morf-prediction` |

### Quality (2)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Quality/Phred_Score_Handling.md` | `phred-score-handling` |
| `docs/algorithms/Quality/Quality_Statistics.md` | `quality-statistics` |

### Repeat_Analysis (5)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Repeat_Analysis/Direct_Repeat_Detection.md` | `direct-repeat-detection` |
| `docs/algorithms/Repeat_Analysis/Inverted_Repeat_Detection.md` | `inverted-repeat-detection` |
| `docs/algorithms/Repeat_Analysis/Microsatellite_Detection.md` | `microsatellite-detection` |
| `docs/algorithms/Repeat_Analysis/Palindrome_Detection.md` | `palindrome-detection` |
| `docs/algorithms/Repeat_Analysis/Tandem_Repeat_Detection.md` | `tandem-repeat-detection` |

### RnaStructure (13)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/RnaStructure/Dot_Bracket_Notation.md` | `dot-bracket-notation` |
| `docs/algorithms/RnaStructure/Hairpin_Energy_Calculation.md` | `hairpin-energy-calculation` |
| `docs/algorithms/RnaStructure/Inverted_Repeats.md` | `inverted-repeats` |
| `docs/algorithms/RnaStructure/Minimum_Free_Energy.md` | `minimum-free-energy` |
| `docs/algorithms/RnaStructure/Pseudoknot_Detection.md` | `pseudoknot-detection` |
| `docs/algorithms/RnaStructure/Pseudoknot_Prediction.md` | `pseudoknot-prediction` |
| `docs/algorithms/RnaStructure/Pseudoknot_Prediction_Recursive.md` | `pseudoknot-prediction-recursive` |
| `docs/algorithms/RnaStructure/RNA_Base_Pairing.md` | `rna-base-pairing` |
| `docs/algorithms/RnaStructure/RNA_Free_Energy.md` | `rna-free-energy` |
| `docs/algorithms/RnaStructure/RNA_Partition_Function.md` | `rna-partition-function` |
| `docs/algorithms/RnaStructure/RNA_Secondary_Structure.md` | `rna-secondary-structure` |
| `docs/algorithms/RnaStructure/RNA_Stemloop.md` | `rna-stemloop` |
| `docs/algorithms/RnaStructure/Turner_McCaskill_Partition_Function.md` | `turner-mccaskill-partition-function` |

### Sequence_Composition (8)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Sequence_Composition/GC_Skew.md` | `gc-skew` |
| `docs/algorithms/Sequence_Composition/Linguistic_Complexity.md` | `linguistic-complexity` |
| `docs/algorithms/Sequence_Composition/RNA_Complement.md` | `rna-complement` |
| `docs/algorithms/Sequence_Composition/Replication_Origin_Prediction.md` | `replication-origin-prediction` |
| `docs/algorithms/Sequence_Composition/Sequence_Composition.md` | `sequence-composition` |
| `docs/algorithms/Sequence_Composition/Sequence_Composition_Statistics.md` | `sequence-composition-statistics` |
| `docs/algorithms/Sequence_Composition/Sequence_Validation.md` | `sequence-validation` |
| `docs/algorithms/Sequence_Composition/Shannon_Entropy.md` | `shannon-entropy` |

### Splicing (3)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Splicing/Acceptor_Site_Detection.md` | `acceptor-site-detection` |
| `docs/algorithms/Splicing/Donor_Site_Detection.md` | `donor-site-detection` |
| `docs/algorithms/Splicing/Gene_Structure_Prediction.md` | `gene-structure-prediction` |

### Statistics (11)

| Algorithm doc | Expected slug |
| --- | --- |
| `docs/algorithms/Statistics/Codon_Frequencies.md` | `codon-frequencies` |
| `docs/algorithms/Statistics/DNA_Thermodynamics.md` | `dna-thermodynamics` |
| `docs/algorithms/Statistics/Dinucleotide_Analysis.md` | `dinucleotide-analysis` |
| `docs/algorithms/Statistics/Entropy_Profile.md` | `entropy-profile` |
| `docs/algorithms/Statistics/GC_Content_Profile.md` | `gc-content-profile` |
| `docs/algorithms/Statistics/Hydrophobicity_Analysis.md` | `hydrophobicity-analysis` |
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
  docs, not algorithm units — ingest only if a navigational need arises.
- The `docs/Evidence/**` campaign (175 of 213 remaining) is the primary driver: each
  Evidence ingest typically creates or extends the concept that also covers the
  matching algorithm doc, clearing a pending row here as a side effect.
