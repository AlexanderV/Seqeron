# Algorithm Docs Reconciliation Backlog

Reconciliation of the 245 canonical algorithm reference docs under `docs/algorithms/` against the Evidence-driven concept pages. `docs/algorithms/**` is NOT excluded from coverage.

- **done-via-concept** — a concept page already synthesizes this algorithm; the doc path has been added to that page's `sources:` frontmatter.
- **pending** — no concept page yet; a real gap folded into the main ingest campaign (expected future page slug shown).
- **meta** — index/map files, not an algorithm (no page expected).

Generated at commit `9ce49bade5c1`.

## (root)

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `CANONICAL_MAP.md` | meta | — |
| `README.md` | meta | — |

## Alignment

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Alignment_Statistics.md` | done-via-concept | [[alignment-statistics]] |
| `Global_Alignment_Needleman_Wunsch.md` | done-via-concept | [[global-alignment-needleman-wunsch]] |
| `Local_Alignment_Smith_Waterman.md` | pending | `local-alignment-smith-waterman` (expected) |
| `Multiple_Sequence_Alignment.md` | done-via-concept | [[multiple-sequence-alignment]] |
| `Semi_Global_Alignment.md` | done-via-concept | [[semi-global-alignment-fitting]] |

## Analysis

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Open_Reading_Frame_Detection.md` | pending | `open-reading-frame-detection` (expected) |
| `Sequence_Similarity.md` | pending | `sequence-similarity` (expected) |

## Annotation

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `GFF3_IO.md` | pending | `gff3-io` (expected) |
| `Gene_Prediction.md` | pending | `gene-prediction` (expected) |
| `ORF_Detection.md` | pending | `orf-detection` (expected) |
| `Promoter_Detection.md` | pending | `promoter-detection` (expected) |
| `Relative_Synonymous_Codon_Usage.md` | done-via-concept | [[relative-synonymous-codon-usage]] |
| `Repetitive_Element_Detection.md` | done-via-concept | [[repetitive-element-detection]] |

## Assembly

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Assembly_Statistics.md` | done-via-concept | [[assembly-statistics]] |
| `Consensus_Computation.md` | done-via-concept | [[consensus-sequence]] |
| `Coverage_Calculation.md` | done-via-concept | [[coverage-depth-calculation]] |
| `De_Bruijn_Graph_Assembly.md` | done-via-concept | [[de-bruijn-graph-assembly]] |
| `Error_Correction.md` | done-via-concept | [[kmer-spectrum-error-correction]] |
| `Overlap_Layout_Consensus.md` | done-via-concept | [[overlap-layout-consensus-assembly]] |
| `Quality_Trimming.md` | done-via-concept | [[quality-trimming-running-sum]] |

## Chromosome_Analysis

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Aneuploidy_Detection.md` | done-via-concept | [[aneuploidy-detection]] |
| `Centromere_Analysis.md` | done-via-concept | [[centromere-analysis]] |
| `Higher_Order_Repeat_Detection.md` | pending | `higher-order-repeat-detection` (expected) |
| `Karyotype_Analysis.md` | done-via-concept | [[karyotype-analysis]] |
| `Synteny_Analysis.md` | done-via-concept | [[synteny-and-rearrangement-detection]] |
| `Telomere_Analysis.md` | done-via-concept | [[telomere-analysis]] |

## Codon

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Codon_Usage_Statistics.md` | pending | `codon-usage-statistics` (expected) |
| `Effective_Number_of_Codons.md` | done-via-concept | [[effective-number-of-codons]] |
| `Relative_Synonymous_Codon_Usage.md` | done-via-concept | [[relative-synonymous-codon-usage]] |

## Codon_Optimization

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `CAI_Calculation.md` | done-via-concept | [[codon-adaptation-index]] |
| `Codon_Usage_Analysis.md` | done-via-concept | [[codon-usage-comparison]] |
| `Rare_Codon_Detection.md` | done-via-concept | [[rare-codon-analysis]] |
| `Sequence_Optimization.md` | done-via-concept | [[codon-optimization]] |

## Comparative_Genomics

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Average_Nucleotide_Identity.md` | done-via-concept | [[average-nucleotide-identity]] |
| `Conserved_Gene_Clusters.md` | done-via-concept | [[conserved-gene-clusters-common-intervals]] |
| `Dot_Plot_Generation.md` | done-via-concept | [[dot-plot-word-match]] |
| `Genome_Comparison.md` | done-via-concept | [[genome-comparison-core-dispensable]] |
| `Genome_Rearrangement_Detection.md` | done-via-concept | [[genome-rearrangement-breakpoint-distance]] |
| `Ortholog_Identification.md` | done-via-concept | [[ortholog-detection-reciprocal-best-hits]] |
| `Reciprocal_Best_Hits.md` | done-via-concept | [[ortholog-detection-reciprocal-best-hits]] |
| `Reversal_Distance.md` | done-via-concept | [[genome-rearrangement-breakpoint-distance]] |
| `Synteny_Block_Detection.md` | done-via-concept | [[synteny-and-rearrangement-detection]] |

## Complexity

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `DUST_Score.md` | pending | `dust-score` (expected) |
| `K-mer_Entropy.md` | pending | `k-mer-entropy` (expected) |
| `Lempel_Ziv_Complexity.md` | pending | `lempel-ziv-complexity` (expected) |
| `Windowed_Complexity.md` | pending | `windowed-complexity` (expected) |

## Epigenetics

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Bisulfite_Sequencing_Analysis.md` | pending | `bisulfite-sequencing-analysis` (expected) |
| `Chromatin_State_Prediction.md` | pending | `chromatin-state-prediction` (expected) |
| `CpG_Site_Detection.md` | pending | `cpg-site-detection` (expected) |
| `Differentially_Methylated_Regions.md` | pending | `differentially-methylated-regions` (expected) |
| `Epigenetic_Age_Estimation.md` | pending | `epigenetic-age-estimation` (expected) |
| `Methylation_Analysis.md` | pending | `methylation-analysis` (expected) |

## Extended_Annotation

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Coding_Potential_Calculation.md` | done-via-concept | [[coding-potential-hexamer-score]] |

## Extended_Assembly

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Contig_Merging.md` | done-via-concept | [[contig-merge-overlap-collapse]] |
| `Scaffolding.md` | done-via-concept | [[scaffolding]] |

## Extended_GC_Skew_Analysis

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `AT_Skew.md` | pending | `at-skew` (expected) |
| `Comprehensive_GC_Analysis.md` | pending | `comprehensive-gc-analysis` (expected) |

## FileIO

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `BED_Parsing.md` | pending | `bed-parsing` (expected) |
| `EMBL_Parsing.md` | pending | `embl-parsing` (expected) |
| `FASTA_Parsing.md` | pending | `fasta-parsing` (expected) |
| `FASTQ_Parsing.md` | pending | `fastq-parsing` (expected) |
| `GFF_Parsing.md` | pending | `gff-parsing` (expected) |
| `GenBank_Parsing.md` | pending | `genbank-parsing` (expected) |
| `VCF_Parsing.md` | pending | `vcf-parsing` (expected) |

## Genomic_Analysis

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Tandem_Repeat_Detection.md` | pending | `tandem-repeat-detection` (expected) |

## K-mer

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Asynchronous_K-mer_Counting.md` | pending | `asynchronous-k-mer-counting` (expected) |
| `Both_Strand_Kmer_Counting.md` | pending | `both-strand-kmer-counting` (expected) |
| `K-mer_Counting.md` | pending | `k-mer-counting` (expected) |
| `K-mer_Euclidean_Distance.md` | pending | `k-mer-euclidean-distance` (expected) |
| `K-mer_Frequency_Analysis.md` | pending | `k-mer-frequency-analysis` (expected) |
| `K-mer_Generation.md` | pending | `k-mer-generation` (expected) |
| `K-mer_Positions.md` | pending | `k-mer-positions` (expected) |
| `K-mer_Search.md` | pending | `k-mer-search` (expected) |
| `K-mer_Statistics.md` | pending | `k-mer-statistics` (expected) |
| `Unique_And_MinCount_Kmers.md` | pending | `unique-and-mincount-kmers` (expected) |

## Metagenomics

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Alpha_Diversity.md` | pending | `alpha-diversity` (expected) |
| `Antibiotic_Resistance_Detection.md` | pending | `antibiotic-resistance-detection` (expected) |
| `Beta_Diversity.md` | pending | `beta-diversity` (expected) |
| `Functional_Prediction.md` | pending | `functional-prediction` (expected) |
| `Genome_Binning.md` | pending | `genome-binning` (expected) |
| `PanGenome_Core_Accessory.md` | pending | `pangenome-core-accessory` (expected) |
| `Pathway_Enrichment_ORA.md` | pending | `pathway-enrichment-ora` (expected) |
| `Significant_Taxa_Detection.md` | pending | `significant-taxa-detection` (expected) |
| `Taxonomic_Classification.md` | pending | `taxonomic-classification` (expected) |
| `Taxonomic_Profile.md` | pending | `taxonomic-profile` (expected) |

## MiRNA

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `MiRNA_Target_Pairing.md` | pending | `mirna-target-pairing` (expected) |
| `Pre_miRNA_Detection.md` | pending | `pre-mirna-detection` (expected) |
| `Seed_Sequence_Analysis.md` | pending | `seed-sequence-analysis` (expected) |
| `Target_Site_Prediction.md` | pending | `target-site-prediction` (expected) |

## MolTools

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `DNA_Dimer_Tm.md` | pending | `dna-dimer-tm` (expected) |
| `DNA_Hairpin_Folding_Tm.md` | pending | `dna-hairpin-folding-tm` (expected) |
| `DNA_Hairpin_Special_Loop_Bonus.md` | pending | `dna-hairpin-special-loop-bonus` (expected) |
| `Guide_RNA_Design.md` | pending | `guide-rna-design` (expected) |
| `Hybridization_Probe_Design.md` | pending | `hybridization-probe-design` (expected) |
| `LNA_Adjusted_Nearest_Neighbor_Tm.md` | pending | `lna-adjusted-nearest-neighbor-tm` (expected) |
| `Melting_Temperature.md` | pending | `melting-temperature` (expected) |
| `NearestNeighbor_Salt_Corrected_Tm.md` | pending | `nearestneighbor-salt-corrected-tm` (expected) |
| `Off_Target_Analysis.md` | pending | `off-target-analysis` (expected) |
| `PAM_Site_Detection.md` | pending | `pam-site-detection` (expected) |
| `Primer3_Penalty_Objective.md` | pending | `primer3-penalty-objective` (expected) |
| `Primer_Design.md` | pending | `primer-design` (expected) |
| `Primer_Structure_Analysis.md` | pending | `primer-structure-analysis` (expected) |
| `Probe_Validation.md` | pending | `probe-validation` (expected) |
| `Restriction_Digest_Simulation.md` | pending | `restriction-digest-simulation` (expected) |
| `Restriction_Enzyme_Filtering.md` | pending | `restriction-enzyme-filtering` (expected) |
| `Restriction_Site_Detection.md` | pending | `restriction-site-detection` (expected) |

## Motif_Analysis

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Known_Motif_Search.md` | pending | `known-motif-search` (expected) |

## Motif_Discovery

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Overrepresented_Kmer_Discovery.md` | pending | `overrepresented-kmer-discovery` (expected) |
| `Regulatory_Elements.md` | pending | `regulatory-elements` (expected) |
| `Shared_Motifs.md` | pending | `shared-motifs` (expected) |

## Oncology

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Allele_Specific_Copy_Number_Derivation.md` | pending | `allele-specific-copy-number-derivation` (expected) |
| `Cancer_Cell_Fraction_Estimation.md` | pending | `cancer-cell-fraction-estimation` (expected) |
| `Cancer_Variant_Annotation.md` | pending | `cancer-variant-annotation` (expected) |
| `Clinical_Actionability_Assessment.md` | pending | `clinical-actionability-assessment` (expected) |
| `Clonal_Hematopoiesis_Filtering.md` | pending | `clonal-hematopoiesis-filtering` (expected) |
| `Clonal_Subclonal_Classification.md` | pending | `clonal-subclonal-classification` (expected) |
| `Complex_Rearrangement_Classification.md` | pending | `complex-rearrangement-classification` (expected) |
| `Copy_Number_Alteration_Classification.md` | pending | `copy-number-alteration-classification` (expected) |
| `CtDNA_Analysis.md` | pending | `ctdna-analysis` (expected) |
| `Driver_Mutation_Detection.md` | pending | `driver-mutation-detection` (expected) |
| `Focal_Amplification_Detection.md` | pending | `focal-amplification-detection` (expected) |
| `Fusion_Breakpoint_Analysis.md` | pending | `fusion-breakpoint-analysis` (expected) |
| `Fusion_Gene_Detection.md` | pending | `fusion-gene-detection` (expected) |
| `HLA_Nomenclature_And_Allele_Specific_LOH.md` | pending | `hla-nomenclature-and-allele-specific-loh` (expected) |
| `HRD_Score.md` | pending | `hrd-score` (expected) |
| `Homozygous_Deletion_Detection.md` | pending | `homozygous-deletion-detection` (expected) |
| `Immune_Infiltration_Estimation.md` | pending | `immune-infiltration-estimation` (expected) |
| `Known_Fusion_Database_Lookup.md` | pending | `known-fusion-database-lookup` (expected) |
| `Loss_Of_Heterozygosity.md` | pending | `loss-of-heterozygosity` (expected) |
| `MHC_Peptide_Binding_Classification.md` | pending | `mhc-peptide-binding-classification` (expected) |
| `MRD_Detection.md` | pending | `mrd-detection` (expected) |
| `Microsatellite_Instability_Detection.md` | pending | `microsatellite-instability-detection` (expected) |
| `Mutational_Process_Classification.md` | pending | `mutational-process-classification` (expected) |
| `Mutational_Signature_Exposure_Bootstrap.md` | pending | `mutational-signature-exposure-bootstrap` (expected) |
| `Mutational_Signature_Extraction_NMF.md` | pending | `mutational-signature-extraction-nmf` (expected) |
| `Mutational_Signature_Fitting.md` | pending | `mutational-signature-fitting` (expected) |
| `Neoantigen_Peptide_Generation.md` | pending | `neoantigen-peptide-generation` (expected) |
| `SBS96_Trinucleotide_Context_Catalog.md` | pending | `sbs96-trinucleotide-context-catalog` (expected) |
| `Sequencing_Artifact_Detection.md` | pending | `sequencing-artifact-detection` (expected) |
| `Somatic_Mutation_Calling.md` | pending | `somatic-mutation-calling` (expected) |
| `Tumor_Gene_Expression_Outlier.md` | pending | `tumor-gene-expression-outlier` (expected) |
| `Tumor_Heterogeneity_Analysis.md` | pending | `tumor-heterogeneity-analysis` (expected) |
| `Tumor_Mutational_Burden.md` | pending | `tumor-mutational-burden` (expected) |
| `Tumor_Phylogeny_Reconstruction.md` | pending | `tumor-phylogeny-reconstruction` (expected) |
| `Tumor_Ploidy_Estimation.md` | pending | `tumor-ploidy-estimation` (expected) |
| `Tumor_Purity_Estimation.md` | pending | `tumor-purity-estimation` (expected) |
| `Variant_Allele_Frequency.md` | pending | `variant-allele-frequency` (expected) |

## PanGenome

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Gene_Clustering.md` | pending | `gene-clustering` (expected) |
| `Pan_Genome_Growth_Model.md` | pending | `pan-genome-growth-model` (expected) |
| `Phylogenetic_Marker_Selection.md` | pending | `phylogenetic-marker-selection` (expected) |

## Pattern_Matching

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Approximate_Matching_Hamming.md` | pending | `approximate-matching-hamming` (expected) |
| `Consensus_From_Alignment.md` | pending | `consensus-from-alignment` (expected) |
| `Edit_Distance.md` | pending | `edit-distance` (expected) |
| `Exact_Pattern_Search.md` | pending | `exact-pattern-search` (expected) |
| `Frequent_Words_With_Mismatches.md` | pending | `frequent-words-with-mismatches` (expected) |
| `IUPAC_Degenerate_Consensus.md` | pending | `iupac-degenerate-consensus` (expected) |
| `IUPAC_Degenerate_Matching.md` | pending | `iupac-degenerate-matching` (expected) |
| `Position_Weight_Matrix.md` | pending | `position-weight-matrix` (expected) |
| `Suffix_Tree.md` | pending | `suffix-tree` (expected) |

## Phylogenetics

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Bootstrap_Analysis.md` | pending | `bootstrap-analysis` (expected) |
| `Distance_Matrix.md` | pending | `distance-matrix` (expected) |
| `Newick_Format.md` | pending | `newick-format` (expected) |
| `Tree_Comparison.md` | pending | `tree-comparison` (expected) |
| `Tree_Construction.md` | pending | `tree-construction` (expected) |
| `Tree_Statistics.md` | pending | `tree-statistics` (expected) |

## Population_Genetics

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Allele_Frequency.md` | pending | `allele-frequency` (expected) |
| `Ancestry_Estimation.md` | pending | `ancestry-estimation` (expected) |
| `Diversity_Statistics.md` | pending | `diversity-statistics` (expected) |
| `F_Statistics.md` | pending | `f-statistics` (expected) |
| `Hardy_Weinberg_Test.md` | pending | `hardy-weinberg-test` (expected) |
| `Integrated_Haplotype_Score.md` | pending | `integrated-haplotype-score` (expected) |
| `Linkage_Disequilibrium.md` | pending | `linkage-disequilibrium` (expected) |
| `Runs_Of_Homozygosity.md` | pending | `runs-of-homozygosity` (expected) |

## ProteinMotif

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Coiled_Coil_Prediction.md` | pending | `coiled-coil-prediction` (expected) |
| `Common_Motif_Finding.md` | pending | `common-motif-finding` (expected) |
| `Domain_Prediction.md` | pending | `domain-prediction` (expected) |
| `Low_Complexity_Region_Detection.md` | pending | `low-complexity-region-detection` (expected) |
| `Motif_Search.md` | pending | `motif-search` (expected) |
| `PROSITE_Pattern_Matching.md` | pending | `prosite-pattern-matching` (expected) |
| `Pattern_Matching_Methods.md` | pending | `pattern-matching-methods` (expected) |
| `Profile_HMM_Domain_Detection.md` | pending | `profile-hmm-domain-detection` (expected) |
| `Signal_Peptide_Prediction.md` | pending | `signal-peptide-prediction` (expected) |
| `Transmembrane_Helix_Prediction.md` | pending | `transmembrane-helix-prediction` (expected) |

## ProteinPred

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Disorder_Prediction.md` | pending | `disorder-prediction` (expected) |
| `Disorder_Propensity.md` | pending | `disorder-propensity` (expected) |
| `Disordered_Region_Detection.md` | pending | `disordered-region-detection` (expected) |
| `Low_Complexity_Region_Detection.md` | done-via-concept | [[protein-low-complexity-seg]] |
| `MoRF_Prediction.md` | pending | `morf-prediction` (expected) |

## Quality

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Phred_Score_Handling.md` | pending | `phred-score-handling` (expected) |
| `Quality_Statistics.md` | pending | `quality-statistics` (expected) |

## Repeat_Analysis

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Direct_Repeat_Detection.md` | pending | `direct-repeat-detection` (expected) |
| `Inverted_Repeat_Detection.md` | pending | `inverted-repeat-detection` (expected) |
| `Microsatellite_Detection.md` | pending | `microsatellite-detection` (expected) |
| `Palindrome_Detection.md` | pending | `palindrome-detection` (expected) |
| `Repeat_Detection.md` | pending | `repeat-detection` (expected) |
| `Tandem_Repeat_Detection.md` | pending | `tandem-repeat-detection` (expected) |

## RnaStructure

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Dot_Bracket_Notation.md` | pending | `dot-bracket-notation` (expected) |
| `Hairpin_Energy_Calculation.md` | pending | `hairpin-energy-calculation` (expected) |
| `Inverted_Repeats.md` | pending | `inverted-repeats` (expected) |
| `Minimum_Free_Energy.md` | pending | `minimum-free-energy` (expected) |
| `Pseudoknot_Detection.md` | pending | `pseudoknot-detection` (expected) |
| `Pseudoknot_Prediction.md` | pending | `pseudoknot-prediction` (expected) |
| `Pseudoknot_Prediction_Recursive.md` | pending | `pseudoknot-prediction-recursive` (expected) |
| `RNA_Base_Pairing.md` | pending | `rna-base-pairing` (expected) |
| `RNA_Free_Energy.md` | pending | `rna-free-energy` (expected) |
| `RNA_Partition_Function.md` | pending | `rna-partition-function` (expected) |
| `RNA_Secondary_Structure.md` | pending | `rna-secondary-structure` (expected) |
| `RNA_Stemloop.md` | pending | `rna-stemloop` (expected) |
| `Turner_McCaskill_Partition_Function.md` | pending | `turner-mccaskill-partition-function` (expected) |

## Sequence_Comparison

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Common_Region_Detection.md` | pending | `common-region-detection` (expected) |

## Sequence_Composition

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `GC_Skew.md` | pending | `gc-skew` (expected) |
| `Linguistic_Complexity.md` | pending | `linguistic-complexity` (expected) |
| `RNA_Complement.md` | pending | `rna-complement` (expected) |
| `Replication_Origin_Prediction.md` | pending | `replication-origin-prediction` (expected) |
| `Sequence_Composition.md` | pending | `sequence-composition` (expected) |
| `Sequence_Composition_Statistics.md` | pending | `sequence-composition-statistics` (expected) |
| `Sequence_Validation.md` | pending | `sequence-validation` (expected) |
| `Shannon_Entropy.md` | pending | `shannon-entropy` (expected) |

## Splicing

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Acceptor_Site_Detection.md` | pending | `acceptor-site-detection` (expected) |
| `Donor_Site_Detection.md` | pending | `donor-site-detection` (expected) |
| `Gene_Structure_Prediction.md` | pending | `gene-structure-prediction` (expected) |

## Statistics

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Codon_Frequencies.md` | pending | `codon-frequencies` (expected) |
| `DNA_Thermodynamics.md` | pending | `dna-thermodynamics` (expected) |
| `Dinucleotide_Analysis.md` | pending | `dinucleotide-analysis` (expected) |
| `Entropy_Profile.md` | pending | `entropy-profile` (expected) |
| `GC_Content_Profile.md` | pending | `gc-content-profile` (expected) |
| `Hydrophobicity_Analysis.md` | pending | `hydrophobicity-analysis` (expected) |
| `Isoelectric_Point.md` | pending | `isoelectric-point` (expected) |
| `Melting_Temperature.md` | pending | `melting-temperature` (expected) |
| `Molecular_Weight_Calculation.md` | pending | `molecular-weight-calculation` (expected) |
| `Secondary_Structure_Prediction.md` | pending | `secondary-structure-prediction` (expected) |
| `Sequence_Summary.md` | pending | `sequence-summary` (expected) |

## StructuralVar

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Breakpoint_Detection.md` | pending | `breakpoint-detection` (expected) |
| `Copy_Number_Variation.md` | pending | `copy-number-variation` (expected) |
| `SV_Detection.md` | pending | `sv-detection` (expected) |

## Transcriptome

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Alternative_Splicing.md` | pending | `alternative-splicing` (expected) |
| `Differential_Expression.md` | pending | `differential-expression` (expected) |
| `Expression_Quantification.md` | pending | `expression-quantification` (expected) |

## Translation

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Codon_Translation.md` | pending | `codon-translation` (expected) |
| `Protein_Translation.md` | pending | `protein-translation` (expected) |
| `Six_Frame_Translation.md` | pending | `six-frame-translation` (expected) |

## Variants

| Doc | Status | Concept page / expected slug |
| --- | --- | --- |
| `Indel_Detection.md` | pending | `indel-detection` (expected) |
| `SNP_Detection.md` | pending | `snp-detection` (expected) |
| `Variant_Annotation.md` | pending | `variant-annotation` (expected) |
| `Variant_Detection.md` | pending | `variant-detection` (expected) |
