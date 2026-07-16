---
type: index
title: "Ingestion backlog — docs/algorithms reconciliation + queued sources"
tags: [meta, coverage]
created: 2026-07-09
updated: 2026-07-16
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

Status at generation: **179** algorithm docs covered-via-concept, **53** pending across 9 domains
(Quality/Quality_Statistics → [[fastq-quality-statistics]] resolved 2026-07-16, closing the Quality
domain (last pending doc) (REUSE: reconciled against the existing Evidence-derived concept, which
already synthesized the mean/median/min/max, population variance/σ (÷N), inclusive Q20/Q30 thresholds,
the Q30 NGS benchmark, the arithmetic-vs-probability mean distinction, and the `5?I`/`5II?`/single-`I`
oracles; enriched only with the spec's genuinely-distinct implementation content: the three
`QualityScoreAnalyzer` entry points, the multi-read `PerPositionMeanQuality` delegate variant, the
O(n log n) sort-dominated / O(n) Q30 complexity, and the sibling `CalculateExpectedErrors` /
`PhredToErrorProbability` error-probability summaries);
Quality/Phred_Score_Handling → [[phred-quality-encoding]] resolved 2026-07-16 (REUSE: reconciled
against the existing Evidence-derived concept, which already synthesized the Q=−10log₁₀p definition,
the Phred+33/Phred+64 offsets and ranges, auto-detection, the score-preserving ±31 re-offset and its
overflow rules; enriched only with the spec's genuinely-distinct implementation content — the three
`QualityScoreAnalyzer` canonical methods (`ParseQualityString`/`ToQualityString`/`ConvertEncoding`) in
Seqeron.Genomics.IO, the O(n) single-pass cost, the `Auto`→`DetectEncoding` (parse) vs `Auto`→Phred+33
(encode) resolution, the INV-03 round-trip invariant, the ArgumentNullException/ArgumentOutOfRangeException
contract, the suffix-tree-not-applicable note, and the legacy `QualityStringToPhred`/`PhredToQualityString`
helpers that lack range validation));
ProteinPred/MoRF_Prediction → [[morf-prediction-dip-in-disorder]] resolved 2026-07-16, closing the
ProteinPred domain (last pending doc) (REUSE: reconciled against the existing Evidence-derived
concept, which already documents the dip-in-disorder criterion, the 10–70 Mohan band, the TOP-IDP
disorder score, α/β/ι sub-types and the paywalled-flank assumption; enriched only with the spec's
genuinely-distinct implementation content: the `PredictMoRFs(string, int, int)` entry point in
`DisorderPredictor.cs`, the O(n·w)/O(n) cost with the 21-residue window, the 0-based-inclusive /
non-overlapping / start-ordered output contract, the explicit `(0.5 − mean d)/0.5` clamped score
formula, the 0.5-MoRF-vs-0.542-TOP-IDP dual-threshold distinction, the suffix-tree-not-applicable
decision and the P/L/I worked oracle (Start 29–End 50, Score ≈ 0.2759 → 0.3996));
ProteinMotif/Signal_Peptide_Prediction → [[protein-domain-and-signal-peptide-prediction]] resolved
2026-07-16 (REUSE: the dedicated von Heijne 1986 signal-peptide spec, unit PROTMOTIF-SP-001, is already
fully owned by this Evidence-derived concept — which documents the current EMBOSS `sigcleave` weight-matrix
method with the P17644 oracle; the spec describes the SAME correct method (no superseded tripartite model —
that was the sibling `Domain_Prediction.md`), so no contradiction; enriched only with the spec's
genuinely-distinct implementation content: the `PredictSignalPeptide`/`BuildWeightMatrix` entry points with
the log-odds matrix built once at static init, the O(n·15)=O(n) time / O(1) space cost, the fixed-width-PWM
suffix-tree-not-applicable decision, and the min-15-residue-window accepted deviation);
ProteinMotif/Profile_HMM_Domain_Detection → [[protein-domain-and-signal-peptide-prediction]] resolved
2026-07-16 (REUSE: the per-algorithm profile-HMM spec is the SAME test unit — PROTMOTIF-DOMAIN-001 — as
`Domain_Prediction`, already owned by this Evidence-derived concept, which had a dedicated `FindDomainsByHmm`
Plan7 section covering Viterbi/Forward, hmmsearch-parity `pre_score`, null2, the Gumbel/exponential E-value
layer, multi-domain envelope decomposition and stochastic-traceback clustering; enriched only with the
spec's genuinely-distinct implementation content: the HMMER3/f parser layout (−ln-to-5-decimals storage,
`*`→−∞, COMPO/BEGIN, the 7-transition node order), the two-row O(n·M)/O(M) Viterbi/Forward DP shape, the
glocal-default-vs-local-multihit-opt-in distinction (spec Deviation #1) with `minBitScore`/`Z` defaults and
the `FindDomainHitsByHmm`/`ScoreDomainHmmEValue`/`FindDomainEnvelopes` method surface, and the INV-HMM-01/03/05 invariants);
ProteinMotif/Pattern_Matching_Methods → [[protein-motif-pattern-search]] resolved 2026-07-16 (REUSE:
this Evidence-derived concept already cites the PROTMOTIF-PATTERN-001 Evidence and is ahead of the spec;
added the spec to `sources:` and enriched only its genuinely-distinct implementation content — `FindDomains`
as the fourth pattern-matching primitive delegating to the same lookahead+IC engine and wrapping hits as
`ProteinDomain`, plus the suffix-tree-evaluated-not-used matcher decision);
ProteinMotif/PROSITE_Pattern_Matching → [[protein-motif-pattern-search]] resolved 2026-07-16 (REUSE:
the dedicated PROSITE spec, unit PROTMOTIF-PROSITE-001, is already fully covered by this Evidence-derived
concept — which is ahead of the spec; enriched only with the spec's explicit scope boundary: syntax
converter + regex-search wrapper with no PROSITE profile/matrix scanning, no external-catalog lookup, no
ScanProsite result metadata, and repository-defined `Score`/`EValue`);
ProteinMotif/Motif_Search → [[protein-motif-pattern-search]] resolved 2026-07-16 (REUSE: the primary
per-algorithm spec is exactly this Evidence-derived concept, unit PROTMOTIF-FIND-001; enriched only with
the `ParseRegexAllowedCounts` regex-walk helper that supplies the per-position allowed-residue count
feeding the information-content score);
ProteinMotif/Domain_Prediction → [[protein-domain-and-signal-peptide-prediction]] resolved 2026-07-16
(REUSE: the primary per-algorithm spec reconciles against the existing Evidence-derived concept, which is
already ahead of the spec — it documents the current von Heijne 1986 weight-matrix signal-peptide method
and the Plan7 profile-HMM engine, whereas the spec still describes the superseded/fabricated tripartite
n/h/c model; enriched with the `FindMotifByPattern` score-delegation / information-content-score provenance
and the O(n·d), d=3 domain-scan cost);
ProteinMotif/Common_Motif_Finding → [[common-protein-motifs]] resolved 2026-07-16 (REUSE: the primary
per-algorithm spec is exactly this Evidence-derived concept; enriched with the O(p·n) `FindCommonMotifs`
complexity, the `FindMotifByPattern` delegation + entry-point decomposition, the suffix-tree-not-applicable
design decision, and the `FindAllKnownMotifs` registry-naming deviation);
ProteinMotif/Coiled_Coil_Prediction → [[coiled-coil-prediction]] resolved 2026-07-16 (REUSE: the primary
per-algorithm spec is exactly this Evidence-derived concept; enriched with the `BestHeptadOccupancy`/`BuildRegion`
helper decomposition, the single-pass array precompute + lazy-yield shape, and the suffix-tree-not-applicable note);
Population_Genetics/Runs_Of_Homozygosity → [[runs-of-homozygosity-inbreeding]] resolved 2026-07-16,
closing the Population_Genetics domain (REUSE: reconciled against the existing Evidence-derived concept;
enriched with the O(n log n)/O(n) `FindROH` cost, the O(m)/O(1) F_ROH cost, the inclusive-run vs
half-open-segment coordinate distinction, and the eager argument-validation / `genomeLength ≤ 0 → 0`
contract);
K-mer, Metagenomics, MolTools, Oncology, PanGenome and Pattern_Matching domains now fully covered;
Phylogenetics/Tree_Construction → [[distance-based-tree-construction]] resolved 2026-07-15 (REUSE: the primary UPGMA/NJ spec is exactly this Evidence-derived hinge concept; enriched with the `distanceMethod=JukesCantor` default, O(n²) space, the cluster-index dictionary detail, and the §5.4 accepted-deviation / §6.2 builder limitations);
Phylogenetics/Newick_Format → [[distance-based-tree-construction]] resolved 2026-07-15 (REUSE: the Newick I/O layer serializes this hinge concept's `PhyloNode` output — a format serializer, not a separate algorithm, so no dedicated concept, per [[phylo-newick-001-evidence]]);
Phylogenetics/Distance_Matrix → [[evolutionary-distance-matrix]] resolved 2026-07-15 (reconciled against the existing Evidence-derived concept);
Phylogenetics/Tree_Comparison → [[tree-comparison-metrics]] resolved 2026-07-15 (reconciled against the existing Evidence-derived concept; enriched with the O(n² log n) impl cost, helper methods, and test-pinned oracles);
Phylogenetics/Tree_Statistics → [[tree-statistics]] resolved 2026-07-15, closing the Phylogenetics domain (REUSE: reconciled against the existing Evidence-derived concept; enriched with the `PhylogeneticAnalyzer` method set, the O(n)/O(h) recurrences, the `yield` pre-order leaf order, the `EmptyTreeHeight` constant, and the topological-vs-Biopython-`depths()` distinction);
Pattern_Matching/Suffix_Tree → [[suffix-tree]] resolved 2026-07-15, closing the Pattern_Matching domain (last pending doc);
PanGenome/Phylogenetic_Marker_Selection → [[phylogenetic-marker-selection]] resolved 2026-07-15, closing the PanGenome domain;
Oncology/Variant_Allele_Frequency → [[variant-allele-frequency-and-binomial-ci]] resolved 2026-07-15, closing the Oncology domain;
Oncology/Tumor_Purity_Estimation → [[tumor-purity-from-mutation-vaf]] resolved 2026-07-15;
Oncology/Tumor_Ploidy_Estimation → [[tumor-ploidy-estimation-and-whole-genome-doubling]] resolved 2026-07-15;
Oncology/Tumor_Phylogeny_Reconstruction → [[tumor-phylogeny-clonal-tree-reconstruction]] resolved 2026-07-15;
Oncology/Tumor_Mutational_Burden → [[tumor-mutational-burden]] resolved 2026-07-15;
Oncology/Tumor_Heterogeneity_Analysis → [[intratumor-heterogeneity-metrics]] resolved 2026-07-15;
Oncology/Tumor_Gene_Expression_Outlier → [[expression-outlier-zscore-signature-score]] resolved 2026-07-15;
Oncology/Sequencing_Artifact_Detection → [[sequencing-artifact-detection]] resolved 2026-07-15;
Oncology/SBS96_Trinucleotide_Context_Catalog → [[sbs96-mutational-signature-catalog]] resolved 2026-07-14;
Oncology/Neoantigen_Peptide_Generation → [[neoantigen-peptide-generation]] resolved 2026-07-14;
Oncology/Mutational_Signature_Fitting → [[mutational-signature-fitting-and-extraction]] resolved 2026-07-14;
Oncology/Mutational_Signature_Extraction_NMF → [[mutational-signature-fitting-and-extraction]] resolved 2026-07-14;
Oncology/Mutational_Signature_Exposure_Bootstrap → [[signature-exposure-bootstrap-confidence-intervals]] resolved 2026-07-14;
Oncology/Mutational_Process_Classification → [[mutational-process-classification]] resolved 2026-07-14;
Oncology/Microsatellite_Instability_Detection → [[microsatellite-instability-detection]] resolved 2026-07-14;
Oncology/MRD_Detection → [[tumor-informed-mrd-detection]] resolved 2026-07-14;
Oncology/MHC_Peptide_Binding_Classification → [[mhc-peptide-binding-prediction]] resolved 2026-07-14;
Oncology/Loss_Of_Heterozygosity → [[loss-of-heterozygosity-detection]] resolved 2026-07-14;
Oncology/Known_Fusion_Database_Lookup → [[gene-fusion-nomenclature-known-fusion-lookup]] resolved 2026-07-14;
Oncology/Immune_Infiltration_Estimation → [[immune-infiltration-deconvolution]] resolved 2026-07-14;
Oncology/Homozygous_Deletion_Detection → [[homozygous-deletion-detection]] resolved 2026-07-14;
Oncology/HRD_Score → [[homologous-recombination-deficiency-score]] resolved 2026-07-14;
Oncology/HLA_Nomenclature_And_Allele_Specific_LOH → [[hla-nomenclature-and-allele-specific-loh]] resolved 2026-07-14;
Oncology/Fusion_Gene_Detection → [[gene-fusion-detection-read-evidence]] resolved 2026-07-14;
Oncology/Fusion_Breakpoint_Analysis → [[fusion-breakpoint-frame-and-protein-prediction]] resolved 2026-07-14;
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
| `docs/algorithms/PanGenome/Gene_Clustering.md` | [[pan-genome-gene-clustering]] |
| `docs/algorithms/PanGenome/Pan_Genome_Growth_Model.md` | [[pan-genome-heaps-law-fit]] |
| `docs/algorithms/PanGenome/Phylogenetic_Marker_Selection.md` | [[phylogenetic-marker-selection]] |
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
| `docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md` | [[approximate-pattern-matching-mismatches]] |
| `docs/algorithms/Pattern_Matching/Consensus_From_Alignment.md` | [[consensus-from-alignment]] |
| `docs/algorithms/Pattern_Matching/Frequent_Words_With_Mismatches.md` | [[approximate-pattern-matching-mismatches]] |
| `docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md` | [[exact-pattern-search]] |
| `docs/algorithms/Pattern_Matching/Edit_Distance.md` | [[edit-distance]] |
| `docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Consensus.md` | [[iupac-degenerate-consensus]] |
| `docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md` | [[iupac-degenerate-matching]] |
| `docs/algorithms/Pattern_Matching/Position_Weight_Matrix.md` | [[position-weight-matrix]] |
| `docs/algorithms/Pattern_Matching/Suffix_Tree.md` | [[suffix-tree]] |
| `docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md` | [[coiled-coil-prediction]] |
| `docs/algorithms/ProteinMotif/Common_Motif_Finding.md` | [[common-protein-motifs]] |
| `docs/algorithms/ProteinMotif/Domain_Prediction.md` | [[protein-domain-and-signal-peptide-prediction]] |
| `docs/algorithms/ProteinMotif/Low_Complexity_Region_Detection.md` | [[protein-low-complexity-seg]] |
| `docs/algorithms/ProteinMotif/Motif_Search.md` | [[protein-motif-pattern-search]] |
| `docs/algorithms/ProteinMotif/PROSITE_Pattern_Matching.md` | [[protein-motif-pattern-search]] |
| `docs/algorithms/ProteinMotif/Pattern_Matching_Methods.md` | [[protein-motif-pattern-search]] |
| `docs/algorithms/ProteinMotif/Profile_HMM_Domain_Detection.md` | [[protein-domain-and-signal-peptide-prediction]] |
| `docs/algorithms/ProteinMotif/Signal_Peptide_Prediction.md` | [[protein-domain-and-signal-peptide-prediction]] |
| `docs/algorithms/ProteinMotif/Transmembrane_Helix_Prediction.md` | [[transmembrane-helix-prediction]] |
| `docs/algorithms/ProteinPred/Disorder_Prediction.md` | [[intrinsic-disorder-prediction-top-idp]] |
| `docs/algorithms/ProteinPred/Disorder_Propensity.md` | [[intrinsic-disorder-prediction-top-idp]] |
| `docs/algorithms/ProteinPred/Disordered_Region_Detection.md` | [[intrinsic-disorder-prediction-top-idp]] |
| `docs/algorithms/ProteinPred/Low_Complexity_Region_Detection.md` | [[protein-low-complexity-seg]] |
| `docs/algorithms/ProteinPred/MoRF_Prediction.md` | [[morf-prediction-dip-in-disorder]] |
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
| `docs/algorithms/Oncology/Fusion_Breakpoint_Analysis.md` | [[fusion-breakpoint-frame-and-protein-prediction]] |
| `docs/algorithms/Oncology/Fusion_Gene_Detection.md` | [[gene-fusion-detection-read-evidence]] |
| `docs/algorithms/Oncology/HLA_Nomenclature_And_Allele_Specific_LOH.md` | [[hla-nomenclature-and-allele-specific-loh]] |
| `docs/algorithms/Oncology/HRD_Score.md` | [[homologous-recombination-deficiency-score]] |
| `docs/algorithms/Oncology/Homozygous_Deletion_Detection.md` | [[homozygous-deletion-detection]] |
| `docs/algorithms/Oncology/Immune_Infiltration_Estimation.md` | [[immune-infiltration-deconvolution]] |
| `docs/algorithms/Oncology/Known_Fusion_Database_Lookup.md` | [[gene-fusion-nomenclature-known-fusion-lookup]] |
| `docs/algorithms/Oncology/Loss_Of_Heterozygosity.md` | [[loss-of-heterozygosity-detection]] |
| `docs/algorithms/Oncology/MHC_Peptide_Binding_Classification.md` | [[mhc-peptide-binding-prediction]] |
| `docs/algorithms/Oncology/MRD_Detection.md` | [[tumor-informed-mrd-detection]] |
| `docs/algorithms/Oncology/Microsatellite_Instability_Detection.md` | [[microsatellite-instability-detection]] |
| `docs/algorithms/Oncology/Mutational_Process_Classification.md` | [[mutational-process-classification]] |
| `docs/algorithms/Oncology/Mutational_Signature_Exposure_Bootstrap.md` | [[signature-exposure-bootstrap-confidence-intervals]] |
| `docs/algorithms/Oncology/Mutational_Signature_Extraction_NMF.md` | [[mutational-signature-fitting-and-extraction]] |
| `docs/algorithms/Oncology/Mutational_Signature_Fitting.md` | [[mutational-signature-fitting-and-extraction]] |
| `docs/algorithms/Oncology/Neoantigen_Peptide_Generation.md` | [[neoantigen-peptide-generation]] |
| `docs/algorithms/Oncology/SBS96_Trinucleotide_Context_Catalog.md` | [[sbs96-mutational-signature-catalog]] |
| `docs/algorithms/Oncology/Sequencing_Artifact_Detection.md` | [[sequencing-artifact-detection]] |
| `docs/algorithms/Oncology/Somatic_Mutation_Calling.md` | [[somatic-variant-calling-tumor-normal]] |
| `docs/algorithms/Oncology/Tumor_Gene_Expression_Outlier.md` | [[expression-outlier-zscore-signature-score]] |
| `docs/algorithms/Oncology/Tumor_Heterogeneity_Analysis.md` | [[intratumor-heterogeneity-metrics]] |
| `docs/algorithms/Oncology/Tumor_Mutational_Burden.md` | [[tumor-mutational-burden]] |
| `docs/algorithms/Oncology/Tumor_Phylogeny_Reconstruction.md` | [[tumor-phylogeny-clonal-tree-reconstruction]] |
| `docs/algorithms/Oncology/Tumor_Ploidy_Estimation.md` | [[tumor-ploidy-estimation-and-whole-genome-doubling]] |
| `docs/algorithms/Oncology/Tumor_Purity_Estimation.md` | [[tumor-purity-from-mutation-vaf]] |
| `docs/algorithms/Oncology/Variant_Allele_Frequency.md` | [[variant-allele-frequency-and-binomial-ci]] |
| `docs/algorithms/Phylogenetics/Bootstrap_Analysis.md` | [[phylogenetic-bootstrap-support]] |
| `docs/algorithms/Phylogenetics/Distance_Matrix.md` | [[evolutionary-distance-matrix]] |
| `docs/algorithms/Phylogenetics/Newick_Format.md` | [[distance-based-tree-construction]] |
| `docs/algorithms/Phylogenetics/Tree_Comparison.md` | [[tree-comparison-metrics]] |
| `docs/algorithms/Phylogenetics/Tree_Construction.md` | [[distance-based-tree-construction]] |
| `docs/algorithms/Phylogenetics/Tree_Statistics.md` | [[tree-statistics]] |
| `docs/algorithms/Population_Genetics/Allele_Frequency.md` | [[allele-genotype-frequencies]] |
| `docs/algorithms/Population_Genetics/Ancestry_Estimation.md` | [[ancestry-estimation-admixture]] |
| `docs/algorithms/Population_Genetics/Diversity_Statistics.md` | [[genetic-diversity-statistics]] |
| `docs/algorithms/Population_Genetics/F_Statistics.md` | [[population-differentiation-fst]] |
| `docs/algorithms/Population_Genetics/Hardy_Weinberg_Test.md` | [[hardy-weinberg-equilibrium-test]] |
| `docs/algorithms/Population_Genetics/Integrated_Haplotype_Score.md` | [[selection-scan-ihs-ehh]] |
| `docs/algorithms/Population_Genetics/Linkage_Disequilibrium.md` | [[linkage-disequilibrium]] |
| `docs/algorithms/Population_Genetics/Runs_Of_Homozygosity.md` | [[runs-of-homozygosity-inbreeding]] |
| `docs/algorithms/Quality/Phred_Score_Handling.md` | [[phred-quality-encoding]] |
| `docs/algorithms/Quality/Quality_Statistics.md` | [[fastq-quality-statistics]] |

## Pending (fold into the ingest campaign)

The per-domain pending tables (53 algorithm docs across 9 domains, no concept page yet) live in **[[backlog-pending]]** to keep this hub under the page-size cap. A pending row is resolved when a concept page lists the algorithm doc in `sources:`, at which point it moves to the *Covered via concept* table above.

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
  `README.md` is **resolved index-only (2026-07-16)** — a coverage exclusion, not a
  wiki page. It is purely navigational (section headers linking to the per-algorithm
  folders, every one of which is already synthesized by a concept in the *Covered via
  concept* table) plus a Canonicalization section whose content is already owned by
  [[canonical-algorithm-map]]. No distinct synthesis to capture, so no dedicated page.
- The `docs/Evidence/**` campaign (175 of 213 remaining) is the primary driver: each
  Evidence ingest typically creates or extends the concept that also covers the
  matching algorithm doc, clearing a pending row here as a side effect.
