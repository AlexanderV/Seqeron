---
type: concept
title: "MCP tool catalog — the thin wrapper surface"
tags: [mcp, reference]
sources:
  - docs/mcp/README.md
  - docs/mcp/MCP_STATUS.md
source_commit: b3f950caf701615bb8a0296df6c5368d26dde7ec
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:three-front-doors
      source: mcp-readme
      evidence: "README frames the MCP tool call as 'just a different front door' onto the same validated algorithm engine (one of the three front doors)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:skill-layer
      source: mcp-readme
      evidence: "README: the skill layer attaches the right MCP server on demand and keeps the 427 tool schemas out of the model's context."
      confidence: high
      status: current
---

# MCP tool catalog — the thin wrapper surface

The **427 MCP tools across 11 servers** are not new algorithms. Each is a thin
[Model Context Protocol](https://modelcontextprotocol.io/) wrapper that validates a JSON-schema
input, delegates to exactly one `Seqeron.Genomics` (or `SuffixTree`) library method, and returns a
structured record. As `docs/mcp/README.md` puts it: *"Same math as the C# API; the tool call is just
a different front door."* This page is the **map from every tool to the concept its wrapped
algorithm belongs to** — it replaces per-tool ingestion. It does not restate schemas or algorithm
behaviour; each tool ships its own reference doc at `docs/mcp/tools/<server>/<tool>.md`, and every
tool has a gold-standard binding + Schema+Binding test per [[mcp-checklist]] / [[mcp-plan]] (campaign
COMPLETE 2026-07-01, see `docs/mcp/MCP_STATUS.md`).

This surface is one of the [[three-front-doors]] (skills, C# API, MCP over one engine) and sits
beneath the [[skill-layer]], which attaches the right server on demand so the 427
schemas stay out of the model's context. The front-door overview is [[mcp-readme]]; the design
standards are [[mcp-plan]] and [[mcp-checklist]]; the reference operating prompt is [[mcp-prompt]].

## Coverage

Deterministic mapping (`.claude/skills/wiki-ingest-doc/scripts/mcp_map.py --all --catalog`) binds **214 / 427**
tools to an existing concept (via authoritative method-ID / algorithm-doc / literal bridges, or
high-confidence name overlap), touching **124 distinct concepts** (the `parse_gff3` / `to_gff3`
pair was re-homed from [[bed-format-parsing]] onto the dedicated [[gff3-io]] concept; the
`GcSkewCalculator` pair `gc_skew` / `windowed_gc_skew` was mapped onto the new [[gc-skew]]
concept; the scalar linguistic-complexity pair `complexity_linguistic` / `linguistic_complexity`
was re-homed from [[windowed-sequence-complexity-profile]] onto the new [[linguistic-complexity]]
concept, leaving the windowed profile with just `windowed_complexity`; the `SequenceExtensions`
predicate pair `is_valid_dna` / `is_valid_rna` was mapped onto the new [[sequence-validation]]
concept). The remaining **213** are listed per server as *unmapped* — a tool with no confidently
matched existing concept, not necessarily a missing algorithm. No new concept page was created; the
unmapped tools are recorded as gaps (see the closing section) so a future ingest can decide whether a
cluster deserves its own concept.

Reading the rows: `[[concept-slug]] — tool_a, tool_b` means those tools all wrap methods that belong
to that concept (one binding per library Method ID). Unmapped tools are grouped by their wrapped C#
class so algorithm families that lack a concept are visible.

## Tool surface by server

### Sequence — 35 tools (17 mapped, 18 unmapped)

`Seqeron.Mcp.Sequence` — DNA/RNA/Protein models, composition, complexity, k-mers, Tm.

- [[crispr-guide-rna-design]] — `iupac_matches`
- [[dust-low-complexity-score]] — `complexity_dust_score`, `complexity_mask_low`
- [[genetic-code-translation]] — `translate_dna`, `translate_rna`
- [[isoelectric-point]] — `isoelectric_point`
- [[iupac-degenerate-consensus]] — `iupac_code`
- [[k-mer-euclidean-distance]] — `kmer_distance`
- [[k-mer-statistics]] — `kmer_entropy`
- [[molecular-weight]] — `molecular_weight_nucleotide`, `molecular_weight_protein`
- [[nucleotide-composition-skew]] — `nucleotide_composition`
- [[sequence-complexity-compression-lempel-ziv]] — `complexity_compression_ratio`
- [[linguistic-complexity]] — `complexity_linguistic`, `linguistic_complexity`
- [[sequence-validation]] — `is_valid_dna`, `is_valid_rna`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _DnaSequence_: `dna_reverse_complement`, `dna_validate`
  - _IupacDnaSequence_: `iupac_match`
  - _KmerAnalyzer_: `kmer_analyze`, `kmer_count`
  - _ProteinSequence_: `protein_validate`
  - _RnaSequence_: `rna_from_dna`, `rna_validate`
  - _SequenceComplexity_: `complexity_kmer_entropy`, `complexity_shannon`
  - _SequenceExtensions_: `complement_base`, `gc_content`
  - _SequenceStatistics_: `amino_acid_composition`, `hydrophobicity`, `melting_temperature`, `shannon_entropy`, `summarize_sequence`, `thermodynamics`

### Parsers — 41 tools (2 mapped, 39 unmapped)

`Seqeron.Mcp.Parsers` — FASTA/FASTQ/GenBank/GFF/VCF/BED/EMBL parsing & utilities.

- [[fastq-quality-statistics]] — `fastq_filter`, `fastq_statistics`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _BedParser_: `bed_filter`, `bed_intersect`, `bed_merge`, `bed_parse`
  - _EmblParser_: `embl_features`, `embl_parse`, `embl_statistics`
  - _FastaParser_: `fasta_format`, `fasta_parse`, `fasta_write`
  - _FastqParser_: `fastq_detect_encoding`, `fastq_encode_quality`, `fastq_error_to_phred`, `fastq_format`, `fastq_parse`, `fastq_phred_to_error`, `fastq_trim_adapter`, `fastq_trim_quality`, `fastq_write`
  - _GenBankParser_: `genbank_extract_sequence`, `genbank_features`, `genbank_parse`, `genbank_parse_location`, `genbank_statistics`
  - _GffParser_: `gff_filter`, `gff_parse`, `gff_statistics`
  - _VcfParser_: `vcf_classify`, `vcf_filter`, `vcf_has_flag`, `vcf_is_het`, `vcf_is_hom_alt`, `vcf_is_hom_ref`, `vcf_is_indel`, `vcf_is_snp`, `vcf_parse`, `vcf_statistics`, `vcf_variant_length`, `vcf_write`

### Alignment — 22 tools (18 mapped, 4 unmapped)

`Seqeron.Mcp.Alignment` — pairwise & multiple sequence alignment, assembly graph.

- [[alignment-statistics]] — `alignment_statistics`, `format_alignment`
- [[approximate-pattern-matching-mismatches]] — `find_with_mismatches`
- [[assembly-statistics]] — `assembly_stats`
- [[consensus-sequence]] — `compute_consensus`
- [[contig-merge-overlap-collapse]] — `merge_contigs`
- [[coverage-depth-calculation]] — `calculate_coverage`
- [[de-bruijn-graph-assembly]] — `assemble_de_bruijn`
- [[global-alignment-needleman-wunsch]] — `global_align`
- [[kmer-spectrum-error-correction]] — `error_correct_reads`
- [[local-alignment-smith-waterman]] — `local_align`
- [[multiple-sequence-alignment]] — `multiple_align`
- [[overlap-layout-consensus-assembly]] — `assemble_olc`, `find_all_overlaps`, `find_overlap`
- [[quality-trimming-running-sum]] — `quality_trim_reads`
- [[scaffolding]] — `scaffold_contigs`
- [[semi-global-alignment-fitting]] — `semi_global_align`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _ApproximateMatcher_: `find_best_match`, `find_with_edits`, `frequent_kmers_with_mismatches`
  - _SequenceAssembler_: `sequence_identity`

### Analysis — 91 tools (62 mapped, 29 unmapped)

`Seqeron.Mcp.Analysis` — k-mer, motif, repeat, complexity, comparative & structural genomics.

- [[algorithm-validation-evidence]] — `parse_dot_bracket`
- [[average-nucleotide-identity]] — `calculate_ani`
- [[base-composition]] — `analyze_gc_content`
- [[both-strand-kmer-counting]] — `count_kmers_both_strands`
- [[k-mer-counting]] — `count_kmers`
- [[codon-usage-comparison]] — `codon_frequencies`
- [[coiled-coil-prediction]] — `predict_coiled_coils`
- [[common-protein-motifs]] — `find_protein_motifs`
- [[conserved-gene-clusters-common-intervals]] — `find_conserved_clusters`
- [[dot-plot-word-match]] — `generate_dot_plot`
- [[dust-low-complexity-score]] — `dust_score`, `mask_low_complexity`
- [[genome-comparison-core-dispensable]] — `compare_genomes`
- [[genome-rearrangement-breakpoint-distance]] — `detect_rearrangements`, `reversal_distance`
- [[hydrophobicity-gravy-and-profile]] — `hydrophobicity_profile`
- [[intrinsic-disorder-prediction-top-idp]] — `disorder_propensity`, `is_disorder_promoting`, `predict_disorder`
- [[k-mer-euclidean-distance]] — `kmer_distance`, `kmer_frequencies`
- [[k-mer-generation]] — `generate_all_kmers`
- [[k-mer-positions]] — `kmer_positions`
- [[kmer-spectrum-error-correction]] — `kmer_spectrum`
- [[known-motif-search]] — `find_exact_motif`, `find_known_motifs`
- [[longest-common-substring]] — `find_common_regions`
- [[longest-repeated-substring]] — `find_repeats`
- [[nucleotide-composition-skew]] — `at_skew`
- [[open-reading-frame-detection]] — `find_open_reading_frames`
- [[ortholog-detection-reciprocal-best-hits]] — `find_orthologs`, `find_reciprocal_best_hits`
- [[overrepresented-kmer-discovery]] — `discover_motifs`
- [[protein-domain-and-signal-peptide-prediction]] — `find_protein_domains`, `predict_signal_peptide`
- [[protein-low-complexity-seg]] — `find_protein_low_complexity_regions`, `predict_low_complexity_seg`
- [[protein-motif-pattern-search]] — `find_motif_by_pattern`, `find_motif_by_prosite`, `prosite_to_regex`
- [[protein-secondary-structure-chou-fasman]] — `predict_chou_fasman`
- [[regulatory-element-detection]] — `create_pwm`, `find_degenerate_motif`, `find_regulatory_elements`, `scan_with_pwm`
- [[repetitive-element-detection]] — `find_tandem_repeats`, `tandem_repeat_summary`
- [[gc-skew]] — `gc_skew`, `windowed_gc_skew`
- [[replication-origin-cumulative-skew]] — `cumulative_gc_skew`, `predict_replication_origin`
- [[rna-dot-bracket-notation]] — `validate_dot_bracket`
- [[rna-free-energy-turner-model]] — `terminal_mismatch_energy`
- [[rna-minimum-free-energy-folding]] — `minimum_free_energy`
- [[rna-pseudoknot-detection]] — `detect_pseudoknots`
- [[rna-stem-loop-enumeration]] — `find_stem_loops`, `hairpin_loop_energy`
- [[shared-motifs]] — `find_shared_motifs`
- [[synteny-and-rearrangement-detection]] — `find_syntenic_blocks`
- [[unique-and-mincount-kmers]] — `kmers_with_min_count`, `unique_kmers`
- [[windowed-sequence-complexity-profile]] — `windowed_complexity`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _DisorderPredictor_: `predict_morfs`
  - _GenomicAnalyzer_: `find_motif`
  - _KmerAnalyzer_: `analyze_kmers`, `find_clumps`, `most_frequent_kmers`
  - _MotifFinder_: `generate_consensus`
  - _ProteinMotifFinder_: `predict_transmembrane_helices`
  - _RepeatFinder_: `find_direct_repeats`, `find_inverted_repeats`, `find_microsatellites`, `find_palindromes`
  - _RnaSecondaryStructure_: `base_pair_type`, `bulge_loop_energy`, `can_pair`, `dangling_end_energy`, `find_rna_inverted_repeats`, `flush_coaxial_stacking`, `internal_loop_energy`, `mismatch_coaxial_stacking`, `multibranch_loop_energy`, `predict_rna_structure`, `rna_complement_base`, `stem_energy`
  - _SequenceComplexity_: `compression_ratio`, `find_low_complexity_regions`
  - _SequenceStatistics_: `dinucleotide_frequencies`, `dinucleotide_ratios`, `entropy_profile`, `gc_content_profile`

### Annotation — 97 tools (40 mapped, 57 unmapped)

`Seqeron.Mcp.Annotation` — genes/ORFs/promoters, variants, epigenetics, miRNA, splicing, SVs, transcriptomics.

- [[alternative-splicing-psi]] — `detect_alternative_splicing`
- [[bisulfite-methylation-calling]] — `methylation_from_bisulfite`, `simulate_bisulfite_conversion`
- [[chromatin-state-prediction]] — `annotate_histone_modifications`, `find_accessible_regions`, `predict_chromatin_state`
- [[coding-potential-hexamer-score]] — `coding_potential`
- [[cpg-island-detection]] — `cpg_observed_expected`, `find_cpg_islands`, `find_cpg_sites`
- [[differential-expression]] — `differential_expression`
- [[differentially-methylated-regions]] — `find_dmrs`
- [[discordant-pair-sv-detection]] — `find_discordant_pairs`
- [[epigenetic-age-horvath-clock]] — `epigenetic_age`
- [[expression-quantification]] — `calculate_tpm`
- [[gff3-io]] — `parse_gff3`, `to_gff3`
- [[gene-structure-prediction-intron-exon]] — `predict_gene_structure`, `predict_introns`
- [[germline-variant-calling-snp-indel]] — `find_indels`, `find_snps`, `titv_ratio`
- [[methylation-context-classification]] — `find_methylation_sites`
- [[mirna-target-site-prediction]] — `analyze_target_context`, `site_accessibility`
- [[pre-mirna-hairpin-detection]] — `find_pre_mirna_hairpins`
- [[prokaryotic-gene-prediction-rbs]] — `predict_genes`
- [[promoter-detection]] — `find_promoter_motifs`
- [[read-depth-cnv-segmentation]] — `segment_copy_number`
- [[regulatory-element-detection]] — `annotate_regulatory_elements`
- [[relative-synonymous-codon-usage]] — `codon_usage`
- [[repetitive-element-detection]] — `find_repetitive_elements`
- [[rna-base-pairing]] — `can_pair`, `is_wobble_pair`
- [[seed-sequence-analysis]] — `compare_seed_regions`, `group_by_seed_family`, `mirna_seed_sequence`
- [[splice-acceptor-site-prediction]] — `find_branch_points`
- [[variant-effect-annotation-vep]] — `annotate_variants`, `predict_variant_effect`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _EpigeneticsAnalyzer_: `methylation_profile`, `predict_imprinted_genes`
  - _GenomeAnnotator_: `find_orfs`, `find_ribosome_binding_sites`, `longest_orfs_per_frame`
  - _MiRnaAnalyzer_: `align_mirna_to_target`, `create_mirna`, `find_mirna_target_sites`, `find_similar_mirnas`, `generate_seed_variants`, `rna_reverse_complement`
  - _SpliceSitePredictor_: `find_acceptor_sites`, `find_donor_sites`, `find_retained_intron_candidates`, `is_within_coding_region`, `maxent_score`
  - _StructuralVariantAnalyzer_: `annotate_svs`, `assemble_breakpoint_sequence`, `cluster_discordant_pairs`, `cluster_split_reads`, `filter_svs`, `find_microhomology`, `find_split_reads`, `genotype_sv`, `identify_cnvs`, `merge_overlapping_svs`
  - _TranscriptomeAnalyzer_: `build_coexpression_network`, `cluster_genes_by_expression`, `detect_differential_splicing`, `detect_isoform_switching`, `enrichment_score`, `find_dominant_isoforms`, `find_skipped_exon_events`, `log2_transform`, `over_representation_analysis`, `pearson_correlation`, `perform_pca`, `quantile_normalize`, `rnaseq_quality_metrics`
  - _VariantAnnotator_: `annotate_variant_on_transcripts`, `calculate_conservation`, `classify_variant`, `find_conserved_elements`, `format_vcf_info`, `impact_level`, `normalize_variant`, `parse_vcf_variant`, `predict_pathogenicity`, `predict_tf_binding_change`
  - _VariantCaller_: `call_variants`, `call_variants_from_alignment`, `classify_mutation`, `find_deletions`, `find_insertions`, `find_snps_direct`, `variant_statistics`, `variants_to_vcf`

### Phylogenetics — 13 tools (8 mapped, 5 unmapped)

`Seqeron.Mcp.Phylogenetics` — distances, tree building, phylogenetic statistics.

- [[evolutionary-distance-matrix]] — `build_tree_from_matrix`, `distance_matrix`
- [[mirna-target-site-prediction]] — `parse_newick`
- [[phylogenetic-bootstrap-support]] — `bootstrap_support`
- [[tree-comparison-metrics]] — `mrca`, `patristic_distance`, `robinson_foulds_distance`
- [[tree-statistics]] — `tree_depth`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _PhylogeneticAnalyzer_: `build_phylogenetic_tree`, `pairwise_distance`, `to_newick`, `tree_leaves`, `tree_length`

### Population — 18 tools (13 mapped, 5 unmapped)

`Seqeron.Mcp.Population` — population genetics (Fst, diversity, LD, selection).

- [[allele-genotype-frequencies]] — `allele_frequencies`, `filter_variants_by_maf`
- [[ancestry-estimation-admixture]] — `estimate_ancestry`
- [[genetic-diversity-statistics]] — `diversity_statistics`, `tajimas_d`
- [[hardy-weinberg-equilibrium-test]] — `hardy_weinberg_test`
- [[linkage-disequilibrium]] — `haplotype_blocks`, `linkage_disequilibrium`
- [[population-differentiation-fst]] — `fst`, `pairwise_fst`
- [[runs-of-homozygosity-inbreeding]] — `inbreeding_from_roh`, `runs_of_homozygosity`
- [[selection-scan-ihs-ehh]] — `integrated_haplotype_score`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _PopulationGeneticsAnalyzer_: `f_statistics`, `minor_allele_frequency`, `nucleotide_diversity`, `scan_selection_signals`, `wattersons_theta`

### Metagenomics — 19 tools (10 mapped, 9 unmapped)

`Seqeron.Mcp.Metagenomics` — taxonomic classification & community profiling.

- [[alpha-diversity]] — `alpha_diversity`
- [[beta-diversity]] — `beta_diversity`
- [[functional-prediction]] — `predict_functions`
- [[metagenomic-binning]] — `bin_contigs`
- [[pan-genome-heaps-law-fit]] — `fit_heaps_law`
- [[phylogenetic-marker-selection]] — `select_phylogenetic_markers`
- [[significant-taxa-detection]] — `differential_abundance`
- [[taxonomic-classification]] — `build_kmer_database`, `classify_reads`
- [[taxonomic-profile]] — `taxonomic_profile`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _MetagenomicsAnalyzer_: `find_resistance_genes`, `functional_diversity`
  - _PanGenomeAnalyzer_: `accessory_genes`, `cluster_genes`, `construct_pangenome`, `core_gene_clusters`, `core_genome_alignment`, `find_genome_specific_genes`, `gene_presence_absence_matrix`

### Chromosome — 32 tools (17 mapped, 15 unmapped)

`Seqeron.Mcp.Chromosome` — chromosome-scale analysis (karyotype, centromere, synteny).

- [[aneuploidy-detection]] — `detect_aneuploidy`, `identify_whole_chromosome_aneuploidy`
- [[assembly-statistics]] — `assembly_statistics`, `au_n`, `find_gaps`, `nx_statistics`
- [[centromere-analysis]] — `analyze_centromere`, `arm_ratio`, `classify_chromosome_by_arm_ratio`
- [[karyotype-analysis]] — `analyze_karyotype`, `detect_ploidy`
- [[repetitive-element-detection]] — `find_repetitive_regions`
- [[scaffolding]] — `analyze_scaffolds`
- [[synteny-and-rearrangement-detection]] — `detect_rearrangements`, `find_synteny_blocks`
- [[telomere-analysis]] — `analyze_telomeres`, `estimate_telomere_length_from_ts_ratio`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _ChromosomeAnalyzer_: `estimate_cell_divisions_from_telomere_length`, `find_heterochromatin_regions`, `predict_g_bands`
  - _GenomeAssemblyAnalyzer_: `assess_completeness`, `compare_assemblies`, `estimate_completeness_from_kmers`, `extract_contigs`, `find_suspicious_regions`, `find_syntenic_blocks_assemblies`, `find_tandem_repeats`, `gap_distribution`, `length_distribution`, `local_quality`, `nx_curve`, `repeat_content`

### MolTools — 47 tools (22 mapped, 25 unmapped)

`Seqeron.Mcp.MolTools` — primer/probe/CRISPR design, codon optimization, restriction.

- [[codon-adaptation-index]] — `cai_from_organism_table`, `codon_adaptation_index`
- [[codon-optimization]] — `optimize_codons`, `reduce_secondary_structure`
- [[codon-usage-comparison]] — `compare_codon_usage`
- [[crispr-guide-rna-design]] — `crispr_system_info`, `design_guide_rnas`, `evaluate_guide_rna`, `find_pam_sites`
- [[effective-number-of-codons]] — `effective_number_of_codons`
- [[primer-dimer-thermodynamics-tm]] — `primer_dimer`, `primer_melting_temperature`, `primer_melting_temperature_salt`
- [[rare-codon-analysis]] — `find_rare_codons`
- [[relative-synonymous-codon-usage]] — `count_codons`, `rscu`
- [[restriction-enzyme-filtering]] — `find_all_restriction_sites`, `find_restriction_sites`, `get_enzyme`
- [[taqman-probe-design-rules]] — `design_antisense_probes`, `design_probes`, `design_tiling_probes`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _CodonOptimizer_: `build_codon_table`, `remove_restriction_sites`
  - _CodonUsageAnalyzer_: `codon_usage_statistics`
  - _CrisprDesigner_: `crispr_specificity_score`, `find_off_targets`
  - _PrimerDesigner_: `design_primers`, `evaluate_primer`, `generate_primer_candidates`, `hairpin_potential`, `longest_dinucleotide_repeat`, `longest_homopolymer`, `three_prime_stability`
  - _ProbeDesigner_: `analyze_oligo`, `design_molecular_beacon`, `oligo_concentration_from_absorbance`, `oligo_extinction_coefficient`, `validate_probe`
  - _RestrictionAnalyzer_: `blunt_cutters`, `compatible_enzymes`, `digest_summary`, `enzymes_by_cut_length`, `enzymes_compatible`, `restriction_digest`, `restriction_map`, `sticky_cutters`

### Core — 12 tools (6 mapped, 6 unmapped)

`SuffixTree.Mcp.Core` — suffix-tree search, edit/Hamming distance, k-mer similarity.

- [[kmer-jaccard-similarity]] — `calculate_similarity`
- [[longest-common-substring]] — `find_longest_common_region`, `suffix_tree_lcs`
- [[longest-repeated-substring]] — `find_longest_repeat`, `suffix_tree_find_all`, `suffix_tree_lrs`

Unmapped (no confidently-matched concept), grouped by wrapped class:
  - _ApproximateMatcher_: `count_approximate_occurrences`, `edit_distance`, `hamming_distance`
  - _SuffixTree_: `suffix_tree_contains`, `suffix_tree_count`, `suffix_tree_stats`

## Gaps & unmapped clusters (follow-ups)

The unmapped tools cluster around a handful of library classes with no synthesizing concept page yet. The largest clusters — each a candidate for a future concept ingest — are:

- `TranscriptomeAnalyzer` (Annotation) — 13 unmapped tools
- `RnaSecondaryStructure` (Analysis) — 12 unmapped tools
- `GenomeAssemblyAnalyzer` (Chromosome) — 12 unmapped tools
- `VcfParser` (Parsers) — 12 unmapped tools
- `SequenceStatistics` (Analysis) — 10 unmapped tools
- `StructuralVariantAnalyzer` (Annotation) — 10 unmapped tools
- `VariantAnnotator` (Annotation) — 10 unmapped tools
- `FastqParser` (Parsers) — 9 unmapped tools
- `VariantCaller` (Annotation) — 8 unmapped tools
- `RestrictionAnalyzer` (MolTools) — 8 unmapped tools
- `PanGenomeAnalyzer` (Metagenomics) — 7 unmapped tools
- `PrimerDesigner` (MolTools) — 7 unmapped tools
- `ApproximateMatcher` (Alignment) — 6 unmapped tools
- `KmerAnalyzer` (Analysis) — 6 unmapped tools
- `MiRnaAnalyzer` (Annotation) — 6 unmapped tools
- `SpliceSitePredictor` (Annotation) — 5 unmapped tools
- `ProbeDesigner` (MolTools) — 5 unmapped tools
- `GenBankParser` (Parsers) — 5 unmapped tools
- `PhylogeneticAnalyzer` (Phylogenetics) — 5 unmapped tools
- `PopulationGeneticsAnalyzer` (Population) — 5 unmapped tools
- `SequenceComplexity` (Analysis) — 4 unmapped tools
- `RepeatFinder` (Analysis) — 4 unmapped tools
- `GenomeAnnotator` (Annotation) — 4 unmapped tools
- `BedParser` (Parsers) — 4 unmapped tools
- `SequenceExtensions` (Sequence) — 4 unmapped tools

Smaller unmapped sets (RNA free-energy loop terms on `RnaSecondaryStructure`, the `Parsers` server's per-format parse/filter/statistics tools, individual popgen/phylo/chromosome utilities) are recorded inline above. Most are recorded as **gaps**, not missing algorithms: many wrap helper/query methods (`suffix_tree_contains`, `edit_distance`, `vcf_is_snp`) or fine-grained energy terms that a concept page already covers narratively even though the deterministic map could not bind them by method ID. Parsers (39/41 unmapped) is the clearest whole-server gap — file-format parsing has no concept pages yet.
