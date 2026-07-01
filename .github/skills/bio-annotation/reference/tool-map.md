# bio-annotation tool map — ~188 tools by family

Human index for the **Annotation** (97) + **Analysis** (91) servers. Grouped by workflow family.
Each row: `tool` · server · one-line purpose · `Method ID`. Open the per-tool doc for the full I/O
schema — **point, don't duplicate**: docs live at `docs/mcp/tools/{annotation,analysis}/<tool>.md`.

Servers: `A` = Annotation, `X` = Analysis. `⚠` = guarded / documented-limited (see
[`pipelines.md`](pipelines.md) → Envelope).

> This is a hand-curated router. The **canonical generated index** (all tools with links) is
> `_generated/tools.md`, produced by `scripts/skills/gen-catalog.py`. If a tool is missing here,
> use `seqeron-discovery` (`python3 scripts/skills/find-tool.py <kw> --server annotation|analysis`).

---

## 1. Structural annotation — ORFs / genes / promoters / RBS / GFF3 / codon

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `find_orfs` | A | ORFs across up to six frames | `GenomeAnnotator.FindOrfs` |
| `longest_orfs_per_frame` | A | Longest ORF per frame | `GenomeAnnotator.FindLongestOrfsPerFrame` |
| `find_open_reading_frames` | X | ORFs in all six frames (analysis path) | `GenomicAnalyzer.FindOpenReadingFrames` |
| `predict_genes` | A | ORF-based gene prediction | `GenomeAnnotator.PredictGenes` |
| `coding_potential` | A | CPAT hexamer coding-potential score | `GenomeAnnotator.CalculateCodingPotential` |
| `find_promoter_motifs` | A | −10 (Pribnow/TATA) + −35 boxes | `GenomeAnnotator.FindPromoterMotifs` |
| `find_ribosome_binding_sites` | A | Shine–Dalgarno upstream of ORFs | `GenomeAnnotator.FindRibosomeBindingSites` |
| `find_repetitive_elements` | A | Tandem + inverted repeats | `GenomeAnnotator.FindRepetitiveElements` |
| `codon_usage` | A | In-frame codon counts | `GenomeAnnotator.GetCodonUsage` |
| `codon_frequencies` | X | Codon usage frequencies | `SequenceStatistics.CalculateCodonFrequencies` |
| `parse_gff3` | A | Parse GFF3 → features | `GenomeAnnotator.ParseGff3` |
| `to_gff3` | A | Serialize gene annotations → GFF3 | `GenomeAnnotator.ToGff3` |

## 2. Variant calling + annotation + classification

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `call_variants` | A | SNPs+indels via global alignment | `VariantCaller.CallVariants` |
| `call_variants_from_alignment` | A | Variants from pre-aligned strings | `VariantCaller.CallVariantsFromAlignment` |
| `find_snps` | A | SNPs only (aligned) | `VariantCaller.FindSnps` |
| `find_snps_direct` | A | SNPs by direct positional compare | `VariantCaller.FindSnpsDirect` |
| `find_indels` | A | Indels (ins+del) | `VariantCaller.FindIndels` |
| `find_insertions` | A | Insertions only | `VariantCaller.FindInsertions` |
| `find_deletions` | A | Deletions only | `VariantCaller.FindDeletions` |
| `annotate_variants` | A | Call + annotate all variants (one shot) | `VariantCaller.AnnotateVariants` |
| `annotate_variant_on_transcripts` | A | VEP-like single-variant annotation | `VariantAnnotator.AnnotateVariant` |
| `classify_variant` | A | SNV/Ins/Del/MNV/Indel/Complex | `VariantAnnotator.ClassifyVariant` |
| `classify_mutation` | A | Transition/transversion/other | `VariantCaller.ClassifyMutation` |
| `normalize_variant` | A | Trim + classify | `VariantAnnotator.NormalizeVariant` |
| `predict_variant_effect` | A | Protein-level effect of a SNP in CDS | `VariantCaller.PredictEffect` |
| `predict_pathogenicity` | A | ACMG-like pathogenicity ⚠clinical | `VariantAnnotator.PredictPathogenicity` |
| `impact_level` | A | Consequence → High/Moderate/Low/Modifier | `VariantAnnotator.GetImpactLevel` |
| `calculate_conservation` | A | PhyloP/PhastCons/GERP-like scores | `VariantAnnotator.CalculateConservation` |
| `find_conserved_elements` | A | Conserved element runs from scores | `VariantAnnotator.FindConservedElements` |
| `annotate_regulatory_elements` | A | Regulatory regions overlapping a variant | `VariantAnnotator.AnnotateRegulatoryElements` |
| `predict_tf_binding_change` | A | TF-binding score change from a SNV | `VariantAnnotator.PredictTfBindingChange` |
| `titv_ratio` | A | Ti/Tv from a variant list | `VariantCaller.CalculateTiTvRatio` |
| `variant_statistics` | A | Totals + Ti/Tv + density | `VariantCaller.CalculateStatistics` |
| `parse_vcf_variant` | A | VCF fields → Variant record | `VariantAnnotator.ParseVcfVariant` |
| `format_vcf_info` | A | Annotation → VCF INFO string | `VariantAnnotator.FormatAsVcfInfo` |
| `variants_to_vcf` | A | Variants → VCF v4.2 lines | `VariantCaller.ToVcfLines` |

## 3. Structural variants / CNV / breakpoints

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `find_discordant_pairs` | A | Discordant read-pair SV signatures | `StructuralVariantAnalyzer.FindDiscordantPairs` |
| `cluster_discordant_pairs` | A | Cluster pairs → SV candidates | `StructuralVariantAnalyzer.ClusterDiscordantPairs` |
| `find_split_reads` | A | Split reads from soft-clipped CIGAR | `StructuralVariantAnalyzer.FindSplitReads` |
| `cluster_split_reads` | A | Cluster split reads → breakpoints | `StructuralVariantAnalyzer.ClusterSplitReads` |
| `assemble_breakpoint_sequence` | A | Assemble junction from clipped frags | `StructuralVariantAnalyzer.AssembleBreakpointSequence` |
| `find_microhomology` | A | Longest microhomology at a junction | `StructuralVariantAnalyzer.FindMicrohomology` |
| `identify_cnvs` | A | CN segments → del/dup SVs | `StructuralVariantAnalyzer.IdentifyCNVs` |
| `segment_copy_number` | A | Segment CN probes → plateaus | `StructuralVariantAnalyzer.SegmentCopyNumber` |
| `genotype_sv` | A | Genotype SV from support counts | `StructuralVariantAnalyzer.GenotypeSV` |
| `filter_svs` | A | Filter SVs by quality/support/length | `StructuralVariantAnalyzer.FilterSVs` |
| `merge_overlapping_svs` | A | Merge overlapping same-type SVs | `StructuralVariantAnalyzer.MergeOverlappingSVs` |
| `annotate_svs` | A | Annotate SVs w/ genes + impact | `StructuralVariantAnalyzer.AnnotateSVs` |
| `detect_rearrangements` | X | Breakpoints of signed gene-order perm | `ComparativeGenomics.DetectRearrangements` |

## 4. Motif discovery & scanning (DNA + protein)

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `discover_motifs` | X | De novo overrepresented k-mer motifs | `MotifFinder.DiscoverMotifs` |
| `find_exact_motif` | X | Exact motif positions (suffix tree) | `MotifFinder.FindExactMotif` |
| `find_motif` | X | Exact motif occurrences | `GenomicAnalyzer.FindMotif` |
| `find_degenerate_motif` | X | IUPAC-degenerate matches | `MotifFinder.FindDegenerateMotif` |
| `find_known_motifs` | X | Search a set of known motifs at once | `GenomicAnalyzer.FindKnownMotifs` |
| `find_shared_motifs` | X | k-mers shared across sequences | `MotifFinder.FindSharedMotifs` |
| `find_regulatory_elements` | X | Built-in regulatory motif scan | `MotifFinder.FindRegulatoryElements` |
| `generate_consensus` | X | IUPAC consensus from aligned seqs | `MotifFinder.GenerateConsensus` |
| `create_pwm` | X | Log-odds PWM from aligned DNA | `MotifFinder.CreatePwm` |
| `scan_with_pwm` | X | Scan a sequence with a PWM | `MotifFinder.ScanWithPwm` |
| `find_protein_motifs` | X | PROSITE-style protein motif catalog | `ProteinMotifFinder.FindCommonMotifs` |
| `find_motif_by_pattern` | X | Regex match in a protein | `ProteinMotifFinder.FindMotifByPattern` |
| `find_motif_by_prosite` | X | PROSITE-pattern match in a protein | `ProteinMotifFinder.FindMotifByProsite` |
| `find_protein_domains` | X | Domains via exact PROSITE patterns | `ProteinMotifFinder.FindDomains` |
| `prosite_to_regex` | X | PROSITE pattern → .NET regex | `ProteinMotifFinder.ConvertPrositeToRegex` |

## 5. Repeat analysis

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `find_tandem_repeats` | X | Consecutive repeating units | `GenomicAnalyzer.FindTandemRepeats` |
| `tandem_repeat_summary` | X | Aggregate microsatellite stats | `RepeatFinder.GetTandemRepeatSummary` |
| `find_microsatellites` | X | STRs / microsatellites | `RepeatFinder.FindMicrosatellites` |
| `find_direct_repeats` | X | Direct repeats (spacer-separated) | `RepeatFinder.FindDirectRepeats` |
| `find_inverted_repeats` | X | Inverted repeats / hairpins (DNA) | `RepeatFinder.FindInvertedRepeats` |
| `find_palindromes` | X | DNA palindromes (restriction sites) | `RepeatFinder.FindPalindromes` |
| `find_repeats` | X | All repeats ≥ minLength | `GenomicAnalyzer.FindRepeats` |
| `find_rna_inverted_repeats` | X | RNA hairpin stems | `RnaSecondaryStructure.FindInvertedRepeats` |

## 6. Complexity / low-complexity masking

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `find_low_complexity_regions` | X | Entropy-thresholded LCRs (DNA) | `SequenceComplexity.FindLowComplexityRegions` |
| `mask_low_complexity` | X | Mask low-complexity windows | `SequenceComplexity.MaskLowComplexity` |
| `dust_score` | X | DUST low-complexity score | `SequenceComplexity.CalculateDustScore` |
| `compression_ratio` | X | LZ-normalized complexity | `SequenceComplexity.EstimateCompressionRatio` |
| `windowed_complexity` | X | Windowed Shannon + linguistic complexity | `SequenceComplexity.CalculateWindowedComplexity` |
| `entropy_profile` | X | Sliding-window Shannon entropy | `SequenceStatistics.CalculateEntropyProfile` |
| `find_protein_low_complexity_regions` | X | Protein LCRs via SEG | `ProteinMotifFinder.FindLowComplexityRegions` |
| `predict_low_complexity_seg` | X | SEG LCRs in a protein | `DisorderPredictor.PredictLowComplexityRegions` |

## 7. k-mer & composition analysis

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `count_kmers` | X | Count every k-mer | `KmerAnalyzer.CountKmers` |
| `count_kmers_both_strands` | X | k-mers over both strands | `KmerAnalyzer.CountKmersBothStrands` |
| `analyze_kmers` | X | Aggregate k-mer composition stats | `KmerAnalyzer.AnalyzeKmers` |
| `kmer_frequencies` | X | Normalized k-mer frequencies | `KmerAnalyzer.GetKmerFrequencies` |
| `kmer_spectrum` | X | Frequency-of-frequencies | `KmerAnalyzer.GetKmerSpectrum` |
| `kmer_positions` | X | 0-based positions of a k-mer | `KmerAnalyzer.FindKmerPositions` |
| `most_frequent_kmers` | X | k-mers tied for max count | `KmerAnalyzer.FindMostFrequentKmers` |
| `kmers_with_min_count` | X | k-mers ≥ minCount | `KmerAnalyzer.FindKmersWithMinCount` |
| `unique_kmers` | X | Singletons (count == 1) | `KmerAnalyzer.FindUniqueKmers` |
| `generate_all_kmers` | X | Enumerate whole k-mer space | `KmerAnalyzer.GenerateAllKmers` |
| `find_clumps` | X | k-mers clumping in a window | `KmerAnalyzer.FindClumps` |
| `kmer_distance` | X | Alignment-free composition distance | `KmerAnalyzer.KmerDistance` |
| `analyze_gc_content` | X | Comprehensive GC analysis | `GcSkewCalculator.AnalyzeGcContent` |
| `gc_content_profile` | X | Windowed GC content | `SequenceStatistics.CalculateGcContentProfile` |
| `gc_skew` / `cumulative_gc_skew` / `windowed_gc_skew` | X | GC skew (whole / cumulative / windowed) | `GcSkewCalculator.CalculateGcSkew` / `CalculateCumulativeGcSkew` / `CalculateWindowedGcSkew` |
| `at_skew` | X | Whole-sequence AT skew | `GcSkewCalculator.CalculateAtSkew` |
| `predict_replication_origin` | X | Origin/terminus from cumulative skew | `GcSkewCalculator.PredictReplicationOrigin` |
| `dinucleotide_frequencies` | X | Adjacent dinucleotide frequencies | `SequenceStatistics.CalculateDinucleotideFrequencies` |
| `dinucleotide_ratios` | X | Observed/expected dinucleotide odds | `SequenceStatistics.CalculateDinucleotideRatios` |
| `find_common_regions` | X | Common substrings of two DNA seqs | `GenomicAnalyzer.FindCommonRegions` |

## 8. Splicing

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `find_donor_sites` | A | 5′ (donor) splice sites | `SpliceSitePredictor.FindDonorSites` |
| `find_acceptor_sites` | A | 3′ (acceptor) splice sites | `SpliceSitePredictor.FindAcceptorSites` |
| `find_branch_points` | A | Intron branch-point candidates | `SpliceSitePredictor.FindBranchPoints` |
| `predict_introns` | A | Pair donor+acceptor → introns | `SpliceSitePredictor.PredictIntrons` |
| `predict_gene_structure` | A | Exon/intron structure | `SpliceSitePredictor.PredictGeneStructure` |
| `detect_alternative_splicing` | A | Candidate AS patterns | `SpliceSitePredictor.DetectAlternativeSplicing` |
| `find_retained_intron_candidates` | A | Likely retained introns | `SpliceSitePredictor.FindRetainedIntronCandidates` |
| `is_within_coding_region` | A | Heuristic coding-region check | `SpliceSitePredictor.IsWithinCodingRegion` |
| `maxent_score` | A | MaxEntScan-like splice-site score | `SpliceSitePredictor.CalculateMaxEntScore` |

## 9. Epigenetics / methylation

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `find_cpg_sites` | A | CpG dinucleotide positions | `EpigeneticsAnalyzer.FindCpGSites` |
| `find_cpg_islands` | A | Gardiner-Garden & Frommer islands | `EpigeneticsAnalyzer.FindCpGIslands` |
| `cpg_observed_expected` | A | CpG O/E ratio | `EpigeneticsAnalyzer.CalculateCpGObservedExpected` |
| `find_methylation_sites` | A | Candidate mC sites + context | `EpigeneticsAnalyzer.FindMethylationSites` |
| `methylation_from_bisulfite` | A | Per-CpG methylation from BS reads | `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite` |
| `simulate_bisulfite_conversion` | A | In-silico bisulfite conversion | `EpigeneticsAnalyzer.SimulateBisulfiteConversion` |
| `methylation_profile` | A | Aggregate global/CpG/CHG/CHH profile | `EpigeneticsAnalyzer.GenerateMethylationProfile` |
| `find_dmrs` | A | Differentially methylated regions | `EpigeneticsAnalyzer.FindDMRs` |
| `epigenetic_age` | A | Horvath-style epigenetic age | `EpigeneticsAnalyzer.CalculateEpigeneticAge` |
| `predict_imprinted_genes` | A | Allele-specific-methylation imprinting | `EpigeneticsAnalyzer.PredictImprintedGenes` |
| `annotate_histone_modifications` | A | Chromatin state from a histone mark | `EpigeneticsAnalyzer.AnnotateHistoneModifications` |
| `predict_chromatin_state` | A | State from histone-mod signals | `EpigeneticsAnalyzer.PredictChromatinState` |
| `find_accessible_regions` | A | ATAC-like accessible peaks | `EpigeneticsAnalyzer.FindAccessibleRegions` |

## 10. miRNA ⚠ (MIRNA-TARGET-001 / MIRNA-CLEAVAGE-001 guarded)

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `create_mirna` | A | Build MiRna record (T→U) | `MiRnaAnalyzer.CreateMiRna` |
| `mirna_seed_sequence` | A | Seed region (pos 2–8) | `MiRnaAnalyzer.GetSeedSequence` |
| `find_mirna_target_sites` | A | miRNA target sites in mRNA ⚠ | `MiRnaAnalyzer.FindTargetSites` |
| `align_mirna_to_target` | A | Duplex alignment + stats | `MiRnaAnalyzer.AlignMiRnaToTarget` |
| `analyze_target_context` | A | AU-content + positional context ⚠ | `MiRnaAnalyzer.AnalyzeTargetContext` |
| `site_accessibility` | A | Target-site accessibility | `MiRnaAnalyzer.CalculateSiteAccessibility` |
| `compare_seed_regions` | A | Compare two miRNA seeds | `MiRnaAnalyzer.CompareSeedRegions` |
| `find_similar_mirnas` | A | miRNAs with similar seeds | `MiRnaAnalyzer.FindSimilarMiRnas` |
| `group_by_seed_family` | A | Group by identical seed | `MiRnaAnalyzer.GroupBySeedFamily` |
| `generate_seed_variants` | A | SNVs of a seed | `MiRnaAnalyzer.GenerateSeedVariants` |
| `find_pre_mirna_hairpins` | A | Pre-miRNA hairpin candidates | `MiRnaAnalyzer.FindPreMiRnaHairpins` |
| `can_pair` (A) | A | RNA base-pair test | `MiRnaAnalyzer.CanPair` |
| `is_wobble_pair` | A | G-U wobble test | `MiRnaAnalyzer.IsWobblePair` |
| `rna_reverse_complement` | A | RNA reverse complement | `MiRnaAnalyzer.GetReverseComplement` |

## 11. Transcriptome / RNA-seq

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `calculate_tpm` | A | TPM + FPKM from counts | `TranscriptomeAnalyzer.CalculateTPM` |
| `differential_expression` | A | Two-group DE analysis | `TranscriptomeAnalyzer.AnalyzeDifferentialExpression` |
| `log2_transform` | A | log2(x+pseudocount) | `TranscriptomeAnalyzer.Log2Transform` |
| `quantile_normalize` | A | Quantile-normalize vectors | `TranscriptomeAnalyzer.QuantileNormalize` |
| `pearson_correlation` | A | Pearson correlation of vectors | `TranscriptomeAnalyzer.CalculatePearsonCorrelation` |
| `perform_pca` | A | Project onto first 2 PCs | `TranscriptomeAnalyzer.PerformPCA` |
| `cluster_genes_by_expression` | A | Cluster by profile correlation | `TranscriptomeAnalyzer.ClusterGenesByExpression` |
| `build_coexpression_network` | A | Co-expression network | `TranscriptomeAnalyzer.BuildCoExpressionNetwork` |
| `enrichment_score` | A | GSEA-like running-sum score | `TranscriptomeAnalyzer.CalculateEnrichmentScore` |
| `over_representation_analysis` | A | Pathway/gene-set ORA | `TranscriptomeAnalyzer.PerformOverRepresentationAnalysis` |
| `rnaseq_quality_metrics` | A | Basic RNA-seq QC metrics | `TranscriptomeAnalyzer.CalculateQualityMetrics` |
| `detect_differential_splicing` | A | Differential splicing (2 conditions) | `TranscriptomeAnalyzer.DetectDifferentialSplicing` |
| `detect_isoform_switching` | A | Isoform-usage switching | `TranscriptomeAnalyzer.DetectIsoformSwitching` |
| `find_dominant_isoforms` | A | Dominant isoform per gene | `TranscriptomeAnalyzer.FindDominantIsoforms` |
| `find_skipped_exon_events` | A | PSI for skipped-exon candidates | `TranscriptomeAnalyzer.FindSkippedExonEvents` |

## 12. RNA secondary structure ⚠ (RNA-STRUCT-001 documented limit)

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `predict_rna_structure` | X | Greedy stem-loop structure | `RnaSecondaryStructure.PredictStructure` |
| `minimum_free_energy` | X | Zuker MFE | `RnaSecondaryStructure.CalculateMinimumFreeEnergy` |
| `find_stem_loops` | X | Hairpin stem-loop candidates | `RnaSecondaryStructure.FindStemLoops` |
| `detect_pseudoknots` | X | Crossing base pairs | `RnaSecondaryStructure.DetectPseudoknots` |
| `parse_dot_bracket` / `validate_dot_bracket` | X | Dot-bracket ↔ base pairs | `RnaSecondaryStructure.ParseDotBracket` / `ValidateDotBracket` |
| `can_pair` (X) / `base_pair_type` / `rna_complement_base` | X | RNA pairing primitives | `RnaSecondaryStructure.CanPair` / `GetBasePairType` / `GetComplement` |
| `*_energy` (hairpin/internal/bulge/multibranch/stem/terminal_mismatch/dangling_end/coaxial) | X | Turner 2004 free-energy terms | `RnaSecondaryStructure.Calculate*` / `Get*` |

## 13. Protein features ⚠ (DISORDER-REGION-001 guarded)

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `predict_disorder` | X | TOP-IDP disorder prediction ⚠ | `DisorderPredictor.PredictDisorder` |
| `disorder_propensity` | X | Per-residue TOP-IDP propensity | `DisorderPredictor.GetDisorderPropensity` |
| `is_disorder_promoting` | X | Dunker disorder-promoting set | `DisorderPredictor.IsDisorderPromoting` |
| `predict_morfs` | X | MoRFs within IDRs | `DisorderPredictor.PredictMoRFs` |
| `predict_signal_peptide` | X | von Heijne cleavage-site | `ProteinMotifFinder.PredictSignalPeptide` |
| `predict_transmembrane_helices` | X | Hydropathy TM-helix prediction | `ProteinMotifFinder.PredictTransmembraneHelices` |
| `predict_coiled_coils` | X | Heptad-repeat coiled coils | `ProteinMotifFinder.PredictCoiledCoils` |
| `predict_chou_fasman` | X | Chou-Fasman helix/sheet/turn | `SequenceStatistics.PredictSecondaryStructure` |
| `hydrophobicity_profile` | X | Kyte-Doolittle hydropathy | `SequenceStatistics.CalculateHydrophobicityProfile` |

## 14. Comparative genomics

| Tool | Srv | Purpose | Method ID |
|---|---|---|---|
| `calculate_ani` | X | Average Nucleotide Identity | `ComparativeGenomics.CalculateANI` |
| `compare_genomes` | X | End-to-end comparative pipeline | `ComparativeGenomics.CompareGenomes` |
| `find_orthologs` | X | Best-hit ortholog pairs | `ComparativeGenomics.FindOrthologs` |
| `find_reciprocal_best_hits` | X | RBH pairs | `ComparativeGenomics.FindReciprocalBestHits` |
| `find_syntenic_blocks` | X | Collinear ortholog runs | `ComparativeGenomics.FindSyntenicBlocks` |
| `find_conserved_clusters` | X | Common-interval gene clusters | `ComparativeGenomics.FindConservedClusters` |
| `generate_dot_plot` | X | Matching k-mer coordinates | `ComparativeGenomics.GenerateDotPlot` |
| `reversal_distance` | X | Lower-bound reversal distance | `ComparativeGenomics.CalculateReversalDistance` |

---

Not exhaustively tabled above (open the per-tool doc or use `seqeron-discovery`): the remaining
Turner-2004 RNA energy terms, and any tool added after this map was written — the **generated**
`_generated/tools.md` is the authoritative complete list.
