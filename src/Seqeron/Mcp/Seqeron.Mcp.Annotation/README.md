# Seqeron.Mcp.Annotation

MCP server — **Gene/ORF/promoter annotation, variant calling & effect, epigenetics, miRNA, splicing, SVs, transcriptomics.**

Exposes **97 tools** wrapping the `Seqeron.Genomics` library. Every tool has an
explicit JSON input/output schema, a Schema+Binding test, and per-tool docs under
[`docs/mcp/tools/annotation/`](../../../../docs/mcp/tools/annotation) — see the
campaign ledger [`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Annotation
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Annotation"]`). See [`docs/mcp/README.md`](../../../../docs/mcp/README.md).

## Tools (97)

| Tool | Description |
|------|-------------|
| `align_mirna_to_target` | Align a miRNA against a target sequence and compute duplex statistics (matches, mismatches, G-U wobbles, gaps, free energy). |
| `analyze_target_context` | Compute AU content and positional context score around a miRNA target site within an mRNA. |
| `annotate_histone_modifications` | Annotate intervals with predicted chromatin state from a histone mark and signal level. |
| `annotate_regulatory_elements` | Find regulatory regions overlapping a variant. |
| `annotate_svs` | Annotate structural variants with overlapping genes and exons and a coarse functional impact (HIGH / MODERATE / MODIFIER / LOW). |
| `annotate_variant_on_transcripts` | VEP-like annotation: predict consequence, impact, codon/AA change, SIFT/PolyPhen for one variant against transcripts. |
| `annotate_variants` | Call and annotate all variants between reference and query (effect + Ti/Tv classification). |
| `assemble_breakpoint_sequence` | Heuristically assemble a breakpoint-junction sequence from split-read clipped fragments (returns null if no reads). |
| `build_coexpression_network` | Build a gene co-expression network by emitting all gene pairs whose Pearson correlation magnitude meets or exceeds the threshold. |
| `calculate_conservation` | Compute PhyloP-, PhastCons-, GERP-like conservation scores from per-position multi-species alleles. |
| `calculate_tpm` | Compute Transcripts-Per-Million (TPM) and FPKM expression values from per-gene raw counts and gene lengths. |
| `call_variants` | Detect SNPs and indels between two DNA sequences using global alignment. |
| `call_variants_from_alignment` | Detect variants from already-aligned reference and query strings (gaps as '-'). |
| `can_pair` | Test whether two RNA bases can pair (Watson-Crick A-U, G-C, or G-U wobble). |
| `classify_mutation` | Classify a SNP as transition, transversion, or other. |
| `classify_variant` | Classify a (ref, alt) allele pair as SNV / Insertion / Deletion / MNV / Indel / Complex. |
| `cluster_discordant_pairs` | Cluster discordant read pairs into structural-variant (SV) candidates. |
| `cluster_genes_by_expression` | Cluster genes by similarity of their expression profiles using a k-means-like procedure on Pearson correlation. |
| `cluster_split_reads` | Cluster split reads into breakpoint candidates. |
| `coding_potential` | Score the coding potential of a DNA sequence using the CPAT hexamer usage-bias log-likelihood (Wang et al. |
| `codon_usage` | Count codon occurrences (in-frame, 5'→3') for a coding sequence. |
| `compare_seed_regions` | Compare seed regions of two miRNAs and report Hamming-style match counts and seed-family identity. |
| `cpg_observed_expected` | Compute CpG observed/expected ratio for a sequence. |
| `create_mirna` | Build a MiRna record (T→U normalized) with seed metadata. |
| `detect_alternative_splicing` | Detect candidate alternative splicing patterns (exon skipping, alternative 5'/3' splice sites). |
| `detect_differential_splicing` | Detect splicing events whose PSI differs between two conditions by more than the given delta-PSI threshold. |
| `detect_isoform_switching` | Detect genes where the dominant transcript isoform switches between two conditions based on changes in usage proportions. |
| `differential_expression` | Run a simple two-group differential-expression analysis using log2 fold change, Welch-style t-test p-values, and Benjamini–Hochberg FDR c… |
| `enrichment_score` | Compute a GSEA-like running-sum enrichment score for a gene set against a ranked gene list. |
| `epigenetic_age` | Estimate epigenetic age (Horvath-clock-style) from methylation values at clock CpGs. |
| `filter_svs` | Filter structural variants by quality, supporting-read count, and length bounds. |
| `find_acceptor_sites` | Find candidate 3' (acceptor) splice sites (canonical AG, optional U12 AC). |
| `find_accessible_regions` | Identify accessible chromatin regions (ATAC-seq-like peaks) from per-position accessibility signal. |
| `find_branch_points` | Find branch-point candidates using a YNYURAC-like position weight matrix. |
| `find_conserved_elements` | Identify conserved genomic elements from per-position conservation scores (runs above threshold). |
| `find_cpg_islands` | Identify CpG islands using Gardiner-Garden & Frommer criteria (length, GC content, CpG O/E). |
| `find_cpg_sites` | Return all CpG dinucleotide start positions in a sequence. |
| `find_deletions` | Detect only deletions between two DNA sequences. |
| `find_discordant_pairs` | Identify discordant read pairs (interchromosomal, abnormal insert size, or abnormal orientation). |
| `find_dmrs` | Identify differentially methylated regions (DMRs) between two methylation samples using a sliding window. |
| `find_dominant_isoforms` | For each gene, identify the most-expressed transcript isoform and report its share of the total per-gene expression. |
| `find_donor_sites` | Find candidate 5' (donor) splice sites (canonical GT/GU, optional GC and U12 AT/AU). |
| `find_indels` | Detect insertions and deletions (indels) between two DNA sequences. |
| `find_insertions` | Detect only insertions between two DNA sequences. |
| `find_methylation_sites` | Find candidate cytosine methylation sites and classify context (CpG/CHG/CHH). |
| `find_microhomology` | Find the longest shared microhomology between the 3' end of a left flank and the 5' end of a right flank (up to maxLength nt). |
| `find_mirna_target_sites` | Find miRNA target sites in an mRNA (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer per Bartel 2009 and TargetScan). |
| `find_orfs` | Find all open reading frames (ORFs) in a DNA sequence across forward and reverse strands. |
| `find_pre_mirna_hairpins` | Identify pre-miRNA hairpin candidates with stem/loop layout and Turner 2004 nearest-neighbor free energy. |
| `find_promoter_motifs` | Find -10 (Pribnow/TATA) and -35 box bacterial promoter motifs. |
| `find_repetitive_elements` | Find tandem and inverted repeats in a DNA sequence. |
| `find_retained_intron_candidates` | Find short, moderately-scored intron candidates likely to be retained in some transcripts. |
| `find_ribosome_binding_sites` | Locate Shine–Dalgarno (ribosome binding site) motifs upstream of ORFs. |
| `find_similar_mirnas` | Find miRNAs in a database whose seed region is within `maxMismatches` Hamming distance of a query miRNA. |
| `find_skipped_exon_events` | Compute Percent Spliced In (PSI = inclusion / (inclusion + skipping)) for candidate skipped-exon events. |
| `find_snps` | Detect only SNPs between two DNA sequences using global alignment. |
| `find_snps_direct` | Detect SNPs by direct positional comparison without alignment. |
| `find_split_reads` | Find split reads from soft-clipped CIGAR alignments. |
| `format_vcf_info` | Format a variant annotation as a VCF INFO field string (GENE/TRANSCRIPT/CONSEQUENCE/IMPACT/HGVSp/HGVSc/SIFT/POLYPHEN). |
| `generate_seed_variants` | Enumerate single-nucleotide variants of a seed sequence (the original plus all single-position substitutions over A/C/G/U). |
| `genotype_sv` | Genotype a structural variant (0/0, 0/1, 1/1, ./.) from reference and alternate supporting-read counts. |
| `group_by_seed_family` | Group miRNAs by identical seed sequence (seed family). |
| `identify_cnvs` | Convert non-baseline copy-number segments into deletion / duplication structural variants. |
| `impact_level` | Map a ConsequenceType to an ImpactLevel (High/Moderate/Low/Modifier). |
| `is_within_coding_region` | Heuristic check whether a sequence position lies in a coding region downstream of an in-frame start codon. |
| `is_wobble_pair` | Test whether two bases form a G-U wobble pair. |
| `log2_transform` | Apply log2(x + pseudocount) to each value (used to stabilize variance of count-like expression data). |
| `longest_orfs_per_frame` | Return the longest ORF in each of the (up to) six reading frames. |
| `maxent_score` | Compute a MaxEntScan-like log-likelihood score for a donor or acceptor motif. |
| `merge_overlapping_svs` | Merge structural variants of the same type whose reciprocal overlap fraction exceeds a threshold. |
| `methylation_from_bisulfite` | Compute per-CpG methylation levels from bisulfite sequencing reads aligned to a reference. |
| `methylation_profile` | Aggregate methylation sites into a global/CpG/CHG/CHH profile with per-position methylation. |
| `mirna_seed_sequence` | Extract the canonical seed region (positions 2-8) from a miRNA sequence. |
| `normalize_variant` | Normalize a variant by left-trimming common prefixes and right-trimming common suffixes; |
| `over_representation_analysis` | Pathway / gene-set over-representation analysis (ORA) using a hypergeometric / Fisher's exact-test approximation. |
| `parse_gff3` | Parse a GFF3 annotation document into structured features. |
| `parse_vcf_variant` | Build a Variant record from VCF fields and classify its type. |
| `pearson_correlation` | Compute the Pearson product-moment correlation coefficient between two expression vectors of equal length. |
| `perform_pca` | Project samples onto the first two principal components computed from the most variable genes (approximate PCA). |
| `predict_chromatin_state` | Predict chromatin state from H3K4me3 / H3K4me1 / H3K27ac / H3K36me3 / H3K27me3 / H3K9me3 signal levels (0..1). |
| `predict_gene_structure` | Predict exon/intron gene structure with greedy non-overlapping intron selection. |
| `predict_genes` | Predict gene annotations from a DNA sequence using ORF-based heuristics. |
| `predict_imprinted_genes` | Predict imprinted genes from allele-specific (maternal vs paternal) methylation differences. |
| `predict_introns` | Predict introns by pairing donor and acceptor splice sites and locating branch points. |
| `predict_pathogenicity` | ACMG-like pathogenicity prediction combining annotation, frequency, conservation, ClinVar, and functional evidence. |
| `predict_tf_binding_change` | Predict transcription-factor binding score changes induced by a SNV against IUPAC motifs. |
| `predict_variant_effect` | Predict the protein-level effect of a variant in a coding sequence (Synonymous/Missense/Nonsense/StopLoss/Frameshift/Unknown). |
| `quantile_normalize` | Quantile-normalize multiple expression vectors of equal length so that each sample shares the same value distribution. |
| `rna_reverse_complement` | Reverse-complement of an RNA sequence (A↔U, G↔C, unknown bases → N). |
| `rnaseq_quality_metrics` | Compute basic RNA-seq QC metrics: mapping rate, exonic rate, rRNA rate, and number of detected genes. |
| `segment_copy_number` | Segment copy-number probe data into log-ratio plateaus (CBS-like). |
| `simulate_bisulfite_conversion` | Simulate bisulfite conversion (unmethylated C → T) given a set of methylated positions to protect. |
| `site_accessibility` | Estimate miRNA target site accessibility from local secondary-structure density (lower structure ⇒ higher accessibility). |
| `titv_ratio` | Compute transition/transversion (Ti/Tv) ratio from a list of variants. |
| `to_gff3` | Serialize gene annotations to GFF3 lines (with header). |
| `variant_statistics` | Compute summary variant statistics between reference and query (totals, Ti/Tv, density). |
| `variants_to_vcf` | Format variants as VCF v4.2 lines (header + records). |
