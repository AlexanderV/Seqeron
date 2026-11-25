# Seqeron.Mcp.Annotation

MCP server — **Gene/ORF/promoter annotation, variant calling & effect, epigenetics, miRNA, splicing, SVs, transcriptomics.**

Exposes **112 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/annotation/`](../../../../docs/mcp/tools/annotation). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Annotation
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Annotation"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (112)

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
| `calculate_randfold_pvalue` | randfold randomisation test (Bonnet, Wuyts, Rouze & Van de Peer 2004, Bioinformatics 20:2911; dinucleotide mode randfold -d): the input is MFE-folded (canonical Turner 2004 engine), `iterations` Altschul-Erickson dinucleotide shuffles are folded, and P = (1 + #{shuffle dG <= input dG}) / (iterations + 1) (\|diff\| <= 1e-7 counts |
| `calculate_tpm` | Compute Transcripts-Per-Million (TPM) and FPKM expression values from per-gene raw counts and gene lengths. |
| `call_variants` | Detect SNPs and indels between two DNA sequences using global alignment. |
| `call_variants_from_alignment` | Detect variants from already-aligned reference and query strings (gaps as '-'). |
| `can_pair` | Test whether two RNA bases can pair (Watson-Crick A-U, G-C, or G-U wobble). |
| `classify_mutation` | Classify a SNP as transition, transversion, or other. |
| `classify_pre_mirna` | Classify a candidate hairpin as natural pre-miRNA vs dinucleotide-shuffled background with the bundled IN-REPO ridge-logistic model over MFE / AMFE / MFEI / GC / paired-fraction features (Bonnet 2004 shuffle convention, Zhang 2006 MFEI). |
| `classify_pre_mirna_triplet_svm` | Published Triplet-SVM pre-miRNA classifier (Xue, Li, He, Liu, Li & Zhang 2005, BMC Bioinformatics 6:310) with the authors' own trained libsvm RBF model: retrieves the stem-loops with >= minBasePairs paired bases per arm (2_get_stemloop.pl), codes each single stem-loop as 32 triplet-element frequencies and returns the libsvm deci |
| `classify_variant` | Classify a (ref, alt) allele pair as SNV / Insertion / Deletion / MNV / Indel / Complex. |
| `cluster_discordant_pairs` | Cluster discordant read pairs into structural-variant (SV) candidates. |
| `cluster_genes_by_expression` | Cluster genes by similarity of their expression profiles using a k-means-like procedure on Pearson correlation. |
| `cluster_split_reads` | Cluster split reads into breakpoint candidates. |
| `coding_potential` | Score the coding potential of a DNA sequence using the CPAT hexamer log-likelihood (Wang et al. 2013); -1 when no hexamer is scorable. |
| `codon_usage` | Count in-frame codons of a coding DNA or RNA sequence (U read as T; triplets with other symbols skipped without shifting the frame). |
| `compare_seed_regions` | Compare seed regions of two miRNAs and report Hamming-style match counts and seed-family identity. |
| `compute_ta_3utr` | TargetScan context++ feature TA_3UTR (Garcia et al. 2011; Agarwal et al. 2015): TA = log10(N), N = total number of non-overlapping canonical 8mer + 7mer-m8 + 7mer-A1 sites of the miRNA seed (nt 2-8) across the supplied 3'UTRs (bare 6mers and offset 6mers are not counted; scan is left-to-right greedy). |
| `cpat_coding_probability` | CPAT (Wang et al. 2013): ORF candidates scored by transcript length, ORF length, Fickett TESTCODE and hexamer score; logistic coding probability. |
| `cpg_observed_expected` | Compute CpG observed/expected ratio for a sequence. |
| `create_mirna` | Build a MiRna record (T→U normalized) with seed metadata. |
| `detect_alternative_splicing` | Detect candidate alternative splicing events (SUPPA/rMATS: SE, A5SS, A3SS, MXE) relative to the predicted gene structure. |
| `detect_differential_splicing` | Detect splicing events whose PSI differs between two conditions by more than the given delta-PSI threshold. |
| `detect_isoform_switching` | Detect genes where the dominant transcript isoform switches between two conditions based on changes in usage proportions. |
| `differential_expression` | Run a simple two-group differential-expression analysis using log2 fold change, Welch-style t-test p-values, and Benjamini–Hochberg FDR c… |
| `enrichment_score` | Compute a GSEA-like running-sum enrichment score for a gene set against a ranked gene list. |
| `epigenetic_age` | Estimate epigenetic age (Horvath-clock-style) from methylation values at clock CpGs. |
| `fickett_scan` | Sliding-window Fickett (1982) TESTCODE with coding calls, as EMBOSS tcode. |
| `fickett_score` | Fickett (1982) TESTCODE coding statistic of a whole DNA sequence (published or CPAT table set). |
| `filter_svs` | Filter structural variants by quality, supporting-read count, and length bounds. |
| `find_acceptor_sites` | Find 3' (acceptor) splice sites: Shapiro & Senapathy 1987 score (default) or MaxEntScan score3; optional U12 acceptors (UUCCUUAAC branch site). |
| `find_accessible_regions` | Identify accessible chromatin regions (ATAC-seq-like peaks) from per-position accessibility signal. |
| `find_branch_points` | Find branch-point candidates: human yUnAy (Gao et al. 2008) or yeast UACUAAC. |
| `find_centered_sites` | Find centered miRNA target sites (Shin et al. 2010, Mol Cell 38:789): sites without perfect seed pairing (miRNA nt 2-7 not all Watson-Crick) that pair the miRNA centre with >= 11 contiguous Watson-Crick pairs (nt 4-14, 5-15 or 4-15; G:U wobble does not count). |
| `find_conserved_elements` | Identify conserved genomic elements from per-position conservation scores (runs above threshold). |
| `find_cpg_islands` | Identify CpG islands using Gardiner-Garden & Frommer criteria (length, GC content, CpG O/E). |
| `find_cpg_sites` | Return all CpG dinucleotide start positions in a sequence. |
| `find_deletions` | Detect only deletions between two DNA sequences. |
| `find_discordant_pairs` | Identify discordant read pairs (interchromosomal, abnormal insert size, or abnormal orientation). |
| `find_dmrs` | Identify differentially methylated regions (DMRs) between two methylation samples using a sliding window. |
| `find_dominant_isoforms` | For each gene, identify the most-expressed transcript isoform and report its share of the total per-gene expression. |
| `find_donor_sites` | Find 5' (donor) splice sites: Shapiro & Senapathy 1987 score (default) or MaxEntScan score5; optional GC and U12 /RTATCCTY donors. |
| `find_indels` | Detect insertions and deletions (indels) between two DNA sequences. |
| `find_insertions` | Detect only insertions between two DNA sequences. |
| `find_methylation_sites` | Find candidate cytosine methylation sites and classify context (CpG/CHG/CHH). |
| `find_microhomology` | Find the longest shared microhomology between the 3' end of a left flank and the 5' end of a right flank (up to maxLength nt). |
| `find_mirna_target_sites` | Find miRNA target sites in an mRNA (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer per Bartel 2009 and TargetScan). |
| `find_orfs` | Find all open reading frames (ORFs) in a DNA sequence across forward and reverse strands (delegates to `Translator`; stop-to-stop or start-to-stop). |
| `find_pre_mirna_hairpins` | Identify pre-miRNA hairpin candidates with stem/loop layout and Turner 2004 nearest-neighbor free energy. |
| `find_promoter_motifs` | Find paired σ70 -35 (TTGACA) / -10 (TATAAT) boxes, <= 2 mismatches each, 15-21 bp spacer (Harley & Reynolds 1987); score = 1 - mismatches/6. |
| `find_repetitive_elements` | Find exact tandem repeats (maximal runs, primitive ACGT unit 1-60 bp) and maximal exact inverted repeats (EMBOSS palindrome model, arm >= 1) via RepeatFinder. |
| `find_retained_intron_candidates` | Predicted introns that keep the reading frame open (length multiple of 3, no in-frame stop; Jacob & Smith 2017); optional Braunschweig 2014 filters. |
| `find_ribosome_binding_sites` | Locate Shine–Dalgarno (ribosome binding site) motifs upstream of ORFs; score = motif length / 6 (default) or Prodigal log-odds trained on the input genome (`scoreMode`). |
| `find_similar_mirnas` | Find miRNAs in a database whose seed region is within `maxMismatches` Hamming distance of a query miRNA. |
| `find_skipped_exon_events` | Compute Percent Spliced In (PSI = inclusion / (inclusion + skipping)) for candidate skipped-exon events. |
| `find_snps` | Detect only SNPs between two DNA sequences using global alignment. |
| `find_snps_direct` | Detect SNPs by direct positional comparison without alignment. |
| `find_split_reads` | Find split reads from soft-clipped CIGAR alignments. |
| `find_targetscan_site_groups` | Port of targetscan_70.pl: finds 8mer / 7mer-m8 / 7mer-A1 / 6mer seed-match sites of a miRNA family in every species of a multi-species 3'UTR alignment (alignedUtrs: species identifier → aligned UTR with '-' gaps; lower-case letters flag ORF overlap; the first 14 nt of each UTR are masked as TargetScan does) and groups sites that |
| `fold_mirna_duplex` | Minimum-free-energy miRNA·target duplex (RNAduplex / duplexfold model): one contiguous duplex with internal loops and bulges up to maxLoopSize per side, scored with Turner 2004 nearest-neighbor parameters (kcal/mol, 37 °C). |
| `format_vcf_info` | Format a variant annotation as a VCF INFO field string (GENE/TRANSCRIPT/CONSEQUENCE/IMPACT/HGVSp/HGVSc/SIFT/POLYPHEN). |
| `generate_seed_variants` | Enumerate single-nucleotide variants of a seed sequence (the original plus all single-position substitutions over A/C/G/U). |
| `genotype_sv` | Genotype a structural variant (0/0, 0/1, 1/1, ./.) from reference and alternate supporting-read counts. |
| `group_by_seed_family` | Group miRNAs by identical seed sequence (seed family). |
| `identify_cnvs` | Convert non-baseline copy-number segments into deletion / duplication structural variants. |
| `impact_level` | Map a ConsequenceType to an ImpactLevel (High/Moderate/Low/Modifier). |
| `is_within_coding_region` | Whether a position lies in an open reading frame opened by an upstream AUG with no in-frame stop before it. |
| `is_wobble_pair` | Test whether two bases form a G-U wobble pair. |
| `log2_transform` | Apply log2(x + pseudocount) to each value (used to stabilize variance of count-like expression data). |
| `longest_orfs_per_frame` | Return the longest ORF in each of the (up to) six reading frames. |
| `maxent_score` | MaxEntScan (Yeo & Burge 2004) splice-site score in bits: Donor = score5 on the 9-nt window, Acceptor = score3 on the 23-nt window. |
| `merge_overlapping_svs` | Merge structural variants of the same type whose reciprocal overlap fraction exceeds a threshold. |
| `methylation_from_bisulfite` | Compute per-CpG methylation levels from bisulfite sequencing reads aligned to a reference. |
| `methylation_profile` | Aggregate methylation sites into a global/CpG/CHG/CHH profile with per-position methylation. |
| `mirna_seed_sequence` | Extract the canonical seed region (positions 2-8) from a miRNA sequence. |
| `normalize_variant` | Normalize a variant by left-trimming common prefixes and right-trimming common suffixes; |
| `over_representation_analysis` | Pathway / gene-set over-representation analysis (ORA) using a hypergeometric / Fisher's exact-test approximation. |
| `parse_gff3` | Parse a GFF3 (SO v1.26) document: percent-decoded fields, multi-valued attributes, stops at `##FASTA`; spec-invalid feature lines are skipped. |
| `parse_vcf_variant` | Build a Variant record from VCF fields and classify its type. |
| `pearson_correlation` | Compute the Pearson product-moment correlation coefficient between two expression vectors of equal length. |
| `perform_pca` | Project samples onto the first two principal components computed from the most variable genes (approximate PCA). |
| `pooled_rscu` | RSCU pooled over several in-frame coding sequences for an NCBI genetic code (CodonW `-totals`; stop codons excluded; unobserved family = 0). |
| `predict_chromatin_state` | Predict chromatin state from H3K4me3 / H3K4me1 / H3K27ac / H3K36me3 / H3K27me3 / H3K9me3 signal levels (0..1). |
| `predict_drosha_dicer_cleavage` | Predict Drosha (basal) and Dicer (apical) cleavage sites of a pri-miRNA hairpin by the published measuring rules applied to its stem (base pairs counted up the unbranched stem from basalJunction; coordinates 0-based, spans inclusive; 2-nt 3' overhang). |
| `predict_gene_structure` | Predict exon/intron structure: greedy intron selection (default, not coding-aware) or `model=splice_dp`, the exact MaxEntScan splice + open-reading-frame DP. |
| `predict_genes` | Predict genes from ORFs; `startSelection` LongestOrf (default), ProdigalScore or ProdigalDynamicProgramming (Prodigal V2.6.3 model trained on the input genome, >= 20000 nt). |
| `predict_imprinted_genes` | Predict imprinted genes from allele-specific (maternal vs paternal) methylation differences. |
| `predict_introns` | Predict introns by pairing donors and acceptors (U2 GT-AG/GC-AG; U12 needs donor + UUCCUUAAC branch site) and locating yUnAy branch points. |
| `predict_pathogenicity` | ACMG-like pathogenicity prediction combining annotation, frequency, conservation, ClinVar, and functional evidence. |
| `predict_tf_binding_change` | Predict transcription-factor binding score changes induced by a SNV against IUPAC motifs. |
| `predict_variant_effect` | Predict the protein-level effect of a variant in a coding sequence (Synonymous/Missense/Nonsense/StopLoss/Frameshift/Unknown). |
| `quantile_normalize` | Quantile-normalize multiple expression vectors of equal length so that each sample shares the same value distribution. |
| `rna_reverse_complement` | Reverse-complement of an RNA sequence (A↔U, G↔C, unknown bases → N). |
| `rnaseq_quality_metrics` | Compute basic RNA-seq QC metrics: mapping rate, exonic rate, rRNA rate, and number of detected genes. |
| `score_target_site_context_plus_plus` | TargetScan 7 context++ score (Agarwal et al. 2015, eLife 4:e05005; port of targetscan_70_context_scores.pl) of an 8mer, 7mer-m8, 7mer-A1 or 6mer site in a 3'UTR (siteStart/siteEnd 0-based inclusive; lengths 8/7/7/6). |
| `score_target_site_grimson` | TargetScan 5.0 / Grimson et al. (2007, Mol Cell 27:91) context score of an 8mer, 7mer-m8 or 7mer-A1 site in a 3'UTR (index 0 = first nt after the stop codon): site type + 3' pairing + local AU + position contributions in log2 fold-change units (negative = repressed). |
| `score_targetscan_conservation` | End-to-end TargetScan 7.0 conservation of a miRNA family's sites in a multi-species 3'UTR alignment: site groups (targetscan_70.pl) → UTR bin 1-10 from the median per-position branch length (targetscan_70_BL_bins.pl) → per-site branch length on the bin tree, PCT with the bundled TargetScan 7.0 sigmoid rows and the 'Conserved' ca |
| `segment_copy_number` | Segment copy-number probe data into log-ratio plateaus (CBS-like). |
| `simulate_bisulfite_conversion` | Simulate bisulfite conversion (unmethylated C → T) given a set of methylated positions to protect. |
| `site_accessibility` | Estimate miRNA target site accessibility from local secondary-structure density (lower structure ⇒ higher accessibility). |
| `titv_ratio` | Compute transition/transversion (Ti/Tv) ratio from a list of variants. |
| `to_gff3` | Serialize gene annotations to spec-valid GFF3 (header, CDS phase per gene, percent-encoding); rejects invalid strands and empty/inverted intervals. |
| `variant_statistics` | Compute summary variant statistics between reference and query (totals, Ti/Tv, density). |
| `variants_to_vcf` | Format variants as VCF v4.2 lines (header + records). |
