# Gotchas

Non-obvious limitations and sharp edges.

Return to the [[index|Wiki Index]].

- [[predict-genes-emits-every-orf]] — `PredictGenes` emits every qualifying ORF as a CDS; it does not resolve overlaps, rank genes, or integrate RBS evidence.
- [[research-grade-limitations]] — beta, not for clinical use; simplified-subset implementations; internal-only validation.
- [[differential-expression-is-welch-t-not-deseq2]] — `differential_expression` is a two-group log2FC + Welch-t + BH screen, not the DESeq2/edgeR NB-GLM; small-N calls are unreliable.
- [[variant-calling-is-alignment-not-pileup]] — `find_snps`/`find_indels` compare one query to a reference (no genotypes/depth); indels aren't left-normalized, so positions drift in repeats.
- [[replication-origin-assumes-single-circular-chromosome]] — `predict_replication_origin` maps min-skew→origin only for a single circular chromosome; meaningless on linear/draft/eukaryotic input.
- [[ani-is-directional-use-reciprocal]] — `calculate_ani` fragments only the query, so ANI(A,B)≠ANI(B,A); use reciprocal ANI for a symmetric species-boundary call.
- [[p-distance-is-uncorrected-and-jc-k2p-saturate]] — p-distance underestimates divergence (no multiple-hit correction); JC69/K2P return +∞ at saturation.
- [[shared-motifs-counts-sequences-not-occurrences]] — `find_shared_motifs` counts how many sequences contain a motif (first occurrence only), not total occurrences.
- [[chao1-degrades-to-observed-richness-on-proportions]] — `alpha_diversity` Chao1 silently returns observed richness for non-integer/proportional abundances.
- [[kmer-counts-are-not-canonical-collapsed]] — `count_kmers` is forward-only and `count_kmers_both_strands` is additive; neither is Jellyfish `-C` canonical collapsing.
- [[crispr-guide-and-offtarget-scores-are-heuristic]] — CRISPR guide (`evaluate_guide_rna`) and off-target (`crispr_specificity_score`) scores are heuristics, not calibrated efficacy/risk.
- [[olc-assembly-is-greedy-not-optimal]] — `assemble_olc` is a greedy, order-dependent, error-free-read heuristic, not an exact shortest-common-superstring.
- [[rbh-orthologs-use-kmer-jaccard-not-blast]] — `find_orthologs` ranks best hits by 5-mer Jaccard, not BLAST bit-score, so RBH sets won't match a BLAST pipeline.
- [[gene-structure-is-splice-motif-not-coding-aware]] — `predict_gene_structure` pairs GT-AG splice motifs; it is not coding/frame-aware and splicing means intron removal.
- [[disorder-prediction-is-top-idp-propensity-not-trained]] — `predict_disorder` is a TOP-IDP propensity window, not a trained IUPred/DISOPRED probability.
- [[coiled-coil-score-is-heptad-occupancy-heuristic]] — `predict_coiled_coils` scores heptad a/d-core hydrophobic occupancy, not a COILS/Marcoil probability.
- [[chou-fasman-returns-propensity-profile-not-state-calls]] — `predict_chou_fasman` returns Pα/Pβ/Pt propensity profiles, not H/E/C state assignments.
- [[methylation-context-classifies-not-calls-methylation]] — `find_methylation_sites` classifies cytosine context (CpG/CHG/CHH) from sequence; it does not call methylation.
- [[chromatin-state-is-mark-presence-rules-not-chromhmm]] — `predict_chromatin_state` uses present-mark combinatorics (magnitudes ignored), not ChromHMM on read counts.
- [[vep-needs-reference-window-for-coding-consequences]] — `annotate_variants` without a reference window returns only the coarse coding term (no missense/synonymous/stop).
- [[read-depth-cnv-misses-copy-neutral-and-nan-to-cn2]] — `segment_copy_number` is depth-only (blind to copy-neutral events); NaN windows silently become neutral CN2.
- [[discordant-pairs-are-signatures-not-breakpoints]] — `find_discordant_pairs` yields SV signatures/intervals and type, not base-pair breakpoints.
- [[find-protein-domains-is-prosite-pattern-not-pfam]] — `find_protein_domains` matches exact PROSITE patterns + 3 bundled profiles, not a full Pfam/HMMER scan.
- [[taxonomic-profile-is-read-tally-not-marker-abundance]] — `taxonomic_profile` is read-tally relative abundance over classified reads, not MetaPhlAn marker abundance.
