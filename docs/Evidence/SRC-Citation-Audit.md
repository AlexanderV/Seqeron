# SRC citation audit — per-occurrence table (2026-10-14)

Census script: every `*.md` / `*.cs` tracked file outside `wiki/`, URLs on `wikipedia.org` / `addgene.org`. Each distinct URL
was fetched 2026-10-14 (Wikipedia: `index.php?action=raw` wikitext, redirects followed; Addgene: page text). Verdicts come
from four NUMBER and two DEFINITION verification passes (claim context vs fetched page), summarised in
`docs/Validation/review-2026-09/SRC.md`. Line numbers are those of the census (HEAD at 2026-10-14 ~23:30 UTC); later edits may shift them.
Two URLs were truncated by the first census regex at `,` / `'` (`N50,_L50,_and_related_statistics`, `Fisher's_exact_test`);
the verifiers checked those occurrences against the real articles (fetched separately, HTTP 200).

Classes: NUMBER = backs a numeric value; DEFINITION = backs prose; BIBLIOGRAPHY = reference-list entry (dead-link check only).
Verdicts: SUPPORTED (page states the number), NO_NUMBER (classifier over-call; checked as a definition), SUPPORTED_DEF,
NOT_STATED / NOT_SUPPORTED_DEF (defect — corrected, see SRC.md F-entries), LINK-OK (HTTP 200, not a disambiguation page).

## Distinct URLs

| URL | HTTP 2026-10-14 | redirect | occurrences |
|---|---|---|---|
| https://blog.addgene.org/plasmids-101-how-to-verify-your-plasmid | 200 | — | 1 |
| https://blog.addgene.org/plasmids-101-restriction-cloning | 200 | — | 1 |
| https://en.wikipedia.org/w/index.php?title=N50,_L50,_and_related_statistics&action=raw | 200 | — | 1 |
| https://en.wikipedia.org/wiki/A_Mathematical_Theory_of_Communication | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Allele_frequency | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Alpha_diversity | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Amino_acid | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Aneuploidy | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Approximate_string_matching | 200 | — | 1 |
| https://en.wikipedia.org/wiki/BED_(file_format) | 200 | — | 4 |
| https://en.wikipedia.org/wiki/BLAST_(biotechnology) | 200 | — | 3 |
| https://en.wikipedia.org/wiki/BamHI | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Base_pair | 200 | — | 6 |
| https://en.wikipedia.org/wiki/Beta_diversity | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Binning_(metagenomics) | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity#Example | 200 | — | 2 |
| https://en.wikipedia.org/wiki/CIBERSORT | 200 | — | 2 |
| https://en.wikipedia.org/wiki/CRISPR | 200 | — | 1 |
| https://en.wikipedia.org/wiki/CRISPR_gene_editing | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Centromere | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Chi-squared_test | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Chromosomal_rearrangement | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Chromosome | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Clustal | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Codon_Adaptation_Index | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Codon_usage_bias | 200 | — | 9 |
| https://en.wikipedia.org/wiki/Coiled_coil | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Comparative_genomics | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Complementarity_(molecular_biology) | 200 | — | 6 |
| https://en.wikipedia.org/wiki/Consensus_sequence | 200 | — | 14 |
| https://en.wikipedia.org/wiki/Copy_number_variation | 200 | — | 2 |
| https://en.wikipedia.org/wiki/CpG_site | 200 | — | 4 |
| https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables | 200 | — | 1 |
| https://en.wikipedia.org/wiki/DNA_melting | 200 | Nucleic acid thermodynamics | 1 |
| https://en.wikipedia.org/wiki/DNA_microarray | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Direct_repeat | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Distance_matrices_in_phylogeny | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Diversity_index | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Diversity_index#Shannon_index | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) | 200 | — | 5 |
| https://en.wikipedia.org/wiki/E-box | 200 | — | 2 |
| https://en.wikipedia.org/wiki/EcoRI | 200 | — | 4 |
| https://en.wikipedia.org/wiki/EcoRV | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Edit_distance | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Entropy_(information_theory) | 200 | — | 19 |
| https://en.wikipedia.org/wiki/Exon | 200 | — | 2 |
| https://en.wikipedia.org/wiki/F-statistics | 200 | — | 4 |
| https://en.wikipedia.org/wiki/FASTA_format | 200 | — | 3 |
| https://en.wikipedia.org/wiki/FASTQ_format | 200 | — | 5 |
| https://en.wikipedia.org/wiki/False_discovery_rate | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Fisher%27s_exact_test | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Fisher's_exact_test | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Fixation_index | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Fluorescence_in_situ_hybridization | 200 | — | 1 |
| https://en.wikipedia.org/wiki/GC-content | 200 | — | 14 |
| https://en.wikipedia.org/wiki/GC_box | 200 | — | 2 |
| https://en.wikipedia.org/wiki/GC_skew | 200 | — | 24 |
| https://en.wikipedia.org/wiki/GenBank | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Gene_prediction | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Gene_structure | 200 | — | 1 |
| https://en.wikipedia.org/wiki/General_feature_format | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Genetic_code | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Genotype_frequency | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Guide_RNA | 200 | — | 3 |
| https://en.wikipedia.org/wiki/H3K27ac | 200 | — | 1 |
| https://en.wikipedia.org/wiki/H3K27me3 | 200 | — | 1 |
| https://en.wikipedia.org/wiki/H3K36me3 | 200 | — | 1 |
| https://en.wikipedia.org/wiki/H3K4me1 | 200 | — | 1 |
| https://en.wikipedia.org/wiki/H3K4me3 | 200 | — | 1 |
| https://en.wikipedia.org/wiki/H3K9me3 | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Hamming_distance | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Haplotype_block | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Hardy-Weinberg_principle | 200 | Hardy–Weinberg principle | 3 |
| https://en.wikipedia.org/wiki/Heptad_repeat | 200 | — | 4 |
| https://en.wikipedia.org/wiki/HindIII | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Hybridization_probe | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Hydrophilicity_plot | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Informative_site | 200 | — | 6 |
| https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | 200 | — | 9 |
| https://en.wikipedia.org/wiki/Intron | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Inverted_repeat | 200 | — | 12 |
| https://en.wikipedia.org/wiki/Jaccard_index | 200 | — | 6 |
| https://en.wikipedia.org/wiki/K-means_clustering | 200 | — | 2 |
| https://en.wikipedia.org/wiki/K-mer | 200 | — | 39 |
| https://en.wikipedia.org/wiki/K-mer#Introduction | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Karyotype | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Kozak_consensus_sequence | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Levenshtein_distance | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | 200 | — | 11 |
| https://en.wikipedia.org/wiki/Linkage_disequilibrium | 200 | — | 2 |
| https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Longest_common_substring | 200 | — | 8 |
| https://en.wikipedia.org/wiki/Longest_repeated_substring_problem | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Metagenomics | 200 | — | 4 |
| https://en.wikipedia.org/wiki/MicroRNA | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Microsatellite | 200 | — | 10 |
| https://en.wikipedia.org/wiki/Minimum_evolution | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Minor_allele_frequency | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Minor_spliceosome | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Models_of_DNA_evolution | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Molecular_recognition_feature | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Most_recent_common_ancestor | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Multiple_sequence_alignment | 200 | — | 12 |
| https://en.wikipedia.org/wiki/Mutational_signatures | 200 | — | 1 |
| https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | 200 | — | 10 |
| https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#Choosing_a_scoring_system | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Neighbor_joining | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Neighbor_joining#Conclusion:_additive_distances | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Neighbor_joining#Example | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Neighbor_joining#First_branch_length_estimation | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Neighbor_joining#First_joining | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Newick_format | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Non-negative_least_squares | 200 | — | 7 |
| https://en.wikipedia.org/wiki/Non-negative_matrix_factorization | 200 | — | 7 |
| https://en.wikipedia.org/wiki/Nuclear_localization_sequence | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Nucleic_acid_notation | 200 | — | 21 |
| https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | 200 | — | 7 |
| https://en.wikipedia.org/wiki/Nucleotide_diversity | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Off-target_genome_editing | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Open_reading_frame | 200 | — | 8 |
| https://en.wikipedia.org/wiki/PROSITE | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Palindromic_sequence | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Pan-genome | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Phred_quality_score | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Phylogenetic_tree | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Ploidy | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Poisson_distribution#Random_variate_generation | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Polypyrimidine_tract | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Position_weight_matrix | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Pribnow_box | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Primer_(molecular_biology) | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Primer_dimer | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Promoter_(genetics) | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Protein_domain | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Protospacer_adjacent_motif | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Pseudoknot | 200 | — | 11 |
| https://en.wikipedia.org/wiki/PstI | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Quantile | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Quantile_normalization | 200 | — | 5 |
| https://en.wikipedia.org/wiki/RNA_splicing | 200 | — | 6 |
| https://en.wikipedia.org/wiki/Reading_frame | 200 | — | 16 |
| https://en.wikipedia.org/wiki/Relative_abundance_distribution | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Repeated_sequence_(DNA) | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Restriction_digest | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Restriction_enzyme | 200 | — | 11 |
| https://en.wikipedia.org/wiki/Restriction_enzyme#Examples | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Restriction_enzyme#Recognition_site | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Restriction_map | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Restriction_site | 200 | — | 6 |
| https://en.wikipedia.org/wiki/Ribosome-binding_site | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Robinson–Foulds_metric | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics) | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Sequence_alignment | 200 | — | 10 |
| https://en.wikipedia.org/wiki/Sequence_logo | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Sequence_motif | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Shine%E2%80%93Dalgarno_sequence | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Shine-Dalgarno_sequence | 200 | Shine–Dalgarno sequence | 2 |
| https://en.wikipedia.org/wiki/Signal_peptide | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | 200 | — | 7 |
| https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Example | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Species_evenness | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Species_richness | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Spliceosome | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Start_codon | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Stem-loop | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Sticky_and_blunt_ends | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Stop_codon | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Student%27s_t-distribution | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Substitution_model | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Suffix_tree | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Synteny | 200 | — | 10 |
| https://en.wikipedia.org/wiki/TATA_box | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Tajima%27s_D | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Tajima%27s_D#Example | 200 | — | 5 |
| https://en.wikipedia.org/wiki/Tandem_repeat | 200 | — | 14 |
| https://en.wikipedia.org/wiki/Telomere | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Tetraloop | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Transition_(genetics) | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Translation_(biology) | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Transversion | 200 | — | 8 |
| https://en.wikipedia.org/wiki/Tree_(abstract_data_type) | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Tree_(graph_theory) | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Tumor_microenvironment | 200 | — | 1 |
| https://en.wikipedia.org/wiki/UPGMA | 200 | — | 4 |
| https://en.wikipedia.org/wiki/UPGMA#The_UPGMA_dendrogram | 200 | — | 1 |
| https://en.wikipedia.org/wiki/UPGMA#Working_example | 200 | — | 3 |
| https://en.wikipedia.org/wiki/Variant_Call_Format | 200 | — | 2 |
| https://en.wikipedia.org/wiki/Watterson_estimator | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Welch%27s_t-test | 200 | — | 4 |
| https://en.wikipedia.org/wiki/Wobble_base_pair | 200 | — | 10 |
| https://en.wikipedia.org/wiki/Zinc_finger | 200 | — | 1 |
| https://en.wikipedia.org/wiki/Zygosity | 200 | — | 1 |
| https://www.addgene.org/guides/crispr/ | 200 | — | 3 |
| https://www.addgene.org/protocols/primer-design/ | 200 | — | 2 |
| https://www.addgene.org/protocols/restriction-digest/ | 200 | — | 2 |

## NUMBER and DEFINITION occurrences

| id | file:line | class | verdict | drives | page evidence / primary | scope |
|---|---|---|---|---|---|---|
| O001 | `docs/Evidence/ALIGN-GLOBAL-001-Evidence.md:12` | NUMBER | SUPPORTED | test | Match: +1 / Mismatch or Indel: −1 ... +−++−−+− −> 1*4 + (−1)*4 = 0 ... the first row being 0, −1, −2, −3, −4, −5, −6, −7 ... F_{0j} = d*j, F_{i0} = d*i · primary: Needleman & Wunsch 1970 J Mol Biol 48:443 (method); linear-gap recurrence as in Durbin et al. 1998 §2.3 — no primary listed in doc, not fetched — — (+1/−1/−1 scheme and GC | in-scope |
| O002 | `docs/Evidence/ALIGN-GLOBAL-001-Evidence.md:60` | DEFINITION | NOT_SUPPORTED_DEF | doc-only | Global alignments ... roughly equal size ... can start and/or end in gaps. A general global alignment technique is the Needleman–Wunsch algorithm | in-scope |
| O003 | `docs/Evidence/ALIGN-MULTI-001-Evidence.md:13` | NUMBER | SUPPORTED | doc-only | Multiple sequence alignment (MSA) is the process or the result of sequence alignment of three or more biological sequences ... BAliBase benchmark found that at least 24% ... | in-scope |
| O004 | `docs/Evidence/ALIGN-MULTI-001-Evidence.md:30` | NUMBER | NO_NUMBER | code | The gap opening penalty and gap extension penalty parameters can be adjusted by the user. ... \|asterisk ... fully conserved residue | in-scope |
| O005 | `docs/Evidence/ALIGN-MULTI-001-Evidence.md:45` | NUMBER | NO_NUMBER | code | the '''consensus sequence''' ... is the calculated sequence of most frequent residues, either [[nucleotide]] or [[amino acid]], found at each position in a [[sequence alignment]] | in-scope |
| O007 | `docs/Evidence/ALIGN-MULTI-001-Evidence.md:206` | DEFINITION | SUPPORTED_DEF | doc-only | optimizes the sum of all of the pairs of characters at each position in the alignment (the so-called sum of pair score) ... iterative methods can return to previously calculated pa | in-scope |
| O013 | `docs/Evidence/ANNOT-REPEAT-001-Evidence.md:13` | NUMBER | NOT_STATED | code | When between 10 and 60 nucleotides are repeated, it is called a minisatellite. Those with fewer are known as microsatellites or short tandem repeats. · primary: 'two or more': Benson G. (1999) Tandem repeats finder, NAR 27(2):573 (PMC148217); '1–6 bp': Tóth G. et al. (2000) Genome Res 10:967 (PMC310925) — Benson abstract: 'A tande | in-scope |
| O014 | `docs/Evidence/ANNOT-REPEAT-001-Evidence.md:26` | NUMBER | NO_NUMBER | code | a single stranded sequence of [[nucleotide]]s followed downstream by its reverse complement ... can be any length including zero. ... When the intervening length is zero, the compo | in-scope |
| O015 | `docs/Evidence/ANNOT-REPEAT-001-Evidence.md:85` | DEFINITION | SUPPORTED_DEF | test | repetitions are directly adjacent to each other, e.g. ATTCG ATTCG ATTCG, in which the sequence ATTCG is repeated three times | in-scope |
| O018 | `docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md:53` | DEFINITION | SUPPORTED_DEF | doc-only | the calculated sequence of most frequent residues, either nucleotide or amino acid, found at each position in a sequence alignment ... Y represents any pyrimidine, and R indicates  | in-scope |
| O020 | `docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md:49` | DEFINITION | SUPPORTED_DEF | doc-only | consisting of sequences separated by gaps of known length ... Bambus uses a greedy algorithm, defined as such because it joins together contigs with the most links first | in-scope |
| O022 | `docs/Evidence/ASSEMBLY-STATS-001-Evidence.md:23` | NUMBER | SUPPORTED | test | Assembly A contains six contigs of lengths 80 kbp, 70 kbp, 50 kbp, 40 kbp, 30 kbp, and 20 kbp. The sum size of assembly A is 290 kbp, the N50 contig length is 70 kbp ... The size o · primary: Miller JR, Koren S, Sutton G (2010) Genomics 95:315-327 (PMC2874646; listed in same Evidence) — §1.2: 'The contig N50 is the length of the smallest contig in the set that  | in-scope |
| O023 | `docs/Evidence/ASSEMBLY-STATS-001-Evidence.md:23` | NUMBER | SUPPORTED | test | Assembly A contains six contigs of lengths 80 kbp, 70 kbp, 50 kbp, 40 kbp, 30 kbp, and 20 kbp. The sum size of assembly A is 290 kbp, the N50 contig length is 70 kbp ... assembly B · primary: Miller, Koren & Sutton (2010) Genomics 95(6):315-327, PMC2874646 (same Evidence doc) — §1.2: 'The contig N50 is the length of the smallest contig in the set that contains  | in-scope |
| O027 | `docs/Evidence/CHROM-CENT-001-Evidence.md:13` | NUMBER | NOT_STATED | code | α-satellite (or alphoid), although a number of other sequence types are found in this region. · primary: Hartley G, O'Neill RJ (2019) Alpha satellite DNA biology (PMC6121732), citing Willard 1985 / Waye & Willard 1987 — already listed in the same Evidence file (line ~117) — v | in-scope |
| O028 | `docs/Evidence/CHROM-CENT-001-Evidence.md:24` | NUMBER | SUPPORTED | doc-only | \| A \| 1–3 \| Large, metacentric or submetacentric ... \| D \| 13–15 \| Medium-sized, acrocentric, with satellite ... \| G \| 21–22, Y \| Very small, acrocentric | in-scope |
| O029 | `docs/Evidence/CHROM-CENT-001-Evidence.md:33` | NUMBER | NO_NUMBER | code | Constitutive heterochromatin, which is never expressed. It is located around the centromere and usually contains repetitive sequences. ... p arms (from the French petit, small) ... · primary: PMC6121732 (Willard 1985; Waye & Willard 1987), same Evidence doc line 117 — 'Alpha satellite DNA is composed of fundamental 171bp monomeric repeat units.' (171 bp) | in-scope |
| O030 | `docs/Evidence/CHROM-KARYO-001-Evidence.md:15` | NUMBER | SUPPORTED | doc-only | contains 22 pairs of autosomal chromosomes and one pair of sex chromosomes ... The most common karyotypes for females contain two X chromosomes and are denoted 46,XX | in-scope |
| O031 | `docs/Evidence/CHROM-KARYO-001-Evidence.md:25` | NUMBER | SUPPORTED | code | monoploid (1 set), diploid (2 sets), triploid (3 sets), tetraploid (4 sets) ... 23 homologous pairs, providing a full complement of 46 chromosomes · primary: standard karyotype facts (ISCN); not re-fetched — - (2n=46, n=23) | in-scope |
| O032 | `docs/Evidence/CHROM-KARYO-001-Evidence.md:35` | NUMBER | NOT_STATED | code | \| 1 \|\| Monosomy ... \| 2 \|\| Disomy ... \| 3 \|\| Trisomy ... \| 4/5 \|\| tetrasomy/pentasomy \|\| ''Tetrasomy'' and ''pentasomy'' are the presence of four or five copies of a  · primary: not found (NCBI Bookshelf Griffiths 'Introduction to Genetic Analysis' blocked by captcha; no OA primary located in 3 fetches) — - (0 copies = nullisomy (unverified)) | in-scope |
| O036 | `docs/Evidence/CHROM-TELO-001-Evidence.md:16` | NUMBER | NO_NUMBER | doc-only | In many species, the sequence repeats are enriched in [[guanine]], e.g. TTAGGG in [[vertebrate]]s,<ref name=pmid2780561>{{cite journal \| vauthors = Meyne J ... | in-scope |
| O037 | `docs/Evidence/CODON-CAI-001-Evidence.md:13` | NUMBER | NO_NUMBER | doc-only | is the most widespread technique for analyzing codon usage bias ... CAI measures the deviation of a given protein coding gene sequence with respect to a reference set of genes | OUT-OF-SCOPE:B02 |
| O038 | `docs/Evidence/CODON-OPT-001-Evidence.md:13` | NUMBER | NO_NUMBER | doc-only | heterologous expression consider global nucleotide content such as local mRNA folding, codon pair bias, a codon ramp, codon harmonization ... presence of rare codons ... can lead t | OUT-OF-SCOPE:B02 |
| O039 | `docs/Evidence/CODON-OPT-001-Evidence.md:21` | NUMBER | NOT_STATED | test | w_i = f_i / max(f_j) ... The CAI of a gene is simply defined as the geometric mean of the weight associated to each codon over the length (L) · primary: no primary fetched (Sharp & Li 1987 PMC340524 has no OA full text); bound follows mathematically from the page's own definition (w_i ≤ 1 ⇒ geometric mean ≤ 1) — - (0 < CAI | OUT-OF-SCOPE:B02 |
| O040 | `docs/Evidence/CODON-RARE-001-Evidence.md:15` | NUMBER | NO_NUMBER | doc-only | Amino acids incorporation at more frequent codons occurs faster than that of rare codons. | OUT-OF-SCOPE:B02 |
| O041 | `docs/Evidence/CODON-STATS-001-Evidence.md:25` | NUMBER | NO_NUMBER | code | the frequency of the most frequent synonymous codon f_j for that amino acid ... CAI= (\Pi_{i=1}^{L} w_{i})^{\frac{1}{L}} | OUT-OF-SCOPE:B02 |
| O042 | `docs/Evidence/CODON-USAGE-001-Evidence.md:17` | NUMBER | SUPPORTED | doc-only | There are 64 different codons (61 codons coding for amino acids and 3 stop codons) but only 20 different translated amino acids. The overabundance in the number of codons allows ma | OUT-OF-SCOPE:B02 |
| O043 | `docs/Evidence/COMPGEN-COMPARE-001-Evidence.md:39` | NUMBER | NOT_STATED | doc-only | conservation of blocks of order within two sets of chromosomes that are being compared with each other. These blocks are referred to as syntenic blocks. · primary: not found: ScienceDirect topic page fetched but the phrase / 50.2% / Acropora not present in retrievable text — - (-) | in-scope |
| O045 | `docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md:25` | NUMBER | NO_NUMBER | doc-only | When the residues of both sequences match at the same location on the plot, a dot is drawn at the corresponding position. ... The main diagonal represents the sequence's alignment  | in-scope |
| O047 | `docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md:47` | DEFINITION | SUPPORTED_DEF | doc-only | Shared synteny (also known as conserved synteny) describes preserved co-localization of genes ... dynamic programming is used to select the best scoring path of shared homologous g | in-scope |
| O049 | `docs/Evidence/DISORDER-MORF-001-Evidence.md:29` | NUMBER | SUPPORTED | code | Molecular recognition features (MoRFs) are small (10-70 residues) intrinsically disordered regions in proteins that undergo a disorder-to-order transition upon binding to their par · primary: Mohan et al. 2006 J Mol Biol 362(5):1043-59 (same doc, line 20) — 'relatively short (10-70 residues), loosely structured protein regions' (10–70 residues) | in-scope |
| O051 | `docs/Evidence/DISORDER-PRED-001-Evidence.md:28` | NUMBER | SUPPORTED | code | Based on DISOPRED2 prediction, long (>30 residue) disordered segments occur in 2.0% of archaean, 4.2% of eubacterial and 33.0% of eukaryotic proteins · primary: Ward JJ et al. (2004) J Mol Biol 337:635-645 (PMID 15019783) — abstract: 'Putative, long (>30 residue) disordered segments are found to occur in 2.0% of archaean, 4.2% of  | in-scope |
| O052 | `docs/Evidence/DISORDER-PRED-001-Evidence.md:44` | NUMBER | SUPPORTED | code | Isoleucine \| 4.5 \| Valine \| 4.2 \| Leucine \| 3.8 \| Phenylalanine \| 2.8 \| Cysteine \| 2.5 \| Methionine \| 1.9 \| Alanine \| 1.8 ... Arginine \| -4.5 · primary: Kyte & Doolittle (1982) J Mol Biol 157:105-132 (paywalled); reference implementation Biopython Bio/SeqUtils/ProtParamData.py `kd` — ProtParamData.py line 16-19 (I 4.5, V 4 | in-scope |
| O053 | `docs/Evidence/DISORDER-PRED-001-Evidence.md:81` | NUMBER | SUPPORTED | code | Histidine's imidazole group has a pKa of 6.0, and is only around 10% protonated at neutral pH. (table: Hydropathy index col. cites Kyte J, Doolittle RF; Arg −4.5, Ile 4.5) · primary: Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132 (cited in DisorderPredictor.cs:1106) — Table 2 hydropathy scale (abstract fetched, scale values not in abstract; not re- | in-scope |
| O057 | `docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md:31` | NUMBER | SUPPORTED | code | The following hydrophilic, charged amino acids A, R, G, Q, S, P, E and K have been characterized as disorder-promoting ... W, C, F, I, Y, V, L, and N ... H, M, T and D are ambiguou · primary: Dunker AK et al. (2001) J Mol Graph Model 19:26-59 (PMID 11381529; cited by the page and in the same Evidence file) — order/disorder-promoting residue classes (not re-fetc | in-scope |
| O059 | `docs/Evidence/DISORDER-REGION-001-Evidence.md:89` | NUMBER | SUPPORTED | code | long (>30 residue) disordered segments occur in 2.0% of archaean, 4.2% of eubacterial and 33.0% of eukaryotic proteins · primary: Ward JJ et al. (2004) J Mol Biol 337:635-645 (PMID 15019783) — abstract: 'Putative, long (>30 residue) disordered segments are found to occur in 2.0% of archaean, 4.2% ... | in-scope |
| O060 | `docs/Evidence/EPIGEN-CHROM-001-Evidence.md:59` | DEFINITION | SUPPORTED_DEF | doc-only | it is highly enriched at active promoters near transcription start sites (TSS) <ref>... Liang G, Lin JC ... | in-scope |
| O061 | `docs/Evidence/EPIGEN-CHROM-001-Evidence.md:71` | DEFINITION | SUPPORTED_DEF | doc-only | H3K4me1 is enriched at active and primed enhancers. ... Nature Genetics \|volume=50 \|issue=1 \|pages=4–5 \|year=2018 \|last1=Rada-Iglesias | in-scope |
| O062 | `docs/Evidence/EPIGEN-CHROM-001-Evidence.md:83` | DEFINITION | SUPPORTED_DEF | doc-only | H3K27ac which is an active enhancer mark ... Creyghton ... Histone H3K27ac separates active from poised enhancers and predicts developmental state ... PNAS ... 2010 ... volume=107  | in-scope |
| O063 | `docs/Evidence/EPIGEN-CHROM-001-Evidence.md:95` | DEFINITION | SUPPORTED_DEF | doc-only | This tri-methylation is associated with the downregulation of nearby genes via the formation of heterochromatic regions. ... Ferrari KJ ... PRC2, mediates the tri-methylation of hi | in-scope |
| O064 | `docs/Evidence/EPIGEN-CHROM-001-Evidence.md:107` | DEFINITION | SUPPORTED_DEF | doc-only | the tri-methylation at the 9th lysine residue of the histone H3 protein and is often associated with heterochromatin. ... Nicetto D, Donahue G ... | in-scope |
| O065 | `docs/Evidence/EPIGEN-CHROM-001-Evidence.md:119` | NUMBER | NO_NUMBER | code | a mark that indicates the tri-[[methylation]] at the 36th [[lysine]] residue of the histone H3 protein and often associated with [[gene\|gene bodies]]. | in-scope |
| O066 | `docs/Evidence/EPIGEN-CPG-001-Evidence.md:13` | NUMBER | SUPPORTED | code | a region with at least 200 bp, a GC percentage greater than 50%, and an observed-to-expected CpG ratio greater than 60% ... greater than 500 bp ... GC content greater than 55% ...  · primary: Gardiner-Garden M, Frommer M (1987) J Mol Biol 196:261-282 (listed in same doc line 28) — paper body (abstract fetched via EuropePMC does not state the thresholds; body no | in-scope |
| O068 | `docs/Evidence/EPIGEN-DMR-001-Evidence.md:64` | NUMBER | SUPPORTED | test | \tfrac{10!~14!~12!~12!}{1!~9!~11!~3!~24!} \approx 0.001346 ... one-tailed test, with p approximately 0.001346076 + 0.000033652 = 0.001379728 · primary: hypergeometric probability (Fisher 1922/1935) — verified by direct computation C(10,1)C(14,11)/C(24,12) — computed 2026-10-14 (0.0013460761879) | in-scope |
| O069 | `docs/Evidence/EPIGEN-DMR-001-Evidence.md:96` | DEFINITION | SUPPORTED | test | (re-fetched real cited page Fisher's_exact_test) ... = 10! 14! 12! 12! / (1! 9! 11! 3! 24!) ≈ 0.001346076 · primary: Fisher RA (1922/1935) hypergeometric formula; 0.001346076 is pure arithmetic of Wikipedia's own 'studying by gender' table — Fisher's exact test, Example section (0.001346 | in-scope |
| O071 | `docs/Evidence/GENOMIC-COMMON-001-Evidence.md:13` | NUMBER | SUPPORTED | test | The strings "ABABC", "BABCA" and "ABCBA" have only one longest common substring, viz. "ABC" of length 3. ... in Θ(n+m) time with the help of a generalized suffix tree · primary: Gusfield (1997) Algorithms on Strings, Trees and Sequences (cited on page, not re-fetched) — generalized suffix tree LCSubstring (Θ(n+m); examples ADA/ANA and ABC are Wiki | in-scope |
| O074 | `docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md:38` | NUMBER | SUPPORTED | test | the common restriction enzyme EcoRI recognizes the palindromic sequence GAATTC and cuts between the G and the A on both the top and bottom strands. · primary: REBASE / NEB EcoRI G^AATTC (not re-fetched) — - (GAATTC, 6 bp) | in-scope |
| O076 | `docs/Evidence/GENOMIC-ORF-001-Evidence.md:27` | NUMBER | SUPPORTED | doc-only | reading frames are defined as spans of DNA sequence between the start and stop codons ... there are six possible frame translations ... some authors say that an ORF should have a m | in-scope |
| O078 | `docs/Evidence/GENOMIC-REPEAT-001-Evidence.md:26` | DEFINITION | SUPPORTED | test | In the figure with the string "ATCGATCGA$", the longest substring that repeats at least twice is "ATCGA". ... linear time Θ(n) ... deepest internal node in the tree with more than  · primary: Worked example is Wikipedia's own illustration (trivially verifiable); algorithm from Gusfield 1997 / CMU notes cited elsewhere in the same Evidence doc — article lead (AT | in-scope |
| O080 | `docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md:13` | NUMBER | NO_NUMBER | code | The Jaccard index measures similarity between finite non-empty sample sets and is defined as the size of the intersection divided by the size of the union ... By definition, 0 \le  · primary: Jaccard P (1901) Bull Soc Vaudoise Sci Nat 37(142):547-579 (cited by page and doc) — definition; not fetched (no numeric constant (range 0..1)) | in-scope |
| O082 | `docs/Evidence/GENOMIC-TANDEM-001-Evidence.md:26` | NUMBER | SUPPORTED | test | e.g. ATTCG ATTCG ATTCG, in which the sequence ATTCG is repeated three times. ... Tandem repeats constitute about 8% ... implicated in more than 50 lethal human diseases ... on the  · primary: MeSH 'Tandem Repeat Sequences' (cited on page for the definition/example); not re-fetched — - (ATTCG ×3, period 5, length 15) | in-scope |
| O084 | `docs/Evidence/KMER-ASYNC-001-Evidence.md:13` | NUMBER | SUPPORTED | test | a sequence of length L will have L - k + 1 [k-mers] ... n^{k} total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA) ... AGAT would have four  · primary: Combinatorial identity (L-k+1 windows; 4^k strings); Compeau & Pevzner (2015) Bioinformatics Algorithms (cited in KmerAnalyzer tests) — not fetched; hand-derivable (L-k+1, | in-scope |
| O085 | `docs/Evidence/KMER-ASYNC-001-Evidence.md:133` | DEFINITION | SUPPORTED | test | the sequence AGAT would have ... three 2-mers (AG, GA, AT) ... a sequence of length L will have L - k + 1 k-mers ... n^{k} possible k-mers ... \|2 \|GT, TA, AG, GA, AG, GC, CT, TG, · primary: Definition-level / arithmetic: counts derive directly from Wikipedia's GTAGAGCTGT k-mer table (k=2: 9 total, 7 distinct; k=3: 8/8; k=4: 7/7); n^k universe size is combinat | in-scope |
| O087 | `docs/Evidence/KMER-GENERATE-001-Evidence.md:13` | NUMBER | SUPPORTED | test | a sequence of length L will have L - k + 1 [k-mers] ... n^{k} total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA) ... AGAT would have four  · primary: Combinatorial identity (L-k+1 windows; 4^k strings); Compeau & Pevzner (2015) Bioinformatics Algorithms (cited in KmerAnalyzer tests) — not fetched; hand-derivable (L-k+1, | in-scope |
| O089 | `docs/Evidence/KMER-POSITIONS-001-Evidence.md:26` | NUMBER | SUPPORTED | test | a sequence of length L will have L - k + 1 [k-mers] ... n^{k} total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA) ... AGAT would have four  · primary: Combinatorial identity (L-k+1 windows; 4^k strings); Compeau & Pevzner (2015) Bioinformatics Algorithms (cited in KmerAnalyzer tests) — not fetched; hand-derivable (L-k+1, | in-scope |
| O090 | `docs/Evidence/KMER-POSITIONS-001-Evidence.md:87` | DEFINITION | SUPPORTED | test | the sequence AGAT would have ... three 2-mers (AG, GA, AT) ... a sequence of length L will have L - k + 1 k-mers ... n^{k} possible k-mers ... \|2 \|GT, TA, AG, GA, AG, GC, CT, TG, · primary: Definition-level / arithmetic: counts derive directly from Wikipedia's GTAGAGCTGT k-mer table (k=2: 9 total, 7 distinct; k=3: 8/8; k=4: 7/7); n^k universe size is combinat | in-scope |
| O092 | `docs/Evidence/KMER-STATS-001-Evidence.md:13` | NUMBER | SUPPORTED | test | a sequence of length L will have L - k + 1 [k-mers] ... n^{k} total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA) ... AGAT would have four  · primary: Combinatorial identity (L-k+1 windows; 4^k strings); Compeau & Pevzner (2015) Bioinformatics Algorithms (cited in KmerAnalyzer tests) — not fetched; hand-derivable (L-k+1, | in-scope |
| O093 | `docs/Evidence/KMER-STATS-001-Evidence.md:94` | DEFINITION | SUPPORTED | test | the sequence AGAT would have ... three 2-mers (AG, GA, AT) ... a sequence of length L will have L - k + 1 k-mers ... n^{k} possible k-mers ... \|2 \|GT, TA, AG, GA, AG, GC, CT, TG, · primary: Definition-level / arithmetic: counts derive directly from Wikipedia's GTAGAGCTGT k-mer table (k=2: 9 total, 7 distinct; k=3: 8/8; k=4: 7/7); n^k universe size is combinat | in-scope |
| O095 | `docs/Evidence/KMER-UNIQUE-001-Evidence.md:13` | NUMBER | SUPPORTED | test | a sequence of length L will have L - k + 1 [k-mers] ... n^{k} total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA) ... AGAT would have four  · primary: Combinatorial identity (L-k+1 windows; 4^k strings); Compeau & Pevzner (2015) Bioinformatics Algorithms (cited in KmerAnalyzer tests) — not fetched; hand-derivable (L-k+1, | in-scope |
| O097 | `docs/Evidence/META-ALPHA-001-Evidence.md:17` | DEFINITION | SUPPORTED_DEF | doc-only | ==Shannon index== ... ==Simpson index== The Simpson index was introduced in 1949 by Edward H. Simpson | in-scope |
| O098 | `docs/Evidence/META-ALPHA-001-Evidence.md:22` | DEFINITION | SUPPORTED_DEF | doc-only | '''alpha diversity''' ('''α-diversity''') is the mean [[species diversity]] in a site at a local scale. | in-scope |
| O099 | `docs/Evidence/META-ALPHA-001-Evidence.md:27` | DEFINITION | SUPPORTED_DEF | doc-only | Species richness is simply a count of species ... species diversity is a scientific assessment of both species richness and ... | in-scope |
| O100 | `docs/Evidence/META-ALPHA-001-Evidence.md:32` | DEFINITION | NOT_SUPPORTED_DEF | code | The Shannon index ... H = - sum_{i=1}^S p_i ln p_i · primary: Pielou EC (1966) J Theor Biol 13:131-144, listed in the same Evidence (item 9, line 50) + scikit-bio 0.7.4 pielou_e cross-check (line 179) — META-ALPHA-001-Evidence.md lin | in-scope |
| O104 | `docs/Evidence/META-BETA-001-Evidence.md:48` | DEFINITION | SUPPORTED | test | This leads to BC_{jk} = 1 - (2 \times 10) / (17 + 16) = 0.39 · primary: Bray JR & Curtis JT (1957) Ecol Monogr 27:325 — formula only; 0.39 is the arithmetic of Wikipedia's own aquarium example — formula BC=1-2C/(S_j+S_k) (13/33 ≈ 0.3939 (deriv | in-scope |
| O105 | `docs/Evidence/META-BIN-001-Evidence.md:8` | NUMBER | NO_NUMBER | doc-only | Binning methods can be based on either compositional sequence features (such as GC-content or tetranucleotide frequencies) or sequence read mapping coverage across samples, or both | in-scope |
| O106 | `docs/Evidence/META-CLASS-001-Evidence.md:22` | NUMBER | NO_NUMBER | doc-only | Binning is the process of associating a particular sequence with an organism. In similarity-based binning, methods such as BLAST are used to rapidly search for phylogenetic markers | in-scope |
| O107 | `docs/Evidence/META-PROF-001-Evidence.md:20` | NUMBER | NO_NUMBER | doc-only | MetaPhlAn, Meteor and AMPHORA are methods based on unique clade-specific markers for estimating organismal relative abundances ... Metagenomes are most commonly compared by analyzi | in-scope |
| O108 | `docs/Evidence/META-PROF-001-Evidence.md:28` | NUMBER | NOT_STATED | test | the relative abundance distribution (RAD) or species abundance distribution (SAD) describes the relationship between the number of species observed in a field study as a function o · primary: Definitional (relative abundance = n_i / Σn); MetaPhlAn (Segata et al. 2012 Nat Methods 9:811) reports per-rank relative abundances in % summing to 100 — not fetched — — ( | in-scope |
| O109 | `docs/Evidence/META-TAXA-001-Evidence.md:13` | NUMBER | SUPPORTED | test | take each tortoise in turn ... U_T = 11 ... we get 5, 5, 5, 5, 5, 0, so U_H = 25 ... sum ... U = 36, which is 6×6 ... 12 + 6 + 5 + 4 + 3 + 2 = 32 · primary: SciPy scipy.stats.mannwhitneyu docs (same Evidence doc §2) / Mann & Whitney 1947 Ann Math Stat 18:50-60 — U statistic + normal approx with tie correction (U1 = R1 − n1(n1+ | in-scope |
| O111 | `docs/Evidence/MIRNA-PAIR-001-Evidence.md:54` | DEFINITION | SUPPORTED_DEF | doc-only | A wobble base pair is a pairing between two nucleotides in RNA molecules that does not follow Watson–Crick base pair rules ... table: G pairs C or U; cites Crick FH, Codon--anticod | in-scope |
| O113 | `docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md:13` | NUMBER | SUPPORTED | code | molecules containing 21–23 nucleotides ... These hairpin loop structures are composed of about 70 nucleotides ... The product resulting has a two-nucleotide overhang at its 3' end · primary: Han J et al. (2006) Cell 125:887 / Lee Y et al. (2003) (2-nt overhang); Park JE et al. (2011) Nature 475:201 (~22 nt); Bartel (2004) Cell 116:281 (~70 nt) — all already ci | in-scope |
| O117 | `docs/Evidence/MOTIF-CONS-001-Evidence.md:13` | DEFINITION | SUPPORTED_DEF | doc-only | the calculated sequence of most frequent residues, either nucleotide or amino acid, found at each position in a sequence alignment ... Y represents any pyrimidine, and R indicates  | in-scope |
| O119 | `docs/Evidence/MOTIF-GENERATE-001-Evidence.md:37` | NUMBER | NO_NUMBER | code | Nomenclature Committee of the International Union of Biochemistry (NC-IUB)\|year=1984\|title=Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences\|journal=Nuclei · primary: NC-IUB (1984) Nucleic Acids Res 13(9):3021-3030 — Table 1 (IUPAC codes (no numeric value)) | in-scope |
| O121 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:13` | NUMBER | SUPPORTED | code | The TATA box is usually located 25-35 base pairs upstream of the transcription start site. ... consensus sequence 5'-TATA(A/T)A(A/T)-3' ... (yeast) 5'-TATA(A/T)A(A/T)(A/G)-3' · primary: Bucher 1990 J Mol Biol 212:563-578 (cited in same entry) — abstract: weight matrices from 502 unrelated Pol II promoters (502; TATA core TATAAA) | in-scope |
| O122 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:25` | NUMBER | SUPPORTED | code | The Pribnow box ... is a sequence of TATAAT of six nucleotides · primary: Harley CB, Reynolds RP (1987) Nucleic Acids Res 15(5):2343-2361 (in same doc; abstract fetched via EuropePMC) — Abstract ('all bases in the -35 (TTGACA) and -10 (TATAAT) h | in-scope |
| O123 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:48` | NUMBER | NO_NUMBER | code | It has a [[consensus sequence]] GGGCGG which is position-dependent and orientation-independent. | in-scope |
| O124 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:59` | NUMBER | SUPPORTED | code | The Kozak sequence was determined by sequencing of 699 vertebrate mRNAs ... defined as 5'-(gcc)gccRccAUGG-3' ... nucleotides at positions +4 (i.e. G in the consensus) and −3 (i.e.  · primary: Kozak M (1987) NAR 15(20):8125-8148, 'An analysis of 5'-noncoding sequences from 699 vertebrate messenger RNAs' (cited) — title/abstract (699 mRNAs; gccRccAUGG) | in-scope |
| O125 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:71` | NUMBER | SUPPORTED | code | The six-base [[consensus sequence]] is '''AGGAGG''' ... ACCUCCU ... (~8 bases upstream of the start codon AUG) · primary: Shine J, Dalgarno L (1975) Nature 254:34-38 / (1974) PNAS 71:1342 — cited in GenomeAnnotator.cs and same Evidence — 16S rRNA 3' ...ACCUCCUUA complementary to AGGAGG (AGGAG | in-scope |
| O126 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:94` | NUMBER | NO_NUMBER | code | Its specific DNA sequence, CANNTG (where N can be any [[nucleotide]]), with a palindromic [[canonical sequence]] of CACGTG | in-scope |
| O133 | `docs/Evidence/ONCO-CCF-001-Evidence.md:57` | NUMBER | NO_NUMBER | doc-only | Assignment step: Assign each observation to the cluster with the nearest mean (centroid): that with the least squared Euclidean distance. ... Update step: Recalculate means (centro | in-scope |
| O135 | `docs/Evidence/ONCO-FUSION-001-Evidence.md:113` | NUMBER | NO_NUMBER | code | There are three reading frames that can be read in this 5′→3′ direction, each beginning from a different nucleotide in a triplet. · primary: Lodish, Molecular Cell Biology 6th ed. p.121 (cited on page; not re-fetched) — - (codon = 3 nt) | in-scope |
| O137 | `docs/Evidence/ONCO-FUSION-003-Evidence.md:56` | NUMBER | NO_NUMBER | code | as a sequence of triplets. Where these triplets equate to amino acids or stop signals during translation · primary: Lodish 6th ed. p.121 / Badger & Olsen 1999 Mol Biol Evol 16(4):512–24 (cited on page) — - (triplet = 3) | in-scope |
| O140 | `docs/Evidence/ONCO-IMMUNE-001-Evidence.md:13` | NUMBER | NOT_STATED | code | (page is a 3 kB stub) 'CIBERSORT, and its successor called CIBERSORTx, are a set of bioinformatics tools used to deconvolute cell type proportions and gene expression profiles from · primary: Newman et al. 2015 Nat Methods 12:453 via Chen et al. 2018 CIBERSORT protocol, PMC5895181 (Evidence §8 of same doc) — Evidence line 142: 'a signature matrix file consistin | in-scope |
| O141 | `docs/Evidence/ONCO-IMMUNE-001-Evidence.md:74` | DEFINITION | SUPPORTED_DEF | doc-only | cytotoxic T cells (CD8+), T helper 1 (Th-1) cells and regulatory T cells (Tregs) ... M1-polarized ... anti-inflammatory M2-polarized macrophages ... A high concentration is general | in-scope |
| O142 | `docs/Evidence/ONCO-IMMUNE-001-Evidence.md:86` | NUMBER | NO_NUMBER | doc-only | This score is a Kolmogorov–Smirnov-like statistic ... ES=P(S, i)=P_{hit}(S, i)-P_{miss}(S,i)=max(\|P_{hit} | in-scope |
| O143 | `docs/Evidence/ONCO-SIG-001-Evidence.md:73` | DEFINITION | SUPPORTED_DEF | doc-only | [[DNA]] \|\| adenine(A), thymine(T), guanine(G), cytosine(C)\|\| A = T, G ≡ C | in-scope |
| O145 | `docs/Evidence/ONCO-SIG-002-Evidence.md:49` | NUMBER | NO_NUMBER | code | The first widely used algorithm for solving this problem is an active-set method published by Lawson and Hanson in their 1974 book Solving Least Squares Problems · primary: Lawson CL, Hanson RJ (1974) Solving Least Squares Problems, Ch. 23 — algorithm NNLS (book; not fetched) (no numeric constant) | in-scope |
| O146 | `docs/Evidence/ONCO-SIG-002-Evidence.md:107` | DEFINITION | SUPPORTED_DEF | doc-only | Current algorithms are sub-optimal in that they only guarantee finding a local minimum, rather than a global minimum of the cost function; ... are matrices of ones when V = WH | in-scope |
| O149 | `docs/Evidence/ONCO-SIG-003-Evidence.md:76` | NUMBER | SUPPORTED | code | R‑7, Excel, Python, SciPy‑(1,1), Julia-(1,1), Maple‑6, NumPy, Guava \| (N − 1)p + 1 \| Linear interpolation of the modes ... R‑7. This is R's and Julia's default method. · primary: Hyndman RJ, Fan Y (1996) Am Stat 50(4):361-365 (cited by page and doc) — Definition 7 (not fetched) (h=(N-1)p+1 (1-based) = p(n-1) 0-based) | in-scope |
| O153 | `docs/Evidence/ONCO-VAF-001-Evidence.md:24` | NUMBER | SUPPORTED | code | a 95% confidence interval requires alpha = 0.05, thereby producing z_{.05} = 1.96 · primary: Standard normal quantile (Wilson 1927 interval uses z = Phi^-1(1-alpha/2)) — Phi^-1(0.975) (1.959964 ~ 1.96) | in-scope |
| O155 | `docs/Evidence/PANGEN-CORE-001-Evidence.md:42` | NUMBER | SUPPORTED | code | If α ≤ 1 then the pangenome is considered open. Otherwise (α > 1) it is considered closed. ... 2000 genomes is composed by 89,000 different gene families ... drops to zero when the · primary: micropan heaps() / Tettelin et al. (2008) — already cited in PanGenomeAnalyzer.cs (HeapsAlphaOpenThreshold) — heaps(): alpha < 1 open, alpha > 1 closed (alpha threshold 1. | in-scope |
| O157 | `docs/Evidence/PANGEN-MARKER-001-Evidence.md:48` | NUMBER | SUPPORTED | code | an informative site is a position in the relevant set of aligned sequences at which there are at least two different character states and each of those states occurs in at least tw · primary: Zvelebil & Baum (2008) Understanding Bioinformatics, Garland (cited by page & doc) — textbook, not fetched — definition of parsimony-informative site (≥2 states, each in ≥ | in-scope |
| O158 | `docs/Evidence/PANGEN-MARKER-001-Evidence.md:80` | DEFINITION | SUPPORTED_DEF | doc-only | a position in the relevant set of aligned sequences at which there are at least two different character states and each of those states occurs in at least two of the sequences. | in-scope |
| O160 | `docs/Evidence/PARSE-BED-001-Evidence.md:50` | DEFINITION | SUPPORTED_DEF | doc-only | These columns must be separated by spaces or tabs, the latter being recommended for reasons of compatib[ility] ... Each row of a file must have the same number of columns ... was p | in-scope |
| O162 | `docs/Evidence/PARSE-FASTA-001-Evidence.md:13` | DEFINITION | SUPPORTED_DEF | doc-only | The description line (defline) or header/identifier line, which begins with [the] greater-than character (">") followed by a description of the sequence (all in a single line) | in-scope |
| O163 | `docs/Evidence/PARSE-FASTA-001-Evidence.md:54` | DEFINITION | SUPPORTED_DEF | doc-only | Table: Weak W A,T; Strong S C,G; Amino M A,C; Keto K G,T; Purine R A,G; Pyrimidine Y C,T; Not A B; Not C D; Not G H; Not T V; Any N. Cites NC-IUB NAR 13(9):3021-3030 and IUPAC-IUB  | in-scope |
| O164 | `docs/Evidence/PARSE-FASTQ-001-Evidence.md:15` | NUMBER | NOT_STATED | code | A FASTQ file has four line-separated fields per sequence ... Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126 ... Starting with Illumina 1.3 and be · primary: FastQC PhredEncoding.java (already listed as source 5 in the same Evidence file) — getFastQEncodingOffset: lowestChar < 33 -> error; Phred64 path requires lowestChar >= IL | in-scope |
| O166 | `docs/Evidence/PARSE-GENBANK-001-Evidence.md:22` | NUMBER | NOT_STATED | doc-only | As of GenBank release 271.0 (April 2026), the database contained 53.90 trillion bases and 6.27 billion sequence records | in-scope |
| O167 | `docs/Evidence/PARSE-GFF-001-Evidence.md:13` | NUMBER | NOT_STATED | doc-only | All GFF formats (GFF2, GFF3 and GTF) are tab delimited with 9 fields per line. They all share the same structure for the first 7 fields, while differing in the content and format o | in-scope |
| O168 | `docs/Evidence/PARSE-VCF-001-Evidence.md:13` | NUMBER | SUPPORTED | code | The format was developed in 2010 for the 1000 Genomes Project ... latest_release_version = 4.5 ... latest_release_date 2024\|10\|09 ... tab separated into 8 mandatory columns ... P · primary: VCFv4.5 specification (samtools/hts-specs VCFv4.5.tex) — l.362 'The mandatory header line names the 8 fixed, mandatory columns'; l.386 'There are 8 fixed fields per record | in-scope |
| O172 | `docs/Evidence/PHYLO-DIST-001-Evidence.md:14` | NUMBER | SUPPORTED | code | \hat{d}=-{3\over4} \ln({1-{4\over3}p}) ... K = - {1\over2}\ln((1-2p-q) \sqrt{1-2q}) ... equilibrium base frequencies are equal (π_A = π_G = π_C = π_T =0.25 · primary: Jukes & Cantor 1969; Kimura 1980 J Mol Evol 16:111-120 (cited on page; not re-fetched) — - (JC69 -3/4 ln(1-4p/3); K80 -1/2 ln((1-2P-Q)√(1-2Q))) | in-scope |
| O173 | `docs/Evidence/PHYLO-DIST-001-Evidence.md:30` | NUMBER | NO_NUMBER | doc-only | When the model is time reversible, this can be performed between any two sequences, even if one is not the ancestor of the other | in-scope |
| O174 | `docs/Evidence/PHYLO-DIST-001-Evidence.md:41` | NUMBER | NO_NUMBER | doc-only | Distance is often defined as the fraction of mismatches at aligned positions, with gaps either ignored or counted as mismatches. | in-scope |
| O175 | `docs/Evidence/PHYLO-NEWICK-001-Evidence.md:16` | NUMBER | SUPPORTED | doc-only | \|title=Interpretation of "Newick's 8:45" Tree Format \|url=https://phylipweb.github.io/phylip/newick_doc.html | in-scope |
| O177 | `docs/Evidence/PHYLO-STATS-001-Evidence.md:13` | NUMBER | NO_NUMBER | code | A ''leaf'' is a vertex with no children. ... The ''depth'' of a vertex is the length of the path to its root ... The root has depth zero ... a tree with only a single vertex ... ha | in-scope |
| O178 | `docs/Evidence/PHYLO-STATS-001-Evidence.md:27` | NUMBER | SUPPORTED | test | an external node (also known as an outer node, leaf node, or terminal node) is any node that does not have child nodes ... the root node has depth zero, leaf nodes have height zero · primary: none in doc besides Wikipedia/Biopython/DendroPy; CS convention — not fetched — — (height −1 for empty tree) | in-scope |
| O179 | `docs/Evidence/PHYLO-STATS-001-Evidence.md:63` | DEFINITION | SUPPORTED_DEF | doc-only | selects the branching pattern with the smallest total branch length ... later independently developed further by Andrey Rzhetsky and Masatoshi Nei (1992) | in-scope |
| O198 | `docs/Evidence/POP-LD-001.md:15` | DEFINITION | SUPPORTED_DEF | doc-only | D = P_{AB} - P_{A}P_{B} ... D' = \frac D {D_\max} ... correlation coefficient between pairs of loci, usually expressed as its square, r^2 | in-scope |
| O199 | `docs/Evidence/POP-LD-001.md:24` | DEFINITION | SUPPORTED_DEF | doc-only | one based on whether a given genomic sequence displays higher linkage disequilibrium than a predetermined threshold ... In 2001, Patil et al. ... percent of haplotypes are represen | in-scope |
| O200 | `docs/Evidence/PROBE-VALID-001-Evidence.md:13` | NUMBER | NO_NUMBER | code | negative scoring matrix cells are set to zero. Traceback procedure starts at the highest scoring matrix cell and proceeds until a cell with score zero is encountered · primary: Smith & Waterman 1981 J Mol Biol 147(1):195-197 (cited) — recurrence with 0 term (n/a) | in-scope |
| O201 | `docs/Evidence/PROBE-VALID-001-Evidence.md:25` | DEFINITION | NOT_SUPPORTED_DEF | doc-only | BLAST2 produces a single alignment with gaps that can include all of the initially found HSP regions. | in-scope |
| O204 | `docs/Evidence/PROTMOTIF-CC-001-Evidence.md:29` | NUMBER | NO_NUMBER | code | the positions in the heptad repeat are usually labeled abcdefg, where a and d are the hydrophobic positions, often being occupied by isoleucine, leucine, or valine. · primary: Mason JM, Arndt KM (2004) ChemBioChem 5(2):170-176 (cited by page and doc) — review; not fetched ({I,L,V} core set (no numeric value)) | in-scope |
| O205 | `docs/Evidence/PROTMOTIF-CC-001-Evidence.md:46` | DEFINITION | SUPPORTED_DEF | doc-only | a repeating pattern of seven amino acids ... H represents hydrophobic residues, C represents, typically, charged residues, and P represents polar ... predominantly leucine in the d | in-scope |
| O208 | `docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md:13` | NUMBER | NOT_STATED | doc-only | Domains vary in length from about 50 amino acids up to 250 amino acids in length. The shortest domains, such as zinc fingers, are stabilized by metal ions or disulfide bridges. | in-scope |
| O209 | `docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md:27` | NUMBER | NOT_STATED | doc-only | is a short peptide (usually 16–30 amino acids long) ... The core of the signal peptide contains a long stretch of hydrophobic amino acids (about 5–16 residues long) | in-scope |
| O210 | `docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md:39` | NUMBER | SUPPORTED | code | X2-Cys-X2,4-Cys-X12-His-X3,4,5-His ... Pfam = PF00096 ... PS00028 · primary: PROSITE PS00028 (prosite.expasy.org/PS00028, cited in same doc) — pattern (C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H) | in-scope |
| O211 | `docs/Evidence/PROTMOTIF-FIND-001-Evidence.md:35` | NUMBER | SUPPORTED | code | This pattern may be written as N{P}[ST]{P} ... The signature of the C2H2-type zinc finger domain is: C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H · primary: PROSITE PS00028 ZINC_FINGER_C2H2_1 / PS00001 ASN_GLYCOSYLATION (same doc lines 20, 30) — pattern entries (C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H; N-{P}-[ST]-{P}) | in-scope |
| O212 | `docs/Evidence/PROTMOTIF-FIND-001-Evidence.md:48` | NUMBER | SUPPORTED | doc-only | PROSITE was created in 1988 by Amos Bairoch ... 'The 20 years of PROSITE' ... 2007 ... volume=36\|pages=D245 ... De Castro E ... 2006 \|volume=34 \|issue=Web Server issue \|pages=W | in-scope |
| O215 | `docs/Evidence/QUALITY-STATS-001-Evidence.md:13` | NUMBER | SUPPORTED | code | Q = -10 log_{10} P ... \| 10 \| 1 in 10 \| 90% \| 20 \| 1 in 100 \| 99% \| 30 \| 1 in 1000 \| 99.9% \| 40 \| 1 in 10,000 \| 99.99% \| 50 \| 1 in 100,000 \| 99.999% \| 60 \| 1 in 1, · primary: Ewing B, Green P (1998) Genome Res 8:186-194 (listed in same doc) — defines Q = -10 log10 P (not fetched) (Q=-10log10P) | in-scope |
| O218 | `docs/Evidence/REP-STR-001-Evidence.md:15` | NUMBER | SUPPORTED | test | generally ten nucleotides or less ... and are typically repeated 5–50 times. For example, the sequence TATATATATA is a dinucleotide microsatellite, and GTCGTCGTCGTCGTC is a trinucl · primary: Worked examples are hand-verifiable test fixtures; GGAT history: Weller P, Jeffreys AJ et al. (1984) EMBO J 3:439-446 (cited by page) — not fetched (TATATATATA, GTCGTCGTCG | in-scope |
| O220 | `docs/Evidence/RESTR-FILTER-001-Evidence.md:13` | NUMBER | NO_NUMBER | test | In a blunt-ended molecule, both strands terminate in a base pair. ... blunt ends are always compatible with each other. ... An overhang is a stretch of unpaired nucleotides ... Lon | OUT-OF-SCOPE:B08 |
| O221 | `docs/Evidence/RESTR-FILTER-001-Evidence.md:26` | NUMBER | SUPPORTED | test | recognition sites that are usually undivided and palindromic and 4–8 nucleotides in length ... It can either cleave at the center of both strands to yield a blunt end, or at a stag · primary: not looked up (out of scope) —  () | OUT-OF-SCOPE:B08 |
| O222 | `docs/Evidence/RESTR-FILTER-001-Evidence.md:41` | NUMBER | NOT_STATED | test | They recognize a specific DNA sequence, usually short (3 to 8 bp), and cut it, producing either blunt or overhung ends · primary: REBASE / per-enzyme recognition sites (EcoRI GAATTC 6 bp, TaqI TCGA 4 bp, NotI GCGGCCGC 8 bp) - also on Wikipedia 'Restriction enzyme#Examples' table — Examples table (4/6 | OUT-OF-SCOPE:B08 |
| O226 | `docs/Evidence/RNA-INVERT-001-Evidence.md:28` | NUMBER | NO_NUMBER | code | For example, 5'---TTACG nnnnnn CGTAA---3' is an inverted repeat sequence. When the intervening length is zero, the composite sequence is a [[palindromic sequence]]. | in-scope |
| O228 | `docs/Evidence/RNA-PAIR-001-Evidence.md:13` | NUMBER | SUPPORTED | doc-only | Crick FH \| title = Codon--anticodon pairing: the wobble hypothesis \| journal = Journal of Molecular Biology \| volume = 19 \| issue = 2 \| pages = 548–555 \| date = August 1966 \ | in-scope |
| O230 | `docs/Evidence/RNA-PAIR-001-Evidence.md:28` | NUMBER | SUPPORTED | doc-only | Top, a '''G•C''' base pair with three [[hydrogen bond]]s. Bottom, an '''A•T''' base pair with two hydrogen bonds. ... Watson–Crick base pairs (guanine–cytosine and adenine-uracil) | in-scope |
| O231 | `docs/Evidence/RNA-PAIR-001-Evidence.md:32` | DEFINITION | SUPPORTED_DEF | doc-only | in RNA molecules ... Watson–Crick base pairs (guanine–cytosine and adenine-uracil) ... a G•C base pair with three hydrogen bonds. Bottom, an A•T base pair with two hydrogen | in-scope |
| O232 | `docs/Evidence/RNA-PAIR-001-Evidence.md:42` | DEFINITION | SUPPORTED_DEF | doc-only | Table: Weak W A,T; Strong S C,G; Amino M A,C; Keto K G,T; Purine R A,G; Pyrimidine Y C,T; Not A B; Not C D; Not G H; Not T V; Any N. Cites NC-IUB NAR 13(9):3021-3030 and IUPAC-IUB  | in-scope |
| O233 | `docs/Evidence/RNA-PAIR-001-Evidence.md:46` | DEFINITION | SUPPORTED_DEF | doc-only | Table: Weak W A,T; Strong S C,G; Amino M A,C; Keto K G,T; Purine R A,G; Pyrimidine Y C,T; Not A B; Not C D; Not G H; Not T V; Any N. Cites NC-IUB NAR 13(9):3021-3030 and IUPAC-IUB  | in-scope |
| O237 | `docs/Evidence/RNA-PKPREDICT-001-Evidence.md:44` | NUMBER | NO_NUMBER | doc-only | containing at least two stem-loop structures in which half of one stem is intercalated between the two halves of another stem ... The simple telomerase P2b-P3 example in the articl | in-scope |
| O239 | `docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md:50` | DEFINITION | SUPPORTED_DEF | doc-only | containing at least two stem-loop structures in which half of one stem is intercalated between the two halves of another stem ... base pairing in pseudoknots is not well nested | in-scope |
| O244 | `docs/Evidence/SEQ-ATSKEW-001-Evidence.md:35` | NUMBER | SUPPORTED | code | GC skew = (G - C)/(G + C) ... AT skew = (A − T)/(A + T) ... skew spectra ranges from −1, which corresponds to G = 0 or A = 0, to +1, which corresponds to T= 0 or C = 0 · primary: Charneski et al. 2011 PLoS Genet / Lobry 1996 Mol Biol Evol 13:660 (same doc); Biopython GC_skew analogue — AT skew (A−T)/(A+T) in [−1,1] ((A−T)/(A+T), ±1) | in-scope |
| O246 | `docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md:39` | DEFINITION | SUPPORTED_DEF | doc-only | '''Codon usage bias''' refers to differences in the frequency of occurrence of synonymous codons in coding DNA. | in-scope |
| O248 | `docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md:70` | NUMBER | NO_NUMBER | doc-only | while u + v <= n do / if S[i + v] = S[u + v] then / v := v + 1 / else ... if i = u then // all the pointers have been treated / C := C + 1 | in-scope |
| O250 | `docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md:71` | NUMBER | SUPPORTED | code | The maximal entropy of an event with ''n'' different outcomes is log_b(n): it is attained by the uniform probability distribution. ... The minimum surprise is when p = 0 ... or p = · primary: Shannon CE (1948) A Mathematical Theory of Communication, BSTJ 27:379-423 — Part I, Sec. 6, property 2: 'For a given n, H is a maximum and equal to log n when all the p_i  | in-scope |
| O252 | `docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md:13` | NUMBER | SUPPORTED | code | the number of nodes in the most balanced tree ... at the tree level i is either 4^i or N-i+1, whichever is smaller ... in the range 0<C<1 · primary: Troyanskaya et al. 2002 / Trifonov 1990 (cited in doc header; not re-fetched) — - (min(4^i, N−i+1)) | in-scope |
| O253 | `docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md:39` | NUMBER | SUPPORTED | code | The maximal entropy of an event with ''n'' different outcomes is log_b(n) ... the summand 0 log_b(0) is taken to be 0 · primary: Shannon CE (1948) A Mathematical Theory of Communication, BSTJ 27:379-423 — Part I, Sec. 6, property 2: 'For a given n, H is a maximum and equal to log n when all the p_i  | in-scope |
| O258 | `docs/Evidence/SEQ-COMPOSITION-001-Evidence.md:27` | NUMBER | SUPPORTED | code | GC skew = (G - C)/(G + C) ... AT skew = (A − T)/(A + T) ... skew spectra ranges from −1, which corresponds to G = 0 or A = 0, to +1, which corresponds to T= 0 or C = 0 ... Molecula · primary: Biopython Bio.SeqUtils.GC_skew source (listed in same Evidence doc) — GC_skew: skew = (g - c) / (g + c) ((G−C)/(G+C), range [−1,1]) | in-scope |
| O260 | `docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md:13` | NUMBER | SUPPORTED | doc-only | \|journal=[[Bell System Technical Journal]] \|volume=27 \|issue=3 \|pages=379–423 \|date=July 1948 \|doi=10.1002/j.1538-7305.1948.tb01338.x | in-scope |
| O261 | `docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md:26` | NUMBER | SUPPORTED | code | Base 2 gives the unit of [[bit]]s (or "shannons") ... The maximal entropy of an event with ''n'' different outcomes is log_b(n): it is attained by the uniform probability distribut · primary: Shannon CE (1948) A Mathematical Theory of Communication, BSTJ 27:379-423 — Part I, Sec. 6, property 2: 'For a given n, H is a maximum and equal to log n when all the p_i  | in-scope |
| O264 | `docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md:13` | NUMBER | SUPPORTED | code | GC-content percentage is calculated as<ref>{{cite book \| author =Madigan, MT. and Martinko JM. ... :<math>\cfrac{G+C}{A+T+G+C} \times 100%</math> · primary: Biopython Bio.SeqUtils.gc_fraction (Cock et al. 2009, Bioinformatics 25:1422) / Madigan & Martinko 2003 (cited by Wikipedia) — gc_fraction: (G+C)/len over A,T,G,C; x100 fo | in-scope |
| O265 | `docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md:24` | NUMBER | SUPPORTED | code | GC skew = (G - C)/(G + C) ... AT skew = (A − T)/(A + T) ... skew spectra ranges from −1, which corresponds to G = 0 or A = 0, to +1, which corresponds to T= 0 or C = 0 · primary: Biopython Bio.SeqUtils.GC_skew source (listed in same Evidence doc) — GC_skew: skew = (g - c) / (g + c) ((G−C)/(G+C), range [−1,1]) | in-scope |
| O269 | `docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md:13` | NUMBER | SUPPORTED | code | <math>\cfrac{G+C}{A+T+G+C} \times 100%</math> ... whereas the AT/GC ratio is calculated as <math>\cfrac{A+T}{G+C}</math> · primary: Biopython Bio.SeqUtils.gc_fraction (Cock et al. 2009, Bioinformatics 25:1422) / Madigan & Martinko 2003 (cited by Wikipedia) — gc_fraction: (G+C)/len over A,T,G,C; x100 fo | in-scope |
| O271 | `docs/Evidence/SEQ-REPLICATION-001-Evidence.md:49` | NUMBER | SUPPORTED | code | With the flipped definition of GC skew, the maximum value of the cumulative skew corresponds to the terminal, and the minimum value corresponds to the origin of replication. · primary: Grigoriev 1998 NAR 26:2286-2290 (same doc, line 30ff) — 'a sum of (G−C)/(G+C) in adjacent windows'; global minimum at origin, maximum near terminus ((G−C)/(G+C); min=origi | in-scope |
| O273 | `docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md:13` | NUMBER | SUPPORTED | code | four out of any six contiguous amino acids were sufficient to nucleate helix, and three out of any contiguous five were sufficient for a sheet ... originally 1.03 was set as the he · primary: Chen H, Gu F, Huang Z (2006) BMC Bioinformatics 7(Suppl 4):S14, Methods rules 1-3 (cited in same Evidence doc line ~98 and SequenceStatistics.cs:2113) — Methods: 4 of 6 he | in-scope |
| O275 | `docs/Evidence/SEQ-STATS-001-Evidence.md:26` | DEFINITION | SUPPORTED_DEF | doc-only | GC skew = (G − C)/(G + C) ... AT skew = (A − T)/(A + T) ... positive GC skew represents richness of G over C and the negative GC skew represents richness of C over G | in-scope |
| O277 | `docs/Evidence/SEQ-SUMMARY-001-Evidence.md:42` | NUMBER | SUPPORTED | code | Base 2 gives the unit of [[bit]]s (or "shannons") ... The maximal entropy of an event with ''n'' different outcomes is log_b(n) · primary: Shannon CE (1948) A Mathematical Theory of Communication, BSTJ 27:379-423 — Part I, Sec. 6, property 2: 'For a given n, H is a maximum and equal to log n when all the p_i  | in-scope |
| O278 | `docs/Evidence/SEQ-SUMMARY-001-Evidence.md:55` | NUMBER | SUPPORTED | code | Complexity (C) of a sequence fragment ... can be directly calculated as the product of vocabulary-usage measures (U_i) ... in the range 0<C<1 · primary: Gabrielian & Bolshoy 1999 (the page attributes the product formula to Gabrielian1999, not Trifonov 1990) — not re-fetched (C = U1·U2·…·Uw, 0<C<1) | in-scope |
| O281 | `docs/Evidence/SEQ-THERMO-001-Evidence.md:59` | NUMBER | NOT_STATED | code | Terminal A/T base pair \| 9.6 \|\| 17.2 \|\| 4.31 ... Terminal G/C base pair \| 0.4 \|\| −11.7 \|\| 4.05 ... T_m = \frac{\Delta H^\circ}{\Delta S^\circ+R\ln([A]_{total} - [B]_{tota · primary: Biopython Bio/SeqUtils/MeltingTemp.py Tm_NN (SantaLucia 1998 PNAS 95:1460 eq. 3; PMC19045 full text not retrievable) — MeltingTemp.py ~l.1051: k = (dnac1 - dnac2/2) ; if s | in-scope |
| O283 | `docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md:13` | NUMBER | SUPPORTED | doc-only | canonical splicing or termed the lariat pathway, which accounts for more than 99% of splicing ... Y-rich-N-C-A-G-[cut]-G (acceptor site) | in-scope |
| O284 | `docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md:28` | NUMBER | SUPPORTED | doc-only | rich with pyrimidine nucleotides, especially uracil, and is usually 15–20 base pairs long, located about 5–40 base pairs before the 3' end of the intron to be spliced. | in-scope |
| O285 | `docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md:43` | NUMBER | NOT_STATED | doc-only | The 3' splice site can be further defined by a variable length of polypyrimidines, called the polypyrimidine tract (PPT), which serves the dual function of recruiting factors to th | in-scope |
| O286 | `docs/Evidence/SPLICE-DONOR-001-Evidence.md:13` | NUMBER | NOT_STATED | code | The consensus sequence for an intron ... is: G-G-[cut]-G-U-R-A-G-U (donor site) ... accounts for more than 99% of splicing · primary: Shapiro & Senapathy 1987 NAR 15:7155 (same doc line 44) and Yeo & Burge 2004 MaxEntScan (line 70) — 'nucleotide frequency distributions at positions -3 to +6 relative to t | in-scope |
| O287 | `docs/Evidence/SPLICE-DONOR-001-Evidence.md:27` | NUMBER | NO_NUMBER | code | Introns typically have a GU nucleotide sequence at the 5' end splice site, and an AG at the 3' end splice site. ... U1 snRNP binding to the 5' end splice site · primary: Canonical GT-AG rule (Breathnach & Chambon 1981; Burge et al. 1999 cited by page); not fetched — - (GT donor (no numeric value)) | in-scope |
| O290 | `docs/Evidence/SPLICE-PREDICT-001-Evidence.md:13` | NUMBER | SUPPORTED | doc-only | the shortest known [[metazoan]] intron length is 30 base pairs (bp) belonging to the human ''MST1L'' ... ≥3.6 [[megabase]] (Mb) intron ... the splicing reaction is likely to be 99. | in-scope |
| O291 | `docs/Evidence/SPLICE-PREDICT-001-Evidence.md:33` | NUMBER | SUPPORTED | doc-only | Across all eukaryotic genes in GenBank, there were (in 2002), on average, 5.48 exons per protein coding gene. The average exon encoded 30-36 amino acids. | in-scope |
| O292 | `docs/Evidence/SPLICE-PREDICT-001-Evidence.md:49` | DEFINITION | NOT_SUPPORTED_DEF | doc-only | Once spliced together, the exons form a single continuous protein-coding regions, and the splice boundaries are not detectable ... Introns are extremely rare in prokaryotes | in-scope |
| O296 | `docs/Evidence/TRANS-DIFF-001-Evidence.md:39` | NUMBER | NO_NUMBER | code | Unlike in Student's t-test, the denominator is not based on a pooled variance estimate. ... ν ≈ (s1²/N1 + s2²/N2)² / (s1⁴/(N1²ν1) + s2⁴/(N2²ν2)) ... s_i ... corrected sample standa · primary: Welch 1947 Biometrika 34:28-35 (cited; not fetched) — — (Welch–Satterthwaite df) | in-scope |
| O297 | `docs/Evidence/TRANS-DIFF-001-Evidence.md:52` | NUMBER | NO_NUMBER | code | For t > 0, F(t) = 1 - 1/2 I_{x(t)}(nu/2, 1/2), where x(t) = nu/(t^2+nu). · primary: Standard identity (Abramowitz & Stegun 1964, 26.7.1 / 26.5.27); not fetched — - (no numeric constant) | in-scope |
| O298 | `docs/Evidence/TRANS-DIFF-001-Evidence.md:64` | NUMBER | NO_NUMBER | code | For a given alpha, find the largest k for which P_(k) <= (k/m) alpha ... Journal of the Royal Statistical Society, Series B \| volume = 57 \| issue = 1 \| pages = 289–300 · primary: Benjamini Y, Hochberg Y (1995) JRSS-B 57:289-300; R stats p.adjust BH (cited in same doc) — step-up procedure / R pmin(1, cummin(n/i * p[o])) (no numeric constant) | in-scope |
| O302 | `docs/Evidence/TRANS-EXPR-001-Evidence.md:52` | NUMBER | SUPPORTED | test | values are tied in rank, they should instead be assigned the mean of the values corresponding to the ranks they would normally represent ... ((4.67 + 5.67)/2 = 5.17) · primary: Bolstad et al. 2003 Bioinformatics 19:185 (method; worked example is Wikipedia's own) — not re-fetched (rank-mean method; 5.17 tie value is Wikipedia's example) | in-scope |
| O304 | `docs/Evidence/TRANS-PROT-001-Evidence.md:12` | NUMBER | NO_NUMBER | doc-only | Termination of elongation depends on the release factor eRF1 that recognizes all three stop codons ... (UAA, UAG, or UGA) ... In eukaryotes and archaea, the amino acid encoded by t | OUT-OF-SCOPE:B02 |
| O305 | `docs/Evidence/TRANS-PROT-001-Evidence.md:28` | NUMBER | NO_NUMBER | code | Any given sequence of DNA can therefore be read in six different ways ... The mRNA is single-stranded and therefore only contains three possible reading frames, of which only one | OUT-OF-SCOPE:B02 |
| O306 | `docs/Evidence/TRANS-PROT-001-Evidence.md:42` | NUMBER | SUPPORTED | doc-only | in a randomly generated DNA sequence with an equal percentage of each nucleotide, a stop-codon would be expected once every 21 codons ... small ORFs usually < 100 codons in length | OUT-OF-SCOPE:B02 |
| O307 | `docs/Evidence/TRANS-SIXFRAME-001-Evidence.md:70` | NUMBER | SUPPORTED | code | Any given sequence of DNA can therefore be read in six different ways: Three reading frames in one direction (starting at different nucleotides) and three in the opposite direction · primary: Lodish 2007 6th ed. p.121 (cited on page; not re-fetched) — - (6 frames) | OUT-OF-SCOPE:B02 |
| O309 | `docs/Evidence/VARIANT-CALL-001-Evidence.md:30` | DEFINITION | SUPPORTED | doc-only | approximately two out of three single nucleotide polymorphisms (SNPs) are transitions (cites Collins & Jukes 1994) ... Although there are twice as many possible transversions, tran | in-scope |
| O310 | `docs/Evidence/VARIANT-CALL-001-Evidence.md:44` | NUMBER | SUPPORTED | code | a single (two ring) purine (A or G) is changed for a (one ring) pyrimidine (T or C), or vice versa. ... each of the 8 nucleotide changes between a purine and a pyrimidine ... The 4 · primary: Futuyma (2013) Evolution 3rd ed. (cited on page; not re-fetched) — - (8 directed = 4 undirected transversion pairs; 4 directed = 2 transition pairs) | in-scope |
| O313 | `docs/Evidence/VARIANT-SNP-001-Evidence.md:28` | DEFINITION | NOT_STATED | doc-only | Although there are two possible transversions but only one possible transition per base, transition mutations are more likely | in-scope |
| O568 | `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs:179` | NUMBER | NO_NUMBER | code | :<math>F_{0j} = d*j</math> :<math>F_{i0} = d*i</math> · primary: Durbin, Eddy, Krogh & Mitchison 1998 Biological Sequence Analysis §2.3 (not fetched) — eq. 2.8 boundary F(i,0) = −i·d (linear gap init) | in-scope |
| O569 | `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs:1512` | NUMBER | NO_NUMBER | code | a heuristic search known as progressive technique ... developed by Da-Fei Feng and Doolittle in 1987 ... clustering method such as neighbor-joining or unweighted pair group method  · primary: Feng & Doolittle 1987 J Mol Evol 25:351-360 (cited in same comment) — progressive method (n/a (method)) | in-scope |
| O570 | `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs:1927` | DEFINITION | SUPPORTED_DEF | doc-only | optimizes the sum of all of the pairs of characters at each position in the alignment (the so-called sum of pair score) ... iterative methods can return to previously calculated pa | in-scope |
| O571 | `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:419` | NUMBER | SUPPORTED | code | Since DNA is interpreted in groups of three nucleotides (codons), a DNA strand has three distinct reading frames · primary: Genetic-code triplet structure (definitional); not fetched — - (3 frames per strand) | in-scope |
| O572 | `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.Fusions.cs:31` | NUMBER | SUPPORTED | code | the sequence of nucleotides in a nucleic acid (DNA or RNA) molecule as a sequence of triplets · primary: Lodish, Molecular Cell Biology 6th ed. (2007) p.121, cited on page — not re-fetched (CodonLength = 3) | in-scope |
| O573 | `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.Fusions.cs:125` | NUMBER | NO_NUMBER | code | There are three reading frames that can be read in this 5′→3′ direction, each beginning from a different nucleotide in a triplet. | in-scope |
| O574 | `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.Signatures.cs:850` | NUMBER | NOT_STATED | code | :Until W and H are stable. · primary: SigProfilerExtractor (AlexandrovLab) sigpro.py defaults — sigpro.py lines 341-343: min_nmf_iterations=10000, max_nmf_iterations=1000000, nmf_test_conv=10000 (10000 (as MIN | in-scope |
| O575 | `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.Signatures.cs:875` | NUMBER | NOT_STATED | code | The NMF page gives the Lee–Seung update rules but says nothing about zero denominators, regularisation or an epsilon. · primary: none states 1e-12; closest reference practice: scikit-learn sklearn/decomposition/_nmf.py — line 32 EPSILON = np.finfo(np.float32).eps; floors WH[WH < EPSILON] = EPSILON ( | in-scope |
| O576 | `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:51` | NUMBER | SUPPORTED | code | so that 1 - alpha/2 = 0.975 and z_{.05} = 1.96 · primary: Standard normal quantile; Wilson E.B. (1927) JASA 22:209-212 — Phi^-1(0.975) (1.959964 ~ 1.96) | in-scope |
| O577 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/ApproximateMatcher_HammingDistance_Tests.cs:12` | DEFINITION | SUPPORTED_DEF | doc-only | the Hamming distance between two strings or vectors of equal length is the number of positions at which the corresponding symbols are different. | in-scope |
| O578 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_GlobalAlign_Tests.cs:8` | DEFINITION | SUPPORTED_DEF | doc-only | Gap penalty, here called d ... sum of the scores of all the pairings is the score of the whole alignment candidate | in-scope |
| O579 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_GlobalAlign_Tests.cs:9` | DEFINITION | SUPPORTED_DEF | doc-only | Global alignments, which attempt to align every residue in every sequence, are most useful when the sequences in the query set are similar and of roughly equal size. (This does not | in-scope |
| O580 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_GlobalAlign_Tests.cs:24` | NUMBER | SUPPORTED | test | * Match: +1 ... * Mismatch or Indel: −1 · primary: none needed — illustrative Wikipedia scoring scheme, not a published parameter set —  () | in-scope |
| O581 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_GlobalAlign_Tests.cs:43` | NUMBER | SUPPORTED | test | GCATG-CG ... For the Example above, the score of the alignment would be 0: · primary: none — Wikipedia's own worked example; score is closed-form arithmetic (4x(+1)+2x(−1)+2x(−1)=0) —  () | in-scope |
| O582 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_GlobalAlign_Tests.cs:68` | NUMBER | NO_NUMBER | test | Each of these scenarios is assigned a score and the sum of the scores of all the pairings is the score of the whole alignment candidate. | in-scope |
| O583 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_GlobalAlign_Tests.cs:69` | DEFINITION | SUPPORTED_DEF | doc-only | Global alignments, which attempt to align every residue in every sequence, are most useful when the sequences in the query set are similar and of roughly equal size. (This does not | in-scope |
| O584 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_LocalAlign_Tests.cs:5` | DEFINITION | SUPPORTED | test | Take the alignment of DNA sequences TGTTACGG and GGTTGACTA as an example ... s = +3 if a_i=b_j, -3 otherwise ... a linear gap penalty of W_1 = 2 · primary: Smith & Waterman (1981) J Mol Biol 147:195 (algorithm); example parameters are Wikipedia's own worked example — Wikipedia §Example is the origin of the dataset (+3/-3/W1=2 | in-scope |
| O585 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_LocalAlign_Tests.cs:20` | NUMBER | SUPPORTED | test | Substitution matrix: s(a_i,b_j) = +3 if a_i=b_j, −3 otherwise ... Gap penalty: W_k=2k (a linear gap penalty of W_1 = 2) · primary: none needed — illustrative Wikipedia scheme —  () | in-scope |
| O586 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_LocalAlign_Tests.cs:32` | NUMBER | SUPPORTED | test | Take the alignment of DNA sequences TGTTACGG and GGTTGACTA ... The alignment result is: G T T - A C / G T T G A C · primary: none — Wikipedia's own example; score is arithmetic on the shown alignment — 5 matches x 3 − 1 gap x 2 (13) | in-scope |
| O587 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlignIterative_Tests.cs:35` | NUMBER | NO_NUMBER | test | The MSA program optimizes the sum of all of the pairs of characters at each position in the alignment (the so-called sum of pair score) · primary: Carrillo & Lipman 1988 / Barton & Sternberg 1987 (cited) — SP objective (n/a (definition)) | in-scope |
| O588 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlignProgressive_Tests.cs:19` | NUMBER | NO_NUMBER | test | guide tree ... clustering method such as neighbor-joining or ... UPGMA · primary: Feng & Doolittle 1987 J Mol Evol 25:351-360 — progressive method (n/a) | in-scope |
| O589 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlign_Tests.cs:15` | NUMBER | NO_NUMBER | test | the sum of all of the pairs of characters at each position in the alignment (the so-called sum of pair score) | in-scope |
| O590 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlign_Tests.cs:16` | DEFINITION | SUPPORTED_DEF | doc-only | A guide tree is calculated from the scores of the sequences in the matrix, then subsequently used to build the multiple se[quence alignment] | in-scope |
| O591 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlign_Tests.cs:17` | DEFINITION | SUPPORTED_DEF | doc-only | the calculated sequence of most frequent residues, either nucleotide or amino acid, found at each position in a sequence alignment ... Y represents any pyrimidine, and R indicates  | in-scope |
| O592 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Tests.cs:287` | NUMBER | SUPPORTED | test | e.g. ATTCG ATTCG ATTCG, in which the sequence ATTCG is repeated three times. · primary: MeSH 'Tandem Repeat Sequences' (cited on page); example is Wikipedia's — - (unit ATTCG, 3 copies, length 15) | in-scope |
| O593 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_AnalyzeKmers_Tests.cs:4` | DEFINITION | SUPPORTED | test | the sequence AGAT would have ... three 2-mers (AG, GA, AT) ... a sequence of length L will have L - k + 1 k-mers ... n^{k} possible k-mers ... \|2 \|GT, TA, AG, GA, AG, GC, CT, TG, · primary: Definition-level / arithmetic: counts derive directly from Wikipedia's GTAGAGCTGT k-mer table (k=2: 9 total, 7 distinct; k=3: 8/8; k=4: 7/7); n^k universe size is combinat | in-scope |
| O594 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountKmersAsync_Tests.cs:4` | DEFINITION | SUPPORTED | test | the sequence AGAT would have ... three 2-mers (AG, GA, AT) ... a sequence of length L will have L - k + 1 k-mers ... n^{k} possible k-mers ... \|2 \|GT, TA, AG, GA, AG, GC, CT, TG, · primary: Definition-level / arithmetic: counts derive directly from Wikipedia's GTAGAGCTGT k-mer table (k=2: 9 total, 7 distinct; k=3: 8/8; k=4: 7/7); n^k universe size is combinat | in-scope |
| O595 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountKmers_Tests.cs:387` | NUMBER | SUPPORTED | test | The sequence ATGG has two 3-mers: ATG and TGG. · primary: Combinatorial identity L-k+1 (4-3+1=2) — hand-derivable (2) | in-scope |
| O596 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountKmers_Tests.cs:406` | NUMBER | SUPPORTED | test | k-mers for GTAGAGCTGT ... \|2 \|GT, TA, AG, GA, AG, GC, CT, TG, GT · primary: Wikipedia worked example (no primary; L−k+1 is combinatorial) — k-mer table row k=2 (9 listed, 7 distinct) | in-scope |
| O597 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountKmers_Tests.cs:430` | NUMBER | SUPPORTED | test | More generally, a sequence of length L will have L − k + 1 k-mers · primary: combinatorial identity (no primary needed) — - (L−k+1) | in-scope |
| O598 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountKmers_Tests.cs:452` | NUMBER | SUPPORTED | test | \|4 \|GTAG, TAGA, AGAG, GAGC, AGCT, GCTG, CTGT · primary: Wikipedia worked example — k-mer table row k=4 (7 distinct 4-mers, each once) | in-scope |
| O599 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_FindUniqueAndMinCount_Tests.cs:6` | NUMBER | NO_NUMBER | test | In bioinformatics, k-mers are substrings of length k contained within a biological sequence. | in-scope |
| O600 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_GenerateAllKmers_Tests.cs:4` | DEFINITION | SUPPORTED | test | the sequence AGAT would have ... three 2-mers (AG, GA, AT) ... a sequence of length L will have L - k + 1 k-mers ... n^{k} possible k-mers ... \|2 \|GT, TA, AG, GA, AG, GC, CT, TG, · primary: Definition-level / arithmetic: counts derive directly from Wikipedia's GTAGAGCTGT k-mer table (k=2: 9 total, 7 distinct; k=3: 8/8; k=4: 7/7); n^k universe size is combinat | in-scope |
| O601 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_PWM_Tests.cs:8` | DEFINITION | SUPPORTED_DEF | doc-only | A position weight matrix (PWM) ... Pseudocounts (or Laplace estimators) are often applied when calculating PPMs ... (cites Nishida et al. 2008) | in-scope |
| O602 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_PWM_Tests.cs:292` | NUMBER | SUPPORTED | test | GAGGTAAAC<br />TCCGTAAGT<br />CAGGTTGGA<br />ACAGTCAGT<br />TAGGTCATT<br />TAGGTACTG<br />ATGGTAACT<br />CAGGTATAC<br />TGTGTGAGT<br />AAGGTAAGT · primary: none — Wikipedia's illustrative dataset; consensus and log-odds are arithmetic on it —  () | in-scope |
| O603 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Microsatellite_Tests.cs:7` | DEFINITION | SUPPORTED_DEF | doc-only | the sequence TATATATATA is a dinucleotide microsatellite, and GTCGTCGTCGTCGTC is a trinucleotide microsatellite ... the first microsatellite was characterised in 1984 ... by Weller | in-scope |
| O604 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Microsatellite_Tests.cs:67` | DEFINITION | SUPPORTED_DEF | doc-only | the sequence TATATATATA is a dinucleotide microsatellite, and GTCGTCGTCGTCGTC is a trinucleotide microsatellite ... the first microsatellite was characterised in 1984 ... by Weller | in-scope |
| O605 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Microsatellite_Tests.cs:92` | DEFINITION | SUPPORTED_DEF | doc-only | the sequence TATATATATA is a dinucleotide microsatellite, and GTCGTCGTCGTCGTC is a trinucleotide microsatellite ... the first microsatellite was characterised in 1984 ... by Weller | in-scope |
| O606 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Microsatellite_Tests.cs:168` | DEFINITION | SUPPORTED_DEF | doc-only | the sequence TATATATATA is a dinucleotide microsatellite, and GTCGTCGTCGTCGTC is a trinucleotide microsatellite ... the first microsatellite was characterised in 1984 ... by Weller | in-scope |
| O607 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Palindrome_Tests.cs:10` | DEFINITION | SUPPORTED_DEF | doc-only | identical to the sequence in the same direction (e.g. 5' to 3') on the complementary strand ... Many restriction endonucleases recognize specific palindromic sequences | in-scope |
| O608 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Palindrome_Tests.cs:11` | DEFINITION | SUPPORTED_DEF | doc-only | Many of them are palindromic, meaning the base sequence reads the same backwards and forwards ... EcoRI digestion produces "sticky" ends | in-scope |
| O609 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_CalculateGcContentProfile_Tests.cs:5` | NUMBER | NO_NUMBER | test | '''GC-content''' ... is the percentage of [[nitrogenous bases]] in a DNA or RNA molecule that are either guanine (G) or cytosine (C) · primary: Biopython Bio.SeqUtils.gc_fraction (Cock et al. 2009, Bioinformatics 25:1422) / Madigan & Martinko 2003 (cited by Wikipedia) — gc_fraction: (G+C)/len over A,T,G,C; x100 fo | in-scope |
| O610 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/GenomeAnnotator_Gene_Tests.cs:462` | NUMBER | SUPPORTED | test | The six-base consensus sequence is AGGAGG ... Shine J, Dalgarno L ... 1975 ... Nature \| volume = 254 \| issue = 5495\| pages = 34–38 · primary: Shine J, Dalgarno L (1975) Nature 254:34-38 (cited in test) — not fetched (AGGAGG) | in-scope |
| O611 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Synteny_Tests.cs:9` | NUMBER | NO_NUMBER | doc-only | These blocks are referred to as syntenic blocks. ... This is typically done using a version of the MCScan algorithm | in-scope |
| O612 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Synteny_Tests.cs:10` | NUMBER | NO_NUMBER | doc-only | In comparative genomics, synteny is the preserved order of genes on chromosomes of related species indicating their descent from a common ancestor | in-scope |
| O613 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Synteny_Tests.cs:11` | NUMBER | NO_NUMBER | doc-only | classes of events, like deletions, duplications, inversions, and translocations | in-scope |
| O614 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Telomere_Tests.cs:8` | DEFINITION | SUPPORTED_DEF | doc-only | the sequence repeats are enriched in guanine, e.g. TTAGGG in vertebrates (cites Meyne et al. 1989) | in-scope |
| O615 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/SequenceExtensions_CalculateGcContent_Tests.cs:7` | NUMBER | SUPPORTED | test | <math>\cfrac{G+C}{A+T+G+C} \times 100%</math> · primary: Biopython Bio.SeqUtils.gc_fraction (Cock et al. 2009, Bioinformatics 25:1422) / Madigan & Martinko 2003 (cited by Wikipedia) — gc_fraction: (G+C)/len over A,T,G,C; x100 fo | in-scope |
| O616 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Metagenomics/MetagenomicsAnalyzer_BetaDiversity_Tests.cs:21` | NUMBER | SUPPORTED | test | \| Goldfish\|\| 6 \|\| 10\|\| 6 ... \| Guppy\|\| 7 \|\| 0\|\| 0 ... \| Rainbow fish\|\| 4 \|\| 6 \|\| 4 ... C_{jk} = 6 + 0 + 4 = 10 ... BC_{jk} = 1 - (2 \times 10) / (17 + 16) = 0. · primary: Bray & Curtis (1957) Ecol Monogr 27(4):325-349 (cited on page) — defines the index BC = 1 - 2C/(S_j+S_k); the aquarium numbers are Wikipedia's own worked example (formula  | in-scope |
| O617 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_GuideRNA_Tests.cs:8` | DEFINITION | SUPPORTED_DEF | doc-only | The PAM sequence (NGG) serves as a binding signal for Cas9 ... 20 bp spacer/targeting sequence that is used to guide Cas9 to target DNA | OUT-OF-SCOPE:B08 |
| O618 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_GuideRNA_Tests.cs:9` | DEFINITION | SUPPORTED_DEF | doc-only | The targeting specificity of CRISPR-Cas9 is determined by the 20-nucleotide (nt) sequence at the 5' end of the gRNA. The desired target sequence must precede the Protospacer Adjace | OUT-OF-SCOPE:B08 |
| O619 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_GuideRNA_Tests.cs:10` | DEFINITION | SUPPORTED_DEF | doc-only | The canonical Cas9 PAM is the sequence 5'-NGG-3', where "N" is any nucleobase followed by two guanine nucleobases | OUT-OF-SCOPE:B08 |
| O620 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Digest_Tests.cs:5` | DEFINITION | SUPPORTED_DEF | doc-only | In a restriction digest, DNA molecules are cleaved at specific regions of 4-12 nucleotides in length (restriction sites) | OUT-OF-SCOPE:B08 |
| O621 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Digest_Tests.cs:6` | DEFINITION | SUPPORTED_DEF | doc-only | Many of them are palindromic, meaning the base sequence reads the same backwards and forwards ... EcoRI digestion produces "sticky" ends | OUT-OF-SCOPE:B08 |
| O622 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Digest_Tests.cs:7` | DEFINITION | SUPPORTED_DEF | doc-only | A restriction map is a map of known restriction sites within a sequence of DNA | OUT-OF-SCOPE:B08 |
| O623 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Digest_Tests.cs:8` | DEFINITION | SUPPORTED_DEF | doc-only | Restriction Digest of Plasmid DNA ... Incubate tube at appropriate temperature (usually 37 °C) for 1 hour | OUT-OF-SCOPE:B08 |
| O624 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Digest_Tests.cs:289` | NUMBER | NO_NUMBER | test | Enzymes that only cut once allow you to more easily and accurately visualize the full size of your construct. · primary: Univ. Illinois MolBio 'Restriction Digestion/Gel Electrophoresis Assignment 1' (cited in same comment block) — quoted: 'If you cut a circle once, you get one linear fragme | OUT-OF-SCOPE:B08 |
| O625 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:5` | DEFINITION | SUPPORTED_DEF | doc-only | Many of them are palindromic, meaning the base sequence reads the same backwards and forwards ... EcoRI digestion produces "sticky" ends | OUT-OF-SCOPE:B08 |
| O626 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:6` | DEFINITION | SUPPORTED_DEF | doc-only | EcoRI recognizes the palindromic sequence GAATTC and cuts between the G and the A on both the top and bottom strands | OUT-OF-SCOPE:B08 |
| O627 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:7` | DEFINITION | SUPPORTED_DEF | doc-only | EcoRI creates sticky ends with 5' end overhangs of AATT. The nucleic acid recognition sequence where the enzyme cuts is G↓AATTC | OUT-OF-SCOPE:B08 |
| O628 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:708` | NUMBER | SUPPORTED | test | EcoRI creates sticky ends with 5' end overhangs of AATT. The nucleic acid recognition sequence where the enzyme cuts is G↓AATTC · primary: REBASE EcoRI entry (standard reference); not fetched — - (G^AATTC) | OUT-OF-SCOPE:B08 |
| O629 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:738` | NUMBER | SUPPORTED | test | BamHI binds at the recognition sequence 5'-GGATCC-3', and cleaves these sequences just after the 5'-guanine on each strand. This cleavage results in sticky ends which are 4 bp long · primary: REBASE BamHI entry (G^GATCC) — not fetched — enzyme record (cut after G, 4-nt 5' overhang GATC) | OUT-OF-SCOPE:B08 |
| O630 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:768` | NUMBER | SUPPORTED | test | cleaves the DNA palindromic sequence AAGCTT ... The cleavage of this sequence between the AA's results in 5' overhangs on the DNA called sticky ends · primary: REBASE HindIII entry; not fetched — - (A^AGCTT) | OUT-OF-SCOPE:B08 |
| O631 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:799` | NUMBER | SUPPORTED | test | PstI cleaves DNA at the recognition sequence 5′-CTGCA/G-3′ generating fragments with 3′-cohesive termini. ... This cleavage yields [[sticky ends]] 4 base pairs long. · primary: Thermo Fisher PstI (ER0611) product page, cited by the page (not fetched; out of scope) —  (CTGCA^G, 4-nt 3' overhang) | OUT-OF-SCOPE:B08 |
| O632 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:824` | DEFINITION | SUPPORTED_DEF | test | It creates blunt ends. The enzyme recognizes the palindromic 6-base DNA sequence 5'-GAT\|ATC-3' and makes a blunt end at the vertical line. · primary: Schildkraut I et al. (1984) Gene 27:327-329 'The cleavage site for the restriction endonuclease EcoRV is 5'-GAT/ATC' — title (cut after position 3 (GAT\|ATC)) | OUT-OF-SCOPE:B08 |
| O633 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:845` | DEFINITION | SUPPORTED_DEF | doc-only | Many of them are palindromic, meaning the base sequence reads the same backwards and forwards. | OUT-OF-SCOPE:B08 |
| O634 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_FindSites_Tests.cs:871` | NUMBER | NOT_STATED | test | EcoRI 5'GAATTC ... 5'---G AATTC---3' 3'---CTTAA G---5' (table also has BamHI, HindIII, TaqI, NotI, Sau3AI, PvuII, SmaI, AluI, HaeIII, EcoRV, KpnI, PstI, SacI, SalI, ScaI, SpeI, Sph · primary: REBASE (rebase.neb.com) enzyme entries; not fetched — BglII A^GATCT, ApaI GGGCC^C (1/5 and 5/1) | OUT-OF-SCOPE:B08 |
| O635 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ClassifySbsContext_Tests.cs:7` | NUMBER | NO_NUMBER | doc-only | Match up between two DNA bases (adenine and thymine ... Complementarity is achieved by distinct interactions between nucleobases: adenine, thymine (uracil in RNA), guanine and cyto | in-scope |
| O636 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_FitSignatures_Tests.cs:6` | NUMBER | NO_NUMBER | test | active-set method published by Lawson and Hanson in their 1974 book Solving Least Squares Problems | in-scope |
| O637 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:365` | NUMBER | SUPPORTED | test | distance matrix computed from the 5S ribosomal RNA sequence alignment of five bacteria: Bacillus subtilis (a), Bacillus stearothermophilus (b), ... \| 17 \|\| 21 \|\| 31 \|\| 23 · primary: Sokal & Michener 1958 (method); matrix from Erdmann/Olsen 1988 5S rRNA data cited by page — not fetched — — (matrix {17,21,31,23,30,34,21,28,39,43}) | in-scope |
| O638 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:381` | NUMBER | SUPPORTED | test | five taxa (a,b,c,d,e) and the following distance matrix D: a 0 5 9 9 8 / b 5 0 10 10 9 / c 9 10 0 8 7 / d 9 10 8 0 3 / e 8 9 7 3 0 ... A distance matrix ... is said to be 'additive · primary: Test fixture (worked example); NJ algorithm Saitou N, Nei M (1987) Mol Biol Evol 4:406-425 — expected tree hand-derivable from the additive matrix (5x5 matrix) | in-scope |
| O639 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:399` | NUMBER | SUPPORTED | test | D_2((a,b),e)=(D_1(a,e) + D_1(b,e))/2=(23+21)/2=22 ... smallest D1 = 17 (a,b) · primary: Sokal & Michener 1958 (method) — worked values are Wikipedia's example — — (17, 22) | in-scope |
| O640 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:426` | NUMBER | SUPPORTED | test | δ(a,u)=δ(b,u)=17/2=8.5 ... δ(a,v)=δ(b,v)=δ(e,v)=22/2=11 ... δ(u,v)=δ(e,v)−δ(a,u)=11−8.5=2.5 · primary: Sokal & Michener 1958 (method) — worked values are Wikipedia's example — — (8.5, 2.5, 11) | in-scope |
| O641 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:486` | NUMBER | SUPPORTED | test | delta(((a,b),e),r)=delta((c,d),r)=33/2=16.5 ... It is ultrametric because all tips (a to e) are equidistant from r · primary: Test fixture (worked example; matrix a-e 17/21/31/23... identical to page); UPGMA: Sokal & Michener (1958) — hand-derivable (16.5) | in-scope |
| O642 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:511` | DEFINITION | SUPPORTED_DEF | doc-only | if we move from any taxon to any other along the branches of the tree, and sum the lengths of the branches traversed ... Given an additive distance matrix as input, neighbor joinin | in-scope |
| O643 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:543` | NUMBER | SUPPORTED | test | In the example above, Q_1(a,b)=-50. This is the smallest value of Q_1, so we join elements a and b. ... both pairs have the minimal Q_2 value of -28 · primary: Saitou & Nei 1987 Mol Biol Evol 4:406 (Q-criterion; Studier & Keppler 1988 form) — not fetched — method (worked values −50/−28 are Wikipedia's example, derivable from its  | in-scope |
| O644 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Phylogenetics/PhylogeneticAnalyzer_TreeConstruction_Tests.cs:570` | NUMBER | SUPPORTED | test | \delta(a,u)=... =\frac{5}{2} + \frac{31-34}{6} = 2 ... \delta(b,u)=d(a,b)-\delta(a,u) \quad = 5-2 = 3 · primary: Saitou & Nei (1987) branch-length formula; the 5-taxon example is Wikipedia's own — - (δ(a,u)=2, δ(b,u)=3 (computed from the formula)) | in-scope |
| O645 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_Diversity_Tests.cs:113` | NUMBER | SUPPORTED | test | Person Y  00000 00000 00000 00000 / Person A  00100 00000 00100 00010 / Person B  00000 00000 00100 00010 ... The average number of polymorphisms is ... = 2 · primary: none needed — data are Wikipedia's example; π = 2/20 = 0.1 follows from the page's average of 2 —  (0.1) | in-scope |
| O646 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_Diversity_Tests.cs:414` | NUMBER | NOT_STATED | test | n=5 individuals and S=4 segregating sites / a1=1/1+1/2+1/3+1/4=2.08 / M=4/2.08=1.92 / d = 2 - 1.92= .08 · primary: Tajima F (1989) Genetics 123:585-595, doi:10.1093/genetics/123.3.585 (listed in docs/Evidence/POP-DIV-001-Evidence.md row 7) — eq. 38 with a1,a2,b1,b2,c1,c2,e1,e2; recompu | in-scope |
| O647 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_Diversity_Tests.cs:445` | DEFINITION | NOT_STATED | test | Since there were n=5 individuals and S=4 segregating sites; a1=2.08; M=4/2.08=1.92; d = 2 - 1.92 = .08 · primary: Tajima F (1989) Genetics 123:585-595, eq. 38 (D = d / sqrt(e1*S + e2*S*(S-1))) — Applied to the Wikipedia dataset (n=5, S=4, d=0.08): a1=2.0833, a2=1.4236, e1=0.00960, e2= | in-scope |
| O648 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_FStatistics_Tests.cs:11` | NUMBER | NO_NUMBER | doc-only | Developed as a special case of Wright's F-statistics, it is one of the most commonly used statistics in population genetics. | in-scope |
| O649 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_FStatistics_Tests.cs:12` | NUMBER | NO_NUMBER | doc-only | Wright 1965 'The Interpretation of Population Structure by F-Statistics with Special Regard to Systems of Mating' Evolution 19(3):395–420 | in-scope |
| O650 | `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs:9` | DEFINITION | SUPPORTED | test | This data is from E. B. Ford (1971) on the scarlet tiger [moth] ... = (1469 - 1467.4)^2/1467.4 + (138 - 141.2)^2/141.2 + (5 - 3.4)^2/3.4 · primary: Ford EB (1971) Ecological Genetics (Ford dataset 1469/138/5) — cited in the test header alongside Wikipedia — scarlet tiger moth genotype counts (1469 / 138 / 5) | in-scope |
| O651 | `tests/Seqeron/Seqeron.Mcp.Population.Tests/DiversityStatisticsTests.cs:23` | NUMBER | NOT_STATED | test | n=5 individuals and S=4 segregating sites / a1=1/1+1/2+1/3+1/4=2.08 / M=4/2.08=1.92 / d = 2 - 1.92= .08 · primary: Tajima F (1989) Genetics 123:585-595, doi:10.1093/genetics/123.3.585 (listed in docs/Evidence/POP-DIV-001-Evidence.md row 7) — eq. 38 with a1,a2,b1,b2,c1,c2,e1,e2; recompu | in-scope |
| O755 | `tests/TestSpecs/SEQ-COMP-001.md:59` | DEFINITION | SUPPORTED_DEF | doc-only | Table: Weak W A,T; Strong S C,G; Amino M A,C; Keto K G,T; Purine R A,G; Pyrimidine Y C,T; Not A B; Not C D; Not G H; Not T V; Any N. Cites NC-IUB NAR 13(9):3021-3030 and IUPAC-IUB  | in-scope |

## BIBLIOGRAPHY occurrences (link check)

| id | file:line | URL | result |
|---|---|---|---|
| O006 | `docs/Evidence/ALIGN-MULTI-001-Evidence.md:140` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O008 | `docs/Evidence/ALIGN-MULTI-001-Evidence.md:265` | https://en.wikipedia.org/wiki/Multiple_sequence_alignment | LINK-OK |
| O009 | `docs/Evidence/ALIGN-SEMI-001-Evidence.md:14` | https://en.wikipedia.org/wiki/Sequence_alignment | LINK-OK |
| O010 | `docs/Evidence/ALIGN-SEMI-001-Evidence.md:15` | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | LINK-OK |
| O011 | `docs/Evidence/ALIGN-SEMI-001-Evidence.md:128` | https://en.wikipedia.org/wiki/Sequence_alignment | LINK-OK |
| O012 | `docs/Evidence/ALIGN-SEMI-001-Evidence.md:129` | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | LINK-OK |
| O016 | `docs/Evidence/ANNOT-REPEAT-001-Evidence.md:139` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O017 | `docs/Evidence/ANNOT-REPEAT-001-Evidence.md:140` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O019 | `docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md:158` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O021 | `docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md:127` | https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics) | LINK-OK |
| O024 | `docs/Evidence/ASSEMBLY-STATS-001-Evidence.md:149` | https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics | LINK-OK |
| O025 | `docs/Evidence/CHROM-ANEU-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Aneuploidy | LINK-OK |
| O026 | `docs/Evidence/CHROM-ANEU-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Copy_number_variation | LINK-OK |
| O033 | `docs/Evidence/CHROM-SYNT-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Synteny | LINK-OK |
| O034 | `docs/Evidence/CHROM-SYNT-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Comparative_genomics | LINK-OK |
| O035 | `docs/Evidence/CHROM-SYNT-001-Evidence.md:18` | https://en.wikipedia.org/wiki/Chromosomal_rearrangement | LINK-OK |
| O044 | `docs/Evidence/COMPGEN-COMPARE-001-Evidence.md:124` | https://en.wikipedia.org/wiki/Synteny | LINK-OK |
| O046 | `docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md:152` | https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) | LINK-OK |
| O048 | `docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md:142` | https://en.wikipedia.org/wiki/Synteny | LINK-OK |
| O050 | `docs/Evidence/DISORDER-MORF-001-Evidence.md:159` | https://en.wikipedia.org/wiki/Molecular_recognition_feature | LINK-OK |
| O054 | `docs/Evidence/DISORDER-PRED-001-Evidence.md:232` | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | LINK-OK |
| O055 | `docs/Evidence/DISORDER-PRED-001-Evidence.md:233` | https://en.wikipedia.org/wiki/Hydrophilicity_plot | LINK-OK |
| O056 | `docs/Evidence/DISORDER-PRED-001-Evidence.md:234` | https://en.wikipedia.org/wiki/Amino_acid | LINK-OK |
| O058 | `docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md:135` | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | LINK-OK |
| O067 | `docs/Evidence/EPIGEN-CPG-001-Evidence.md:175` | https://en.wikipedia.org/wiki/CpG_site | LINK-OK |
| O070 | `docs/Evidence/EPIGEN-DMR-001-Evidence.md:149` | https://en.wikipedia.org/wiki/Fisher's_exact_test | LINK-OK |
| O072 | `docs/Evidence/GENOMIC-COMMON-001-Evidence.md:100` | https://en.wikipedia.org/wiki/Longest_common_substring | LINK-OK |
| O073 | `docs/Evidence/GENOMIC-COMMON-001-Evidence.md:101` | https://en.wikipedia.org/wiki/Longest_common_substring | LINK-OK |
| O075 | `docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md:123` | https://en.wikipedia.org/wiki/Restriction_site | LINK-OK |
| O077 | `docs/Evidence/GENOMIC-ORF-001-Evidence.md:130` | https://en.wikipedia.org/wiki/Open_reading_frame | LINK-OK |
| O079 | `docs/Evidence/GENOMIC-REPEAT-001-Evidence.md:141` | https://en.wikipedia.org/wiki/Longest_repeated_substring_problem | LINK-OK |
| O081 | `docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md:109` | https://en.wikipedia.org/wiki/Jaccard_index | LINK-OK |
| O083 | `docs/Evidence/GENOMIC-TANDEM-001-Evidence.md:98` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O086 | `docs/Evidence/KMER-ASYNC-001-Evidence.md:192` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O088 | `docs/Evidence/KMER-GENERATE-001-Evidence.md:166` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O091 | `docs/Evidence/KMER-POSITIONS-001-Evidence.md:131` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O094 | `docs/Evidence/KMER-STATS-001-Evidence.md:158` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O096 | `docs/Evidence/KMER-UNIQUE-001-Evidence.md:163` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O101 | `docs/Evidence/META-BETA-001-Evidence.md:15` | https://en.wikipedia.org/wiki/Beta_diversity | LINK-OK |
| O102 | `docs/Evidence/META-BETA-001-Evidence.md:20` | https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity | LINK-OK |
| O103 | `docs/Evidence/META-BETA-001-Evidence.md:25` | https://en.wikipedia.org/wiki/Jaccard_index | LINK-OK |
| O110 | `docs/Evidence/META-TAXA-001-Evidence.md:139` | https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test | LINK-OK |
| O112 | `docs/Evidence/MIRNA-PAIR-001-Evidence.md:178` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O114 | `docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md:199` | https://en.wikipedia.org/wiki/MicroRNA | LINK-OK |
| O115 | `docs/Evidence/MIRNA-SEED-001-Evidence.md:7` | https://en.wikipedia.org/wiki/MicroRNA | LINK-OK |
| O116 | `docs/Evidence/MIRNA-TARGET-001-Evidence.md:12` | https://en.wikipedia.org/wiki/MicroRNA | LINK-OK |
| O118 | `docs/Evidence/MOTIF-CONS-001-Evidence.md:171` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O120 | `docs/Evidence/MOTIF-GENERATE-001-Evidence.md:146` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O127 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:253` | https://en.wikipedia.org/wiki/TATA_box | LINK-OK |
| O128 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:253` | https://en.wikipedia.org/wiki/Pribnow_box | LINK-OK |
| O129 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:253` | https://en.wikipedia.org/wiki/GC_box | LINK-OK |
| O130 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:253` | https://en.wikipedia.org/wiki/E-box | LINK-OK |
| O131 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:253` | https://en.wikipedia.org/wiki/Shine%E2%80%93Dalgarno_sequence | LINK-OK |
| O132 | `docs/Evidence/MOTIF-REGULATORY-001-Evidence.md:253` | https://en.wikipedia.org/wiki/Kozak_consensus_sequence | LINK-OK |
| O134 | `docs/Evidence/ONCO-CCF-001-Evidence.md:194` | https://en.wikipedia.org/wiki/K-means_clustering | LINK-OK |
| O136 | `docs/Evidence/ONCO-FUSION-001-Evidence.md:248` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O138 | `docs/Evidence/ONCO-FUSION-003-Evidence.md:138` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O139 | `docs/Evidence/ONCO-HETERO-001-Evidence.md:167` | https://en.wikipedia.org/wiki/Diversity_index#Shannon_index | LINK-OK |
| O144 | `docs/Evidence/ONCO-SIG-001-Evidence.md:178` | https://en.wikipedia.org/wiki/Complementarity_(molecular_biology) | LINK-OK |
| O147 | `docs/Evidence/ONCO-SIG-002-Evidence.md:318` | https://en.wikipedia.org/wiki/Non-negative_least_squares | LINK-OK |
| O148 | `docs/Evidence/ONCO-SIG-002-Evidence.md:324` | https://en.wikipedia.org/wiki/Non-negative_matrix_factorization | LINK-OK |
| O150 | `docs/Evidence/ONCO-SIG-003-Evidence.md:181` | https://en.wikipedia.org/wiki/Quantile | LINK-OK |
| O151 | `docs/Evidence/ONCO-SIG-003-Evidence.md:182` | https://en.wikipedia.org/wiki/Non-negative_least_squares | LINK-OK |
| O152 | `docs/Evidence/ONCO-SIG-003-Evidence.md:183` | https://en.wikipedia.org/wiki/Poisson_distribution#Random_variate_generation | LINK-OK |
| O154 | `docs/Evidence/ONCO-VAF-001-Evidence.md:139` | https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval | LINK-OK |
| O156 | `docs/Evidence/PANGEN-CORE-001-Evidence.md:192` | https://en.wikipedia.org/wiki/Pan-genome | LINK-OK |
| O159 | `docs/Evidence/PANGEN-MARKER-001-Evidence.md:140` | https://en.wikipedia.org/wiki/Informative_site | LINK-OK |
| O161 | `docs/Evidence/PARSE-BED-001-Evidence.md:174` | https://en.wikipedia.org/wiki/BED_(file_format) | LINK-OK |
| O165 | `docs/Evidence/PARSE-FASTQ-001-Evidence.md:159` | https://en.wikipedia.org/wiki/FASTQ_format | LINK-OK |
| O169 | `docs/Evidence/PHYLO-COMP-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Robinson–Foulds_metric | LINK-OK |
| O170 | `docs/Evidence/PHYLO-COMP-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Most_recent_common_ancestor | LINK-OK |
| O171 | `docs/Evidence/PHYLO-COMP-001-Evidence.md:18` | https://en.wikipedia.org/wiki/Phylogenetic_tree | LINK-OK |
| O176 | `docs/Evidence/PHYLO-NEWICK-001-Evidence.md:196` | https://en.wikipedia.org/wiki/Newick_format | LINK-OK |
| O180 | `docs/Evidence/PHYLO-STATS-001-Evidence.md:144` | https://en.wikipedia.org/wiki/Tree_(graph_theory) | LINK-OK |
| O181 | `docs/Evidence/PHYLO-STATS-001-Evidence.md:145` | https://en.wikipedia.org/wiki/Tree_(abstract_data_type) | LINK-OK |
| O182 | `docs/Evidence/PHYLO-STATS-001-Evidence.md:148` | https://en.wikipedia.org/wiki/Minimum_evolution | LINK-OK |
| O183 | `docs/Evidence/PHYLO-TREE-001-Evidence.md:16` | https://en.wikipedia.org/wiki/UPGMA | LINK-OK |
| O184 | `docs/Evidence/PHYLO-TREE-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Neighbor_joining | LINK-OK |
| O185 | `docs/Evidence/PHYLO-TREE-001-Evidence.md:18` | https://en.wikipedia.org/wiki/Phylogenetic_tree | LINK-OK |
| O186 | `docs/Evidence/POP-DIV-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Nucleotide_diversity | LINK-OK |
| O187 | `docs/Evidence/POP-DIV-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Watterson_estimator | LINK-OK |
| O188 | `docs/Evidence/POP-DIV-001-Evidence.md:18` | https://en.wikipedia.org/wiki/Tajima%27s_D | LINK-OK |
| O189 | `docs/Evidence/POP-DIV-001-Evidence.md:19` | https://en.wikipedia.org/wiki/Zygosity | LINK-OK |
| O190 | `docs/Evidence/POP-DIV-001-Evidence.md:119` | https://en.wikipedia.org/wiki/Tajima%27s_D#Example | LINK-OK |
| O191 | `docs/Evidence/POP-FREQ-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Allele_frequency | LINK-OK |
| O192 | `docs/Evidence/POP-FREQ-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Minor_allele_frequency | LINK-OK |
| O193 | `docs/Evidence/POP-FREQ-001-Evidence.md:18` | https://en.wikipedia.org/wiki/Genotype_frequency | LINK-OK |
| O194 | `docs/Evidence/POP-FST-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Fixation_index | LINK-OK |
| O195 | `docs/Evidence/POP-FST-001-Evidence.md:17` | https://en.wikipedia.org/wiki/F-statistics | LINK-OK |
| O196 | `docs/Evidence/POP-HW-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Hardy-Weinberg_principle | LINK-OK |
| O197 | `docs/Evidence/POP-HW-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Chi-squared_test | LINK-OK |
| O202 | `docs/Evidence/PROBE-VALID-001-Evidence.md:231` | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | LINK-OK |
| O203 | `docs/Evidence/PROBE-VALID-001-Evidence.md:232` | https://en.wikipedia.org/wiki/BLAST_(biotechnology) | LINK-OK |
| O206 | `docs/Evidence/PROTMOTIF-CC-001-Evidence.md:130` | https://en.wikipedia.org/wiki/Coiled_coil | LINK-OK |
| O207 | `docs/Evidence/PROTMOTIF-CC-001-Evidence.md:131` | https://en.wikipedia.org/wiki/Heptad_repeat | LINK-OK |
| O213 | `docs/Evidence/PROTMOTIF-FIND-001-Evidence.md:189` | https://en.wikipedia.org/wiki/Sequence_motif | LINK-OK |
| O214 | `docs/Evidence/PROTMOTIF-FIND-001-Evidence.md:190` | https://en.wikipedia.org/wiki/PROSITE | LINK-OK |
| O216 | `docs/Evidence/QUALITY-STATS-001-Evidence.md:144` | https://en.wikipedia.org/wiki/Phred_quality_score | LINK-OK |
| O217 | `docs/Evidence/QUALITY-STATS-001-Evidence.md:145` | https://en.wikipedia.org/wiki/Phred_quality_score | LINK-OK |
| O219 | `docs/Evidence/REP-STR-001-Evidence.md:382` | https://en.wikipedia.org/wiki/Microsatellite | LINK-OK |
| O223 | `docs/Evidence/RESTR-FILTER-001-Evidence.md:135` | https://en.wikipedia.org/wiki/Sticky_and_blunt_ends | LINK-OK |
| O224 | `docs/Evidence/RESTR-FILTER-001-Evidence.md:136` | https://en.wikipedia.org/wiki/Restriction_enzyme | LINK-OK |
| O225 | `docs/Evidence/RESTR-FILTER-001-Evidence.md:137` | https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites | LINK-OK |
| O227 | `docs/Evidence/RNA-INVERT-001-Evidence.md:119` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O229 | `docs/Evidence/RNA-PAIR-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O234 | `docs/Evidence/RNA-PAIR-001-Evidence.md:151` | https://en.wikipedia.org/wiki/Base_pair | LINK-OK |
| O235 | `docs/Evidence/RNA-PAIR-001-Evidence.md:152` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O236 | `docs/Evidence/RNA-PAIR-001-Evidence.md:153` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O238 | `docs/Evidence/RNA-PKPREDICT-001-Evidence.md:140` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O240 | `docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md:132` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O241 | `docs/Evidence/RNA-STEMLOOP-001-Evidence.md:23` | https://en.wikipedia.org/wiki/Stem-loop | LINK-OK |
| O242 | `docs/Evidence/RNA-STEMLOOP-001-Evidence.md:24` | https://en.wikipedia.org/wiki/Tetraloop | LINK-OK |
| O243 | `docs/Evidence/RNA-STEMLOOP-001-Evidence.md:25` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O245 | `docs/Evidence/SEQ-ATSKEW-001-Evidence.md:113` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O247 | `docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md:141` | https://en.wikipedia.org/wiki/Codon_usage_bias | LINK-OK |
| O249 | `docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md:212` | https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity | LINK-OK |
| O251 | `docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md:148` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O254 | `docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md:184` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O255 | `docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md:186` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O256 | `docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md:187` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O257 | `docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md:188` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O259 | `docs/Evidence/SEQ-COMPOSITION-001-Evidence.md:110` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O262 | `docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md:123` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O263 | `docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md:124` | https://en.wikipedia.org/wiki/A_Mathematical_Theory_of_Communication | LINK-OK |
| O266 | `docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md:140` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O267 | `docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md:141` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O268 | `docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md:142` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O270 | `docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md:157` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O272 | `docs/Evidence/SEQ-REPLICATION-001-Evidence.md:191` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O274 | `docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md:285` | https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method | LINK-OK |
| O276 | `docs/Evidence/SEQ-STATS-001-Evidence.md:100` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O279 | `docs/Evidence/SEQ-SUMMARY-001-Evidence.md:132` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O280 | `docs/Evidence/SEQ-SUMMARY-001-Evidence.md:133` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O282 | `docs/Evidence/SEQ-THERMO-001-Evidence.md:160` | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | LINK-OK |
| O288 | `docs/Evidence/SPLICE-DONOR-001-Evidence.md:210` | https://en.wikipedia.org/wiki/RNA_splicing | LINK-OK |
| O289 | `docs/Evidence/SPLICE-DONOR-001-Evidence.md:211` | https://en.wikipedia.org/wiki/Spliceosome | LINK-OK |
| O293 | `docs/Evidence/TRANS-CODON-001-Evidence.md:15` | https://en.wikipedia.org/wiki/Genetic_code | LINK-OK |
| O294 | `docs/Evidence/TRANS-CODON-001-Evidence.md:16` | https://en.wikipedia.org/wiki/Start_codon | LINK-OK |
| O295 | `docs/Evidence/TRANS-CODON-001-Evidence.md:17` | https://en.wikipedia.org/wiki/Stop_codon | LINK-OK |
| O299 | `docs/Evidence/TRANS-DIFF-001-Evidence.md:170` | https://en.wikipedia.org/wiki/False_discovery_rate | LINK-OK |
| O300 | `docs/Evidence/TRANS-DIFF-001-Evidence.md:171` | https://en.wikipedia.org/wiki/Welch%27s_t-test | LINK-OK |
| O301 | `docs/Evidence/TRANS-DIFF-001-Evidence.md:172` | https://en.wikipedia.org/wiki/Student%27s_t-distribution | LINK-OK |
| O303 | `docs/Evidence/TRANS-EXPR-001-Evidence.md:149` | https://en.wikipedia.org/wiki/Quantile_normalization | LINK-OK |
| O308 | `docs/Evidence/TRANS-SIXFRAME-001-Evidence.md:276` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O311 | `docs/Evidence/VARIANT-CALL-001-Evidence.md:199` | https://en.wikipedia.org/wiki/Transition_(genetics) | LINK-OK |
| O312 | `docs/Evidence/VARIANT-CALL-001-Evidence.md:200` | https://en.wikipedia.org/wiki/Transversion | LINK-OK |
| O314 | `docs/Evidence/VARIANT-SNP-001-Evidence.md:134` | https://en.wikipedia.org/wiki/Transversion | LINK-OK |
| O315 | `docs/Validation/reports/ALIGN-LOCAL-001.md:16` | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | LINK-OK |
| O316 | `docs/Validation/reports/ALIGN-MULTI-001.md:28` | https://en.wikipedia.org/wiki/UPGMA | LINK-OK |
| O317 | `docs/Validation/reports/ALIGN-MULTI-001.md:30` | https://en.wikipedia.org/wiki/Multiple_sequence_alignment | LINK-OK |
| O318 | `docs/Validation/reports/ALIGN-MULTI-001.md:154` | https://en.wikipedia.org/wiki/UPGMA | LINK-OK |
| O319 | `docs/Validation/reports/ANNOT-REPEAT-001.md:11` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O320 | `docs/Validation/reports/ANNOT-REPEAT-001.md:13` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O321 | `docs/Validation/reports/ASSEMBLY-STATS-001.md:26` | https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics | LINK-OK |
| O322 | `docs/Validation/reports/CODON-USAGE-001.md:167` | https://en.wikipedia.org/wiki/Codon_usage_bias | LINK-OK |
| O323 | `docs/Validation/reports/COMPGEN-DOTPLOT-001.md:17` | https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) | LINK-OK |
| O324 | `docs/Validation/reports/DISORDER-MORF-001.md:19` | https://en.wikipedia.org/wiki/Molecular_recognition_feature | LINK-OK |
| O325 | `docs/Validation/reports/DISORDER-PROPENSITY-001.md:29` | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | LINK-OK |
| O326 | `docs/Validation/reports/EPIGEN-DMR-001.md:16` | https://en.wikipedia.org/wiki/Fisher%27s_exact_test | LINK-OK |
| O327 | `docs/Validation/reports/GENOMIC-COMMON-001.md:12` | https://en.wikipedia.org/wiki/Longest_common_substring | LINK-OK |
| O328 | `docs/Validation/reports/GENOMIC-MOTIFS-001.md:28` | https://en.wikipedia.org/wiki/EcoRI | LINK-OK |
| O329 | `docs/Validation/reports/GENOMIC-ORF-001.md:14` | https://en.wikipedia.org/wiki/Open_reading_frame | LINK-OK |
| O330 | `docs/Validation/reports/GENOMIC-SIMILARITY-001.md:12` | https://en.wikipedia.org/wiki/Jaccard_index | LINK-OK |
| O331 | `docs/Validation/reports/GENOMIC-TANDEM-001.md:23` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O332 | `docs/Validation/reports/GENOMIC-TANDEM-001.md:32` | https://en.wikipedia.org/wiki/Microsatellite | LINK-OK |
| O333 | `docs/Validation/reports/KMER-ASYNC-001.md:14` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O334 | `docs/Validation/reports/KMER-DIST-001.md:17` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O335 | `docs/Validation/reports/KMER-GENERATE-001.md:14` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O336 | `docs/Validation/reports/KMER-POSITIONS-001.md:16` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O337 | `docs/Validation/reports/KMER-STATS-001.md:13` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O338 | `docs/Validation/reports/KMER-UNIQUE-001.md:12` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O339 | `docs/Validation/reports/MOTIF-CONS-001.md:22` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O340 | `docs/Validation/reports/MOTIF-GENERATE-001.md:15` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O341 | `docs/Validation/reports/ONCO-FUSION-001.md:42` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O342 | `docs/Validation/reports/ONCO-FUSION-003.md:14` | https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables | LINK-OK |
| O343 | `docs/Validation/reports/ONCO-SIG-001.md:18` | https://en.wikipedia.org/wiki/Mutational_signatures | LINK-OK |
| O344 | `docs/Validation/reports/ONCO-SIG-002.md:25` | https://en.wikipedia.org/wiki/Non-negative_matrix_factorization | LINK-OK |
| O345 | `docs/Validation/reports/PANGEN-MARKER-001.md:13` | https://en.wikipedia.org/wiki/Informative_site | LINK-OK |
| O346 | `docs/Validation/reports/PARSE-BED-001.md:21` | https://en.wikipedia.org/wiki/BED_(file_format) | LINK-OK |
| O347 | `docs/Validation/reports/PARSE-FASTA-001.md:11` | https://en.wikipedia.org/wiki/FASTA_format | LINK-OK |
| O348 | `docs/Validation/reports/PARSE-FASTA-001.md:147` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O349 | `docs/Validation/reports/PARSE-FASTQ-001.md:13` | https://en.wikipedia.org/wiki/FASTQ_format | LINK-OK |
| O350 | `docs/Validation/reports/PAT-EXACT-001.md:15` | https://en.wikipedia.org/wiki/Suffix_tree | LINK-OK |
| O351 | `docs/Validation/reports/PAT-IUPAC-001.md:11` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O352 | `docs/Validation/reports/PAT-PWM-001.md:16` | https://en.wikipedia.org/wiki/Position_weight_matrix | LINK-OK |
| O353 | `docs/Validation/reports/PHYLO-NEWICK-001.md:14` | https://en.wikipedia.org/wiki/Newick_format | LINK-OK |
| O354 | `docs/Validation/reports/PRIMER-DESIGN-001.md:19` | https://en.wikipedia.org/wiki/Primer_(molecular_biology) | LINK-OK |
| O355 | `docs/Validation/reports/QUALITY-PHRED-001.md:33` | https://en.wikipedia.org/wiki/FASTQ_format | LINK-OK |
| O356 | `docs/Validation/reports/RNA-INVERT-001.md:12` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O357 | `docs/Validation/reports/RNA-PAIR-001.md:12` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O358 | `docs/Validation/reports/RNA-PAIR-001.md:18` | https://en.wikipedia.org/wiki/Base_pair | LINK-OK |
| O359 | `docs/Validation/reports/RNA-PAIR-001.md:23` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O360 | `docs/Validation/reports/RNA-PSEUDOKNOT-001.md:11` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O361 | `docs/Validation/reports/SEQ-ATSKEW-001.md:21` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O362 | `docs/Validation/reports/SEQ-COMP-001.md:14` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O363 | `docs/Validation/reports/SEQ-COMPLEX-WINDOW-001.md:20` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O364 | `docs/Validation/reports/SEQ-COMPLEX-WINDOW-001.md:22` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O365 | `docs/Validation/reports/SEQ-GC-ANALYSIS-001.md:20` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O366 | `docs/Validation/reports/SEQ-GC-ANALYSIS-001.md:21` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O367 | `docs/Validation/reports/SEQ-GC-PROFILE-001.md:12` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O368 | `docs/Validation/reports/SEQ-REPLICATION-001.md:16` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O369 | `docs/Validation/reports/SEQ-REVCOMP-001.md:12` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O370 | `docs/Validation/reports/SEQ-SUMMARY-001.md:37` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O371 | `docs/Validation/reports/TRANS-EXPR-001.md:16` | https://en.wikipedia.org/wiki/Quantile_normalization | LINK-OK |
| O372 | `docs/Validation/reports/VARIANT-SNP-001.md:17` | https://en.wikipedia.org/wiki/Hamming_distance | LINK-OK |
| O373 | `docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md:156` | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | LINK-OK |
| O374 | `docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md:157` | https://en.wikipedia.org/wiki/Sequence_alignment | LINK-OK |
| O375 | `docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md:170` | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | LINK-OK |
| O376 | `docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md:171` | https://en.wikipedia.org/wiki/Sequence_alignment | LINK-OK |
| O377 | `docs/algorithms/Alignment/Multiple_Sequence_Alignment.md:186` | https://en.wikipedia.org/wiki/Multiple_sequence_alignment | LINK-OK |
| O378 | `docs/algorithms/Alignment/Multiple_Sequence_Alignment.md:187` | https://en.wikipedia.org/wiki/Clustal | LINK-OK |
| O379 | `docs/algorithms/Alignment/Multiple_Sequence_Alignment.md:188` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O380 | `docs/algorithms/Alignment/Semi_Global_Alignment.md:159` | https://en.wikipedia.org/wiki/Sequence_alignment | LINK-OK |
| O381 | `docs/algorithms/Alignment/Semi_Global_Alignment.md:160` | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | LINK-OK |
| O382 | `docs/algorithms/Analysis/Open_Reading_Frame_Detection.md:184` | https://en.wikipedia.org/wiki/Open_reading_frame | LINK-OK |
| O383 | `docs/algorithms/Analysis/Sequence_Similarity.md:185` | https://en.wikipedia.org/wiki/Jaccard_index | LINK-OK |
| O384 | `docs/algorithms/Annotation/GFF3_IO.md:182` | https://en.wikipedia.org/wiki/General_feature_format | LINK-OK |
| O385 | `docs/algorithms/Annotation/Gene_Prediction.md:176` | https://en.wikipedia.org/wiki/Gene_prediction | LINK-OK |
| O386 | `docs/algorithms/Annotation/Gene_Prediction.md:177` | https://en.wikipedia.org/wiki/Shine-Dalgarno_sequence | LINK-OK |
| O387 | `docs/algorithms/Annotation/Gene_Prediction.md:178` | https://en.wikipedia.org/wiki/Ribosome-binding_site | LINK-OK |
| O388 | `docs/algorithms/Annotation/ORF_Detection.md:173` | https://en.wikipedia.org/wiki/Open_reading_frame | LINK-OK |
| O389 | `docs/algorithms/Annotation/Promoter_Detection.md:154` | https://en.wikipedia.org/wiki/Promoter_(genetics) | LINK-OK |
| O390 | `docs/algorithms/Annotation/Promoter_Detection.md:155` | https://en.wikipedia.org/wiki/Pribnow_box | LINK-OK |
| O391 | `docs/algorithms/Annotation/Promoter_Detection.md:156` | https://en.wikipedia.org/wiki/TATA_box | LINK-OK |
| O392 | `docs/algorithms/Annotation/Repetitive_Element_Detection.md:166` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O393 | `docs/algorithms/Annotation/Repetitive_Element_Detection.md:167` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O394 | `docs/algorithms/Assembly/Assembly_Statistics.md:176` | https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics | LINK-OK |
| O395 | `docs/algorithms/Assembly/Consensus_Computation.md:175` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O396 | `docs/algorithms/Chromosome_Analysis/Aneuploidy_Detection.md:199` | https://en.wikipedia.org/wiki/Aneuploidy | LINK-OK |
| O397 | `docs/algorithms/Chromosome_Analysis/Aneuploidy_Detection.md:200` | https://en.wikipedia.org/wiki/Copy_number_variation | LINK-OK |
| O398 | `docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md:185` | https://en.wikipedia.org/wiki/Centromere | LINK-OK |
| O399 | `docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md:188` | https://en.wikipedia.org/wiki/Karyotype | LINK-OK |
| O400 | `docs/algorithms/Chromosome_Analysis/Karyotype_Analysis.md:177` | https://en.wikipedia.org/wiki/Karyotype | LINK-OK |
| O401 | `docs/algorithms/Chromosome_Analysis/Karyotype_Analysis.md:178` | https://en.wikipedia.org/wiki/Ploidy | LINK-OK |
| O402 | `docs/algorithms/Chromosome_Analysis/Karyotype_Analysis.md:179` | https://en.wikipedia.org/wiki/Aneuploidy | LINK-OK |
| O403 | `docs/algorithms/Chromosome_Analysis/Telomere_Analysis.md:196` | https://en.wikipedia.org/wiki/Telomere | LINK-OK |
| O404 | `docs/algorithms/Codon_Optimization/Codon_Usage_Analysis.md:174` | https://en.wikipedia.org/wiki/Codon_usage_bias | LINK-OK |
| O405 | `docs/algorithms/Comparative_Genomics/Dot_Plot_Generation.md:163` | https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) | LINK-OK |
| O406 | `docs/algorithms/Comparative_Genomics/Genome_Comparison.md:176` | https://en.wikipedia.org/wiki/Synteny | LINK-OK |
| O407 | `docs/algorithms/Comparative_Genomics/Synteny_Block_Detection.md:226` | https://en.wikipedia.org/wiki/Synteny | LINK-OK |
| O408 | `docs/algorithms/Complexity/K-mer_Entropy.md:170` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O409 | `docs/algorithms/Complexity/Lempel_Ziv_Complexity.md:162` | https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity | LINK-OK |
| O410 | `docs/algorithms/Complexity/Windowed_Complexity.md:205` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O411 | `docs/algorithms/Epigenetics/CpG_Site_Detection.md:148` | https://en.wikipedia.org/wiki/CpG_site | LINK-OK |
| O412 | `docs/algorithms/Epigenetics/Differentially_Methylated_Regions.md:193` | https://en.wikipedia.org/wiki/Fisher's_exact_test | LINK-OK |
| O413 | `docs/algorithms/Extended_Assembly/Scaffolding.md:197` | https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics) | LINK-OK |
| O414 | `docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md:165` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O415 | `docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md:164` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O416 | `docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md:165` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O417 | `docs/algorithms/FileIO/BED_Parsing.md:197` | https://en.wikipedia.org/wiki/BED_(file_format) | LINK-OK |
| O418 | `docs/algorithms/FileIO/FASTA_Parsing.md:208` | https://en.wikipedia.org/wiki/FASTA_format | LINK-OK |
| O419 | `docs/algorithms/FileIO/FASTQ_Parsing.md:209` | https://en.wikipedia.org/wiki/FASTQ_format | LINK-OK |
| O420 | `docs/algorithms/FileIO/GFF_Parsing.md:195` | https://en.wikipedia.org/wiki/General_feature_format | LINK-OK |
| O421 | `docs/algorithms/FileIO/VCF_Parsing.md:205` | https://en.wikipedia.org/wiki/Variant_Call_Format | LINK-OK |
| O422 | `docs/algorithms/Genomic_Analysis/Tandem_Repeat_Detection.md:178` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O423 | `docs/algorithms/K-mer/Asynchronous_K-mer_Counting.md:224` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O424 | `docs/algorithms/K-mer/K-mer_Counting.md:175` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O425 | `docs/algorithms/K-mer/K-mer_Frequency_Analysis.md:191` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O426 | `docs/algorithms/K-mer/K-mer_Frequency_Analysis.md:192` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O427 | `docs/algorithms/K-mer/K-mer_Generation.md:133` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O428 | `docs/algorithms/K-mer/K-mer_Positions.md:156` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O429 | `docs/algorithms/K-mer/K-mer_Search.md:194` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O430 | `docs/algorithms/K-mer/K-mer_Statistics.md:178` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O431 | `docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md:198` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O432 | `docs/algorithms/Metagenomics/Alpha_Diversity.md:189` | https://en.wikipedia.org/wiki/Diversity_index | LINK-OK |
| O433 | `docs/algorithms/Metagenomics/Alpha_Diversity.md:190` | https://en.wikipedia.org/wiki/Alpha_diversity | LINK-OK |
| O434 | `docs/algorithms/Metagenomics/Genome_Binning.md:206` | https://en.wikipedia.org/wiki/Binning_(metagenomics) | LINK-OK |
| O435 | `docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md:194` | https://en.wikipedia.org/wiki/Pan-genome | LINK-OK |
| O436 | `docs/algorithms/Metagenomics/Significant_Taxa_Detection.md:182` | https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test | LINK-OK |
| O437 | `docs/algorithms/Metagenomics/Taxonomic_Classification.md:195` | https://en.wikipedia.org/wiki/Metagenomics | LINK-OK |
| O438 | `docs/algorithms/Metagenomics/Taxonomic_Profile.md:160` | https://en.wikipedia.org/wiki/Metagenomics | LINK-OK |
| O439 | `docs/algorithms/MiRNA/MiRNA_Target_Pairing.md:180` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O440 | `docs/algorithms/MolTools/Guide_RNA_Design.md:217` | https://www.addgene.org/guides/crispr/ | LINK-OK |
| O441 | `docs/algorithms/MolTools/Guide_RNA_Design.md:218` | https://en.wikipedia.org/wiki/Guide_RNA | LINK-OK |
| O442 | `docs/algorithms/MolTools/Guide_RNA_Design.md:219` | https://en.wikipedia.org/wiki/CRISPR_gene_editing | LINK-OK |
| O443 | `docs/algorithms/MolTools/Guide_RNA_Design.md:220` | https://en.wikipedia.org/wiki/Protospacer_adjacent_motif | LINK-OK |
| O444 | `docs/algorithms/MolTools/Hybridization_Probe_Design.md:228` | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | LINK-OK |
| O445 | `docs/algorithms/MolTools/Hybridization_Probe_Design.md:229` | https://en.wikipedia.org/wiki/Hybridization_probe | LINK-OK |
| O446 | `docs/algorithms/MolTools/Hybridization_Probe_Design.md:230` | https://en.wikipedia.org/wiki/Fluorescence_in_situ_hybridization | LINK-OK |
| O447 | `docs/algorithms/MolTools/Hybridization_Probe_Design.md:231` | https://en.wikipedia.org/wiki/DNA_microarray | LINK-OK |
| O448 | `docs/algorithms/MolTools/Melting_Temperature.md:231` | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | LINK-OK |
| O449 | `docs/algorithms/MolTools/Melting_Temperature.md:232` | https://en.wikipedia.org/wiki/DNA_melting | LINK-OK |
| O450 | `docs/algorithms/MolTools/Off_Target_Analysis.md:215` | https://en.wikipedia.org/wiki/Off-target_genome_editing | LINK-OK |
| O451 | `docs/algorithms/MolTools/PAM_Site_Detection.md:173` | https://en.wikipedia.org/wiki/Protospacer_adjacent_motif | LINK-OK |
| O452 | `docs/algorithms/MolTools/PAM_Site_Detection.md:174` | https://en.wikipedia.org/wiki/CRISPR | LINK-OK |
| O453 | `docs/algorithms/MolTools/Primer_Design.md:401` | https://en.wikipedia.org/wiki/Primer_(molecular_biology) | LINK-OK |
| O454 | `docs/algorithms/MolTools/Primer_Design.md:402` | https://www.addgene.org/protocols/primer-design/ | LINK-OK |
| O455 | `docs/algorithms/MolTools/Probe_Validation.md:215` | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | LINK-OK |
| O456 | `docs/algorithms/MolTools/Probe_Validation.md:216` | https://en.wikipedia.org/wiki/BLAST_(biotechnology) | LINK-OK |
| O457 | `docs/algorithms/MolTools/Probe_Validation.md:217` | https://en.wikipedia.org/wiki/Hybridization_probe | LINK-OK |
| O458 | `docs/algorithms/MolTools/Probe_Validation.md:218` | https://en.wikipedia.org/wiki/DNA_microarray | LINK-OK |
| O459 | `docs/algorithms/MolTools/Probe_Validation.md:219` | https://en.wikipedia.org/wiki/Off-target_genome_editing | LINK-OK |
| O460 | `docs/algorithms/MolTools/Restriction_Digest_Simulation.md:186` | https://en.wikipedia.org/wiki/Restriction_digest | LINK-OK |
| O461 | `docs/algorithms/MolTools/Restriction_Digest_Simulation.md:187` | https://en.wikipedia.org/wiki/Restriction_enzyme | LINK-OK |
| O462 | `docs/algorithms/MolTools/Restriction_Digest_Simulation.md:188` | https://en.wikipedia.org/wiki/Restriction_map | LINK-OK |
| O463 | `docs/algorithms/MolTools/Restriction_Digest_Simulation.md:189` | https://www.addgene.org/protocols/restriction-digest/ | LINK-OK |
| O464 | `docs/algorithms/MolTools/Restriction_Digest_Simulation.md:194` | https://blog.addgene.org/plasmids-101-restriction-cloning | LINK-OK |
| O465 | `docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md:326` | https://en.wikipedia.org/wiki/Sticky_and_blunt_ends | LINK-OK |
| O466 | `docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md:327` | https://en.wikipedia.org/wiki/Restriction_enzyme | LINK-OK |
| O467 | `docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md:328` | https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites | LINK-OK |
| O468 | `docs/algorithms/MolTools/Restriction_Site_Detection.md:209` | https://en.wikipedia.org/wiki/Restriction_enzyme | LINK-OK |
| O469 | `docs/algorithms/MolTools/Restriction_Site_Detection.md:210` | https://en.wikipedia.org/wiki/Restriction_site | LINK-OK |
| O470 | `docs/algorithms/MolTools/Restriction_Site_Detection.md:211` | https://en.wikipedia.org/wiki/EcoRI | LINK-OK |
| O471 | `docs/algorithms/Motif_Analysis/Known_Motif_Search.md:189` | https://en.wikipedia.org/wiki/Restriction_site | LINK-OK |
| O472 | `docs/algorithms/Oncology/Fusion_Breakpoint_Analysis.md:178` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O473 | `docs/algorithms/Oncology/Fusion_Gene_Detection.md:218` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O474 | `docs/algorithms/Oncology/Mutational_Signature_Exposure_Bootstrap.md:227` | https://en.wikipedia.org/wiki/Non-negative_least_squares | LINK-OK |
| O475 | `docs/algorithms/Oncology/Mutational_Signature_Exposure_Bootstrap.md:228` | https://en.wikipedia.org/wiki/Poisson_distribution#Random_variate_generation | LINK-OK |
| O476 | `docs/algorithms/Oncology/Mutational_Signature_Extraction_NMF.md:315` | https://en.wikipedia.org/wiki/Non-negative_matrix_factorization | LINK-OK |
| O477 | `docs/algorithms/Oncology/Mutational_Signature_Fitting.md:227` | https://en.wikipedia.org/wiki/Non-negative_least_squares | LINK-OK |
| O478 | `docs/algorithms/Oncology/SBS96_Trinucleotide_Context_Catalog.md:226` | https://en.wikipedia.org/wiki/Complementarity_(molecular_biology) | LINK-OK |
| O479 | `docs/algorithms/Oncology/Tumor_Heterogeneity_Analysis.md:199` | https://en.wikipedia.org/wiki/Diversity_index#Shannon_index | LINK-OK |
| O480 | `docs/algorithms/Oncology/Variant_Allele_Frequency.md:186` | https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval | LINK-OK |
| O481 | `docs/algorithms/PanGenome/Phylogenetic_Marker_Selection.md:160` | https://en.wikipedia.org/wiki/Informative_site | LINK-OK |
| O482 | `docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md:147` | https://en.wikipedia.org/wiki/Hamming_distance | LINK-OK |
| O483 | `docs/algorithms/Pattern_Matching/Consensus_From_Alignment.md:182` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O484 | `docs/algorithms/Pattern_Matching/Edit_Distance.md:212` | https://en.wikipedia.org/wiki/Levenshtein_distance | LINK-OK |
| O485 | `docs/algorithms/Pattern_Matching/Edit_Distance.md:213` | https://en.wikipedia.org/wiki/Edit_distance | LINK-OK |
| O486 | `docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md:155` | https://en.wikipedia.org/wiki/Suffix_tree | LINK-OK |
| O487 | `docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Consensus.md:188` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O488 | `docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md:168` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O489 | `docs/algorithms/Pattern_Matching/Position_Weight_Matrix.md:206` | https://en.wikipedia.org/wiki/Position_weight_matrix | LINK-OK |
| O490 | `docs/algorithms/Phylogenetics/Distance_Matrix.md:320` | https://en.wikipedia.org/wiki/Models_of_DNA_evolution | LINK-OK |
| O491 | `docs/algorithms/Phylogenetics/Distance_Matrix.md:321` | https://en.wikipedia.org/wiki/Distance_matrices_in_phylogeny | LINK-OK |
| O492 | `docs/algorithms/Phylogenetics/Newick_Format.md:143` | https://en.wikipedia.org/wiki/Newick_format | LINK-OK |
| O493 | `docs/algorithms/Phylogenetics/Tree_Comparison.md:285` | https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric | LINK-OK |
| O494 | `docs/algorithms/Phylogenetics/Tree_Comparison.md:286` | https://en.wikipedia.org/wiki/Most_recent_common_ancestor | LINK-OK |
| O495 | `docs/algorithms/Phylogenetics/Tree_Comparison.md:287` | https://en.wikipedia.org/wiki/Phylogenetic_tree | LINK-OK |
| O496 | `docs/algorithms/Phylogenetics/Tree_Construction.md:258` | https://en.wikipedia.org/wiki/UPGMA | LINK-OK |
| O497 | `docs/algorithms/Phylogenetics/Tree_Construction.md:259` | https://en.wikipedia.org/wiki/Neighbor_joining | LINK-OK |
| O498 | `docs/algorithms/Phylogenetics/Tree_Statistics.md:151` | https://en.wikipedia.org/wiki/Tree_(graph_theory) | LINK-OK |
| O499 | `docs/algorithms/Phylogenetics/Tree_Statistics.md:152` | https://en.wikipedia.org/wiki/Tree_(abstract_data_type) | LINK-OK |
| O500 | `docs/algorithms/Phylogenetics/Tree_Statistics.md:155` | https://en.wikipedia.org/wiki/Minimum_evolution | LINK-OK |
| O501 | `docs/algorithms/Population_Genetics/Allele_Frequency.md:171` | https://en.wikipedia.org/wiki/Allele_frequency | LINK-OK |
| O502 | `docs/algorithms/Population_Genetics/Allele_Frequency.md:172` | https://en.wikipedia.org/wiki/Minor_allele_frequency | LINK-OK |
| O503 | `docs/algorithms/Population_Genetics/Allele_Frequency.md:173` | https://en.wikipedia.org/wiki/Genotype_frequency | LINK-OK |
| O504 | `docs/algorithms/Population_Genetics/Diversity_Statistics.md:240` | https://en.wikipedia.org/wiki/Tajima%27s_D | LINK-OK |
| O505 | `docs/algorithms/Population_Genetics/F_Statistics.md:293` | https://en.wikipedia.org/wiki/Fixation_index | LINK-OK |
| O506 | `docs/algorithms/Population_Genetics/F_Statistics.md:294` | https://en.wikipedia.org/wiki/F-statistics | LINK-OK |
| O507 | `docs/algorithms/Population_Genetics/Hardy_Weinberg_Test.md:211` | https://en.wikipedia.org/wiki/Hardy-Weinberg_principle | LINK-OK |
| O508 | `docs/algorithms/Population_Genetics/Linkage_Disequilibrium.md:179` | https://en.wikipedia.org/wiki/Linkage_disequilibrium | LINK-OK |
| O509 | `docs/algorithms/Population_Genetics/Linkage_Disequilibrium.md:180` | https://en.wikipedia.org/wiki/Haplotype_block | LINK-OK |
| O510 | `docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md:172` | https://en.wikipedia.org/wiki/Coiled_coil | LINK-OK |
| O511 | `docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md:173` | https://en.wikipedia.org/wiki/Heptad_repeat | LINK-OK |
| O512 | `docs/algorithms/ProteinPred/MoRF_Prediction.md:208` | https://en.wikipedia.org/wiki/Molecular_recognition_feature | LINK-OK |
| O513 | `docs/algorithms/Repeat_Analysis/Direct_Repeat_Detection.md:220` | https://en.wikipedia.org/wiki/Direct_repeat | LINK-OK |
| O514 | `docs/algorithms/Repeat_Analysis/Direct_Repeat_Detection.md:221` | https://en.wikipedia.org/wiki/Repeated_sequence_(DNA) | LINK-OK |
| O515 | `docs/algorithms/Repeat_Analysis/Inverted_Repeat_Detection.md:203` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O516 | `docs/algorithms/Repeat_Analysis/Microsatellite_Detection.md:303` | https://en.wikipedia.org/wiki/Microsatellite | LINK-OK |
| O517 | `docs/algorithms/Repeat_Analysis/Microsatellite_Detection.md:304` | https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder | LINK-OK |
| O518 | `docs/algorithms/Repeat_Analysis/Palindrome_Detection.md:192` | https://en.wikipedia.org/wiki/Palindromic_sequence | LINK-OK |
| O519 | `docs/algorithms/Repeat_Analysis/Palindrome_Detection.md:193` | https://en.wikipedia.org/wiki/Restriction_enzyme | LINK-OK |
| O520 | `docs/algorithms/Repeat_Analysis/Repeat_Detection.md:177` | https://en.wikipedia.org/wiki/Longest_repeated_substring_problem | LINK-OK |
| O521 | `docs/algorithms/Repeat_Analysis/Tandem_Repeat_Detection.md:181` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O522 | `docs/algorithms/Repeat_Analysis/Tandem_Repeat_Detection.md:182` | https://en.wikipedia.org/wiki/Microsatellite | LINK-OK |
| O523 | `docs/algorithms/RnaStructure/Inverted_Repeats.md:163` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O524 | `docs/algorithms/RnaStructure/Pseudoknot_Detection.md:179` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O525 | `docs/algorithms/RnaStructure/Pseudoknot_Prediction.md:231` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O526 | `docs/algorithms/RnaStructure/RNA_Base_Pairing.md:152` | https://en.wikipedia.org/wiki/Base_pair | LINK-OK |
| O527 | `docs/algorithms/RnaStructure/RNA_Base_Pairing.md:153` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O528 | `docs/algorithms/RnaStructure/RNA_Base_Pairing.md:154` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O529 | `docs/algorithms/RnaStructure/RNA_Stemloop.md:148` | https://en.wikipedia.org/wiki/Stem-loop | LINK-OK |
| O530 | `docs/algorithms/RnaStructure/RNA_Stemloop.md:149` | https://en.wikipedia.org/wiki/Tetraloop | LINK-OK |
| O531 | `docs/algorithms/RnaStructure/RNA_Stemloop.md:150` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O532 | `docs/algorithms/Sequence_Comparison/Common_Region_Detection.md:147` | https://en.wikipedia.org/wiki/Longest_common_substring | LINK-OK |
| O533 | `docs/algorithms/Sequence_Comparison/Common_Region_Detection.md:149` | https://en.wikipedia.org/wiki/Longest_common_substring | LINK-OK |
| O534 | `docs/algorithms/Sequence_Composition/Replication_Origin_Prediction.md:210` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O535 | `docs/algorithms/Sequence_Composition/Sequence_Composition.md:168` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O536 | `docs/algorithms/Sequence_Composition/Sequence_Composition_Statistics.md:147` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O537 | `docs/algorithms/Sequence_Composition/Sequence_Validation.md:194` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O538 | `docs/algorithms/Sequence_Composition/Shannon_Entropy.md:151` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O539 | `docs/algorithms/Sequence_Composition/Shannon_Entropy.md:152` | https://en.wikipedia.org/wiki/Sequence_logo | LINK-OK |
| O540 | `docs/algorithms/Sequence_Composition/Shannon_Entropy.md:153` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O541 | `docs/algorithms/Splicing/Donor_Site_Detection.md:167` | https://en.wikipedia.org/wiki/RNA_splicing | LINK-OK |
| O542 | `docs/algorithms/Splicing/Donor_Site_Detection.md:168` | https://en.wikipedia.org/wiki/Spliceosome | LINK-OK |
| O543 | `docs/algorithms/Statistics/Codon_Frequencies.md:163` | https://en.wikipedia.org/wiki/Codon_usage_bias | LINK-OK |
| O544 | `docs/algorithms/Statistics/DNA_Thermodynamics.md:262` | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | LINK-OK |
| O545 | `docs/algorithms/Statistics/Entropy_Profile.md:143` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O546 | `docs/algorithms/Statistics/GC_Content_Profile.md:213` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O547 | `docs/algorithms/Statistics/Secondary_Structure_Prediction.md:303` | https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method | LINK-OK |
| O548 | `docs/algorithms/Statistics/Sequence_Summary.md:193` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O549 | `docs/algorithms/Statistics/Sequence_Summary.md:194` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O550 | `docs/algorithms/Transcriptome/Differential_Expression.md:203` | https://en.wikipedia.org/wiki/False_discovery_rate | LINK-OK |
| O551 | `docs/algorithms/Transcriptome/Differential_Expression.md:204` | https://en.wikipedia.org/wiki/Welch%27s_t-test | LINK-OK |
| O552 | `docs/algorithms/Transcriptome/Differential_Expression.md:205` | https://en.wikipedia.org/wiki/Student%27s_t-distribution | LINK-OK |
| O553 | `docs/algorithms/Transcriptome/Expression_Quantification.md:184` | https://en.wikipedia.org/wiki/Quantile_normalization | LINK-OK |
| O554 | `docs/algorithms/Translation/Codon_Translation.md:184` | https://en.wikipedia.org/wiki/Genetic_code | LINK-OK |
| O555 | `docs/algorithms/Translation/Codon_Translation.md:185` | https://en.wikipedia.org/wiki/Start_codon | LINK-OK |
| O556 | `docs/algorithms/Translation/Codon_Translation.md:186` | https://en.wikipedia.org/wiki/Stop_codon | LINK-OK |
| O557 | `docs/algorithms/Translation/Six_Frame_Translation.md:288` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O558 | `docs/algorithms/Variants/SNP_Detection.md:174` | https://en.wikipedia.org/wiki/Transversion | LINK-OK |
| O559 | `docs/algorithms/Variants/Variant_Detection.md:209` | https://en.wikipedia.org/wiki/Transition_(genetics) | LINK-OK |
| O560 | `docs/algorithms/Variants/Variant_Detection.md:210` | https://en.wikipedia.org/wiki/Transversion | LINK-OK |
| O561 | `docs/mcp/tools/population/allele_frequencies.md:100` | https://en.wikipedia.org/wiki/Allele_frequency | LINK-OK |
| O562 | `docs/mcp/tools/population/allele_frequencies.md:100` | https://en.wikipedia.org/wiki/Genotype_frequency | LINK-OK |
| O563 | `docs/mcp/tools/population/f_statistics.md:115` | https://en.wikipedia.org/wiki/F-statistics | LINK-OK |
| O564 | `docs/mcp/tools/population/filter_variants_by_maf.md:93` | https://en.wikipedia.org/wiki/Minor_allele_frequency | LINK-OK |
| O565 | `docs/mcp/tools/population/fst.md:99` | https://en.wikipedia.org/wiki/Fixation_index | LINK-OK |
| O566 | `docs/mcp/tools/population/hardy_weinberg_test.md:112` | https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle | LINK-OK |
| O567 | `docs/mcp/tools/population/minor_allele_frequency.md:89` | https://en.wikipedia.org/wiki/Minor_allele_frequency | LINK-OK |
| O652 | `tests/TestSpecs/ALIGN-GLOBAL-001.md:18` | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | LINK-OK |
| O653 | `tests/TestSpecs/ALIGN-GLOBAL-001.md:19` | https://en.wikipedia.org/wiki/Sequence_alignment | LINK-OK |
| O654 | `tests/TestSpecs/ALIGN-LOCAL-001.md:18` | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | LINK-OK |
| O655 | `tests/TestSpecs/ALIGN-MULTI-001.md:17` | https://en.wikipedia.org/wiki/Multiple_sequence_alignment | LINK-OK |
| O656 | `tests/TestSpecs/ALIGN-MULTI-001.md:18` | https://en.wikipedia.org/wiki/Clustal | LINK-OK |
| O657 | `tests/TestSpecs/ALIGN-MULTI-001.md:19` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O658 | `tests/TestSpecs/ALIGN-MULTI-001.md:202` | https://en.wikipedia.org/wiki/Multiple_sequence_alignment | LINK-OK |
| O659 | `tests/TestSpecs/ALIGN-SEMI-001.md:18` | https://en.wikipedia.org/wiki/Sequence_alignment | LINK-OK |
| O660 | `tests/TestSpecs/ALIGN-SEMI-001.md:19` | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | LINK-OK |
| O661 | `tests/TestSpecs/ANNOT-GFF-001.md:25` | https://en.wikipedia.org/wiki/General_feature_format | LINK-OK |
| O662 | `tests/TestSpecs/ANNOT-PROM-001.md:21` | https://en.wikipedia.org/wiki/Promoter_(genetics) | LINK-OK |
| O663 | `tests/TestSpecs/ANNOT-PROM-001.md:22` | https://en.wikipedia.org/wiki/Pribnow_box | LINK-OK |
| O664 | `tests/TestSpecs/ANNOT-PROM-001.md:23` | https://en.wikipedia.org/wiki/TATA_box | LINK-OK |
| O665 | `tests/TestSpecs/ANNOT-REPEAT-001.md:18` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O666 | `tests/TestSpecs/ANNOT-REPEAT-001.md:19` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O667 | `tests/TestSpecs/ASSEMBLY-CONSENSUS-001.md:20` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O668 | `tests/TestSpecs/ASSEMBLY-SCAFFOLD-001.md:21` | https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics) | LINK-OK |
| O669 | `tests/TestSpecs/ASSEMBLY-STATS-001.md:19` | https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics | LINK-OK |
| O670 | `tests/TestSpecs/COMPGEN-COMPARE-001.md:21` | https://en.wikipedia.org/wiki/Synteny | LINK-OK |
| O671 | `tests/TestSpecs/COMPGEN-DOTPLOT-001.md:21` | https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) | LINK-OK |
| O672 | `tests/TestSpecs/COMPGEN-SYNTENY-001.md:20` | https://en.wikipedia.org/wiki/Synteny | LINK-OK |
| O673 | `tests/TestSpecs/CRISPR-GUIDE-001.md:17` | https://www.addgene.org/guides/crispr/ | LINK-OK |
| O674 | `tests/TestSpecs/CRISPR-GUIDE-001.md:18` | https://en.wikipedia.org/wiki/Guide_RNA | LINK-OK |
| O675 | `tests/TestSpecs/CRISPR-GUIDE-001.md:19` | https://en.wikipedia.org/wiki/Protospacer_adjacent_motif | LINK-OK |
| O676 | `tests/TestSpecs/CRISPR-OFF-001.md:28` | https://en.wikipedia.org/wiki/Off-target_genome_editing | LINK-OK |
| O677 | `tests/TestSpecs/DISORDER-MORF-001.md:21` | https://en.wikipedia.org/wiki/Molecular_recognition_feature | LINK-OK |
| O678 | `tests/TestSpecs/DISORDER-PRED-001.md:19` | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | LINK-OK |
| O679 | `tests/TestSpecs/DISORDER-PRED-001.md:20` | https://en.wikipedia.org/wiki/Hydrophilicity_plot | LINK-OK |
| O680 | `tests/TestSpecs/DISORDER-PRED-001.md:21` | https://en.wikipedia.org/wiki/Amino_acid | LINK-OK |
| O681 | `tests/TestSpecs/DISORDER-PROPENSITY-001.md:20` | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | LINK-OK |
| O682 | `tests/TestSpecs/DISORDER-REGION-001.md:22` | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | LINK-OK |
| O683 | `tests/TestSpecs/EPIGEN-CPG-001.md:21` | https://en.wikipedia.org/wiki/CpG_site | LINK-OK |
| O684 | `tests/TestSpecs/EPIGEN-DMR-001.md:22` | https://en.wikipedia.org/wiki/Fisher's_exact_test | LINK-OK |
| O685 | `tests/TestSpecs/GENOMIC-COMMON-001.md:18` | https://en.wikipedia.org/wiki/Longest_common_substring | LINK-OK |
| O686 | `tests/TestSpecs/GENOMIC-COMMON-001.md:19` | https://en.wikipedia.org/wiki/Longest_common_substring | LINK-OK |
| O687 | `tests/TestSpecs/GENOMIC-MOTIFS-001.md:20` | https://en.wikipedia.org/wiki/Restriction_site | LINK-OK |
| O688 | `tests/TestSpecs/GENOMIC-ORF-001.md:19` | https://en.wikipedia.org/wiki/Open_reading_frame | LINK-OK |
| O689 | `tests/TestSpecs/GENOMIC-REPEAT-001.md:19` | https://en.wikipedia.org/wiki/Longest_repeated_substring_problem | LINK-OK |
| O690 | `tests/TestSpecs/GENOMIC-SIMILARITY-001.md:18` | https://en.wikipedia.org/wiki/Jaccard_index | LINK-OK |
| O691 | `tests/TestSpecs/GENOMIC-TANDEM-001.md:19` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O692 | `tests/TestSpecs/KMER-ASYNC-001.md:18` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O693 | `tests/TestSpecs/KMER-COUNT-001.md:29` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O694 | `tests/TestSpecs/KMER-GENERATE-001.md:18` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O695 | `tests/TestSpecs/KMER-POSITIONS-001.md:19` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O696 | `tests/TestSpecs/KMER-STATS-001.md:18` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O697 | `tests/TestSpecs/KMER-UNIQUE-001.md:18` | https://en.wikipedia.org/wiki/K-mer | LINK-OK |
| O698 | `tests/TestSpecs/META-TAXA-001.md:18` | https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test | LINK-OK |
| O699 | `tests/TestSpecs/MIRNA-PAIR-001.md:21` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O700 | `tests/TestSpecs/MIRNA-PRECURSOR-001.md:22` | https://en.wikipedia.org/wiki/MicroRNA | LINK-OK |
| O701 | `tests/TestSpecs/MOTIF-CONS-001.md:18` | https://en.wikipedia.org/wiki/Consensus_sequence | LINK-OK |
| O702 | `tests/TestSpecs/MOTIF-GENERATE-001.md:20` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O703 | `tests/TestSpecs/ONCO-FUSION-001.md:23` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O704 | `tests/TestSpecs/ONCO-FUSION-003.md:20` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O705 | `tests/TestSpecs/ONCO-HETERO-001.md:22` | https://en.wikipedia.org/wiki/Diversity_index#Shannon_index | LINK-OK |
| O706 | `tests/TestSpecs/ONCO-IMMUNE-001.md:24` | https://en.wikipedia.org/wiki/CIBERSORT | LINK-OK |
| O707 | `tests/TestSpecs/ONCO-SIG-001.md:21` | https://en.wikipedia.org/wiki/Complementarity_(molecular_biology) | LINK-OK |
| O708 | `tests/TestSpecs/ONCO-SIG-002.md:20` | https://en.wikipedia.org/wiki/Non-negative_least_squares | LINK-OK |
| O709 | `tests/TestSpecs/ONCO-SIG-002.md:207` | https://en.wikipedia.org/wiki/Non-negative_matrix_factorization | LINK-OK |
| O710 | `tests/TestSpecs/ONCO-SIG-003.md:24` | https://en.wikipedia.org/wiki/Poisson_distribution#Random_variate_generation | LINK-OK |
| O711 | `tests/TestSpecs/ONCO-VAF-001.md:18` | https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval | LINK-OK |
| O712 | `tests/TestSpecs/PANGEN-CORE-001.md:23` | https://en.wikipedia.org/wiki/Pan-genome | LINK-OK |
| O713 | `tests/TestSpecs/PANGEN-MARKER-001.md:21` | https://en.wikipedia.org/wiki/Informative_site | LINK-OK |
| O714 | `tests/TestSpecs/PAT-APPROX-001.md:18` | https://en.wikipedia.org/wiki/Hamming_distance | LINK-OK |
| O715 | `tests/TestSpecs/PAT-APPROX-001.md:19` | https://en.wikipedia.org/wiki/Approximate_string_matching | LINK-OK |
| O716 | `tests/TestSpecs/PAT-EXACT-001.md:18` | https://en.wikipedia.org/wiki/Suffix_tree | LINK-OK |
| O717 | `tests/TestSpecs/PAT-IUPAC-001.md:18` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O718 | `tests/TestSpecs/PAT-PWM-001.md:18` | https://en.wikipedia.org/wiki/Position_weight_matrix | LINK-OK |
| O719 | `tests/TestSpecs/PHYLO-NEWICK-001.md:193` | https://en.wikipedia.org/wiki/Newick_format | LINK-OK |
| O720 | `tests/TestSpecs/PHYLO-STATS-001.md:18` | https://en.wikipedia.org/wiki/Tree_(graph_theory) | LINK-OK |
| O721 | `tests/TestSpecs/PHYLO-STATS-001.md:19` | https://en.wikipedia.org/wiki/Tree_(abstract_data_type) | LINK-OK |
| O722 | `tests/TestSpecs/PHYLO-STATS-001.md:22` | https://en.wikipedia.org/wiki/Minimum_evolution | LINK-OK |
| O723 | `tests/TestSpecs/PRIMER-DESIGN-001.md:19` | https://en.wikipedia.org/wiki/Primer_(molecular_biology) | LINK-OK |
| O724 | `tests/TestSpecs/PRIMER-DESIGN-001.md:20` | https://www.addgene.org/protocols/primer-design/ | LINK-OK |
| O725 | `tests/TestSpecs/PRIMER-STRUCT-001.md:33` | https://en.wikipedia.org/wiki/Primer_(molecular_biology) | LINK-OK |
| O726 | `tests/TestSpecs/PRIMER-STRUCT-001.md:34` | https://en.wikipedia.org/wiki/Primer_dimer | LINK-OK |
| O727 | `tests/TestSpecs/PRIMER-STRUCT-001.md:35` | https://en.wikipedia.org/wiki/Stem-loop | LINK-OK |
| O728 | `tests/TestSpecs/PRIMER-STRUCT-001.md:36` | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | LINK-OK |
| O729 | `tests/TestSpecs/PROTMOTIF-CC-001.md:20` | https://en.wikipedia.org/wiki/Coiled_coil | LINK-OK |
| O730 | `tests/TestSpecs/PROTMOTIF-CC-001.md:21` | https://en.wikipedia.org/wiki/Heptad_repeat | LINK-OK |
| O731 | `tests/TestSpecs/PROTMOTIF-FIND-001.md:22` | https://en.wikipedia.org/wiki/Sequence_motif | LINK-OK |
| O732 | `tests/TestSpecs/PROTMOTIF-FIND-001.md:24` | https://en.wikipedia.org/wiki/Nuclear_localization_sequence | LINK-OK |
| O733 | `tests/TestSpecs/QUALITY-STATS-001.md:18` | https://en.wikipedia.org/wiki/Phred_quality_score | LINK-OK |
| O734 | `tests/TestSpecs/REP-DIRECT-001.md:36` | https://en.wikipedia.org/wiki/Direct_repeat | LINK-OK |
| O735 | `tests/TestSpecs/REP-DIRECT-001.md:37` | https://en.wikipedia.org/wiki/Repeated_sequence_(DNA) | LINK-OK |
| O736 | `tests/TestSpecs/REP-INV-001.md:31` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O737 | `tests/TestSpecs/REP-INV-001.md:32` | https://en.wikipedia.org/wiki/Stem-loop | LINK-OK |
| O738 | `tests/TestSpecs/REP-INV-001.md:33` | https://en.wikipedia.org/wiki/Palindromic_sequence | LINK-OK |
| O739 | `tests/TestSpecs/REP-PALIN-001.md:30` | https://en.wikipedia.org/wiki/Palindromic_sequence | LINK-OK |
| O740 | `tests/TestSpecs/REP-PALIN-001.md:31` | https://en.wikipedia.org/wiki/Restriction_enzyme | LINK-OK |
| O741 | `tests/TestSpecs/REP-TANDEM-001.md:33` | https://en.wikipedia.org/wiki/Tandem_repeat | LINK-OK |
| O742 | `tests/TestSpecs/REP-TANDEM-001.md:34` | https://en.wikipedia.org/wiki/Microsatellite | LINK-OK |
| O743 | `tests/TestSpecs/RESTR-FILTER-001.md:18` | https://en.wikipedia.org/wiki/Sticky_and_blunt_ends | LINK-OK |
| O744 | `tests/TestSpecs/RESTR-FILTER-001.md:19` | https://en.wikipedia.org/wiki/Restriction_enzyme | LINK-OK |
| O745 | `tests/TestSpecs/RESTR-FILTER-001.md:20` | https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites | LINK-OK |
| O746 | `tests/TestSpecs/RNA-INVERT-001.md:19` | https://en.wikipedia.org/wiki/Inverted_repeat | LINK-OK |
| O747 | `tests/TestSpecs/RNA-PAIR-001.md:19` | https://en.wikipedia.org/wiki/Base_pair | LINK-OK |
| O748 | `tests/TestSpecs/RNA-PAIR-001.md:20` | https://en.wikipedia.org/wiki/Wobble_base_pair | LINK-OK |
| O749 | `tests/TestSpecs/RNA-PAIR-001.md:21` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O750 | `tests/TestSpecs/RNA-PKPREDICT-001.md:20` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O751 | `tests/TestSpecs/RNA-PSEUDOKNOT-001.md:21` | https://en.wikipedia.org/wiki/Pseudoknot | LINK-OK |
| O752 | `tests/TestSpecs/SEQ-ATSKEW-001.md:20` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O753 | `tests/TestSpecs/SEQ-CODON-FREQ-001.md:21` | https://en.wikipedia.org/wiki/Codon_usage_bias | LINK-OK |
| O754 | `tests/TestSpecs/SEQ-COMP-001.md:46` | https://en.wikipedia.org/wiki/Complementarity_(molecular_biology) | LINK-OK |
| O756 | `tests/TestSpecs/SEQ-COMPLEX-COMPRESS-001.md:19` | https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity | LINK-OK |
| O757 | `tests/TestSpecs/SEQ-COMPLEX-KMER-001.md:20` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O758 | `tests/TestSpecs/SEQ-COMPLEX-WINDOW-001.md:18` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O759 | `tests/TestSpecs/SEQ-COMPLEX-WINDOW-001.md:21` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O760 | `tests/TestSpecs/SEQ-COMPOSITION-001.md:19` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O761 | `tests/TestSpecs/SEQ-ENTROPY-PROFILE-001.md:19` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O762 | `tests/TestSpecs/SEQ-GC-001.md:18` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O763 | `tests/TestSpecs/SEQ-GC-ANALYSIS-001.md:18` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O764 | `tests/TestSpecs/SEQ-GC-ANALYSIS-001.md:19` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O765 | `tests/TestSpecs/SEQ-GC-PROFILE-001.md:18` | https://en.wikipedia.org/wiki/GC-content | LINK-OK |
| O766 | `tests/TestSpecs/SEQ-GCSKEW-001.md:33` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O767 | `tests/TestSpecs/SEQ-REPLICATION-001.md:21` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O768 | `tests/TestSpecs/SEQ-SECSTRUCT-001.md:20` | https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method | LINK-OK |
| O769 | `tests/TestSpecs/SEQ-STATS-001.md:20` | https://en.wikipedia.org/wiki/GC_skew | LINK-OK |
| O770 | `tests/TestSpecs/SEQ-SUMMARY-001.md:20` | https://en.wikipedia.org/wiki/Entropy_(information_theory) | LINK-OK |
| O771 | `tests/TestSpecs/SEQ-SUMMARY-001.md:21` | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | LINK-OK |
| O772 | `tests/TestSpecs/SEQ-THERMO-001.md:22` | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | LINK-OK |
| O773 | `tests/TestSpecs/SEQ-VALID-001.md:18` | https://en.wikipedia.org/wiki/Nucleic_acid_notation | LINK-OK |
| O774 | `tests/TestSpecs/SPLICE-ACCEPTOR-001.md:27` | https://en.wikipedia.org/wiki/RNA_splicing | LINK-OK |
| O775 | `tests/TestSpecs/SPLICE-ACCEPTOR-001.md:28` | https://en.wikipedia.org/wiki/Polypyrimidine_tract | LINK-OK |
| O776 | `tests/TestSpecs/SPLICE-ACCEPTOR-001.md:29` | https://en.wikipedia.org/wiki/Minor_spliceosome | LINK-OK |
| O777 | `tests/TestSpecs/SPLICE-DONOR-001.md:21` | https://en.wikipedia.org/wiki/RNA_splicing | LINK-OK |
| O778 | `tests/TestSpecs/SPLICE-DONOR-001.md:22` | https://en.wikipedia.org/wiki/Spliceosome | LINK-OK |
| O779 | `tests/TestSpecs/SPLICE-PREDICT-001.md:22` | https://en.wikipedia.org/wiki/Intron | LINK-OK |
| O780 | `tests/TestSpecs/SPLICE-PREDICT-001.md:23` | https://en.wikipedia.org/wiki/Exon | LINK-OK |
| O781 | `tests/TestSpecs/TRANS-DIFF-001.md:19` | https://en.wikipedia.org/wiki/False_discovery_rate | LINK-OK |
| O782 | `tests/TestSpecs/TRANS-DIFF-001.md:20` | https://en.wikipedia.org/wiki/Welch%27s_t-test | LINK-OK |
| O783 | `tests/TestSpecs/TRANS-DIFF-001.md:21` | https://en.wikipedia.org/wiki/Student%27s_t-distribution | LINK-OK |
| O784 | `tests/TestSpecs/TRANS-EXPR-001.md:21` | https://en.wikipedia.org/wiki/Quantile_normalization | LINK-OK |
| O785 | `tests/TestSpecs/TRANS-SIXFRAME-001.md:22` | https://en.wikipedia.org/wiki/Reading_frame | LINK-OK |
| O786 | `tests/TestSpecs/VARIANT-CALL-001.md:22` | https://en.wikipedia.org/wiki/Transition_(genetics) | LINK-OK |
| O787 | `tests/TestSpecs/VARIANT-CALL-001.md:23` | https://en.wikipedia.org/wiki/Transversion | LINK-OK |
| O788 | `tests/TestSpecs/VARIANT-SNP-001.md:19` | https://en.wikipedia.org/wiki/Transversion | LINK-OK |
