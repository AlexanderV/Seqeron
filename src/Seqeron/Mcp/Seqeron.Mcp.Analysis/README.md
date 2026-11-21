# Seqeron.Mcp.Analysis

MCP server — **K-mer, motif, repeat, complexity, RNA-structure and comparative-genomics analysis.**

Exposes **173 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/analysis/`](../../../../docs/mcp/tools/analysis). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Analysis
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Analysis"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (132)

| Tool | Description |
|------|-------------|
| `alignment_distance_matrix` | Pairwise distance matrix of a multiple alignment as EMBOSS 6.6.0 distmat: matches = columns where both rows carry the same non-gap symbol, g = column… |
| `alphabet_pwm_score_pvalue` | Exact p-value of a score / exact score threshold of a p-value for a PWM over any alphabet (protein, RNA, …) built from aligned instances — the TFM-Pvalue engine of pwm_score_pvalue with K rows… |
| `alphabet_pwm_score_thresholds` | Score thresholds (FPR / FNR / balanced / patser) of an any-alphabet PWM from its discretised score distribution — Biopython pssm.distribution… |
| `analyze_gc_content` | Comprehensive GC report: overall GC content, GC/AT skew, content/skew variances, and windowed GC profiles (`fraction=true` for [0,1]; `ambiguity=remove\|ignore\|weighted` for Biopython `gc_fraction` IUPAC handling). |
| `analyze_kmers` | Aggregate k-mer statistics (Jellyfish stats fields): total, distinct (uniqueKmers/distinctKmers), singleton (count-1) k-mers, min/max/mean count, and Shannon entropy; optional Jellyfish -L/-U count filters and canonical (count -C) / acgtOnly modes. |
| `at_skew` | Whole-sequence AT skew = (A - T) / (A + T). |
| `base_pair_type` | Classification of an RNA base-pair candidate: WatsonCrick, Wobble, or null if bases cannot pair. |
| `bulge_loop_energy` | Free energy of an RNA bulge loop (special-C bonus, n=1 stacking, degeneracy entropy from numStates). |
| `calculate_ani` | Average Nucleotide Identity (ANI) between two genome sequences (fragment-and-match). |
| `can_pair` | Whether two RNA bases form a Watson-Crick (A-U, G-C) or wobble (G-U) pair. |
| `chou_fasman_assign` | Discrete Chou-Fasman (1978) per-residue H/E/T/C secondary-structure assignment; `ruleSet` Standard (default) or Graded1978. |
| `codon_frequencies` | Codon usage frequencies in the specified reading frame (0..2). |
| `compare_genomes` | End-to-end comparative pipeline: RBH orthologs + synteny + rearrangements + summary stats. |
| `compression_ratio` | Estimate sequence repetitiveness as the normalized Lempel-Ziv complexity c/(n/log_b(n)). |
| `consambig` | EMBOSS consambig: per-column union of the residue codes written as one IUPAC / amino-acid letter (lower case at gaps); 1 400 random alignments equal consambig 6.6.0 (B09 A3-13d). |
| `count_kmers` | Counts every k-mer (substring of length k) occurrence in a sequence; optional `canonical` (jellyfish -C) and `acgtOnly` modes. |
| `count_kmers_both_strands` | k-mer counts on the forward strand combined with counts on the reverse-complement strand (kPAL balance); optional `acgtOnly` (kPAL: N/IUPAC accepted, non-ACGT windows skipped). |
| `create_alphabet_pwm` | Build a log-odds position weight matrix over an arbitrary alphabet (protein, RNA, gapped DNA, …) from aligned instances — Biopython motifs.create(instances, alphabet).counts.normalize(pseudocounts).log_odds(background)… |
| `create_pwm` | Build a log-odds Position Weight Matrix (4×L; rows A,C,G,T) from aligned, equal-length DNA sequences (Biopython counts.normalize(pseudocounts).log_odds(background)): scalar pseudocount, per-base pseudocounts, or JASPAR pseudocounts… |
| `cumulative_gc_skew` | Cumulative GC skew along the sequence — minimum approximates origin and maximum approximates terminus of replication. |
| `dangling_end_energy` | 5' or 3' dangling-end stacking energy for an RNA helix end. |
| `detect_pseudoknots` | Identifies pseudoknots as crossing base pairs: i < i' < j < j' for pairs (i,j) and (i',j'). |
| `detect_rearrangements` | Inversions, deletions, insertions inferred from gene-order disagreements between two genomes. |
| `dinucleotide_frequencies` | Frequency of each adjacent dinucleotide over the alphabet A/T/G/C/U. |
| `dinucleotide_ratios` | Karlin dinucleotide odds ratio ρ (single strand) or strand-symmetrized ρ* (`strandSymmetric=true`), optional Karlin over/under-representation classes (`includeRepresentation`). Not the CpG-island O/E. |
| `discover_motifs` | Overrepresented k-mers (de novo motif discovery) in a DNA sequence. |
| `disorder_propensity` | Returns the TOP-IDP propensity value for a single amino acid (Campen 2008). |
| `dust_score` | DUST low-complexity score (BLAST-style, triplet-based) for a DNA sequence. |
| `dyad_analysis` | RSAT dyad-analysis: over-represented spaced dyads M1 n{s} M2 (any / dr / ir / rep) with monad-based expected frequencies, z-scores and binomial occ_P / occ_E / occ_sig… |
| `entropy_profile` | Shannon entropy in sliding windows along the sequence; optional k-mer (block) entropy (`k`) and trailing partial windows (`includePartialWindow`). |
| `find_all_longest_repeats` | Every longest repeated substring of a DNA sequence (all length ties) with all positions; `bothStrands` gives the longest repeat classes over both strands (MUMmer repeat-match); IUPAC input accepted (B09 A2-36 / A2-33). |
| `find_approximate_direct_repeats` | Finds **all maximal k-mismatch (Hamming) direct repeats** (Kurtz et al. 2000 REPuter; Vmatch -h k): two copies of equal length length at 0-… |
| `find_approximate_tandem_repeats` | Finds approximate (imperfect / interrupted) tandem repeats with the **Tandem Repeats Finder** model (Benson 1999; TRF 4.10.0 trf File Match… |
| `find_clump_windows` | (L, t)-clump k-mers with the maximal runs of window starts in which they form clumps (same k-mer set as `find_clumps`). |
| `find_clumps` | Finds k-mers that occur at least minOccurrences times within any sliding window of size windowSize. |
| `find_common_regions` | Right-maximal common regions of two DNA sequences: for each start in sequence2 the single longest substring (length >= minLength) that also occurs in sequence1, distinct strings once, with the first occurrence in sequence1. Not every common substring. `mode` (RightMaximal / MaximalExactMatches / MaximalUniqueReference / MaximalUnique = mummer -maxmatch / -mumreference / -mum) and `strands` (mummer -r / -b) (B09 A2-30); IUPAC input accepted. |
| `find_conserved_clusters` | Gene clusters preserved across multiple genomes. |
| `find_degenerate_motif` | Motif search with IUPAC ambiguity codes (N, R, Y, S, W, K, M, B, D, H, V). |
| `find_degenerate_repeats` | Finds **all maximal degenerate repeats** with at most maxDifferences differences (Kurtz et al. 2000 REPuter; Vmatch -e k / -h k, -p, -allma… |
| `find_direct_repeats` | Identical sequences appearing twice with a spacer between them. |
| `find_exact_motif` | Exact-match motif positions in a DNA sequence via suffix tree. |
| `find_inverted_repeats` | Sequences whose two arms are reverse-complement of each other (hairpin candidates). |
| `find_inverted_repeats_scored` | Finds inverted repeats with EMBOSS einverted (Durbin & Thierry-Mieg 1993; EMBOSS 6.6.0): a local alignment of the sequence against its reve… |
| `find_known_motifs` | Search for a user-provided set of motifs simultaneously (exact by default; optional IUPAC codes, reverse strand and a mismatch budget). Only motifs with at least one hit are returned; positions ascending, 0-based, forward-strand window starts. mismatches = k > 0 allows up to k substituted positions per window (Hamming distance, EMBOSS fuzznuc -pmismatch; no gaps) and adds matchMismatches / reverseMismatches (the distance of each reported hit, parallel to the position arrays). IUPAC sequences accepted (B09 A2-35). |
| `mash_bounds` | Mash error-bounds table (`mash bounds`): error of the Mash / screen distance for sketch sizes 100 - 1 000 000 and distances 0.05 - 0.4. |
| `mash_dist` | Mash distances of query files to reference files (`mash dist`: -t, -v, -d, -i, -n). |
| `mash_triangle` | Mash lower-triangular distance matrix or edge list (`mash triangle`: -E, -v, -d, -i, -n). |
| `mash_screen` | Mash containment of query sketches in a mixture of reads / contigs (`mash screen`: identity, shared-hashes, median-multiplicity, p-value, query-ID, comment; -w, -i, -v). |
| `mash_sketch` | Mash sketch file (`mash sketch -o x.msh`: Cap'n Proto `MinHash.capnp`, byte for byte; genomes and read sets, -k -s -S -n -i -r -m -c -g -I -C -l): `msh:` + base64, accepted by the other Mash tools. |
| `mash_info` | `mash info` of a sketch file: header and table, -H, -t, -c (count histogram), -d (JSON). |
| `diffseq` | EMBOSS diffseq: the regions of difference between two nearly identical sequences (word matches), the 7-tag report in every -rformat and the GFF3 feature files. |
| `wordmatch` | EMBOSS wordmatch: exact matches of a pattern set against large sequences (Rabin-Karp), alignments, statistics log and GFF3 feature files. |
| `mash_paste` | `mash paste`: one sketch file of several (first file's parameters; incompatible files skipped, larger sketches reduced). |
| `emboss_report_file` | The whole report file of fuzznuc / fuzzpro / fuzztran / etandem / equicktandem in every -rformat (simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq), with -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq, several sequences per file, -sbegin / -send / -sreverse and the -origfile of the tandem programs. |
| `megamerger` | EMBOSS megamerger: merge two large overlapping DNA sequences by minimal non-overlapping word-match anchors (merged sequence + report); 600 random pairs equal megamerger 6.6.0 (B09 A3-13b). |
| `nucmer_clusters` | MUMmer nucmer anchor clustering (`mummer -mum|-mumreference|-maxmatch | mgaps`: unique / maximal matches joined by diagonal difference and separation, cut into maximum-score chains); 700 related pairs equal the installed nucmer 3.23 pipeline (B09 A3-9a). |
| `nucmer` | MUMmer 3.23 nucmer: the whole pipeline (anchors, mgaps clusters, postnuc banded Smith-Waterman extension), the .delta file byte for byte (434 runs, all option families). |
| `delta_filter` | MUMmer 3.23 delta-filter: keeps the alignments of a delta file that survive -i -l -u -q -r -g -m -1 (rand() seeded with 1, as the program). |
| `show_coords` | MUMmer 3.23 show-coords: the alignment coordinates of a delta file (table, tab-delimited, btab; -b -c -d -H -I -k -l -L -o -q -r -T -B -g -G -w). |
| `show_aligns` | MUMmer 3.23 show-aligns: the alignments of one reference / query pair drawn as three rows (^ under a difference, . for an indel). |
| `show_snps` | MUMmer 3.23 show-snps: substitutions and indels of a nucmer delta file with buffer, distance, conflict counts and context. |
| `show_diff` | MUMmer 3.23 show-diff: breaks, gaps, jumps, inversions, duplications and sequence jumps read off the many-to-many mapping of a delta file. |
| `show_tiling` | MUMmer 3.23 show-tiling: tiles the aligned query contigs on the reference (path, alignments, XML, contig list, unusable contigs, pseudo molecule). |
| `dnadiff` | MUMmer 3.23 dnadiff: nucmer --maxmatch + delta-filter + show-coords / show-snps / show-diff and the report of two sequence sets. |
| `mummer` | MUMmer mummer: maximal exact matches (-mum / -mumreference / -maxmatch, -n -l -b -r -s -c -F -L), byte for byte. |
| `promer` | MUMmer promer: six-frame amino-acid alignments (BLOSUM 45 / 62 / 80) of two DNA files, delta file byte for byte. |
| `mummer_gaps` | MUMmer gaps: longest consistent set of matches and the others (mgaps format), byte for byte. |
| `delta2blocks` | MUMmer delta2blocks: alignments of a nucmer delta file as CIGAR lines (M / X / I / D). |
| `mummer_annotate` | MUMmer mummer-annotate: edit-distance alignment of the gaps of a gaps file with the error counts. |
| `run_mummer1` | MUMmer run-mummer1: mummer -mum -l 20 -r | gaps | mummer-annotate script, files .out .gaps .align .errorsgaps. |
| `run_mummer3` | MUMmer run-mummer3: mummer -mumreference -b -l 20 | mgaps | combineMUMs script, files .out .gaps .align .errorsgaps. |
| `repeat_match` | MUMmer repeat-match: maximal exact repeats (-n -f -t -E -V) in the order of the program's suffix tree, byte for byte. |
| `exact_tandems` | MUMmer exact-tandems: exact tandem repeats (repeat-match -t | sort | awk), the intended pipeline and the script as installed. |
| `mapview` | MUMmer mapview: FIG maps of show-coords -rl / mgaps alignments with genes (-d -f -m -n -p -v -x1 -x2 -g -I -Ir), byte for byte. |
| `mummerplot` | MUMmer mummerplot: gnuplot script and data files of delta / cluster / mummer / tiling alignments (every option), byte for byte. |
| `delta2maf` | MUMmer delta2maf: nucmer alignments as MAF (a score / s lines), byte for byte. |
| `nucmer2xfig` | MUMmer nucmer2xfig: xfig plot of a show-coords -c -l file, byte for byte. |
| `combine_mums` | MUMmer combineMUMs: extension and alignment of the mgaps clusters (-D -e -n -N -q -r -S -t -v -x -W), byte for byte. |
| `est2genome` | EMBOSS 6.6.0 est2genome: the spliced alignment of ESTs to a genomic sequence (exons, GT/AG and CT/AC introns, both strands, linear-space recursion; 330 runs equal). |
| `seqmatchall` | EMBOSS seqmatchall: every maximal exact word match (>= `wordSize`) of every pair of a sequence set, nucleotide or protein, including the program's N / X word-table scan; verified against seqmatchall 6.6.0 on 330 sets (B09 A3-13). |
| `find_longdust_regions` | Finds low-complexity regions with **longdust** (Li & Li 2025; port of lh3/longdust 1.4-r97 ld_dust1/ld_dust2, MIT) — the k-mer generalisati… |
| `find_longest_common_region_mismatches` | Longest common substring with at most k mismatches (k-mismatch LCS) of two DNA sequences, optionally against the reverse complement (B09 A2-30). |
| `find_low_complexity_intervals` | Returns the **SDUST** low-complexity intervals (Morgulis et al. 2006; symmetric DUST) as 0-based half-open [start, end) pairs in ascending… |
| `find_low_complexity_regions` | Entropy-thresholded contiguous low-complexity DNA regions. |
| `find_microsatellites` | Short Tandem Repeats (STRs): 1-6 bp motif units repeated consecutively. |
| `find_motif` | All exact occurrences of a motif (case-insensitive, overlapping, ascending 0-based positions, forward strand) in a DNA sequence. IUPAC sequences accepted (non-A/C/G/T never matches) (B09 A2-35). |
| `find_motif_by_pattern` | Find all (overlapping) regex pattern matches in a protein sequence. |
| `find_motif_by_prosite` | Convert a PROSITE pattern to regex internally, then scan a protein sequence. |
| `find_open_reading_frames` | ORFs in all 6 frames (3 forward + 3 reverse-complement) from a start codon (default ATG) to the first in-frame stop of the chosen genetic code (default Standard TAA/TAG/TGA); nested ORFs, optional alternative starts / partial ORFs / NCBI ORFfinder 0.4.3 options measured on the binary (`ignoreNested` -n, `regionStart/End` -b/-e, `strand`, `orfFinderMinLength` -ml, tables 27/28/31) and length convention; coordinates in the scanned strand and on the input; IUPAC input (N runs) accepted — an ambiguous codon is a stop only when every concrete codon is, never a start. |
| `find_orthologs` | Best-hit ortholog pairs between two genomes by k-mer similarity (one-directional). |
| `find_palindromes` | Sequences identical to their reverse complement (restriction-site candidates); only A/C/G/T pair, windows with N / IUPAC codes are never reported. |
| `find_pattern_hits` | EMBOSS fuzznuc pattern language: IUPAC letters, classes `[ACG]`, exclusions `{AG}`, repeats `N(3)`, ranges `A(2,4)`, anchors `<` `>`; `maxMismatches` (-pmismatch), `bothStrands` (-complement), one pattern or a pattern file; verified against fuzznuc 6.6.0 on 3300 random patterns (B09 A2-32). |
| `fuzzpro` | EMBOSS fuzzpro: exact protein pattern search (letters, `x`, `[ILV]`, `{P}`, `x(3)`, `x(2,4)`, `<` / `>`; B / Z / J ambiguity letters); verified against fuzzpro 6.6.0 on 900 random patterns (B09 A3-12). |
| `fuzztran` | EMBOSS fuzztran: protein pattern search after translation (`frame` 1 2 3 -1 -2 -3 F R 6, EMBOSS `table`; transeq frame conventions); forward-strand location (and the program's own coordinates for the forward frames); verified against fuzztran 6.6.0 on 760 random cases (B09 A3-12). |
| `find_promoter_elements_by_matrix` | Scan a DNA sequence with the Bucher (1990) promoter weight matrices from JASPAR POLII (TATA box POL012.1, cap/Inr POL002.1, CCAAT box POL00… |
| `find_protein_domains` | Detect common protein domains with EXACT PROSITE patterns: zinc finger C2H2 (PS00028), WD-repeats (PS00678), kinase ATP-binding / Walker… |
| `find_protein_low_complexity_regions` | Low-complexity regions in a protein via the SEG algorithm (Wootton & Federhen 1993): sliding-window Shannon entropy in bits/residue, two-… |
| `find_protein_motifs` | Scan a protein sequence against the built-in PROSITE-style motif catalog (N-glycosylation, kinase phosphorylation sites, ATP/GTP P-loop,… |
| `find_reciprocal_best_hits` | Reciprocal best hits (RBH) for stricter ortholog identification. |
| `find_regulatory_elements` | Scan for built-in regulatory motifs (TATA, CAAT, GC-box, Kozak, Shine-Dalgarno, poly(A), E-box, AP-1, NF-κB, CREB). |
| `find_regulatory_elements_both_strands` | Strand-annotated scan of the built-in regulatory library: every element on the given strand, plus the orientation-independent elements (CAA… |
| `find_repeats` | Repeated substrings of a DNA sequence with all 0-based positions (overlaps counted): every distinct substring of length >= minLength occurring at least twice (default), or only Gusfield maximal repeats (`maximalOnly`). `bothStrands` reports repeat classes over both strands (B09 A2-33); IUPAC input accepted (B09 A2-35). |
| `find_reverse_complement_repeats` | Finds **maximal exact reverse-complement repeat pairs**: copy 1 at firstPosition and copy 2 at secondPosition (both 0-based, firstPosition… |
| `find_rna_inverted_repeats` | Finds antiparallel complementary regions (potential RNA hairpin stems). |
| `find_shared_motifs` | k-mers present in at least minSequences of the input DNA sequences. |
| `find_sigma70_promoters` | Pair bacterial σ70 −35 and −10 boxes (Harley & Reynolds 1987: TTGACA / TATAAT, spacer 15–21 bp, 17 ± 1 in 92 % of promoters)… |
| `find_stem_loops` | Enumerates hairpin stem-loop candidates with stem, loop, and Turner 2004 free energy. |
| `find_supermaximal_repeats` | Finds **supermaximal repeats** (Gusfield 1997 §7.12.1; Vmatch -supermax): maximal repeats that never occur as a substring of another maxima… |
| `find_syntenic_blocks` | Collinear runs of orthologous genes between two genomes. |
| `find_tandem_repeats` | Exact tandem repeats (e.g., ATGATGATG): maximal runs of a primitive unit (minUnitLength..maxUnitLength) repeated >= minRepetitions times; IUPAC input accepted (units must be A/C/G/T only); approximate repeats: `find_approximate_tandem_repeats`. `canonicalRotation` reports the least rotation of each unit (B09 A2-34). |
| `find_tandem_repeats_emboss` | EMBOSS etandem / equicktandem ports (diagonal-score models: running consensus score / shifted-diagonal score), `program`, `minRepeat`, `maxRepeat`, `threshold`, `mismatch`, `uniform`; verified against EMBOSS 6.6.0 on 1700 + 1300 sequences (B09 A3-8). |
| `flush_coaxial_stacking` | Coaxial stacking energy for two RNA helices with no intervening unpaired bases. |
| `gc_content_profile` | GC content in sliding windows along the sequence (`fraction=true` for [0,1]; `ambiguityMode=remove\|ignore\|weighted` for Biopython gc_fraction; `includePartialWindow`). |
| `gc_skew` | Whole-sequence GC skew = (G - C) / (G + C). |
| `generate_all_kmers` | Enumerate the entire k-mer space for an alphabet (default "ACGT"); result size alphabet.Length^k, at most 1,048,576 (4^10) k-mers per call. |
| `generate_cavener_consensus` | Degenerate IUPAC consensus of aligned equal-length DNA sequences by the Cavener (1987) rules (TRANSFAC / Biopython degenerate_consensus): s… |
| `generate_consensus` | IUPAC consensus sequence from aligned equal-length DNA sequences: bases above the per-position inclusion threshold (default >25%) form the NC-IUB symbol; if none passes, the tied most frequent bases. |
| `generate_decipher_consensus` | Consensus with Bioconductor DECIPHER ConsensusSequence semantics (DNA, RNA or protein): at each position the least frequent characters are dropped while they represent less than 'threshold'… |
| `generate_dot_plot` | Coordinates of matching k-mers between two sequences for dot-plot visualization. |
| `generate_dumb_consensus` | Majority-threshold consensus with Biopython SummaryInfo.dumb_consensus semantics: per column the most frequent non-gap residue is emitted i… |
| `generate_emboss_consensus` | Scoring-matrix plurality consensus of an alignment, identical to EMBOSS 6.6.0 'cons' (EDNAFULL for nucleotides, EBLOSUM62 for proteins; N/X where no residue reaches the plurality; residueType 'auto', padRaggedRows)… |
| `hairpin_loop_energy` | Free energy of an RNA hairpin loop (Turner 2004 with special tri/tetra/hexaloops, terminal mismatch, all-C and special-GU adjustments). |
| `hydrophobicity_profile` | Sliding-window Kyte-Doolittle hydropathy (ProtScale / Biopython protein_scale): optional `edgeWeight`, `convention` (Standard\|Biopython), `strict`, `normalizeScale`, `scale`, `weightModel` (Linear\|Exponential). |
| `internal_loop_energy` | Free energy of a generic RNA internal loop (Turner 2004; |
| `is_disorder_promoting` | Whether an amino acid is in Dunker's disorder-promoting set {A, R, G, Q, S, P, E, K} (Dunker 2001). |
| `kmer_d2_statistics` | Background-adjusted D2* / D2S statistics, d2* / d2S dissimilarities, Markov orders used and per-sequence BIC values; optional `markovOrder`, `bothStrands` (CAFE -R). |
| `kmer_distance` | Euclidean distance between k-mer frequency vectors of two sequences; optional `metric` (squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd, jensen_shannon, euclidean_counts, ev = spaced -d EV with pattern 1^k), `markovOrder` (background of d2star/d2shepherd) and `bothStrands` (CAFE -R d2star/d2shepherd; spaced reverse-complement mode for ev). |
| `kmer_frequencies` | Normalized k-mer counts (each value in [0,1], summing to 1); optional `canonical` (jellyfish count -C) / `acgtOnly` (kPAL/Jellyfish). |
| `kmer_jaccard` | k-mer Jaccard index of the distinct k-mer sets (exact, or Mash MinHash estimate with `sketchSize` > 0, then also Mash x/s and p-value), the Mash distance and the exact containment indices; optional `canonical` / `acgtOnly` (Mash/sourmash k-mers); optional `scaled` (1..4294967295, with `canonical` = true) = sourmash FracMinHash estimates (Jaccard, `contained_by`, `maxContainment`), plus `angularSimilarity` and the weighted containments (sourmash `contained_by_weighted`) with `trackAbundance`. |
| `kmer_positions` | Zero-based positions of all (overlapping) occurrences of a k-mer. |
| `kmer_spectrum` | Frequency-of-frequencies: for each occurrence count, how many distinct k-mers reach that count; optional `canonical`/`acgtOnly` and jellyfish histo `low`/`high`/`increment`/`full` (adds `histogram`). |
| `kmers_with_min_count` | k-mers occurring at least minCount (and, if given, at most maxCount) times — jellyfish dump -L/-U — sorted by count descending, ties by k-mer; optional `canonical`/`acgtOnly`. |
| `lempel_ziv_complexity` | Returns the raw **Lempel–Ziv (1976) complexity** c — the number of components of the exhaustive history of the sequence (Lempel & Ziv 1976;… |
| `longdust_score` | Computes the **longdust** score of a whole sequence x (Li & Li 2025, arXiv:2509.07357; lh3/longdust): S_L(x) = Σ_t log c_x(t)! − f(ℓ(x)/4^k… |
| `mask_approximate_tandem_repeats` | Returns TRF's masked sequence (trf ... -m; TRF README: "every location that occurred in a tandem repeat changed to the letter 'N'"): every… |
| `mask_low_complexity` | Mask low-complexity windows (DUST-driven) of a DNA sequence with a chosen character. |
| `minimum_free_energy` | Zuker-style minimum free energy with Turner 2004 parameters (O(n³)). |
| `mismatch_coaxial_stacking` | Mismatch-mediated coaxial stacking energy: terminal mismatch + base + WC/GU bonus. |
| `most_frequent_kmers` | Returns all k-mers tied for the maximum occurrence count; optional `canonical` (jellyfish count -C) / `acgtOnly`. |
| `multibranch_loop_energy` | Free energy of an RNA multibranch loop (Turner 2004 affine model: offset + asymmetry + helix term + stacking + strain). |
| `multiple_alignment_statistics` | Statistics of a gapped multiple alignment: per-sequence comparison with a reference row or the EMBOSS cons consensus (EMBOSS infoalign: seqLength, al… |
| `oligo_analysis` | RSAT oligo-analysis over-representation of the k-mers of one DNA sequence (or a set with extraSequences): occurrences, expected frequency under a background model (equipr… |
| `parse_dot_bracket` | Parses dot-bracket notation into a list of base-pair coordinates. |
| `predict_chou_fasman` | Propensity profile: per-window mean Chou-Fasman helix/sheet/turn parameters (window 1 = per residue); not an assignment (see `chou_fasman_assign`). |
| `predict_coiled_coils` | Heptad-repeat-based coiled-coil prediction. |
| `predict_disorder` | TOP-IDP disorder prediction (Campen 2008): per-residue scores plus contiguous IDR regions with confidence and subtype classification. |
| `predict_low_complexity_seg` | SEG algorithm (Wootton & Federhen 1993/1996) for low-complexity protein regions: trigger window K1 + extension K2. |
| `predict_morfs` | Predicts Molecular Recognition Features within IDRs by hydropathy enrichment (heuristic, Mohan 2006-inspired). |
| `predict_replication_origin` | Predicts replication origin and terminus from cumulative GC skew extrema (optional `circular`; `windowSize` = Grigoriev windowed diagram; `skewIndexWindow` adds the SkewIT Skew Index). |
| `predict_rna_structure` | Greedy non-overlapping stem-loop selection — produces dot-bracket notation, base pairs, stems, pseudoknots, and total MFE. |
| `predict_sigma70_promoters` | σ70 promoter prediction with the Promoter Calculator v1.0 (La Fleur, Hossain & Salis 2022, Nat Commun 13:5159; port of the authors' reference code)… |
| `predict_signal_peptide` | von Heijne (1986) weight-matrix signal-peptide cleavage-site prediction (EMBOSS sigcleave). |
| `predict_transmembrane_helices` | Hydropathy-based transmembrane helix prediction (Kyte-Doolittle, ≥15 aa). |
| `prosite_to_regex` | Translate a PROSITE pattern string to a .NET regex string. |
| `pwm_score_pvalue` | Exact p-value of a PWM score, P(S >= score) for a random background word (i.i.d. or order-m Markov table; TFM-Pvalue search options, exhaustive mode), or the exact score threshold of a p-value (smallest word score t with P(S >= t) <= pValue) — Touzet & Varré 2007 TFM-Pvalue… |
| `pwm_score_thresholds` | Score thresholds of a PWM from its discretised score distribution (Biopython Bio.motifs.thresholds.ScoreDistribution): background false-pos… |
| `reversal_distance` | Lower-bound reversal distance via breakpoint count for two equal-length permutations. |
| `rna_complement_base` | Returns the RNA complement of a single base (Biopython `complement_rna`): A↔U, G↔C, T→A, IUPAC ambiguity codes complemented. |
| `scan_with_alphabet_pwm` | Build a PWM over an arbitrary alphabet from aligned instances (as create_alphabet_pwm) and score every window of a sequence (Biopython pssm.calculate rule)… |
| `scan_with_pwm` | Scan a DNA sequence with a 4×L Position Weight Matrix; |
| `scan_with_pwm_both_strands` | Scan both strands of a DNA sequence with a 4×L PWM (rows A,C,G,T), as Biopython pssm.search(both=True): minus strand scored with the revers… |
| `shared_motifs_significance` | RSAT oligo-analysis matching-sequence statistics: k-mers (or reverse-complement pairs) present in at least minSequences of the input DNA se… |
| `showalign` | EMBOSS showalign: pretty text display of an alignment against the consensus or a reference row (-show A/I/N/S/D, -order, -similarcase, -refseq, width / margin / numbers / ruler / uppercase ranges, -html / -highlight); 600 + 300 random cases equal showalign 6.6.0 (B09 A3-13c). |
| `spaced_word_distance` | Multiple-pattern spaced-word distance (Leimeister et al. 2014): mean over equal-weight patterns of the per-pattern word-vector distance (euclidean, jensen_shannon, euclidean_counts, …), or `ev` = spaced -d EV evolutionary distance (Morgenstern et al. 2015). Optional `acgtOnly` (spaced's reader and N-word rule; non-letters deleted) and `bothStrands` (spaced default reverse-complement mode). |
| `standardize_repeat_motif` | Standardizes a repeat unit two ways. canonicalClass — the MISA repeat-type class "considering sequence complementary" (misa.pl .statistics,… |
| `stem_energy` | Free energy of an RNA stem (Turner 2004 nearest-neighbor stacking + AU/GU terminal penalties). |
| `tandem_repeat_bernoulli_statistics` | Estimates the Tandem Repeats Finder Bernoulli-model parameters of a tandem-repeat tract (Benson 1999): **PM** (match probability) and **PI*… |
| `tandem_repeat_summary` | Aggregate statistics across all microsatellites in a DNA sequence. |
| `terminal_mismatch_energy` | Closing-pair × first-mismatch terminal stacking energy (Turner 2004). |
| `unique_kmers` | k-mers that occur exactly once in the sequence (Jellyfish "Unique", count == 1; not the distinct k-mers), in lexicographic order; optional `canonical`/`acgtOnly`. |
| `validate_dot_bracket` | Validates that all bracket symbols in a dot-bracket string are balanced. |
| `windowed_complexity` | Sliding-window Shannon entropy + linguistic complexity for a DNA sequence. |
| `windowed_gc_skew` | Sliding-window GC skew along a sequence; |
