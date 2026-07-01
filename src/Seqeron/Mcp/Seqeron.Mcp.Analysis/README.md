# Seqeron.Mcp.Analysis

MCP server — **K-mer, motif, repeat, complexity, RNA-structure and comparative-genomics analysis.**

Exposes **91 tools** wrapping the `Seqeron.Genomics` library. Every tool has an
explicit JSON input/output schema, a Schema+Binding test, and per-tool docs under
[`docs/mcp/tools/analysis/`](../../../../docs/mcp/tools/analysis) — see the
campaign ledger [`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Analysis
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Analysis"]`). See [`docs/mcp/README.md`](../../../../docs/mcp/README.md).

## Tools (91)

| Tool | Description |
|------|-------------|
| `analyze_gc_content` | Comprehensive GC report: overall GC content, GC/AT skew, content/skew variances, and windowed GC profiles. |
| `analyze_kmers` | Aggregate k-mer statistics: total, unique, min/max/avg count, and Shannon entropy. |
| `at_skew` | Whole-sequence AT skew = (A - T) / (A + T). |
| `base_pair_type` | Classification of an RNA base-pair candidate: WatsonCrick, Wobble, or null if bases cannot pair. |
| `bulge_loop_energy` | Free energy of an RNA bulge loop (special-C bonus, n=1 stacking, degeneracy entropy from numStates). |
| `calculate_ani` | Average Nucleotide Identity (ANI) between two genome sequences (fragment-and-match). |
| `can_pair` | Whether two RNA bases form a Watson-Crick (A-U, G-C) or wobble (G-U) pair. |
| `codon_frequencies` | Codon usage frequencies in the specified reading frame (0..2). |
| `compare_genomes` | End-to-end comparative pipeline: RBH orthologs + synteny + rearrangements + summary stats. |
| `compression_ratio` | Estimate sequence repetitiveness as the normalized Lempel-Ziv complexity c/(n/log_b(n)). |
| `count_kmers` | Counts every k-mer (substring of length k) occurrence in a sequence. |
| `count_kmers_both_strands` | k-mer counts on the forward strand combined with counts on the reverse-complement strand. |
| `create_pwm` | Build a log-odds Position Weight Matrix (4×L; |
| `cumulative_gc_skew` | Cumulative GC skew along the sequence — minimum approximates origin and maximum approximates terminus of replication. |
| `dangling_end_energy` | 5' or 3' dangling-end stacking energy for an RNA helix end. |
| `detect_pseudoknots` | Identifies pseudoknots as crossing base pairs: i < i' < j < j' for pairs (i,j) and (i',j'). |
| `detect_rearrangements` | Inversions, deletions, insertions inferred from gene-order disagreements between two genomes. |
| `dinucleotide_frequencies` | Frequency of each adjacent dinucleotide over the alphabet A/T/G/C/U. |
| `dinucleotide_ratios` | Observed/expected ratios for each dinucleotide (e.g., CpG ratio for CpG-island detection). |
| `discover_motifs` | Overrepresented k-mers (de novo motif discovery) in a DNA sequence. |
| `disorder_propensity` | Returns the TOP-IDP propensity value for a single amino acid (Campen 2008). |
| `dust_score` | DUST low-complexity score (BLAST-style, triplet-based) for a DNA sequence. |
| `entropy_profile` | Shannon entropy in sliding windows along the sequence. |
| `find_clumps` | Finds k-mers that occur at least minOccurrences times within any sliding window of size windowSize. |
| `find_common_regions` | All common regions between two DNA sequences with length >= minLength. |
| `find_conserved_clusters` | Gene clusters preserved across multiple genomes. |
| `find_degenerate_motif` | Motif search with IUPAC ambiguity codes (N, R, Y, S, W, K, M, B, D, H, V). |
| `find_direct_repeats` | Identical sequences appearing twice with a spacer between them. |
| `find_exact_motif` | Exact-match motif positions in a DNA sequence via suffix tree. |
| `find_inverted_repeats` | Sequences whose two arms are reverse-complement of each other (hairpin candidates). |
| `find_known_motifs` | Search for a user-provided set of motifs simultaneously. |
| `find_low_complexity_regions` | Entropy-thresholded contiguous low-complexity DNA regions. |
| `find_microsatellites` | Short Tandem Repeats (STRs): 1-6 bp motif units repeated consecutively. |
| `find_motif` | All exact occurrences of a motif (case-insensitive) in a DNA sequence via suffix tree. |
| `find_motif_by_pattern` | Find all (overlapping) regex pattern matches in a protein sequence. |
| `find_motif_by_prosite` | Convert a PROSITE pattern to regex internally, then scan a protein sequence. |
| `find_open_reading_frames` | ORFs in all 6 frames (3 forward + 3 reverse-complement) starting ATG and ending TAA/TAG/TGA. |
| `find_orthologs` | Best-hit ortholog pairs between two genomes by k-mer similarity (one-directional). |
| `find_palindromes` | Sequences identical to their reverse complement (restriction-site candidates). |
| `find_protein_domains` | Detect common protein domains with EXACT PROSITE patterns: zinc finger C2H2 (PS00028), WD-repeats (PS00678), kinase ATP-binding / Walker… |
| `find_protein_low_complexity_regions` | Low-complexity regions in a protein via the SEG algorithm (Wootton & Federhen 1993): sliding-window Shannon entropy in bits/residue, two-… |
| `find_protein_motifs` | Scan a protein sequence against the built-in PROSITE-style motif catalog (N-glycosylation, kinase phosphorylation sites, ATP/GTP P-loop,… |
| `find_reciprocal_best_hits` | Reciprocal best hits (RBH) for stricter ortholog identification. |
| `find_regulatory_elements` | Scan for built-in regulatory motifs (TATA, CAAT, GC-box, Kozak, Shine-Dalgarno, poly(A), E-box, AP-1, NF-κB, CREB). |
| `find_repeats` | All repeated substrings of length >= minLength in a DNA sequence, with their positions. |
| `find_rna_inverted_repeats` | Finds antiparallel complementary regions (potential RNA hairpin stems). |
| `find_shared_motifs` | k-mers present in at least minSequences of the input DNA sequences. |
| `find_stem_loops` | Enumerates hairpin stem-loop candidates with stem, loop, and Turner 2004 free energy. |
| `find_syntenic_blocks` | Collinear runs of orthologous genes between two genomes. |
| `find_tandem_repeats` | Consecutive repeating units (e.g., ATGATGATG) of unit-length >= minUnitLength repeated >= minRepetitions times. |
| `flush_coaxial_stacking` | Coaxial stacking energy for two RNA helices with no intervening unpaired bases. |
| `gc_content_profile` | GC content in sliding windows along the sequence. |
| `gc_skew` | Whole-sequence GC skew = (G - C) / (G + C). |
| `generate_all_kmers` | Enumerate the entire k-mer space for an alphabet (default "ACGT\ |
| `generate_consensus` | IUPAC consensus sequence from aligned equal-length DNA sequences (>25% per position threshold). |
| `generate_dot_plot` | Coordinates of matching k-mers between two sequences for dot-plot visualization. |
| `hairpin_loop_energy` | Free energy of an RNA hairpin loop (Turner 2004 with special tri/tetra/hexaloops, terminal mismatch, all-C and special-GU adjustments). |
| `hydrophobicity_profile` | Sliding-window Kyte-Doolittle hydropathy values for a protein sequence. |
| `internal_loop_energy` | Free energy of a generic RNA internal loop (Turner 2004; |
| `is_disorder_promoting` | Whether an amino acid is in Dunker's disorder-promoting set {A, R, G, Q, S, P, E, K} (Dunker 2001). |
| `kmer_distance` | Euclidean distance between k-mer frequency vectors of two sequences. |
| `kmer_frequencies` | Normalized k-mer counts (each value in [0,1], summing to 1). |
| `kmer_positions` | Zero-based positions of all (overlapping) occurrences of a k-mer. |
| `kmer_spectrum` | Frequency-of-frequencies: for each occurrence count, how many distinct k-mers reach that count. |
| `kmers_with_min_count` | k-mers occurring at least minCount times, sorted descending by count. |
| `mask_low_complexity` | Mask low-complexity windows (DUST-driven) of a DNA sequence with a chosen character. |
| `minimum_free_energy` | Zuker-style minimum free energy with Turner 2004 parameters (O(n³)). |
| `mismatch_coaxial_stacking` | Mismatch-mediated coaxial stacking energy: terminal mismatch + base + WC/GU bonus. |
| `most_frequent_kmers` | Returns all k-mers tied for the maximum occurrence count. |
| `multibranch_loop_energy` | Free energy of an RNA multibranch loop (Turner 2004 affine model: offset + asymmetry + helix term + stacking + strain). |
| `parse_dot_bracket` | Parses dot-bracket notation into a list of base-pair coordinates. |
| `predict_chou_fasman` | Per-window helix/sheet/turn propensities for a protein sequence (Chou-Fasman parameters). |
| `predict_coiled_coils` | Heptad-repeat-based coiled-coil prediction. |
| `predict_disorder` | TOP-IDP disorder prediction (Campen 2008): per-residue scores plus contiguous IDR regions with confidence and subtype classification. |
| `predict_low_complexity_seg` | SEG algorithm (Wootton & Federhen 1993/1996) for low-complexity protein regions: trigger window K1 + extension K2. |
| `predict_morfs` | Predicts Molecular Recognition Features within IDRs by hydropathy enrichment (heuristic, Mohan 2006-inspired). |
| `predict_replication_origin` | Predicts replication origin and terminus from cumulative GC skew extrema. |
| `predict_rna_structure` | Greedy non-overlapping stem-loop selection — produces dot-bracket notation, base pairs, stems, pseudoknots, and total MFE. |
| `predict_signal_peptide` | von Heijne (1986) weight-matrix signal-peptide cleavage-site prediction (EMBOSS sigcleave). |
| `predict_transmembrane_helices` | Hydropathy-based transmembrane helix prediction (Kyte-Doolittle, ≥15 aa). |
| `prosite_to_regex` | Translate a PROSITE pattern string to a .NET regex string. |
| `reversal_distance` | Lower-bound reversal distance via breakpoint count for two equal-length permutations. |
| `rna_complement_base` | Returns the RNA complement (A↔U, G↔C) for a single base. |
| `scan_with_pwm` | Scan a DNA sequence with a 4×L Position Weight Matrix; |
| `stem_energy` | Free energy of an RNA stem (Turner 2004 nearest-neighbor stacking + AU/GU terminal penalties). |
| `tandem_repeat_summary` | Aggregate statistics across all microsatellites in a DNA sequence. |
| `terminal_mismatch_energy` | Closing-pair × first-mismatch terminal stacking energy (Turner 2004). |
| `unique_kmers` | k-mers that occur exactly once in the sequence. |
| `validate_dot_bracket` | Validates that all bracket symbols in a dot-bracket string are balanced. |
| `windowed_complexity` | Sliding-window Shannon entropy + linguistic complexity for a DNA sequence. |
| `windowed_gc_skew` | Sliding-window GC skew along a sequence; |
