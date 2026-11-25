# MCP Completion — Source of Truth & Campaign Ledger

> Reality-based ledger of the 663 MCP tools across 21 servers. Campaign COMPLETE (2026-07-01):
> every tool has a gold-standard binding (`[McpServerTool(Name/Title/ReadOnly)]` + `[Description]`,
> structured record, input validation, real `Seqeron.Genomics` delegation), a Schema+Binding NUnit
> test, and `{tool}.md` + `{tool}.mcp.json` docs. Supersedes `docs/mcp-plan.md` / `docs/mcp-checklist.md` (v4).

## Roll-up

| Metric | Done | Total |
|---|---:|---:|
| Gold-standard binding | 663 | 663 |
| Schema+Binding tests | 663 | 663 |
| Docs (.md + .mcp.json) | 663 | 663 |

## Per-server summary

| Server | Project | Tools | B | T | D | Status |
|---|---|---:|:--:|:--:|:--:|---|
| Core | SuffixTree.Mcp.Core | 18 | 18/18 | 18/18 | 18/18 | ✅ done |
| Sequence | Seqeron.Mcp.Sequence | 35 | 35/35 | 35/35 | 35/35 | ✅ done |
| Parsers | Seqeron.Mcp.Parsers | 41 | 41/41 | 41/41 | 41/41 | ✅ done |
| Alignment | Seqeron.Mcp.Alignment | 59 | 59/59 | 59/59 | 59/59 | ✅ done |
| Analysis | Seqeron.Mcp.Analysis | 157 | 157/157 | 157/157 | 157/157 | ✅ done |
| Annotation | Seqeron.Mcp.Annotation | 97 | 97/97 | 97/97 | 97/97 | ✅ done |
| Chromosome | Seqeron.Mcp.Chromosome | 32 | 32/32 | 32/32 | 32/32 | ✅ done |
| Metagenomics | Seqeron.Mcp.Metagenomics | 19 | 19/19 | 19/19 | 19/19 | ✅ done |
| MolTools | Seqeron.Mcp.MolTools | 48 | 48/48 | 48/48 | 48/48 | ✅ done |
| Phylogenetics | Seqeron.Mcp.Phylogenetics | 13 | 13/13 | 13/13 | 13/13 | ✅ done |
| Population | Seqeron.Mcp.Population | 18 | 18/18 | 18/18 | 18/18 | ✅ done |
| OligoThermodynamics | Seqeron.Mcp.OligoThermodynamics | 11 | 11/11 | 11/11 | 11/11 | ✅ done |
| OverrepresentedWords | Seqeron.Mcp.OverrepresentedWords | 3 | 3/3 | 3/3 | 3/3 | ✅ done |
| EmbossOpenReadingFrames | Seqeron.Mcp.EmbossOpenReadingFrames | 4 | 4/4 | 4/4 | 4/4 | ✅ done |
| EmbossMotifs | Seqeron.Mcp.EmbossMotifs | 12 | 12/12 | 12/12 | 12/12 | ✅ done |
| EmbossPairwiseAlignment | Seqeron.Mcp.EmbossPairwiseAlignment | 16 | 16/16 | 16/16 | 16/16 | ✅ done |
| EmbossRepeats | Seqeron.Mcp.EmbossRepeats | 3 | 3/3 | 3/3 | 3/3 | ✅ done |
| EmbossFeatureTools | Seqeron.Mcp.EmbossFeatureTools | 9 | 9/9 | 9/9 | 9/9 | ✅ done |
| EmbossSequenceFormats | Seqeron.Mcp.EmbossSequenceFormats | 42 | 42/42 | 42/42 | 42/42 | ✅ done |
| AlignmentConsensus | Seqeron.Mcp.AlignmentConsensus | 1 | 1/1 | 1/1 | 1/1 | ✅ done |
| MultipleAlignment | Seqeron.Mcp.MultipleAlignment | 25 | 25/25 | 25/25 | 25/25 | ✅ done |

## Per-tool ledger

### Core (SuffixTree.Mcp.Core) — 18 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `calculate_similarity` | ☑ | ☑ | ☑ |
| 2 | `count_approximate_occurrences` | ☑ | ☑ | ☑ |
| 3 | `edit_distance` | ☑ | ☑ | ☑ |
| 4 | `find_longest_common_region` | ☑ | ☑ | ☑ |
| 5 | `find_longest_repeat` | ☑ | ☑ | ☑ |
| 6 | `hamming_distance` | ☑ | ☑ | ☑ |
| 7 | `suffix_tree_all_lcs` | ☑ | ☑ | ☑ |
| 8 | `suffix_tree_all_lrs` | ☑ | ☑ | ☑ |
| 9 | `suffix_tree_contains` | ☑ | ☑ | ☑ |
| 10 | `suffix_tree_count` | ☑ | ☑ | ☑ |
| 11 | `suffix_tree_find_all` | ☑ | ☑ | ☑ |
| 12 | `suffix_tree_find_mems` | ☑ | ☑ | ☑ |
| 13 | `suffix_tree_find_mums` | ☑ | ☑ | ☑ |
| 14 | `suffix_tree_k_common_substrings` | ☑ | ☑ | ☑ |
| 15 | `suffix_tree_lcs` | ☑ | ☑ | ☑ |
| 16 | `suffix_tree_lrs` | ☑ | ☑ | ☑ |
| 17 | `suffix_tree_maximal_repeats` | ☑ | ☑ | ☑ |
| 18 | `suffix_tree_stats` | ☑ | ☑ | ☑ |

### Sequence (Seqeron.Mcp.Sequence) — 35 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `amino_acid_composition` | ☑ | ☑ | ☑ |
| 2 | `complement_base` | ☑ | ☑ | ☑ |
| 3 | `complexity_compression_ratio` | ☑ | ☑ | ☑ |
| 4 | `complexity_dust_score` | ☑ | ☑ | ☑ |
| 5 | `complexity_kmer_entropy` | ☑ | ☑ | ☑ |
| 6 | `complexity_linguistic` | ☑ | ☑ | ☑ |
| 7 | `complexity_mask_low` | ☑ | ☑ | ☑ |
| 8 | `complexity_shannon` | ☑ | ☑ | ☑ |
| 9 | `dna_reverse_complement` | ☑ | ☑ | ☑ |
| 10 | `dna_validate` | ☑ | ☑ | ☑ |
| 11 | `gc_content` | ☑ | ☑ | ☑ |
| 12 | `hydrophobicity` | ☑ | ☑ | ☑ |
| 13 | `is_valid_dna` | ☑ | ☑ | ☑ |
| 14 | `is_valid_rna` | ☑ | ☑ | ☑ |
| 15 | `isoelectric_point` | ☑ | ☑ | ☑ |
| 16 | `iupac_code` | ☑ | ☑ | ☑ |
| 17 | `iupac_match` | ☑ | ☑ | ☑ |
| 18 | `iupac_matches` | ☑ | ☑ | ☑ |
| 19 | `kmer_analyze` | ☑ | ☑ | ☑ |
| 20 | `kmer_count` | ☑ | ☑ | ☑ |
| 21 | `kmer_distance` | ☑ | ☑ | ☑ |
| 22 | `kmer_entropy` | ☑ | ☑ | ☑ |
| 23 | `linguistic_complexity` | ☑ | ☑ | ☑ |
| 24 | `melting_temperature` | ☑ | ☑ | ☑ |
| 25 | `molecular_weight_nucleotide` | ☑ | ☑ | ☑ |
| 26 | `molecular_weight_protein` | ☑ | ☑ | ☑ |
| 27 | `nucleotide_composition` | ☑ | ☑ | ☑ |
| 28 | `protein_validate` | ☑ | ☑ | ☑ |
| 29 | `rna_from_dna` | ☑ | ☑ | ☑ |
| 30 | `rna_validate` | ☑ | ☑ | ☑ |
| 31 | `shannon_entropy` | ☑ | ☑ | ☑ |
| 32 | `summarize_sequence` | ☑ | ☑ | ☑ |
| 33 | `thermodynamics` | ☑ | ☑ | ☑ |
| 34 | `translate_dna` | ☑ | ☑ | ☑ |
| 35 | `translate_rna` | ☑ | ☑ | ☑ |

### Parsers (Seqeron.Mcp.Parsers) — 41 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `bed_filter` | ☑ | ☑ | ☑ |
| 2 | `bed_intersect` | ☑ | ☑ | ☑ |
| 3 | `bed_merge` | ☑ | ☑ | ☑ |
| 4 | `bed_parse` | ☑ | ☑ | ☑ |
| 5 | `embl_features` | ☑ | ☑ | ☑ |
| 6 | `embl_parse` | ☑ | ☑ | ☑ |
| 7 | `embl_statistics` | ☑ | ☑ | ☑ |
| 8 | `fasta_format` | ☑ | ☑ | ☑ |
| 9 | `fasta_parse` | ☑ | ☑ | ☑ |
| 10 | `fasta_write` | ☑ | ☑ | ☑ |
| 11 | `fastq_detect_encoding` | ☑ | ☑ | ☑ |
| 12 | `fastq_encode_quality` | ☑ | ☑ | ☑ |
| 13 | `fastq_error_to_phred` | ☑ | ☑ | ☑ |
| 14 | `fastq_filter` | ☑ | ☑ | ☑ |
| 15 | `fastq_format` | ☑ | ☑ | ☑ |
| 16 | `fastq_parse` | ☑ | ☑ | ☑ |
| 17 | `fastq_phred_to_error` | ☑ | ☑ | ☑ |
| 18 | `fastq_statistics` | ☑ | ☑ | ☑ |
| 19 | `fastq_trim_adapter` | ☑ | ☑ | ☑ |
| 20 | `fastq_trim_quality` | ☑ | ☑ | ☑ |
| 21 | `fastq_write` | ☑ | ☑ | ☑ |
| 22 | `genbank_extract_sequence` | ☑ | ☑ | ☑ |
| 23 | `genbank_features` | ☑ | ☑ | ☑ |
| 24 | `genbank_parse` | ☑ | ☑ | ☑ |
| 25 | `genbank_parse_location` | ☑ | ☑ | ☑ |
| 26 | `genbank_statistics` | ☑ | ☑ | ☑ |
| 27 | `gff_filter` | ☑ | ☑ | ☑ |
| 28 | `gff_parse` | ☑ | ☑ | ☑ |
| 29 | `gff_statistics` | ☑ | ☑ | ☑ |
| 30 | `vcf_classify` | ☑ | ☑ | ☑ |
| 31 | `vcf_filter` | ☑ | ☑ | ☑ |
| 32 | `vcf_has_flag` | ☑ | ☑ | ☑ |
| 33 | `vcf_is_het` | ☑ | ☑ | ☑ |
| 34 | `vcf_is_hom_alt` | ☑ | ☑ | ☑ |
| 35 | `vcf_is_hom_ref` | ☑ | ☑ | ☑ |
| 36 | `vcf_is_indel` | ☑ | ☑ | ☑ |
| 37 | `vcf_is_snp` | ☑ | ☑ | ☑ |
| 38 | `vcf_parse` | ☑ | ☑ | ☑ |
| 39 | `vcf_statistics` | ☑ | ☑ | ☑ |
| 40 | `vcf_variant_length` | ☑ | ☑ | ☑ |
| 41 | `vcf_write` | ☑ | ☑ | ☑ |

### Alignment (Seqeron.Mcp.Alignment) — 59 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `alignment_statistics` | ☑ | ☑ | ☑ |
| 2 | `assemble_de_bruijn` | ☑ | ☑ | ☑ |
| 3 | `assemble_olc` | ☑ | ☑ | ☑ |
| 4 | `assembly_stats` | ☑ | ☑ | ☑ |
| 5 | `calculate_coverage` | ☑ | ☑ | ☑ |
| 6 | `compute_consensus` | ☑ | ☑ | ☑ |
| 7 | `damerau_alignment` | ☑ | ☑ | ☑ |
| 8 | `damerau_levenshtein_distance` | ☑ | ☑ | ☑ |
| 9 | `edit_alignment` | ☑ | ☑ | ☑ |
| 10 | `error_correct_reads` | ☑ | ☑ | ☑ |
| 11 | `find_all_overlaps` | ☑ | ☑ | ☑ |
| 12 | `find_best_match` | ☑ | ☑ | ☑ |
| 13 | `find_edit_end_positions` | ☑ | ☑ | ☑ |
| 14 | `find_overlap` | ☑ | ☑ | ☑ |
| 15 | `find_with_edits` | ☑ | ☑ | ☑ |
| 16 | `find_with_mismatches` | ☑ | ☑ | ☑ |
| 17 | `format_alignment` | ☑ | ☑ | ☑ |
| 18 | `frequent_kmers_with_mismatches` | ☑ | ☑ | ☑ |
| 19 | `frequent_kmers_with_mismatches_and_revcomp` | ☑ | ☑ | ☑ |
| 20 | `general_gap_align` | ☑ | ☑ | ☑ |
| 21 | `global_align` | ☑ | ☑ | ☑ |
| 22 | `global_align_affine` | ☑ | ☑ | ☑ |
| 23 | `global_align_co_optimal` | ☑ | ☑ | ☑ |
| 24 | `local_align` | ☑ | ☑ | ☑ |
| 25 | `local_align_affine` | ☑ | ☑ | ☑ |
| 26 | `local_align_co_optimal` | ☑ | ☑ | ☑ |
| 27 | `local_align_top_n` | ☑ | ☑ | ☑ |
| 28 | `merge_contigs` | ☑ | ☑ | ☑ |
| 29 | `multiple_align` | ☑ | ☑ | ☑ |
| 30 | `multiple_align_anchored` | ☑ | ☑ | ☑ |
| 31 | `multiple_align_exact` | ☑ | ☑ | ☑ |
| 32 | `multiple_align_muscle` | ☑ | ☑ | ☑ |
| 33 | `muscle_msa_command` | ☑ | ☑ | ☑ |
| 34 | `quality_trim_reads` | ☑ | ☑ | ☑ |
| 35 | `scaffold_contigs` | ☑ | ☑ | ☑ |
| 36 | `semi_global_align` | ☑ | ☑ | ☑ |
| 37 | `semi_global_align_affine` | ☑ | ☑ | ☑ |
| 38 | `sequence_identity` | ☑ | ☑ | ☑ |
| 39 | `sum_of_pairs_score` | ☑ | ☑ | ☑ |
| 40 | `multiple_align_mafft` | ☑ | ☑ | ☑ |
| 41 | `dual_affine_align` | ☑ | ☑ | ☑ |
| 42 | `ksw2_extend` | ☑ | ☑ | ☑ |
| 43 | `align_gap_profile` | ☑ | ☑ | ☑ |
| 44 | `minimap2_chains` | ☑ | ☑ | ☑ |
| 45 | `minimap2_map` | ☑ | ☑ | ☑ |
| 46 | `write_coordinate_alignment` | ☑ | ☑ | ☑ |
| 47 | `parse_coordinate_alignment` | ☑ | ☑ | ☑ |
| 48 | `spoa_command` | ☑ | ☑ | ☑ |
| 49 | `write_big_alignment` | ☑ | ☑ | ☑ |
| 50 | `parse_big_alignment` | ☑ | ☑ | ☑ |
| 51 | `local_align_zscore` | ☑ | ☑ | ☑ |
| 52 | `pairwise_aligner` | ☑ | ☑ | ☑ |
| 53 | `alignment_counts` | ☑ | ☑ | ☑ |
| 54 | `alignment_operation` | ☑ | ☑ | ☑ |
| 55 | `parasail_score` | ☑ | ☑ | ☑ |
| 56 | `banded_align` | ☑ | ☑ | ☑ |
| 57 | `codon_aligner` | ☑ | ☑ | ☑ |
| 58 | `multiple_seq_alignment` | ☑ | ☑ | ☑ |
| 59 | `mafft_command` | ☑ | ☑ | ☑ |

### OligoThermodynamics (Seqeron.Mcp.OligoThermodynamics) — 11 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `analyze_oligo` | ☑ | ☑ | ☑ |
| 2 | `design_antisense_probes` | ☑ | ☑ | ☑ |
| 3 | `design_molecular_beacon` | ☑ | ☑ | ☑ |
| 4 | `design_primers` | ☑ | ☑ | ☑ |
| 5 | `design_probes` | ☑ | ☑ | ☑ |
| 6 | `design_probes_primer3` | ☑ | ☑ | ☑ |
| 7 | `design_tiling_probes` | ☑ | ☑ | ☑ |
| 8 | `evaluate_primer` | ☑ | ☑ | ☑ |
| 9 | `generate_primer_candidates` | ☑ | ☑ | ☑ |
| 10 | `primer_dimer` | ☑ | ☑ | ☑ |
| 11 | `validate_probe` | ☑ | ☑ | ☑ |

### OverrepresentedWords (Seqeron.Mcp.OverrepresentedWords) — 3 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `dyad_analysis` | ☑ | ☑ | ☑ |
| 2 | `oligo_analysis` | ☑ | ☑ | ☑ |
| 3 | `shared_motifs_significance` | ☑ | ☑ | ☑ |

### EmbossSequenceFormats (Seqeron.Mcp.EmbossSequenceFormats) — 42 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `cutseq` | ☑ | ☑ | ☑ |
| 2 | `dbtell` | ☑ | ☑ | ☑ |
| 3 | `degapseq` | ☑ | ☑ | ☑ |
| 4 | `descseq` | ☑ | ☑ | ☑ |
| 5 | `drfinddata` | ☑ | ☑ | ☑ |
| 6 | `entret` | ☑ | ☑ | ☑ |
| 7 | `extractseq` | ☑ | ☑ | ☑ |
| 8 | `infoseq` | ☑ | ☑ | ☑ |
| 9 | `listor` | ☑ | ☑ | ☑ |
| 10 | `makenucseq` | ☑ | ☑ | ☑ |
| 11 | `makeprotseq` | ☑ | ☑ | ☑ |
| 12 | `maskambignuc` | ☑ | ☑ | ☑ |
| 13 | `maskambigprot` | ☑ | ☑ | ☑ |
| 14 | `maskseq` | ☑ | ☑ | ☑ |
| 15 | `msbar` | ☑ | ☑ | ☑ |
| 16 | `newseq` | ☑ | ☑ | ☑ |
| 17 | `noreturn` | ☑ | ☑ | ☑ |
| 18 | `nospace` | ☑ | ☑ | ☑ |
| 19 | `notab` | ☑ | ☑ | ☑ |
| 20 | `notseq` | ☑ | ☑ | ☑ |
| 21 | `nthseq` | ☑ | ☑ | ☑ |
| 22 | `nthseqset` | ☑ | ☑ | ☑ |
| 23 | `pasteseq` | ☑ | ☑ | ☑ |
| 24 | `prettyseq` | ☑ | ☑ | ☑ |
| 25 | `revseq` | ☑ | ☑ | ☑ |
| 26 | `seqcount` | ☑ | ☑ | ☑ |
| 27 | `seqret` | ☑ | ☑ | ☑ |
| 28 | `seqretsplit` | ☑ | ☑ | ☑ |
| 29 | `seqxref` | ☑ | ☑ | ☑ |
| 30 | `showdb` | ☑ | ☑ | ☑ |
| 31 | `showorf` | ☑ | ☑ | ☑ |
| 32 | `shuffleseq` | ☑ | ☑ | ☑ |
| 33 | `skipseq` | ☑ | ☑ | ☑ |
| 34 | `splitter` | ☑ | ☑ | ☑ |
| 35 | `transeq` | ☑ | ☑ | ☑ |
| 36 | `trimest` | ☑ | ☑ | ☑ |
| 37 | `trimseq` | ☑ | ☑ | ☑ |
| 38 | `trimspace` | ☑ | ☑ | ☑ |
| 39 | `union` | ☑ | ☑ | ☑ |
| 40 | `wordcount` | ☑ | ☑ | ☑ |
| 41 | `write_alignment` | ☑ | ☑ | ☑ |
| 42 | `yank` | ☑ | ☑ | ☑ |

### EmbossFeatureTools (Seqeron.Mcp.EmbossFeatureTools) — 9 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `coderet` | ☑ | ☑ | ☑ |
| 2 | `extractfeat` | ☑ | ☑ | ☑ |
| 3 | `featcopy` | ☑ | ☑ | ☑ |
| 4 | `featmerge` | ☑ | ☑ | ☑ |
| 5 | `featreport` | ☑ | ☑ | ☑ |
| 6 | `feattext` | ☑ | ☑ | ☑ |
| 7 | `maskfeat` | ☑ | ☑ | ☑ |
| 8 | `showfeat` | ☑ | ☑ | ☑ |
| 9 | `twofeat` | ☑ | ☑ | ☑ |

### EmbossRepeats (Seqeron.Mcp.EmbossRepeats) — 3 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `find_tandem_repeats_emboss` | ☑ | ☑ | ☑ |
| 2 | `seqmatchall` | ☑ | ☑ | ☑ |
| 3 | `wordmatch` | ☑ | ☑ | ☑ |

### EmbossPairwiseAlignment (Seqeron.Mcp.EmbossPairwiseAlignment) — 16 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `alignment_distance_matrix` | ☑ | ☑ | ☑ |
| 2 | `diffseq` | ☑ | ☑ | ☑ |
| 3 | `dotmatcher` | ☑ | ☑ | ☑ |
| 4 | `est2genome` | ☑ | ☑ | ☑ |
| 5 | `matcher_align` | ☑ | ☑ | ☑ |
| 6 | `megamerger` | ☑ | ☑ | ☑ |
| 7 | `merger` | ☑ | ☑ | ☑ |
| 8 | `multiple_alignment_statistics` | ☑ | ☑ | ☑ |
| 9 | `needle_align` | ☑ | ☑ | ☑ |
| 10 | `needleall_align` | ☑ | ☑ | ☑ |
| 11 | `plotcon` | ☑ | ☑ | ☑ |
| 12 | `showalign` | ☑ | ☑ | ☑ |
| 13 | `stretcher_align` | ☑ | ☑ | ☑ |
| 14 | `supermatcher_align` | ☑ | ☑ | ☑ |
| 15 | `water_align` | ☑ | ☑ | ☑ |
| 16 | `wordfinder` | ☑ | ☑ | ☑ |

### EmbossMotifs (Seqeron.Mcp.EmbossMotifs) — 12 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `consambig` | ☑ | ☑ | ☑ |
| 2 | `emboss_report_file` | ☑ | ☑ | ☑ |
| 3 | `find_known_motifs` | ☑ | ☑ | ☑ |
| 4 | `find_motif` | ☑ | ☑ | ☑ |
| 5 | `find_pattern_hits` | ☑ | ☑ | ☑ |
| 6 | `fuzzpro` | ☑ | ☑ | ☑ |
| 7 | `fuzztran` | ☑ | ☑ | ☑ |
| 8 | `generate_emboss_consensus` | ☑ | ☑ | ☑ |
| 9 | `marscan` | ☑ | ☑ | ☑ |
| 10 | `patmatdb` | ☑ | ☑ | ☑ |
| 11 | `patmatmotifs` | ☑ | ☑ | ☑ |
| 12 | `prosextract` | ☑ | ☑ | ☑ |

### EmbossOpenReadingFrames (Seqeron.Mcp.EmbossOpenReadingFrames) — 4 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `checktrans` | ☑ | ☑ | ☑ |
| 2 | `getorf` | ☑ | ☑ | ☑ |
| 3 | `plotorf` | ☑ | ☑ | ☑ |
| 4 | `tcode` | ☑ | ☑ | ☑ |

### AlignmentConsensus (Seqeron.Mcp.AlignmentConsensus) — 1 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `generate_decipher_consensus` | ☑ | ☑ | ☑ |

### MultipleAlignment (Seqeron.Mcp.MultipleAlignment) — 25 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `multiple_align_clustalw` | ☑ | ☑ | ☑ |
| 2 | `multiple_align_consistency` | ☑ | ☑ | ☑ |
| 3 | `multiple_align_iterative` | ☑ | ☑ | ☑ |
| 4 | `multiple_align_progressive` | ☑ | ☑ | ☑ |
| 5 | `multiple_align_tcoffee` | ☑ | ☑ | ☑ |
| 6 | `clustalw_pairwise_align` | ☑ | ☑ | ☑ |
| 7 | `clustalw_quick_pairwise_align` | ☑ | ☑ | ☑ |
| 8 | `tcoffee_slow_pair_align` | ☑ | ☑ | ☑ |
| 9 | `tcoffee_lalign_top_alignments` | ☑ | ☑ | ☑ |
| 10 | `tcoffee_primary_library` | ☑ | ☑ | ☑ |
| 11 | `tcoffee_guide_tree` | ☑ | ☑ | ☑ |
| 12 | `tcoffee_ktup_similarity` | ☑ | ☑ | ☑ |
| 13 | `clustalw_profile_align` | ☑ | ☑ | ☑ |
| 14 | `clustalw_sequences_to_profile` | ☑ | ☑ | ☑ |
| 15 | `clustalw_command` | ☑ | ☑ | ☑ |
| 16 | `multiple_align_kalign` | ☑ | ☑ | ☑ |
| 17 | `kalign_command` | ☑ | ☑ | ☑ |
| 18 | `multiple_align_probcons` | ☑ | ☑ | ☑ |
| 19 | `probcons_command` | ☑ | ☑ | ☑ |
| 20 | `multiple_align_poa` | ☑ | ☑ | ☑ |
| 21 | `poa_command` | ☑ | ☑ | ☑ |
| 22 | `clustalw_guide_tree` | ☑ | ☑ | ☑ |
| 23 | `tcoffee_blast_msa_library` | ☑ | ☑ | ☑ |
| 24 | `tcoffee_tblastx_library` | ☑ | ☑ | ☑ |
| 25 | `tcoffee_evaluate` | ☑ | ☑ | ☑ |

### Analysis (Seqeron.Mcp.Analysis) — 157 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `alphabet_pwm_score_pvalue` | ☑ | ☑ | ☑ |
| 2 | `alphabet_pwm_score_thresholds` | ☑ | ☑ | ☑ |
| 3 | `analyze_gc_content` | ☑ | ☑ | ☑ |
| 4 | `analyze_kmers` | ☑ | ☑ | ☑ |
| 5 | `at_skew` | ☑ | ☑ | ☑ |
| 6 | `base_pair_type` | ☑ | ☑ | ☑ |
| 7 | `bulge_loop_energy` | ☑ | ☑ | ☑ |
| 8 | `calculate_ani` | ☑ | ☑ | ☑ |
| 9 | `can_pair` | ☑ | ☑ | ☑ |
| 10 | `chou_fasman_assign` | ☑ | ☑ | ☑ |
| 11 | `codon_frequencies` | ☑ | ☑ | ☑ |
| 12 | `compare_genomes` | ☑ | ☑ | ☑ |
| 13 | `compression_ratio` | ☑ | ☑ | ☑ |
| 14 | `count_kmers` | ☑ | ☑ | ☑ |
| 15 | `count_kmers_both_strands` | ☑ | ☑ | ☑ |
| 16 | `create_alphabet_pwm` | ☑ | ☑ | ☑ |
| 17 | `create_pwm` | ☑ | ☑ | ☑ |
| 18 | `cumulative_gc_skew` | ☑ | ☑ | ☑ |
| 19 | `dangling_end_energy` | ☑ | ☑ | ☑ |
| 20 | `detect_pseudoknots` | ☑ | ☑ | ☑ |
| 21 | `detect_rearrangements` | ☑ | ☑ | ☑ |
| 22 | `dinucleotide_frequencies` | ☑ | ☑ | ☑ |
| 23 | `dinucleotide_ratios` | ☑ | ☑ | ☑ |
| 24 | `discover_motifs` | ☑ | ☑ | ☑ |
| 25 | `disorder_propensity` | ☑ | ☑ | ☑ |
| 26 | `dust_score` | ☑ | ☑ | ☑ |
| 27 | `entropy_profile` | ☑ | ☑ | ☑ |
| 28 | `find_approximate_direct_repeats` | ☑ | ☑ | ☑ |
| 29 | `find_approximate_tandem_repeats` | ☑ | ☑ | ☑ |
| 30 | `find_clump_windows` | ☑ | ☑ | ☑ |
| 31 | `find_clumps` | ☑ | ☑ | ☑ |
| 32 | `find_common_regions` | ☑ | ☑ | ☑ |
| 33 | `find_conserved_clusters` | ☑ | ☑ | ☑ |
| 34 | `find_degenerate_motif` | ☑ | ☑ | ☑ |
| 35 | `find_degenerate_repeats` | ☑ | ☑ | ☑ |
| 36 | `find_direct_repeats` | ☑ | ☑ | ☑ |
| 37 | `find_exact_motif` | ☑ | ☑ | ☑ |
| 38 | `find_inverted_repeats` | ☑ | ☑ | ☑ |
| 39 | `find_inverted_repeats_scored` | ☑ | ☑ | ☑ |
| 40 | `find_longdust_regions` | ☑ | ☑ | ☑ |
| 41 | `find_low_complexity_intervals` | ☑ | ☑ | ☑ |
| 42 | `find_low_complexity_regions` | ☑ | ☑ | ☑ |
| 43 | `find_microsatellites` | ☑ | ☑ | ☑ |
| 44 | `find_motif_by_pattern` | ☑ | ☑ | ☑ |
| 45 | `find_motif_by_prosite` | ☑ | ☑ | ☑ |
| 46 | `find_open_reading_frames` | ☑ | ☑ | ☑ |
| 47 | `find_orthologs` | ☑ | ☑ | ☑ |
| 48 | `find_palindromes` | ☑ | ☑ | ☑ |
| 49 | `find_promoter_elements_by_matrix` | ☑ | ☑ | ☑ |
| 50 | `find_protein_domains` | ☑ | ☑ | ☑ |
| 51 | `find_protein_low_complexity_regions` | ☑ | ☑ | ☑ |
| 52 | `find_protein_motifs` | ☑ | ☑ | ☑ |
| 53 | `find_reciprocal_best_hits` | ☑ | ☑ | ☑ |
| 54 | `find_regulatory_elements` | ☑ | ☑ | ☑ |
| 55 | `find_regulatory_elements_both_strands` | ☑ | ☑ | ☑ |
| 56 | `find_repeats` | ☑ | ☑ | ☑ |
| 57 | `find_reverse_complement_repeats` | ☑ | ☑ | ☑ |
| 58 | `find_rna_inverted_repeats` | ☑ | ☑ | ☑ |
| 59 | `find_shared_motifs` | ☑ | ☑ | ☑ |
| 60 | `find_sigma70_promoters` | ☑ | ☑ | ☑ |
| 61 | `find_stem_loops` | ☑ | ☑ | ☑ |
| 62 | `find_supermaximal_repeats` | ☑ | ☑ | ☑ |
| 63 | `find_syntenic_blocks` | ☑ | ☑ | ☑ |
| 64 | `find_tandem_repeats` | ☑ | ☑ | ☑ |
| 65 | `flush_coaxial_stacking` | ☑ | ☑ | ☑ |
| 66 | `gc_content_profile` | ☑ | ☑ | ☑ |
| 67 | `gc_skew` | ☑ | ☑ | ☑ |
| 68 | `generate_all_kmers` | ☑ | ☑ | ☑ |
| 69 | `generate_cavener_consensus` | ☑ | ☑ | ☑ |
| 70 | `generate_consensus` | ☑ | ☑ | ☑ |
| 71 | `generate_dot_plot` | ☑ | ☑ | ☑ |
| 72 | `generate_dumb_consensus` | ☑ | ☑ | ☑ |
| 73 | `hairpin_loop_energy` | ☑ | ☑ | ☑ |
| 74 | `hydrophobicity_profile` | ☑ | ☑ | ☑ |
| 75 | `internal_loop_energy` | ☑ | ☑ | ☑ |
| 76 | `is_disorder_promoting` | ☑ | ☑ | ☑ |
| 77 | `kmer_d2_statistics` | ☑ | ☑ | ☑ |
| 78 | `kmer_distance` | ☑ | ☑ | ☑ |
| 79 | `kmer_frequencies` | ☑ | ☑ | ☑ |
| 80 | `kmer_jaccard` | ☑ | ☑ | ☑ |
| 81 | `kmer_positions` | ☑ | ☑ | ☑ |
| 82 | `kmer_spectrum` | ☑ | ☑ | ☑ |
| 83 | `kmers_with_min_count` | ☑ | ☑ | ☑ |
| 84 | `lempel_ziv_complexity` | ☑ | ☑ | ☑ |
| 85 | `longdust_score` | ☑ | ☑ | ☑ |
| 86 | `mask_approximate_tandem_repeats` | ☑ | ☑ | ☑ |
| 87 | `mask_low_complexity` | ☑ | ☑ | ☑ |
| 88 | `minimum_free_energy` | ☑ | ☑ | ☑ |
| 89 | `mismatch_coaxial_stacking` | ☑ | ☑ | ☑ |
| 90 | `most_frequent_kmers` | ☑ | ☑ | ☑ |
| 91 | `multibranch_loop_energy` | ☑ | ☑ | ☑ |
| 92 | `parse_dot_bracket` | ☑ | ☑ | ☑ |
| 93 | `predict_chou_fasman` | ☑ | ☑ | ☑ |
| 94 | `predict_coiled_coils` | ☑ | ☑ | ☑ |
| 95 | `predict_disorder` | ☑ | ☑ | ☑ |
| 96 | `predict_low_complexity_seg` | ☑ | ☑ | ☑ |
| 97 | `predict_morfs` | ☑ | ☑ | ☑ |
| 98 | `predict_replication_origin` | ☑ | ☑ | ☑ |
| 99 | `predict_rna_structure` | ☑ | ☑ | ☑ |
| 100 | `predict_sigma70_promoters` | ☑ | ☑ | ☑ |
| 101 | `predict_signal_peptide` | ☑ | ☑ | ☑ |
| 102 | `predict_transmembrane_helices` | ☑ | ☑ | ☑ |
| 103 | `prosite_to_regex` | ☑ | ☑ | ☑ |
| 104 | `pwm_score_pvalue` | ☑ | ☑ | ☑ |
| 105 | `pwm_score_thresholds` | ☑ | ☑ | ☑ |
| 106 | `reversal_distance` | ☑ | ☑ | ☑ |
| 107 | `rna_complement_base` | ☑ | ☑ | ☑ |
| 108 | `scan_with_alphabet_pwm` | ☑ | ☑ | ☑ |
| 109 | `scan_with_pwm` | ☑ | ☑ | ☑ |
| 110 | `scan_with_pwm_both_strands` | ☑ | ☑ | ☑ |
| 111 | `spaced_word_distance` | ☑ | ☑ | ☑ |
| 112 | `standardize_repeat_motif` | ☑ | ☑ | ☑ |
| 113 | `stem_energy` | ☑ | ☑ | ☑ |
| 114 | `tandem_repeat_bernoulli_statistics` | ☑ | ☑ | ☑ |
| 115 | `tandem_repeat_summary` | ☑ | ☑ | ☑ |
| 116 | `terminal_mismatch_energy` | ☑ | ☑ | ☑ |
| 117 | `unique_kmers` | ☑ | ☑ | ☑ |
| 118 | `validate_dot_bracket` | ☑ | ☑ | ☑ |
| 119 | `windowed_complexity` | ☑ | ☑ | ☑ |
| 120 | `windowed_gc_skew` | ☑ | ☑ | ☑ |
| 121 | `find_all_longest_repeats` | ☑ | ☑ | ☑ |
| 122 | `find_longest_common_region_mismatches` | ☑ | ☑ | ☑ |
| 123 | `nucmer_clusters` | ☑ | ☑ | ☑ |
| 124 | `mash_bounds` | ☑ | ☑ | ☑ |
| 125 | `mash_dist` | ☑ | ☑ | ☑ |
| 126 | `mash_triangle` | ☑ | ☑ | ☑ |
| 127 | `mash_screen` | ☑ | ☑ | ☑ |
| 128 | `mash_sketch` | ☑ | ☑ | ☑ |
| 129 | `mash_info` | ☑ | ☑ | ☑ |
| 130 | `mash_paste` | ☑ | ☑ | ☑ |
| 131 | `nucmer` | ☑ | ☑ | ☑ |
| 132 | `delta_filter` | ☑ | ☑ | ☑ |
| 133 | `show_coords` | ☑ | ☑ | ☑ |
| 134 | `show_aligns` | ☑ | ☑ | ☑ |
| 135 | `show_snps` | ☑ | ☑ | ☑ |
| 136 | `show_diff` | ☑ | ☑ | ☑ |
| 137 | `show_tiling` | ☑ | ☑ | ☑ |
| 138 | `dnadiff` | ☑ | ☑ | ☑ |
| 139 | `mummer` | ☑ | ☑ | ☑ |
| 140 | `promer` | ☑ | ☑ | ☑ |
| 141 | `mummer_gaps` | ☑ | ☑ | ☑ |
| 142 | `delta2blocks` | ☑ | ☑ | ☑ |
| 143 | `mummer_annotate` | ☑ | ☑ | ☑ |
| 144 | `run_mummer1` | ☑ | ☑ | ☑ |
| 145 | `delta2maf` | ☑ | ☑ | ☑ |
| 146 | `nucmer2xfig` | ☑ | ☑ | ☑ |
| 147 | `combine_mums` | ☑ | ☑ | ☑ |
| 148 | `run_mummer3` | ☑ | ☑ | ☑ |
| 149 | `repeat_match` | ☑ | ☑ | ☑ |
| 150 | `exact_tandems` | ☑ | ☑ | ☑ |
| 151 | `mapview` | ☑ | ☑ | ☑ |
| 152 | `mummerplot` | ☑ | ☑ | ☑ |
| 153 | `mash_taxscreen` | ☑ | ☑ | ☑ |
| 154 | `find_longest_common_region_multi` | ☑ | ☑ | ☑ |
| 155 | `find_repeat_pairs` | ☑ | ☑ | ☑ |
| 156 | `alignment_summary_info` | ☑ | ☑ | ☑ |
| 157 | `msa_conservation` | ☑ | ☑ | ☑ |

### Annotation (Seqeron.Mcp.Annotation) — 97 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `align_mirna_to_target` | ☑ | ☑ | ☑ |
| 2 | `analyze_target_context` | ☑ | ☑ | ☑ |
| 3 | `annotate_histone_modifications` | ☑ | ☑ | ☑ |
| 4 | `annotate_regulatory_elements` | ☑ | ☑ | ☑ |
| 5 | `annotate_svs` | ☑ | ☑ | ☑ |
| 6 | `annotate_variant_on_transcripts` | ☑ | ☑ | ☑ |
| 7 | `annotate_variants` | ☑ | ☑ | ☑ |
| 8 | `assemble_breakpoint_sequence` | ☑ | ☑ | ☑ |
| 9 | `build_coexpression_network` | ☑ | ☑ | ☑ |
| 10 | `calculate_conservation` | ☑ | ☑ | ☑ |
| 11 | `calculate_tpm` | ☑ | ☑ | ☑ |
| 12 | `call_variants` | ☑ | ☑ | ☑ |
| 13 | `call_variants_from_alignment` | ☑ | ☑ | ☑ |
| 14 | `can_pair` | ☑ | ☑ | ☑ |
| 15 | `classify_mutation` | ☑ | ☑ | ☑ |
| 16 | `classify_variant` | ☑ | ☑ | ☑ |
| 17 | `cluster_discordant_pairs` | ☑ | ☑ | ☑ |
| 18 | `cluster_genes_by_expression` | ☑ | ☑ | ☑ |
| 19 | `cluster_split_reads` | ☑ | ☑ | ☑ |
| 20 | `coding_potential` | ☑ | ☑ | ☑ |
| 21 | `codon_usage` | ☑ | ☑ | ☑ |
| 22 | `compare_seed_regions` | ☑ | ☑ | ☑ |
| 23 | `cpg_observed_expected` | ☑ | ☑ | ☑ |
| 24 | `create_mirna` | ☑ | ☑ | ☑ |
| 25 | `detect_alternative_splicing` | ☑ | ☑ | ☑ |
| 26 | `detect_differential_splicing` | ☑ | ☑ | ☑ |
| 27 | `detect_isoform_switching` | ☑ | ☑ | ☑ |
| 28 | `differential_expression` | ☑ | ☑ | ☑ |
| 29 | `enrichment_score` | ☑ | ☑ | ☑ |
| 30 | `epigenetic_age` | ☑ | ☑ | ☑ |
| 31 | `filter_svs` | ☑ | ☑ | ☑ |
| 32 | `find_acceptor_sites` | ☑ | ☑ | ☑ |
| 33 | `find_accessible_regions` | ☑ | ☑ | ☑ |
| 34 | `find_branch_points` | ☑ | ☑ | ☑ |
| 35 | `find_conserved_elements` | ☑ | ☑ | ☑ |
| 36 | `find_cpg_islands` | ☑ | ☑ | ☑ |
| 37 | `find_cpg_sites` | ☑ | ☑ | ☑ |
| 38 | `find_deletions` | ☑ | ☑ | ☑ |
| 39 | `find_discordant_pairs` | ☑ | ☑ | ☑ |
| 40 | `find_dmrs` | ☑ | ☑ | ☑ |
| 41 | `find_dominant_isoforms` | ☑ | ☑ | ☑ |
| 42 | `find_donor_sites` | ☑ | ☑ | ☑ |
| 43 | `find_indels` | ☑ | ☑ | ☑ |
| 44 | `find_insertions` | ☑ | ☑ | ☑ |
| 45 | `find_methylation_sites` | ☑ | ☑ | ☑ |
| 46 | `find_microhomology` | ☑ | ☑ | ☑ |
| 47 | `find_mirna_target_sites` | ☑ | ☑ | ☑ |
| 48 | `find_orfs` | ☑ | ☑ | ☑ |
| 49 | `find_pre_mirna_hairpins` | ☑ | ☑ | ☑ |
| 50 | `find_promoter_motifs` | ☑ | ☑ | ☑ |
| 51 | `find_repetitive_elements` | ☑ | ☑ | ☑ |
| 52 | `find_retained_intron_candidates` | ☑ | ☑ | ☑ |
| 53 | `find_ribosome_binding_sites` | ☑ | ☑ | ☑ |
| 54 | `find_similar_mirnas` | ☑ | ☑ | ☑ |
| 55 | `find_skipped_exon_events` | ☑ | ☑ | ☑ |
| 56 | `find_snps` | ☑ | ☑ | ☑ |
| 57 | `find_snps_direct` | ☑ | ☑ | ☑ |
| 58 | `find_split_reads` | ☑ | ☑ | ☑ |
| 59 | `format_vcf_info` | ☑ | ☑ | ☑ |
| 60 | `generate_seed_variants` | ☑ | ☑ | ☑ |
| 61 | `genotype_sv` | ☑ | ☑ | ☑ |
| 62 | `group_by_seed_family` | ☑ | ☑ | ☑ |
| 63 | `identify_cnvs` | ☑ | ☑ | ☑ |
| 64 | `impact_level` | ☑ | ☑ | ☑ |
| 65 | `is_within_coding_region` | ☑ | ☑ | ☑ |
| 66 | `is_wobble_pair` | ☑ | ☑ | ☑ |
| 67 | `log2_transform` | ☑ | ☑ | ☑ |
| 68 | `longest_orfs_per_frame` | ☑ | ☑ | ☑ |
| 69 | `maxent_score` | ☑ | ☑ | ☑ |
| 70 | `merge_overlapping_svs` | ☑ | ☑ | ☑ |
| 71 | `methylation_from_bisulfite` | ☑ | ☑ | ☑ |
| 72 | `methylation_profile` | ☑ | ☑ | ☑ |
| 73 | `mirna_seed_sequence` | ☑ | ☑ | ☑ |
| 74 | `normalize_variant` | ☑ | ☑ | ☑ |
| 75 | `over_representation_analysis` | ☑ | ☑ | ☑ |
| 76 | `parse_gff3` | ☑ | ☑ | ☑ |
| 77 | `parse_vcf_variant` | ☑ | ☑ | ☑ |
| 78 | `pearson_correlation` | ☑ | ☑ | ☑ |
| 79 | `perform_pca` | ☑ | ☑ | ☑ |
| 80 | `predict_chromatin_state` | ☑ | ☑ | ☑ |
| 81 | `predict_gene_structure` | ☑ | ☑ | ☑ |
| 82 | `predict_genes` | ☑ | ☑ | ☑ |
| 83 | `predict_imprinted_genes` | ☑ | ☑ | ☑ |
| 84 | `predict_introns` | ☑ | ☑ | ☑ |
| 85 | `predict_pathogenicity` | ☑ | ☑ | ☑ |
| 86 | `predict_tf_binding_change` | ☑ | ☑ | ☑ |
| 87 | `predict_variant_effect` | ☑ | ☑ | ☑ |
| 88 | `quantile_normalize` | ☑ | ☑ | ☑ |
| 89 | `rna_reverse_complement` | ☑ | ☑ | ☑ |
| 90 | `rnaseq_quality_metrics` | ☑ | ☑ | ☑ |
| 91 | `segment_copy_number` | ☑ | ☑ | ☑ |
| 92 | `simulate_bisulfite_conversion` | ☑ | ☑ | ☑ |
| 93 | `site_accessibility` | ☑ | ☑ | ☑ |
| 94 | `titv_ratio` | ☑ | ☑ | ☑ |
| 95 | `to_gff3` | ☑ | ☑ | ☑ |
| 96 | `variant_statistics` | ☑ | ☑ | ☑ |
| 97 | `variants_to_vcf` | ☑ | ☑ | ☑ |

### Chromosome (Seqeron.Mcp.Chromosome) — 32 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `analyze_centromere` | ☑ | ☑ | ☑ |
| 2 | `analyze_karyotype` | ☑ | ☑ | ☑ |
| 3 | `analyze_scaffolds` | ☑ | ☑ | ☑ |
| 4 | `analyze_telomeres` | ☑ | ☑ | ☑ |
| 5 | `arm_ratio` | ☑ | ☑ | ☑ |
| 6 | `assembly_statistics` | ☑ | ☑ | ☑ |
| 7 | `assess_completeness` | ☑ | ☑ | ☑ |
| 8 | `au_n` | ☑ | ☑ | ☑ |
| 9 | `classify_chromosome_by_arm_ratio` | ☑ | ☑ | ☑ |
| 10 | `compare_assemblies` | ☑ | ☑ | ☑ |
| 11 | `detect_aneuploidy` | ☑ | ☑ | ☑ |
| 12 | `detect_ploidy` | ☑ | ☑ | ☑ |
| 13 | `detect_rearrangements` | ☑ | ☑ | ☑ |
| 14 | `estimate_cell_divisions_from_telomere_length` | ☑ | ☑ | ☑ |
| 15 | `estimate_completeness_from_kmers` | ☑ | ☑ | ☑ |
| 16 | `estimate_telomere_length_from_ts_ratio` | ☑ | ☑ | ☑ |
| 17 | `extract_contigs` | ☑ | ☑ | ☑ |
| 18 | `find_gaps` | ☑ | ☑ | ☑ |
| 19 | `find_heterochromatin_regions` | ☑ | ☑ | ☑ |
| 20 | `find_repetitive_regions` | ☑ | ☑ | ☑ |
| 21 | `find_suspicious_regions` | ☑ | ☑ | ☑ |
| 22 | `find_syntenic_blocks_assemblies` | ☑ | ☑ | ☑ |
| 23 | `find_synteny_blocks` | ☑ | ☑ | ☑ |
| 24 | `find_tandem_repeats` | ☑ | ☑ | ☑ |
| 25 | `gap_distribution` | ☑ | ☑ | ☑ |
| 26 | `identify_whole_chromosome_aneuploidy` | ☑ | ☑ | ☑ |
| 27 | `length_distribution` | ☑ | ☑ | ☑ |
| 28 | `local_quality` | ☑ | ☑ | ☑ |
| 29 | `nx_curve` | ☑ | ☑ | ☑ |
| 30 | `nx_statistics` | ☑ | ☑ | ☑ |
| 31 | `predict_g_bands` | ☑ | ☑ | ☑ |
| 32 | `repeat_content` | ☑ | ☑ | ☑ |

### Metagenomics (Seqeron.Mcp.Metagenomics) — 19 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `accessory_genes` | ☑ | ☑ | ☑ |
| 2 | `alpha_diversity` | ☑ | ☑ | ☑ |
| 3 | `beta_diversity` | ☑ | ☑ | ☑ |
| 4 | `bin_contigs` | ☑ | ☑ | ☑ |
| 5 | `build_kmer_database` | ☑ | ☑ | ☑ |
| 6 | `classify_reads` | ☑ | ☑ | ☑ |
| 7 | `cluster_genes` | ☑ | ☑ | ☑ |
| 8 | `construct_pangenome` | ☑ | ☑ | ☑ |
| 9 | `core_gene_clusters` | ☑ | ☑ | ☑ |
| 10 | `core_genome_alignment` | ☑ | ☑ | ☑ |
| 11 | `differential_abundance` | ☑ | ☑ | ☑ |
| 12 | `find_genome_specific_genes` | ☑ | ☑ | ☑ |
| 13 | `find_resistance_genes` | ☑ | ☑ | ☑ |
| 14 | `fit_heaps_law` | ☑ | ☑ | ☑ |
| 15 | `functional_diversity` | ☑ | ☑ | ☑ |
| 16 | `gene_presence_absence_matrix` | ☑ | ☑ | ☑ |
| 17 | `predict_functions` | ☑ | ☑ | ☑ |
| 18 | `select_phylogenetic_markers` | ☑ | ☑ | ☑ |
| 19 | `taxonomic_profile` | ☑ | ☑ | ☑ |

### MolTools (Seqeron.Mcp.MolTools) — 48 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `blunt_cutters` | ☑ | ☑ | ☑ |
| 2 | `build_codon_table` | ☑ | ☑ | ☑ |
| 3 | `cai_from_organism_table` | ☑ | ☑ | ☑ |
| 4 | `calculate_on_target_doench2014` | ☑ | ☑ | ☑ |
| 5 | `calculate_on_target_rule_set2` | ☑ | ☑ | ☑ |
| 6 | `codon_adaptation_index` | ☑ | ☑ | ☑ |
| 7 | `codon_minmax_profile` | ☑ | ☑ | ☑ |
| 8 | `codon_usage_statistics` | ☑ | ☑ | ☑ |
| 9 | `compare_codon_usage` | ☑ | ☑ | ☑ |
| 10 | `compatible_enzymes` | ☑ | ☑ | ☑ |
| 11 | `count_codons` | ☑ | ☑ | ☑ |
| 12 | `crispr_cfd_specificity_score` | ☑ | ☑ | ☑ |
| 13 | `crispr_mit_specificity_score` | ☑ | ☑ | ☑ |
| 14 | `crispr_specificity_score` | ☑ | ☑ | ☑ |
| 15 | `crispr_system_info` | ☑ | ☑ | ☑ |
| 16 | `design_guide_rnas` | ☑ | ☑ | ☑ |
| 17 | `digest_summary` | ☑ | ☑ | ☑ |
| 18 | `effective_number_of_codons` | ☑ | ☑ | ☑ |
| 19 | `enzymes_by_cut_length` | ☑ | ☑ | ☑ |
| 20 | `enzymes_compatible` | ☑ | ☑ | ☑ |
| 21 | `evaluate_guide_rna` | ☑ | ☑ | ☑ |
| 22 | `find_all_restriction_sites` | ☑ | ☑ | ☑ |
| 23 | `find_conserved_rare_codon_clusters` | ☑ | ☑ | ☑ |
| 24 | `find_off_targets` | ☑ | ☑ | ☑ |
| 25 | `find_pam_sites` | ☑ | ☑ | ☑ |
| 26 | `find_rare_codon_clusters` | ☑ | ☑ | ☑ |
| 27 | `find_rare_codons` | ☑ | ☑ | ☑ |
| 28 | `find_restriction_sites` | ☑ | ☑ | ☑ |
| 29 | `get_enzyme` | ☑ | ☑ | ☑ |
| 30 | `hairpin_potential` | ☑ | ☑ | ☑ |
| 31 | `longest_dinucleotide_repeat` | ☑ | ☑ | ☑ |
| 32 | `longest_homopolymer` | ☑ | ☑ | ☑ |
| 33 | `oligo_concentration_from_absorbance` | ☑ | ☑ | ☑ |
| 34 | `oligo_extinction_coefficient` | ☑ | ☑ | ☑ |
| 35 | `optimize_codons` | ☑ | ☑ | ☑ |
| 36 | `partial_digest` | ☑ | ☑ | ☑ |
| 37 | `primer_melting_temperature` | ☑ | ☑ | ☑ |
| 38 | `primer_melting_temperature_salt` | ☑ | ☑ | ☑ |
| 39 | `reduce_secondary_structure` | ☑ | ☑ | ☑ |
| 40 | `remove_restriction_sites` | ☑ | ☑ | ☑ |
| 41 | `restriction_digest` | ☑ | ☑ | ☑ |
| 42 | `restriction_map` | ☑ | ☑ | ☑ |
| 43 | `rscu` | ☑ | ☑ | ☑ |
| 44 | `sticky_cutters` | ☑ | ☑ | ☑ |
| 45 | `three_prime_stability` | ☑ | ☑ | ☑ |
| 46 | `find_pattern_occurrences` | ☑ | ☑ | ☑ |
| 47 | `find_sites_custom_enzyme` | ☑ | ☑ | ☑ |
| 48 | `restriction_digest_custom_enzyme` | ☑ | ☑ | ☑ |

### Phylogenetics (Seqeron.Mcp.Phylogenetics) — 13 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `bootstrap_support` | ☑ | ☑ | ☑ |
| 2 | `build_phylogenetic_tree` | ☑ | ☑ | ☑ |
| 3 | `build_tree_from_matrix` | ☑ | ☑ | ☑ |
| 4 | `distance_matrix` | ☑ | ☑ | ☑ |
| 5 | `mrca` | ☑ | ☑ | ☑ |
| 6 | `pairwise_distance` | ☑ | ☑ | ☑ |
| 7 | `parse_newick` | ☑ | ☑ | ☑ |
| 8 | `patristic_distance` | ☑ | ☑ | ☑ |
| 9 | `robinson_foulds_distance` | ☑ | ☑ | ☑ |
| 10 | `to_newick` | ☑ | ☑ | ☑ |
| 11 | `tree_depth` | ☑ | ☑ | ☑ |
| 12 | `tree_leaves` | ☑ | ☑ | ☑ |
| 13 | `tree_length` | ☑ | ☑ | ☑ |

### Population (Seqeron.Mcp.Population) — 18 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `allele_frequencies` | ☑ | ☑ | ☑ |
| 2 | `diversity_statistics` | ☑ | ☑ | ☑ |
| 3 | `estimate_ancestry` | ☑ | ☑ | ☑ |
| 4 | `f_statistics` | ☑ | ☑ | ☑ |
| 5 | `filter_variants_by_maf` | ☑ | ☑ | ☑ |
| 6 | `fst` | ☑ | ☑ | ☑ |
| 7 | `haplotype_blocks` | ☑ | ☑ | ☑ |
| 8 | `hardy_weinberg_test` | ☑ | ☑ | ☑ |
| 9 | `inbreeding_from_roh` | ☑ | ☑ | ☑ |
| 10 | `integrated_haplotype_score` | ☑ | ☑ | ☑ |
| 11 | `linkage_disequilibrium` | ☑ | ☑ | ☑ |
| 12 | `minor_allele_frequency` | ☑ | ☑ | ☑ |
| 13 | `nucleotide_diversity` | ☑ | ☑ | ☑ |
| 14 | `pairwise_fst` | ☑ | ☑ | ☑ |
| 15 | `runs_of_homozygosity` | ☑ | ☑ | ☑ |
| 16 | `scan_selection_signals` | ☑ | ☑ | ☑ |
| 17 | `tajimas_d` | ☑ | ☑ | ☑ |
| 18 | `wattersons_theta` | ☑ | ☑ | ☑ |
