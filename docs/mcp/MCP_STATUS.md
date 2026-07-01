# MCP Completion — Source of Truth & Campaign Ledger

> Auto-generated from the actual `[McpServerTool]` bindings on disk (11 servers) and the
> presence of per-tool docs/tests. This SUPERSEDES the stale `docs/mcp-plan.md` and
> `docs/mcp-checklist.md` (v4), which described a 12-server/241-tool design that was never built.

## Definition of Done (per tool)

A tool is ☑ **Done** only when ALL hold:
1. **Binding** — `[McpServerTool(Name=…, Title=…, ReadOnly=…)]` explicit; structured `record` result;
   validates input; calls the real `Seqeron.Genomics` method (no re-implemented logic); error mapping.
2. **Tests** — ≥2 NUnit tests (`{Tool}_Schema_ValidatesCorrectly` + `{Tool}_Binding_InvokesSuccessfully`),
   one file per tool, evidence-based (would fail on a wrong impl), full suite green.
3. **Docs** — `{tool}.md` + `{tool}.mcp.json` under `docs/mcp/tools/<server>/`.

Server is Done when all its tools are Done + `Program.cs` wired + builds + integration test + `README.md`.

## Status legend
- **B** = Binding gold-standard · **T** = Tests · **D** = Docs (.md + .mcp.json)

## Roll-up

| Metric | Done | Total |
|---|---:|---:|
| Tools with gold-standard binding | 381 | 427 |
| Tools with tests (server-level) | 89 | 427 |
| Tools with full docs | 90 | 427 |

## Per-server summary

| Server | Project | Tools | B | T | D | Status |
|---|---|---:|:--:|:--:|:--:|---|
| Core | SuffixTree.Mcp.Core | 12 | 12/12 | 12/12 | 12/12 | ✅ done |
| Sequence | Seqeron.Mcp.Sequence | 35 | 35/35 | 35/35 | 35/35 | ✅ done |
| Parsers | Seqeron.Mcp.Parsers | 41 | 41/41 | 41/41 | 41/41 | ✅ done |
| Alignment | Seqeron.Mcp.Alignment | 22 | 22/22 | 0/22 | 0/22 | 🔧 needs work |
| Analysis | Seqeron.Mcp.Analysis | 91 | 91/91 | 0/91 | 1/91 | 🔧 needs work |
| Annotation | Seqeron.Mcp.Annotation | 97 | 97/97 | 0/97 | 0/97 | 🔧 needs work |
| Chromosome | Seqeron.Mcp.Chromosome | 32 | 32/32 | 0/32 | 0/32 | 🔧 needs work |
| Metagenomics | Seqeron.Mcp.Metagenomics | 19 | 19/19 | 0/19 | 0/19 | 🔧 needs work |
| MolTools | Seqeron.Mcp.MolTools | 47 | 1/47 | 1/47 | 1/47 | 🔧 needs work |
| Phylogenetics | Seqeron.Mcp.Phylogenetics | 13 | 13/13 | 0/13 | 0/13 | 🔧 needs work |
| Population | Seqeron.Mcp.Population | 18 | 18/18 | 0/18 | 0/18 | 🔧 needs work |

## Per-tool ledger

### Core (SuffixTree.Mcp.Core) — 12 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `calculate_similarity` | ☑ | ☑ | ☑ |
| 2 | `count_approximate_occurrences` | ☑ | ☑ | ☑ |
| 3 | `edit_distance` | ☑ | ☑ | ☑ |
| 4 | `find_longest_common_region` | ☑ | ☑ | ☑ |
| 5 | `find_longest_repeat` | ☑ | ☑ | ☑ |
| 6 | `hamming_distance` | ☑ | ☑ | ☑ |
| 7 | `suffix_tree_contains` | ☑ | ☑ | ☑ |
| 8 | `suffix_tree_count` | ☑ | ☑ | ☑ |
| 9 | `suffix_tree_find_all` | ☑ | ☑ | ☑ |
| 10 | `suffix_tree_lcs` | ☑ | ☑ | ☑ |
| 11 | `suffix_tree_lrs` | ☑ | ☑ | ☑ |
| 12 | `suffix_tree_stats` | ☑ | ☑ | ☑ |

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

### Alignment (Seqeron.Mcp.Alignment) — 22 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `alignment_statistics` | ☑ | ☐ | ☐ |
| 2 | `assemble_de_bruijn` | ☑ | ☐ | ☐ |
| 3 | `assemble_olc` | ☑ | ☐ | ☐ |
| 4 | `assembly_stats` | ☑ | ☐ | ☐ |
| 5 | `calculate_coverage` | ☑ | ☐ | ☐ |
| 6 | `compute_consensus` | ☑ | ☐ | ☐ |
| 7 | `error_correct_reads` | ☑ | ☐ | ☐ |
| 8 | `find_all_overlaps` | ☑ | ☐ | ☐ |
| 9 | `find_best_match` | ☑ | ☐ | ☐ |
| 10 | `find_overlap` | ☑ | ☐ | ☐ |
| 11 | `find_with_edits` | ☑ | ☐ | ☐ |
| 12 | `find_with_mismatches` | ☑ | ☐ | ☐ |
| 13 | `format_alignment` | ☑ | ☐ | ☐ |
| 14 | `frequent_kmers_with_mismatches` | ☑ | ☐ | ☐ |
| 15 | `global_align` | ☑ | ☐ | ☐ |
| 16 | `local_align` | ☑ | ☐ | ☐ |
| 17 | `merge_contigs` | ☑ | ☐ | ☐ |
| 18 | `multiple_align` | ☑ | ☐ | ☐ |
| 19 | `quality_trim_reads` | ☑ | ☐ | ☐ |
| 20 | `scaffold_contigs` | ☑ | ☐ | ☐ |
| 21 | `semi_global_align` | ☑ | ☐ | ☐ |
| 22 | `sequence_identity` | ☑ | ☐ | ☐ |

### Analysis (Seqeron.Mcp.Analysis) — 91 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `analyze_gc_content` | ☑ | ☐ | ☐ |
| 2 | `analyze_kmers` | ☑ | ☐ | ☐ |
| 3 | `at_skew` | ☑ | ☐ | ☐ |
| 4 | `base_pair_type` | ☑ | ☐ | ☐ |
| 5 | `bulge_loop_energy` | ☑ | ☐ | ☐ |
| 6 | `calculate_ani` | ☑ | ☐ | ☐ |
| 7 | `can_pair` | ☑ | ☐ | ☐ |
| 8 | `codon_frequencies` | ☑ | ☐ | ☐ |
| 9 | `compare_genomes` | ☑ | ☐ | ☐ |
| 10 | `compression_ratio` | ☑ | ☐ | ☐ |
| 11 | `count_kmers` | ☑ | ☐ | ☐ |
| 12 | `count_kmers_both_strands` | ☑ | ☐ | ☐ |
| 13 | `create_pwm` | ☑ | ☐ | ☐ |
| 14 | `cumulative_gc_skew` | ☑ | ☐ | ☐ |
| 15 | `dangling_end_energy` | ☑ | ☐ | ☐ |
| 16 | `detect_pseudoknots` | ☑ | ☐ | ☐ |
| 17 | `detect_rearrangements` | ☑ | ☐ | ☐ |
| 18 | `dinucleotide_frequencies` | ☑ | ☐ | ☐ |
| 19 | `dinucleotide_ratios` | ☑ | ☐ | ☐ |
| 20 | `discover_motifs` | ☑ | ☐ | ☐ |
| 21 | `disorder_propensity` | ☑ | ☐ | ☐ |
| 22 | `dust_score` | ☑ | ☐ | ☐ |
| 23 | `entropy_profile` | ☑ | ☐ | ☐ |
| 24 | `find_clumps` | ☑ | ☐ | ☐ |
| 25 | `find_common_regions` | ☑ | ☐ | ☐ |
| 26 | `find_conserved_clusters` | ☑ | ☐ | ☐ |
| 27 | `find_degenerate_motif` | ☑ | ☐ | ☐ |
| 28 | `find_direct_repeats` | ☑ | ☐ | ☐ |
| 29 | `find_exact_motif` | ☑ | ☐ | ☐ |
| 30 | `find_inverted_repeats` | ☑ | ☐ | ☐ |
| 31 | `find_known_motifs` | ☑ | ☐ | ☐ |
| 32 | `find_low_complexity_regions` | ☑ | ☐ | ☐ |
| 33 | `find_microsatellites` | ☑ | ☐ | ☐ |
| 34 | `find_motif` | ☑ | ☐ | ☐ |
| 35 | `find_motif_by_pattern` | ☑ | ☐ | ☐ |
| 36 | `find_motif_by_prosite` | ☑ | ☐ | ☐ |
| 37 | `find_open_reading_frames` | ☑ | ☐ | ☐ |
| 38 | `find_orthologs` | ☑ | ☐ | ☐ |
| 39 | `find_palindromes` | ☑ | ☐ | ☐ |
| 40 | `find_protein_domains` | ☑ | ☐ | ☐ |
| 41 | `find_protein_low_complexity_regions` | ☑ | ☐ | ☐ |
| 42 | `find_protein_motifs` | ☑ | ☐ | ☐ |
| 43 | `find_reciprocal_best_hits` | ☑ | ☐ | ☐ |
| 44 | `find_regulatory_elements` | ☑ | ☐ | ☐ |
| 45 | `find_repeats` | ☑ | ☐ | ☐ |
| 46 | `find_rna_inverted_repeats` | ☑ | ☐ | ☐ |
| 47 | `find_shared_motifs` | ☑ | ☐ | ☐ |
| 48 | `find_stem_loops` | ☑ | ☐ | ☐ |
| 49 | `find_syntenic_blocks` | ☑ | ☐ | ☐ |
| 50 | `find_tandem_repeats` | ☑ | ☐ | ☐ |
| 51 | `flush_coaxial_stacking` | ☑ | ☐ | ☐ |
| 52 | `gc_content_profile` | ☑ | ☐ | ☐ |
| 53 | `gc_skew` | ☑ | ☐ | ☐ |
| 54 | `generate_all_kmers` | ☑ | ☐ | ☐ |
| 55 | `generate_consensus` | ☑ | ☐ | ☐ |
| 56 | `generate_dot_plot` | ☑ | ☐ | ☐ |
| 57 | `hairpin_loop_energy` | ☑ | ☐ | ☐ |
| 58 | `hydrophobicity_profile` | ☑ | ☐ | ☐ |
| 59 | `internal_loop_energy` | ☑ | ☐ | ☐ |
| 60 | `is_disorder_promoting` | ☑ | ☐ | ☐ |
| 61 | `kmer_distance` | ☑ | ☐ | ☑ |
| 62 | `kmer_frequencies` | ☑ | ☐ | ☐ |
| 63 | `kmer_positions` | ☑ | ☐ | ☐ |
| 64 | `kmer_spectrum` | ☑ | ☐ | ☐ |
| 65 | `kmers_with_min_count` | ☑ | ☐ | ☐ |
| 66 | `mask_low_complexity` | ☑ | ☐ | ☐ |
| 67 | `minimum_free_energy` | ☑ | ☐ | ☐ |
| 68 | `mismatch_coaxial_stacking` | ☑ | ☐ | ☐ |
| 69 | `most_frequent_kmers` | ☑ | ☐ | ☐ |
| 70 | `multibranch_loop_energy` | ☑ | ☐ | ☐ |
| 71 | `parse_dot_bracket` | ☑ | ☐ | ☐ |
| 72 | `predict_chou_fasman` | ☑ | ☐ | ☐ |
| 73 | `predict_coiled_coils` | ☑ | ☐ | ☐ |
| 74 | `predict_disorder` | ☑ | ☐ | ☐ |
| 75 | `predict_low_complexity_seg` | ☑ | ☐ | ☐ |
| 76 | `predict_morfs` | ☑ | ☐ | ☐ |
| 77 | `predict_replication_origin` | ☑ | ☐ | ☐ |
| 78 | `predict_rna_structure` | ☑ | ☐ | ☐ |
| 79 | `predict_signal_peptide` | ☑ | ☐ | ☐ |
| 80 | `predict_transmembrane_helices` | ☑ | ☐ | ☐ |
| 81 | `prosite_to_regex` | ☑ | ☐ | ☐ |
| 82 | `reversal_distance` | ☑ | ☐ | ☐ |
| 83 | `rna_complement_base` | ☑ | ☐ | ☐ |
| 84 | `scan_with_pwm` | ☑ | ☐ | ☐ |
| 85 | `stem_energy` | ☑ | ☐ | ☐ |
| 86 | `tandem_repeat_summary` | ☑ | ☐ | ☐ |
| 87 | `terminal_mismatch_energy` | ☑ | ☐ | ☐ |
| 88 | `unique_kmers` | ☑ | ☐ | ☐ |
| 89 | `validate_dot_bracket` | ☑ | ☐ | ☐ |
| 90 | `windowed_complexity` | ☑ | ☐ | ☐ |
| 91 | `windowed_gc_skew` | ☑ | ☐ | ☐ |

### Annotation (Seqeron.Mcp.Annotation) — 97 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `align_mirna_to_target` | ☑ | ☐ | ☐ |
| 2 | `analyze_target_context` | ☑ | ☐ | ☐ |
| 3 | `annotate_histone_modifications` | ☑ | ☐ | ☐ |
| 4 | `annotate_regulatory_elements` | ☑ | ☐ | ☐ |
| 5 | `annotate_svs` | ☑ | ☐ | ☐ |
| 6 | `annotate_variant_on_transcripts` | ☑ | ☐ | ☐ |
| 7 | `annotate_variants` | ☑ | ☐ | ☐ |
| 8 | `assemble_breakpoint_sequence` | ☑ | ☐ | ☐ |
| 9 | `build_coexpression_network` | ☑ | ☐ | ☐ |
| 10 | `calculate_conservation` | ☑ | ☐ | ☐ |
| 11 | `calculate_tpm` | ☑ | ☐ | ☐ |
| 12 | `call_variants` | ☑ | ☐ | ☐ |
| 13 | `call_variants_from_alignment` | ☑ | ☐ | ☐ |
| 14 | `can_pair` | ☑ | ☐ | ☐ |
| 15 | `classify_mutation` | ☑ | ☐ | ☐ |
| 16 | `classify_variant` | ☑ | ☐ | ☐ |
| 17 | `cluster_discordant_pairs` | ☑ | ☐ | ☐ |
| 18 | `cluster_genes_by_expression` | ☑ | ☐ | ☐ |
| 19 | `cluster_split_reads` | ☑ | ☐ | ☐ |
| 20 | `coding_potential` | ☑ | ☐ | ☐ |
| 21 | `codon_usage` | ☑ | ☐ | ☐ |
| 22 | `compare_seed_regions` | ☑ | ☐ | ☐ |
| 23 | `cpg_observed_expected` | ☑ | ☐ | ☐ |
| 24 | `create_mirna` | ☑ | ☐ | ☐ |
| 25 | `detect_alternative_splicing` | ☑ | ☐ | ☐ |
| 26 | `detect_differential_splicing` | ☑ | ☐ | ☐ |
| 27 | `detect_isoform_switching` | ☑ | ☐ | ☐ |
| 28 | `differential_expression` | ☑ | ☐ | ☐ |
| 29 | `enrichment_score` | ☑ | ☐ | ☐ |
| 30 | `epigenetic_age` | ☑ | ☐ | ☐ |
| 31 | `filter_svs` | ☑ | ☐ | ☐ |
| 32 | `find_acceptor_sites` | ☑ | ☐ | ☐ |
| 33 | `find_accessible_regions` | ☑ | ☐ | ☐ |
| 34 | `find_branch_points` | ☑ | ☐ | ☐ |
| 35 | `find_conserved_elements` | ☑ | ☐ | ☐ |
| 36 | `find_cpg_islands` | ☑ | ☐ | ☐ |
| 37 | `find_cpg_sites` | ☑ | ☐ | ☐ |
| 38 | `find_deletions` | ☑ | ☐ | ☐ |
| 39 | `find_discordant_pairs` | ☑ | ☐ | ☐ |
| 40 | `find_dmrs` | ☑ | ☐ | ☐ |
| 41 | `find_dominant_isoforms` | ☑ | ☐ | ☐ |
| 42 | `find_donor_sites` | ☑ | ☐ | ☐ |
| 43 | `find_indels` | ☑ | ☐ | ☐ |
| 44 | `find_insertions` | ☑ | ☐ | ☐ |
| 45 | `find_methylation_sites` | ☑ | ☐ | ☐ |
| 46 | `find_microhomology` | ☑ | ☐ | ☐ |
| 47 | `find_mirna_target_sites` | ☑ | ☐ | ☐ |
| 48 | `find_orfs` | ☑ | ☐ | ☐ |
| 49 | `find_pre_mirna_hairpins` | ☑ | ☐ | ☐ |
| 50 | `find_promoter_motifs` | ☑ | ☐ | ☐ |
| 51 | `find_repetitive_elements` | ☑ | ☐ | ☐ |
| 52 | `find_retained_intron_candidates` | ☑ | ☐ | ☐ |
| 53 | `find_ribosome_binding_sites` | ☑ | ☐ | ☐ |
| 54 | `find_similar_mirnas` | ☑ | ☐ | ☐ |
| 55 | `find_skipped_exon_events` | ☑ | ☐ | ☐ |
| 56 | `find_snps` | ☑ | ☐ | ☐ |
| 57 | `find_snps_direct` | ☑ | ☐ | ☐ |
| 58 | `find_split_reads` | ☑ | ☐ | ☐ |
| 59 | `format_vcf_info` | ☑ | ☐ | ☐ |
| 60 | `generate_seed_variants` | ☑ | ☐ | ☐ |
| 61 | `genotype_sv` | ☑ | ☐ | ☐ |
| 62 | `group_by_seed_family` | ☑ | ☐ | ☐ |
| 63 | `identify_cnvs` | ☑ | ☐ | ☐ |
| 64 | `impact_level` | ☑ | ☐ | ☐ |
| 65 | `is_within_coding_region` | ☑ | ☐ | ☐ |
| 66 | `is_wobble_pair` | ☑ | ☐ | ☐ |
| 67 | `log2_transform` | ☑ | ☐ | ☐ |
| 68 | `longest_orfs_per_frame` | ☑ | ☐ | ☐ |
| 69 | `maxent_score` | ☑ | ☐ | ☐ |
| 70 | `merge_overlapping_svs` | ☑ | ☐ | ☐ |
| 71 | `methylation_from_bisulfite` | ☑ | ☐ | ☐ |
| 72 | `methylation_profile` | ☑ | ☐ | ☐ |
| 73 | `mirna_seed_sequence` | ☑ | ☐ | ☐ |
| 74 | `normalize_variant` | ☑ | ☐ | ☐ |
| 75 | `over_representation_analysis` | ☑ | ☐ | ☐ |
| 76 | `parse_gff3` | ☑ | ☐ | ☐ |
| 77 | `parse_vcf_variant` | ☑ | ☐ | ☐ |
| 78 | `pearson_correlation` | ☑ | ☐ | ☐ |
| 79 | `perform_pca` | ☑ | ☐ | ☐ |
| 80 | `predict_chromatin_state` | ☑ | ☐ | ☐ |
| 81 | `predict_gene_structure` | ☑ | ☐ | ☐ |
| 82 | `predict_genes` | ☑ | ☐ | ☐ |
| 83 | `predict_imprinted_genes` | ☑ | ☐ | ☐ |
| 84 | `predict_introns` | ☑ | ☐ | ☐ |
| 85 | `predict_pathogenicity` | ☑ | ☐ | ☐ |
| 86 | `predict_tf_binding_change` | ☑ | ☐ | ☐ |
| 87 | `predict_variant_effect` | ☑ | ☐ | ☐ |
| 88 | `quantile_normalize` | ☑ | ☐ | ☐ |
| 89 | `rna_reverse_complement` | ☑ | ☐ | ☐ |
| 90 | `rnaseq_quality_metrics` | ☑ | ☐ | ☐ |
| 91 | `segment_copy_number` | ☑ | ☐ | ☐ |
| 92 | `simulate_bisulfite_conversion` | ☑ | ☐ | ☐ |
| 93 | `site_accessibility` | ☑ | ☐ | ☐ |
| 94 | `titv_ratio` | ☑ | ☐ | ☐ |
| 95 | `to_gff3` | ☑ | ☐ | ☐ |
| 96 | `variant_statistics` | ☑ | ☐ | ☐ |
| 97 | `variants_to_vcf` | ☑ | ☐ | ☐ |

### Chromosome (Seqeron.Mcp.Chromosome) — 32 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `analyze_centromere` | ☑ | ☐ | ☐ |
| 2 | `analyze_karyotype` | ☑ | ☐ | ☐ |
| 3 | `analyze_scaffolds` | ☑ | ☐ | ☐ |
| 4 | `analyze_telomeres` | ☑ | ☐ | ☐ |
| 5 | `arm_ratio` | ☑ | ☐ | ☐ |
| 6 | `assembly_statistics` | ☑ | ☐ | ☐ |
| 7 | `assess_completeness` | ☑ | ☐ | ☐ |
| 8 | `au_n` | ☑ | ☐ | ☐ |
| 9 | `classify_chromosome_by_arm_ratio` | ☑ | ☐ | ☐ |
| 10 | `compare_assemblies` | ☑ | ☐ | ☐ |
| 11 | `detect_aneuploidy` | ☑ | ☐ | ☐ |
| 12 | `detect_ploidy` | ☑ | ☐ | ☐ |
| 13 | `detect_rearrangements` | ☑ | ☐ | ☐ |
| 14 | `estimate_cell_divisions_from_telomere_length` | ☑ | ☐ | ☐ |
| 15 | `estimate_completeness_from_kmers` | ☑ | ☐ | ☐ |
| 16 | `estimate_telomere_length_from_ts_ratio` | ☑ | ☐ | ☐ |
| 17 | `extract_contigs` | ☑ | ☐ | ☐ |
| 18 | `find_gaps` | ☑ | ☐ | ☐ |
| 19 | `find_heterochromatin_regions` | ☑ | ☐ | ☐ |
| 20 | `find_repetitive_regions` | ☑ | ☐ | ☐ |
| 21 | `find_suspicious_regions` | ☑ | ☐ | ☐ |
| 22 | `find_syntenic_blocks_assemblies` | ☑ | ☐ | ☐ |
| 23 | `find_synteny_blocks` | ☑ | ☐ | ☐ |
| 24 | `find_tandem_repeats` | ☑ | ☐ | ☐ |
| 25 | `gap_distribution` | ☑ | ☐ | ☐ |
| 26 | `identify_whole_chromosome_aneuploidy` | ☑ | ☐ | ☐ |
| 27 | `length_distribution` | ☑ | ☐ | ☐ |
| 28 | `local_quality` | ☑ | ☐ | ☐ |
| 29 | `nx_curve` | ☑ | ☐ | ☐ |
| 30 | `nx_statistics` | ☑ | ☐ | ☐ |
| 31 | `predict_g_bands` | ☑ | ☐ | ☐ |
| 32 | `repeat_content` | ☑ | ☐ | ☐ |

### Metagenomics (Seqeron.Mcp.Metagenomics) — 19 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `accessory_genes` | ☑ | ☐ | ☐ |
| 2 | `alpha_diversity` | ☑ | ☐ | ☐ |
| 3 | `beta_diversity` | ☑ | ☐ | ☐ |
| 4 | `bin_contigs` | ☑ | ☐ | ☐ |
| 5 | `build_kmer_database` | ☑ | ☐ | ☐ |
| 6 | `classify_reads` | ☑ | ☐ | ☐ |
| 7 | `cluster_genes` | ☑ | ☐ | ☐ |
| 8 | `construct_pangenome` | ☑ | ☐ | ☐ |
| 9 | `core_gene_clusters` | ☑ | ☐ | ☐ |
| 10 | `core_genome_alignment` | ☑ | ☐ | ☐ |
| 11 | `differential_abundance` | ☑ | ☐ | ☐ |
| 12 | `find_genome_specific_genes` | ☑ | ☐ | ☐ |
| 13 | `find_resistance_genes` | ☑ | ☐ | ☐ |
| 14 | `fit_heaps_law` | ☑ | ☐ | ☐ |
| 15 | `functional_diversity` | ☑ | ☐ | ☐ |
| 16 | `gene_presence_absence_matrix` | ☑ | ☐ | ☐ |
| 17 | `predict_functions` | ☑ | ☐ | ☐ |
| 18 | `select_phylogenetic_markers` | ☑ | ☐ | ☐ |
| 19 | `taxonomic_profile` | ☑ | ☐ | ☐ |

### MolTools (Seqeron.Mcp.MolTools) — 47 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `analyze_oligo` | ☑ | ☑ | ☑ |
| 2 | `blunt_cutters` ⚠️(bare) | ☐ | ☐ | ☐ |
| 3 | `build_codon_table` ⚠️(bare) | ☐ | ☐ | ☐ |
| 4 | `cai_from_organism_table` ⚠️(bare) | ☐ | ☐ | ☐ |
| 5 | `codon_adaptation_index` ⚠️(bare) | ☐ | ☐ | ☐ |
| 6 | `codon_usage_statistics` ⚠️(bare) | ☐ | ☐ | ☐ |
| 7 | `compare_codon_usage` ⚠️(bare) | ☐ | ☐ | ☐ |
| 8 | `compatible_enzymes` ⚠️(bare) | ☐ | ☐ | ☐ |
| 9 | `count_codons` ⚠️(bare) | ☐ | ☐ | ☐ |
| 10 | `crispr_specificity_score` ⚠️(bare) | ☐ | ☐ | ☐ |
| 11 | `crispr_system_info` ⚠️(bare) | ☐ | ☐ | ☐ |
| 12 | `design_antisense_probes` ⚠️(bare) | ☐ | ☐ | ☐ |
| 13 | `design_guide_rnas` ⚠️(bare) | ☐ | ☐ | ☐ |
| 14 | `design_molecular_beacon` ⚠️(bare) | ☐ | ☐ | ☐ |
| 15 | `design_primers` ⚠️(bare) | ☐ | ☐ | ☐ |
| 16 | `design_probes` ⚠️(bare) | ☐ | ☐ | ☐ |
| 17 | `design_tiling_probes` ⚠️(bare) | ☐ | ☐ | ☐ |
| 18 | `digest_summary` ⚠️(bare) | ☐ | ☐ | ☐ |
| 19 | `effective_number_of_codons` ⚠️(bare) | ☐ | ☐ | ☐ |
| 20 | `enzymes_by_cut_length` ⚠️(bare) | ☐ | ☐ | ☐ |
| 21 | `enzymes_compatible` ⚠️(bare) | ☐ | ☐ | ☐ |
| 22 | `evaluate_guide_rna` ⚠️(bare) | ☐ | ☐ | ☐ |
| 23 | `evaluate_primer` ⚠️(bare) | ☐ | ☐ | ☐ |
| 24 | `find_all_restriction_sites` ⚠️(bare) | ☐ | ☐ | ☐ |
| 25 | `find_off_targets` ⚠️(bare) | ☐ | ☐ | ☐ |
| 26 | `find_pam_sites` ⚠️(bare) | ☐ | ☐ | ☐ |
| 27 | `find_rare_codons` ⚠️(bare) | ☐ | ☐ | ☐ |
| 28 | `find_restriction_sites` ⚠️(bare) | ☐ | ☐ | ☐ |
| 29 | `generate_primer_candidates` ⚠️(bare) | ☐ | ☐ | ☐ |
| 30 | `get_enzyme` ⚠️(bare) | ☐ | ☐ | ☐ |
| 31 | `hairpin_potential` ⚠️(bare) | ☐ | ☐ | ☐ |
| 32 | `longest_dinucleotide_repeat` ⚠️(bare) | ☐ | ☐ | ☐ |
| 33 | `longest_homopolymer` ⚠️(bare) | ☐ | ☐ | ☐ |
| 34 | `oligo_concentration_from_absorbance` ⚠️(bare) | ☐ | ☐ | ☐ |
| 35 | `oligo_extinction_coefficient` ⚠️(bare) | ☐ | ☐ | ☐ |
| 36 | `optimize_codons` ⚠️(bare) | ☐ | ☐ | ☐ |
| 37 | `primer_dimer` ⚠️(bare) | ☐ | ☐ | ☐ |
| 38 | `primer_melting_temperature` ⚠️(bare) | ☐ | ☐ | ☐ |
| 39 | `primer_melting_temperature_salt` ⚠️(bare) | ☐ | ☐ | ☐ |
| 40 | `reduce_secondary_structure` ⚠️(bare) | ☐ | ☐ | ☐ |
| 41 | `remove_restriction_sites` ⚠️(bare) | ☐ | ☐ | ☐ |
| 42 | `restriction_digest` ⚠️(bare) | ☐ | ☐ | ☐ |
| 43 | `restriction_map` ⚠️(bare) | ☐ | ☐ | ☐ |
| 44 | `rscu` ⚠️(bare) | ☐ | ☐ | ☐ |
| 45 | `sticky_cutters` ⚠️(bare) | ☐ | ☐ | ☐ |
| 46 | `three_prime_stability` ⚠️(bare) | ☐ | ☐ | ☐ |
| 47 | `validate_probe` ⚠️(bare) | ☐ | ☐ | ☐ |

### Phylogenetics (Seqeron.Mcp.Phylogenetics) — 13 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `bootstrap_support` | ☑ | ☐ | ☐ |
| 2 | `build_phylogenetic_tree` | ☑ | ☐ | ☐ |
| 3 | `build_tree_from_matrix` | ☑ | ☐ | ☐ |
| 4 | `distance_matrix` | ☑ | ☐ | ☐ |
| 5 | `mrca` | ☑ | ☐ | ☐ |
| 6 | `pairwise_distance` | ☑ | ☐ | ☐ |
| 7 | `parse_newick` | ☑ | ☐ | ☐ |
| 8 | `patristic_distance` | ☑ | ☐ | ☐ |
| 9 | `robinson_foulds_distance` | ☑ | ☐ | ☐ |
| 10 | `to_newick` | ☑ | ☐ | ☐ |
| 11 | `tree_depth` | ☑ | ☐ | ☐ |
| 12 | `tree_leaves` | ☑ | ☐ | ☐ |
| 13 | `tree_length` | ☑ | ☐ | ☐ |

### Population (Seqeron.Mcp.Population) — 18 tools

| # | tool | B | T | D |
|---:|---|:--:|:--:|:--:|
| 1 | `allele_frequencies` | ☑ | ☐ | ☐ |
| 2 | `diversity_statistics` | ☑ | ☐ | ☐ |
| 3 | `estimate_ancestry` | ☑ | ☐ | ☐ |
| 4 | `f_statistics` | ☑ | ☐ | ☐ |
| 5 | `filter_variants_by_maf` | ☑ | ☐ | ☐ |
| 6 | `fst` | ☑ | ☐ | ☐ |
| 7 | `haplotype_blocks` | ☑ | ☐ | ☐ |
| 8 | `hardy_weinberg_test` | ☑ | ☐ | ☐ |
| 9 | `inbreeding_from_roh` | ☑ | ☐ | ☐ |
| 10 | `integrated_haplotype_score` | ☑ | ☐ | ☐ |
| 11 | `linkage_disequilibrium` | ☑ | ☐ | ☐ |
| 12 | `minor_allele_frequency` | ☑ | ☐ | ☐ |
| 13 | `nucleotide_diversity` | ☑ | ☐ | ☐ |
| 14 | `pairwise_fst` | ☑ | ☐ | ☐ |
| 15 | `runs_of_homozygosity` | ☑ | ☐ | ☐ |
| 16 | `scan_selection_signals` | ☑ | ☐ | ☐ |
| 17 | `tajimas_d` | ☑ | ☐ | ☐ |
| 18 | `wattersons_theta` | ☑ | ☐ | ☐ |
