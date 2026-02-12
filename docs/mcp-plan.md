# MCP Implementation Plan v4

## 1. Summary

### Migration v3 → v4

| Aspect | v3 (Single Server) | v4 (Multi-Server) |
|--------|-------------------|-------------------|
| Architecture | 1 monolithic server | **12 focused servers** |
| Total tools | 254 | **241** (after exclusions) |
| Context per server | ~24% | **2-8%** |
| User choice | All or nothing | **Pick needed servers** |

### Design Rationale
Based on [MCP best practices](https://www.klavis.ai/blog/less-is-more-mcp-design-patterns-for-ai-agents):
- Tool definitions consume 5-7% of context window before any prompt
- Multi-server approach keeps each server's context footprint minimal
- Users connect only servers relevant to their workflow

---

## 2. Inputs & Inventory Source

### Base Documents
- **MCP-Server-Plan.md** — architectural plan, server definitions, tool distribution
- **MCP-Methods-Audit.md** — complete audit of 277 public static methods
- **Tool Inventory v4** — 241 tools distributed across 12 servers

### HasDocs Determination Rules
A tool has `HasDocs=true` if its underlying method has:
1. XML documentation (`<summary>`, `<param>`, `<returns>` tags) in source, OR
2. Algorithm documentation in `docs/algorithms/*.md`, OR
3. Entry in official library documentation

### DocRef Format
- XML: `{ClassName}.cs:L{line}#xml`
- Markdown: `docs/algorithms/{path}/{file}.md`

---

## 3. Server Architecture

### 12 MCP Servers

| # | Server | Tools | Focus | Context |
|---|--------|-------|-------|---------|
| 1 | **Core** | 12 | Suffix tree operations, similarity | ~3% |
| 2 | **Sequence** | 35 | DNA/RNA/Protein analysis, k-mers | ~7% |
| 3 | **Parsers** | 45 | File format I/O (FASTA, FASTQ, VCF, GFF, BED) | ~8% |
| 4 | **Alignment** | 15 | Sequence alignment, motifs, repeats | ~4% |
| 5 | **Variants** | 22 | SNP/indel calling, annotation, quality | ~5% |
| 6 | **MolBio** | 38 | Primers, probes, restriction, CRISPR, codons | ~7% |
| 7 | **Assembly** | 13 | Genome assembly, pan-genome | ~3% |
| 8 | **Phylogenetics** | 10 | Tree building, comparative genomics | ~3% |
| 9 | **Population** | 15 | Population genetics, metagenomics | ~4% |
| 10 | **Epigenetics** | 14 | Methylation, miRNA analysis | ~4% |
| 11 | **Structure** | 14 | RNA/Protein structure prediction | ~4% |
| 12 | **Annotation** | 8 | Gene finding, chromosome analysis | ~2% |
| | **TOTAL** | **241** | | |

### Separation Criteria
1. **Domain cohesion** — tools in same server share conceptual domain
2. **Context budget** — each server ≤8% context window
3. **Independent deployment** — servers can be used standalone
4. **Minimal cross-dependencies** — avoid tools requiring multiple servers

---

## 4. Tool Inventory

### Inventory Source: v4 Tool List (241 tools)

Legend:
- ✓ = HasDocs=true (XML or Markdown documentation exists)
- ○ = HasDocs=false (documentation TBD)

---

### Server 1: Core (12 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `suffix_tree_contains` | Check pattern existence | ✓ | SuffixTree.cs#xml | SuffixTree.Contains | 1.0.0 | stable |
| 2 | `suffix_tree_count` | Count pattern occurrences | ✓ | SuffixTree.cs#xml | SuffixTree.CountOccurrences | 1.0.0 | stable |
| 3 | `suffix_tree_find_all` | Find all occurrence positions | ✓ | SuffixTree.cs#xml | SuffixTree.FindAllOccurrences | 1.0.0 | stable |
| 4 | `suffix_tree_lrs` | Longest repeated substring | ✓ | SuffixTree.cs#xml | SuffixTree.LongestRepeatedSubstring | 1.0.0 | stable |
| 5 | `suffix_tree_lcs` | Longest common substring | ✓ | SuffixTree.cs#xml | SuffixTree.LongestCommonSubstring | 1.0.0 | stable |
| 6 | `suffix_tree_stats` | Tree statistics | ✓ | SuffixTree.cs#xml | SuffixTree.Properties | 1.0.0 | stable |
| 7 | `find_longest_repeat` | Find longest repeat in sequence | ✓ | GenomicAnalyzer.cs:L20#xml | GenomicAnalyzer.FindLongestRepeat | 1.0.0 | stable |
| 8 | `find_longest_common_region` | Find common region between sequences | ✓ | GenomicAnalyzer.cs:L178#xml | GenomicAnalyzer.FindLongestCommonRegion | 1.0.0 | stable |
| 9 | `calculate_similarity` | K-mer based similarity score | ✓ | GenomicAnalyzer.cs:L238#xml | GenomicAnalyzer.CalculateSimilarity | 1.0.0 | stable |
| 10 | `hamming_distance` | Hamming distance between sequences | ✓ | docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md | ApproximateMatcher.HammingDistance | 1.0.0 | stable |
| 11 | `edit_distance` | Levenshtein edit distance | ✓ | docs/algorithms/Pattern_Matching/Edit_Distance.md | ApproximateMatcher.EditDistance | 1.0.0 | stable |
| 12 | `count_approximate_occurrences` | Count fuzzy matches | ✓ | ApproximateMatcher.cs:L283#xml | ApproximateMatcher.CountApproximateOccurrences | 1.0.0 | stable |

---

### Server 2: Sequence (35 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `dna_validate` | Validate DNA sequence | ✓ | docs/algorithms/Sequence_Composition/Sequence_Validation.md | DnaSequence.TryCreate | 1.0.0 | stable |
| 2 | `dna_reverse_complement` | Get reverse complement | ✓ | DnaSequence.cs:L149#xml | DnaSequence.GetReverseComplementString | 1.0.0 | stable |
| 3 | `rna_validate` | Validate RNA sequence | ✓ | docs/algorithms/Sequence_Composition/Sequence_Validation.md | RnaSequence.TryCreate | 1.0.0 | stable |
| 4 | `rna_from_dna` | Transcribe DNA to RNA | ✓ | RnaSequence.cs:L147#xml | RnaSequence.FromDna | 1.0.0 | stable |
| 5 | `protein_validate` | Validate protein sequence | ✓ | ProteinSequence.cs:L357#xml | ProteinSequence.TryCreate | 1.0.0 | stable |
| 6 | `nucleotide_composition` | Calculate nucleotide frequencies | ✓ | SequenceStatistics.cs:L48#xml | SequenceStatistics.CalculateNucleotideComposition | 1.0.0 | stable |
| 7 | `amino_acid_composition` | Calculate amino acid frequencies | ✓ | SequenceStatistics.cs:L98#xml | SequenceStatistics.CalculateAminoAcidComposition | 1.0.0 | stable |
| 8 | `molecular_weight_protein` | Calculate protein MW | ✓ | SequenceStatistics.cs:L159#xml | SequenceStatistics.CalculateMolecularWeight | 1.0.0 | stable |
| 9 | `molecular_weight_nucleotide` | Calculate nucleotide MW | ✓ | SequenceStatistics.cs:L180#xml | SequenceStatistics.CalculateNucleotideMolecularWeight | 1.0.0 | stable |
| 10 | `isoelectric_point` | Calculate theoretical pI | ✓ | SequenceStatistics.cs:L228#xml | SequenceStatistics.CalculateIsoelectricPoint | 1.0.0 | stable |
| 11 | `hydrophobicity` | Calculate GRAVY score | ✓ | SequenceStatistics.cs:L306#xml | SequenceStatistics.CalculateHydrophobicity | 1.0.0 | stable |
| 12 | `thermodynamics` | Calculate thermodynamic properties | ✓ | SequenceStatistics.cs:L381#xml | SequenceStatistics.CalculateThermodynamics | 1.0.0 | stable |
| 13 | `melting_temperature` | Calculate Tm | ✓ | docs/algorithms/Molecular_Tools/Melting_Temperature.md | SequenceStatistics.CalculateMeltingTemperature | 1.0.0 | stable |
| 14 | `shannon_entropy` | Calculate Shannon entropy | ✓ | docs/algorithms/Sequence_Composition/Shannon_Entropy.md | SequenceStatistics.CalculateShannonEntropy | 1.0.0 | stable |
| 15 | `linguistic_complexity` | Calculate linguistic complexity | ✓ | docs/algorithms/Sequence_Composition/Linguistic_Complexity.md | SequenceStatistics.CalculateLinguisticComplexity | 1.0.0 | stable |
| 16 | `summarize_sequence` | Comprehensive sequence summary | ✓ | SequenceStatistics.cs:L775#xml | SequenceStatistics.SummarizeNucleotideSequence | 1.0.0 | stable |
| 17 | `gc_content` | Calculate GC percentage | ✓ | SequenceExtensions.cs:L41#xml | SequenceExtensions.CalculateGcContentFast | 1.0.0 | stable |
| 18 | `complement_base` | Get complement of single base | ✓ | SequenceExtensions.cs:L83#xml | SequenceExtensions.GetComplementBase | 1.0.0 | stable |
| 19 | `is_valid_dna` | Check if valid DNA | ✓ | SequenceExtensions.cs:L210#xml | SequenceExtensions.IsValidDna | 1.0.0 | stable |
| 20 | `is_valid_rna` | Check if valid RNA | ✓ | SequenceExtensions.cs:L225#xml | SequenceExtensions.IsValidRna | 1.0.0 | stable |
| 21 | `complexity_linguistic` | DNA linguistic complexity | ✓ | docs/algorithms/Sequence_Composition/Linguistic_Complexity.md | SequenceComplexity.CalculateLinguisticComplexity | 1.0.0 | stable |
| 22 | `complexity_shannon` | DNA Shannon entropy | ✓ | docs/algorithms/Sequence_Composition/Shannon_Entropy.md | SequenceComplexity.CalculateShannonEntropy | 1.0.0 | stable |
| 23 | `complexity_kmer_entropy` | K-mer based entropy | ✓ | SequenceComplexity.cs:L128#xml | SequenceComplexity.CalculateKmerEntropy | 1.0.0 | stable |
| 24 | `complexity_dust_score` | DUST low-complexity score | ✓ | SequenceComplexity.cs:L296#xml | SequenceComplexity.CalculateDustScore | 1.0.0 | stable |
| 25 | `complexity_mask_low` | Mask low-complexity regions | ✓ | SequenceComplexity.cs:L346#xml | SequenceComplexity.MaskLowComplexity | 1.0.0 | stable |
| 26 | `complexity_compression_ratio` | Estimate compression ratio | ✓ | SequenceComplexity.cs:L391#xml | SequenceComplexity.EstimateCompressionRatio | 1.0.0 | stable |
| 27 | `kmer_count` | Count k-mer frequencies | ✓ | docs/algorithms/K-mer_Analysis/K-mer_Counting.md | KmerAnalyzer.CountKmers | 1.0.0 | stable |
| 28 | `kmer_distance` | K-mer based distance | ✓ | KmerAnalyzer.cs:L165#xml | KmerAnalyzer.KmerDistance | 1.0.0 | stable |
| 29 | `kmer_entropy` | K-mer entropy calculation | ✓ | docs/algorithms/K-mer_Analysis/K-mer_Frequency_Analysis.md | KmerAnalyzer.CalculateKmerEntropy | 1.0.0 | stable |
| 30 | `kmer_analyze` | Comprehensive k-mer analysis | ✓ | docs/algorithms/K-mer_Analysis/K-mer_Search.md | KmerAnalyzer.AnalyzeKmers | 1.0.0 | stable |
| 31 | `iupac_code` | Get IUPAC ambiguity code | ✓ | docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md | ISequence.GetIupacCode | 1.0.0 | stable |
| 32 | `iupac_match` | Check if codes match | ✓ | docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md | ISequence.CodesMatch | 1.0.0 | stable |
| 33 | `iupac_matches` | Match nucleotide to IUPAC | ✓ | IupacHelper.cs:L14#xml | IupacHelper.MatchesIupac | 1.0.0 | stable |
| 34 | `translate_dna` | Translate DNA to protein | ✓ | Translator.cs:L20#xml | Translator.Translate(DnaSequence) | 1.0.0 | stable |
| 35 | `translate_rna` | Translate RNA to protein | ✓ | Translator.cs:L37#xml | Translator.Translate(RnaSequence) | 1.0.0 | stable |

---

### Server 3: Parsers (45 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `fasta_parse` | Parse FASTA format | ✓ | FastaParser.cs#xml | FastaParser.Parse | 1.0.0 | stable |
| 2 | `fasta_format` | Format to FASTA string | ✓ | FastaParser.cs:L78#xml | FastaParser.ToFasta | 1.0.0 | stable |
| 3 | `fasta_write` | Write FASTA file | ✓ | FastaParser.cs:L98#xml | FastaParser.WriteFile | 1.0.0 | stable |
| 4 | `fastq_parse` | Parse FASTQ format | ✓ | FastqParser.cs#xml | FastqParser.Parse | 1.0.0 | stable |
| 5 | `fastq_detect_encoding` | Detect quality encoding | ✓ | FastqParser.cs:L148#xml | FastqParser.DetectEncoding | 1.0.0 | stable |
| 6 | `fastq_encode_quality` | Encode quality scores | ✓ | FastqParser.cs:L189#xml | FastqParser.EncodeQualityScores | 1.0.0 | stable |
| 7 | `fastq_phred_to_error` | Phred to error probability | ✓ | FastqParser.cs:L206#xml | FastqParser.PhredToErrorProbability | 1.0.0 | stable |
| 8 | `fastq_error_to_phred` | Error probability to Phred | ✓ | FastqParser.cs:L214#xml | FastqParser.ErrorProbabilityToPhred | 1.0.0 | stable |
| 9 | `fastq_trim_quality` | Trim by quality score | ✓ | FastqParser.cs:L264#xml | FastqParser.TrimByQuality | 1.0.0 | stable |
| 10 | `fastq_trim_adapter` | Trim adapter sequences | ✓ | FastqParser.cs:L297#xml | FastqParser.TrimAdapter | 1.0.0 | stable |
| 11 | `fastq_statistics` | Calculate FASTQ statistics | ✓ | FastqParser.cs:L344#xml | FastqParser.CalculateStatistics | 1.0.0 | stable |
| 12 | `fastq_write` | Write FASTQ file | ✓ | FastqParser.cs:L437#xml | FastqParser.WriteToFile | 1.0.0 | stable |
| 13 | `fastq_format` | Format to FASTQ string | ✓ | FastqParser.cs:L474#xml | FastqParser.ToFastqString | 1.0.0 | stable |
| 14 | `genbank_parse` | Parse GenBank format | ✓ | GenBankParser.cs#xml | GenBankParser.Parse | 1.0.0 | stable |
| 15 | `genbank_parse_location` | Parse feature location | ✓ | GenBankParser.cs:L471#xml | GenBankParser.ParseLocation | 1.0.0 | stable |
| 16 | `genbank_extract_sequence` | Extract sequence by location | ✓ | GenBankParser.cs:L556#xml | GenBankParser.ExtractSequence | 1.0.0 | stable |
| 17 | `embl_parse` | Parse EMBL format | ✓ | EmblParser.cs#xml | EmblParser.Parse | 1.0.0 | stable |
| 18 | `embl_parse_location` | Parse EMBL location | ✓ | EmblParser.cs:L471#xml | EmblParser.ParseLocation | 1.0.0 | stable |
| 19 | `embl_extract_sequence` | Extract EMBL sequence | ✓ | EmblParser.cs:L600#xml | EmblParser.ExtractSequence | 1.0.0 | stable |
| 20 | `gff_parse` | Parse GFF/GTF format | ✓ | docs/algorithms/Annotation/GFF3_IO.md | GffParser.Parse | 1.0.0 | stable |
| 21 | `gff_statistics` | Calculate GFF statistics | ✓ | GffParser.cs:L377#xml | GffParser.CalculateStatistics | 1.0.0 | stable |
| 22 | `gff_write` | Write GFF file | ✓ | docs/algorithms/Annotation/GFF3_IO.md | GffParser.WriteToFile | 1.0.0 | stable |
| 23 | `gff_extract_sequence` | Extract sequence from GFF | ✓ | GffParser.cs:L488#xml | GffParser.ExtractSequence | 1.0.0 | stable |
| 24 | `bed_parse` | Parse BED format | ✓ | BedParser.cs#xml | BedParser.Parse | 1.0.0 | stable |
| 25 | `bed_block_length` | Calculate total block length | ✓ | BedParser.cs:L469#xml | BedParser.GetTotalBlockLength | 1.0.0 | stable |
| 26 | `bed_statistics` | Calculate BED statistics | ✓ | BedParser.cs:L510#xml | BedParser.CalculateStatistics | 1.0.0 | stable |
| 27 | `bed_write` | Write BED file | ✓ | BedParser.cs:L549#xml | BedParser.WriteToFile | 1.0.0 | stable |
| 28 | `bed_extract_sequence` | Extract sequence from BED | ✓ | BedParser.cs:L655#xml | BedParser.ExtractSequence | 1.0.0 | stable |
| 29 | `vcf_parse` | Parse VCF format | ✓ | VcfParser.cs#xml | VcfParser.Parse | 1.0.0 | stable |
| 30 | `vcf_classify` | Classify variant type | ✓ | VcfParser.cs:L376#xml | VcfParser.ClassifyVariant | 1.0.0 | stable |
| 31 | `vcf_is_snp` | Check if variant is SNP | ✓ | VcfParser.cs:L407#xml | VcfParser.IsSNP | 1.0.0 | stable |
| 32 | `vcf_is_indel` | Check if variant is indel | ✓ | VcfParser.cs:L412#xml | VcfParser.IsIndel | 1.0.0 | stable |
| 33 | `vcf_variant_length` | Get variant length | ✓ | VcfParser.cs:L421#xml | VcfParser.GetVariantLength | 1.0.0 | stable |
| 34 | `vcf_is_hom_ref` | Check homozygous reference | ✓ | VcfParser.cs:L529#xml | VcfParser.IsHomRef | 1.0.0 | stable |
| 35 | `vcf_is_hom_alt` | Check homozygous alternate | ✓ | VcfParser.cs:L539#xml | VcfParser.IsHomAlt | 1.0.0 | stable |
| 36 | `vcf_is_het` | Check heterozygous | ✓ | VcfParser.cs:L554#xml | VcfParser.IsHet | 1.0.0 | stable |
| 37 | `vcf_statistics` | Calculate VCF statistics | ✓ | VcfParser.cs:L631#xml | VcfParser.CalculateStatistics | 1.0.0 | stable |
| 38 | `vcf_write` | Write VCF file | ✓ | VcfParser.cs:L702#xml | VcfParser.WriteToFile | 1.0.0 | stable |
| 39 | `vcf_has_flag` | Check INFO flag | ✓ | VcfParser.cs:L835#xml | VcfParser.HasInfoFlag | 1.0.0 | stable |
| 40 | `to_genbank` | Convert to GenBank format | ✓ | SequenceIO.cs:L328#xml | SequenceIO.ToGenBank | 1.0.0 | stable |
| 41 | `to_bed` | Convert to BED format | ✓ | SequenceIO.cs:L545#xml | SequenceIO.ToBed | 1.0.0 | stable |
| 42 | `to_gff` | Convert to GFF format | ✓ | SequenceIO.cs:L662#xml | SequenceIO.ToGff | 1.0.0 | stable |
| 43 | `to_phylip` | Convert to PHYLIP format | ✓ | SequenceIO.cs:L913#xml | SequenceIO.ToPhylip | 1.0.0 | stable |
| 44 | `quality_char_to_phred` | Convert char to Phred | ✓ | QualityScoreAnalyzer.cs:L66#xml | QualityScoreAnalyzer.CharToPhred | 1.0.0 | stable |
| 45 | `quality_phred_to_char` | Convert Phred to char | ✓ | QualityScoreAnalyzer.cs:L75#xml | QualityScoreAnalyzer.PhredToChar | 1.0.0 | stable |

---

### Server 4: Alignment (15 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `align_global` | Needleman-Wunsch alignment | ✓ | SequenceAligner.cs:L55#xml | SequenceAligner.GlobalAlign | 1.0.0 | stable |
| 2 | `align_local` | Smith-Waterman alignment | ✓ | SequenceAligner.cs:L164#xml | SequenceAligner.LocalAlign | 1.0.0 | stable |
| 3 | `align_semi_global` | Semi-global alignment | ✓ | SequenceAligner.cs:L286#xml | SequenceAligner.SemiGlobalAlign | 1.0.0 | stable |
| 4 | `align_multiple` | Multiple sequence alignment | ✓ | SequenceAligner.cs:L495#xml | SequenceAligner.MultipleAlign | 1.0.0 | stable |
| 5 | `align_statistics` | Calculate alignment statistics | ✓ | SequenceAligner.cs:L406#xml | SequenceAligner.CalculateStatistics | 1.0.0 | stable |
| 6 | `align_format` | Format alignment for display | ✓ | SequenceAligner.cs:L446#xml | SequenceAligner.FormatAlignment | 1.0.0 | stable |
| 8 | `motif_find_exact` | Find exact motif matches | ✓ | docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md | MotifFinder.FindExactMotif | 1.0.0 | stable |
| 9 | `motif_find_degenerate` | Find IUPAC degenerate motifs | ✓ | docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md | MotifFinder.FindDegenerateMotif | 1.0.0 | stable |
| 10 | `motif_create_pwm` | Create position weight matrix | ✓ | docs/algorithms/Pattern_Matching/Position_Weight_Matrix.md | MotifFinder.CreatePwm | 1.0.0 | stable |
| 11 | `motif_scan_pwm` | Scan sequence with PWM | ✓ | docs/algorithms/Pattern_Matching/Position_Weight_Matrix.md | MotifFinder.ScanWithPwm | 1.0.0 | stable |
| 12 | `motif_consensus` | Generate consensus sequence | ✓ | MotifFinder.cs:L257#xml | MotifFinder.GenerateConsensus | 1.0.0 | stable |
| 13 | `repeat_find` | Find tandem repeats | ✓ | docs/algorithms/Repeat_Analysis/Tandem_Repeat_Detection.md | RepeatFinder.FindRepeats | 1.0.0 | stable |
| 14 | `repeat_summary` | Get repeat summary statistics | ✓ | RepeatFinder.cs:L358#xml | RepeatFinder.GetTandemRepeatSummary | 1.0.0 | stable |
| 15 | `find_with_mismatches` | Find patterns with mismatches | ✓ | docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md | ApproximateMatcher.FindWithMismatches | 1.0.0 | stable |

---

### Server 5: Variants (22 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `variant_call` | Call variants from sequences | ✓ | VariantCaller.cs#xml | VariantCaller.CallVariants | 1.0.0 | stable |
| 2 | `variant_classify_mutation` | Classify mutation type | ✓ | VariantCaller.cs:L178#xml | VariantCaller.ClassifyMutation | 1.0.0 | stable |
| 3 | `variant_titv_ratio` | Calculate Ti/Tv ratio | ✓ | VariantCaller.cs:L197#xml | VariantCaller.CalculateTiTvRatio | 1.0.0 | stable |
| 4 | `variant_statistics` | Calculate variant statistics | ✓ | VariantCaller.cs:L223#xml | VariantCaller.CalculateStatistics | 1.0.0 | stable |
| 5 | `variant_predict_effect` | Predict variant effect | ✓ | VariantCaller.cs:L263#xml | VariantCaller.PredictEffect | 1.0.0 | stable |
| 6 | `annotator_classify` | Classify variant type | ✓ | VariantAnnotator.cs:L194#xml | VariantAnnotator.ClassifyVariant | 1.0.0 | stable |
| 7 | `annotator_normalize` | Normalize variant representation | ✓ | VariantAnnotator.cs:L226#xml | VariantAnnotator.NormalizeVariant | 1.0.0 | stable |
| 8 | `annotator_impact_level` | Get impact level | ✓ | VariantAnnotator.cs:L507#xml | VariantAnnotator.GetImpactLevel | 1.0.0 | stable |
| 9 | `annotator_pathogenicity` | Predict pathogenicity | ✓ | VariantAnnotator.cs:L728#xml | VariantAnnotator.PredictPathogenicity | 1.0.0 | stable |
| 10 | `annotator_parse_vcf` | Parse VCF variant | ✓ | VariantAnnotator.cs:L1193#xml | VariantAnnotator.ParseVcfVariant | 1.0.0 | stable |
| 11 | `annotator_format_vcf` | Format as VCF INFO | ✓ | VariantAnnotator.cs:L1208#xml | VariantAnnotator.FormatAsVcfInfo | 1.0.0 | stable |
| 12 | `quality_statistics` | Quality score statistics | ✓ | QualityScoreAnalyzer.cs:L146#xml | QualityScoreAnalyzer.CalculateStatistics | 1.0.0 | stable |
| 13 | `quality_trim` | Quality-based trimming | ✓ | QualityScoreAnalyzer.cs:L254#xml | QualityScoreAnalyzer.QualityTrim | 1.0.0 | stable |
| 14 | `quality_sliding_trim` | Sliding window trimming | ✓ | QualityScoreAnalyzer.cs:L308#xml | QualityScoreAnalyzer.SlidingWindowTrim | 1.0.0 | stable |
| 15 | `quality_expected_errors` | Calculate expected errors | ✓ | QualityScoreAnalyzer.cs:L402#xml | QualityScoreAnalyzer.CalculateExpectedErrors | 1.0.0 | stable |
| 16 | `quality_mask_low` | Mask low quality bases | ✓ | QualityScoreAnalyzer.cs:L413#xml | QualityScoreAnalyzer.MaskLowQualityBases | 1.0.0 | stable |
| 17 | `quality_consensus` | Calculate consensus quality | ✓ | QualityScoreAnalyzer.cs:L569#xml | QualityScoreAnalyzer.CalculateConsensusQuality | 1.0.0 | stable |
| 18 | `quality_detect_encoding` | Detect quality encoding | ✓ | QualityScoreAnalyzer.cs:L107#xml | QualityScoreAnalyzer.DetectEncoding | 1.0.0 | stable |
| 19 | `quality_phred_to_error` | Phred to error probability | ✓ | QualityScoreAnalyzer.cs:L128#xml | QualityScoreAnalyzer.PhredToErrorProbability | 1.0.0 | stable |
| 20 | `quality_error_to_phred` | Error to Phred score | ✓ | QualityScoreAnalyzer.cs:L136#xml | QualityScoreAnalyzer.ErrorProbabilityToPhred | 1.0.0 | stable |
| 21 | `quality_phred_string` | Convert to Phred string | ✓ | QualityScoreAnalyzer.cs:L99#xml | QualityScoreAnalyzer.PhredToQualityString | 1.0.0 | stable |
| 22 | `sam_decode_flags` | Decode SAM flags | ✓ | SequenceIO.cs:L749#xml | SequenceIO.SAM_FLAGS | 1.0.0 | stable |

---

### Server 6: MolBio (38 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `primer_design` | Design PCR primers | ✓ | docs/algorithms/MolTools/Primer_Design.md | PrimerDesigner.DesignPrimers | 1.0.0 | stable |
| 2 | `primer_evaluate` | Evaluate primer quality | ✓ | PrimerDesigner.cs:L119#xml | PrimerDesigner.EvaluatePrimer | 1.0.0 | stable |
| 3 | `primer_tm` | Calculate melting temperature | ✓ | docs/algorithms/Molecular_Tools/Melting_Temperature.md | PrimerDesigner.CalculateMeltingTemperature | 1.0.0 | stable |
| 4 | `primer_tm_salt` | Tm with salt correction | ✓ | PrimerDesigner.cs:L219#xml | PrimerDesigner.CalculateMeltingTemperatureWithSalt | 1.0.0 | stable |
| 5 | `primer_gc` | Calculate GC content | ✓ | PrimerDesigner.cs:L229#xml | PrimerDesigner.CalculateGcContent | 1.0.0 | stable |
| 6 | `primer_homopolymer` | Find longest homopolymer | ✓ | docs/algorithms/Molecular_Tools/Primer_Structure_Analysis.md | PrimerDesigner.FindLongestHomopolymer | 1.0.0 | stable |
| 7 | `primer_dinucleotide_repeat` | Find dinucleotide repeats | ✓ | docs/algorithms/Molecular_Tools/Primer_Structure_Analysis.md | PrimerDesigner.FindLongestDinucleotideRepeat | 1.0.0 | stable |
| 8 | `primer_hairpin` | Check hairpin potential | ✓ | docs/algorithms/Molecular_Tools/Primer_Structure_Analysis.md | PrimerDesigner.HasHairpinPotential | 1.0.0 | stable |
| 9 | `primer_dimer` | Check primer dimer | ✓ | docs/algorithms/Molecular_Tools/Primer_Structure_Analysis.md | PrimerDesigner.HasPrimerDimer | 1.0.0 | stable |
| 10 | `primer_3prime_stability` | Calculate 3' stability | ✓ | PrimerDesigner.cs:L414#xml | PrimerDesigner.Calculate3PrimeStability | 1.0.0 | stable |
| 11 | `probe_design_tiling` | Design tiling probes | ✓ | docs/algorithms/MolTools/Hybridization_Probe_Design.md | ProbeDesigner.DesignTilingProbes | 1.0.0 | stable |
| 12 | `probe_validate` | Validate probe | ✓ | docs/algorithms/MolTools/Probe_Validation.md | ProbeDesigner.ValidateProbe | 1.0.0 | stable |
| 13 | `probe_specificity` | Check probe specificity | ✓ | ProbeDesigner.cs:L604#xml | ProbeDesigner.CheckSpecificity | 1.0.0 | stable |
| 14 | `probe_molecular_weight` | Calculate probe MW | ✓ | ProbeDesigner.cs:L647#xml | ProbeDesigner.CalculateMolecularWeight | 1.0.0 | stable |
| 15 | `probe_extinction` | Calculate extinction coefficient | ✓ | ProbeDesigner.cs:L673#xml | ProbeDesigner.CalculateExtinctionCoefficient | 1.0.0 | stable |
| 16 | `probe_concentration` | Calculate concentration | ✓ | ProbeDesigner.cs:L699#xml | ProbeDesigner.CalculateConcentration | 1.0.0 | stable |
| 17 | `restriction_digest_summary` | Get digest summary | ✓ | docs/algorithms/MolTools/Restriction_Digest_Simulation.md | RestrictionAnalyzer.GetDigestSummary | 1.0.0 | stable |
| 18 | `restriction_create_map` | Create restriction map | ✓ | docs/algorithms/MolTools/Restriction_Site_Detection.md | RestrictionAnalyzer.CreateMap | 1.0.0 | stable |
| 19 | `restriction_compatible` | Check enzyme compatibility | ✓ | RestrictionAnalyzer.cs:L375#xml | RestrictionAnalyzer.AreCompatible | 1.0.0 | stable |
| 20 | `crispr_system` | Get CRISPR system info | ✓ | CrisprDesigner.cs:L18#xml | CrisprDesigner.GetSystem | 1.0.0 | stable |
| 21 | `crispr_evaluate_guide` | Evaluate guide RNA | ✓ | docs/algorithms/MolTools/Guide_RNA_Design.md | CrisprDesigner.EvaluateGuideRna | 1.0.0 | stable |
| 22 | `crispr_specificity` | Calculate specificity score | ✓ | docs/algorithms/MolTools/Off_Target_Analysis.md | CrisprDesigner.CalculateSpecificityScore | 1.0.0 | stable |
| 23 | `codon_cai` | Calculate CAI (DnaSequence) | ✓ | CodonUsageAnalyzer.cs:L120#xml | CodonUsageAnalyzer.CalculateCai | 1.0.0 | stable |
| 24 | `codon_enc` | Calculate ENC | ✓ | CodonUsageAnalyzer.cs:L324#xml | CodonUsageAnalyzer.CalculateEnc | 1.0.0 | stable |
| 25 | `codon_statistics` | Get codon usage statistics | ✓ | CodonUsageAnalyzer.cs:L410#xml | CodonUsageAnalyzer.GetStatistics | 1.0.0 | stable |
| 26 | `codon_optimize` | Optimize codon usage | ✓ | CodonOptimizer.cs:L213#xml | CodonOptimizer.OptimizeSequence | 1.0.0 | stable |
| 27 | `codon_cai_optimizer` | Calculate CAI (optimizer) | ✓ | CodonOptimizer.cs:L388#xml | CodonOptimizer.CalculateCAI | 1.0.0 | stable |
| 28 | `codon_remove_restriction` | Remove restriction sites | ✓ | CodonOptimizer.cs:L436#xml | CodonOptimizer.RemoveRestrictionSites | 1.0.0 | stable |
| 29 | `codon_reduce_structure` | Reduce secondary structure | ✓ | CodonOptimizer.cs:L489#xml | CodonOptimizer.ReduceSecondaryStructure | 1.0.0 | stable |
| 30 | `codon_compare_usage` | Compare codon usage | ✓ | CodonOptimizer.cs:L616#xml | CodonOptimizer.CompareCodonUsage | 1.0.0 | stable |
| 31 | `codon_table_from_sequence` | Create codon table | ✓ | CodonOptimizer.cs:L673#xml | CodonOptimizer.CreateCodonTableFromSequence | 1.0.0 | stable |
| 32 | `genetic_code_get` | Get genetic code table | ✓ | GeneticCode.cs:L134#xml | GeneticCode.GetByTableNumber | 1.0.0 | stable |
| 33 | `thermo_wallace_tm` | Wallace rule Tm | ✓ | ThermoConstants.cs:L87#xml | ThermoConstants.CalculateWallaceTm | 1.0.0 | stable |
| 34 | `thermo_marmur_doty_tm` | Marmur-Doty Tm | ✓ | ThermoConstants.cs:L96#xml | ThermoConstants.CalculateMarmurDotyTm | 1.0.0 | stable |
| 35 | `thermo_salt_adjusted_tm` | Salt-adjusted Tm | ✓ | ThermoConstants.cs:L109#xml | ThermoConstants.CalculateSaltAdjustedTm | 1.0.0 | stable |
| 36 | `thermo_salt_correction` | Salt correction factor | ✓ | ThermoConstants.cs:L121#xml | ThermoConstants.CalculateSaltCorrection | 1.0.0 | stable |
| 37 | `report_create` | Create report builder | ✓ | ReportGenerator.cs:L77#xml | ReportGenerator.CreateBuilder | 1.0.0 | stable |
| 38 | `report_sequence_analysis` | Create sequence analysis report | ✓ | ReportGenerator.cs:L817#xml | ReportGenerator.CreateSequenceAnalysisReport | 1.0.0 | stable |

---

### Server 7: Assembly (13 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `assemble_olc` | OLC assembly | ✓ | SequenceAssembler.cs:L48#xml | SequenceAssembler.AssembleOLC | 1.0.0 | stable |
| 2 | `assemble_de_bruijn` | De Bruijn graph assembly | ✓ | SequenceAssembler.cs:L75#xml | SequenceAssembler.AssembleDeBruijn | 1.0.0 | stable |
| 3 | `assemble_identity` | Calculate sequence identity | ✓ | SequenceAssembler.cs:L175#xml | SequenceAssembler.CalculateIdentity | 1.0.0 | stable |
| 4 | `assemble_stats` | Calculate assembly stats | ✓ | SequenceAssembler.cs:L388#xml | SequenceAssembler.CalculateStats | 1.0.0 | stable |
| 5 | `assemble_merge` | Merge contigs | ✓ | SequenceAssembler.cs:L423#xml | SequenceAssembler.MergeContigs | 1.0.0 | stable |
| 6 | `assemble_consensus` | Compute consensus | ✓ | SequenceAssembler.cs:L534#xml | SequenceAssembler.ComputeConsensus | 1.0.0 | stable |
| 7 | `assembly_statistics` | Assembly quality statistics | ✓ | GenomeAssemblyAnalyzer.cs:L121#xml | GenomeAssemblyAnalyzer.CalculateStatistics | 1.0.0 | stable |
| 8 | `assembly_nx` | Calculate Nx statistics | ✓ | GenomeAssemblyAnalyzer.cs:L222#xml | GenomeAssemblyAnalyzer.CalculateNx | 1.0.0 | stable |
| 9 | `assembly_aun` | Calculate auN statistic | ✓ | GenomeAssemblyAnalyzer.cs:L270#xml | GenomeAssemblyAnalyzer.CalculateAuN | 1.0.0 | stable |
| 10 | `assembly_completeness` | Assess assembly completeness | ✓ | GenomeAssemblyAnalyzer.cs:L469#xml | GenomeAssemblyAnalyzer.AssessCompleteness | 1.0.0 | stable |
| 11 | `assembly_compare` | Compare assemblies | ✓ | GenomeAssemblyAnalyzer.cs:L791#xml | GenomeAssemblyAnalyzer.CompareAssemblies | 1.0.0 | stable |
| 12 | `pangenome_construct` | Construct pan-genome | ✓ | PanGenomeAnalyzer.cs:L83#xml | PanGenomeAnalyzer.ConstructPanGenome | 1.0.0 | stable |
| 13 | `pangenome_heaps_law` | Fit Heap's law | ✓ | PanGenomeAnalyzer.cs:L378#xml | PanGenomeAnalyzer.FitHeapsLaw | 1.0.0 | stable |

---

### Server 8: Phylogenetics (10 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `tree_build` | Build phylogenetic tree | ✓ | PhylogeneticAnalyzer.cs:L77#xml | PhylogeneticAnalyzer.BuildTree | 1.0.0 | stable |
| 2 | `tree_pairwise_distance` | Calculate pairwise distance | ✓ | PhylogeneticAnalyzer.cs:L134#xml | PhylogeneticAnalyzer.CalculatePairwiseDistance | 1.0.0 | stable |
| 3 | `tree_to_newick` | Export to Newick format | ✓ | PhylogeneticAnalyzer.cs:L430#xml | PhylogeneticAnalyzer.ToNewick | 1.0.0 | stable |
| 4 | `tree_parse_newick` | Parse Newick string | ✓ | PhylogeneticAnalyzer.cs:L469#xml | PhylogeneticAnalyzer.ParseNewick | 1.0.0 | stable |
| 5 | `tree_length` | Calculate tree length | ✓ | PhylogeneticAnalyzer.cs:L589#xml | PhylogeneticAnalyzer.CalculateTreeLength | 1.0.0 | stable |
| 6 | `tree_depth` | Get tree depth | ✓ | PhylogeneticAnalyzer.cs:L603#xml | PhylogeneticAnalyzer.GetTreeDepth | 1.0.0 | stable |
| 7 | `tree_robinson_foulds` | Robinson-Foulds distance | ✓ | PhylogeneticAnalyzer.cs:L616#xml | PhylogeneticAnalyzer.RobinsonFouldsDistance | 1.0.0 | stable |
| 8 | `tree_patristic_distance` | Patristic distance | ✓ | PhylogeneticAnalyzer.cs:L682#xml | PhylogeneticAnalyzer.PatristicDistance | 1.0.0 | stable |
| 9 | `comparative_genomes` | Compare genomes | ✓ | ComparativeGenomics.cs:L424#xml | ComparativeGenomics.CompareGenomes | 1.0.0 | stable |
| 10 | `comparative_ani` | Calculate ANI | ✓ | ComparativeGenomics.cs:L585#xml | ComparativeGenomics.CalculateANI | 1.0.0 | stable |

---

### Server 9: Population (15 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `pop_maf` | Calculate minor allele frequency | ✓ | PopulationGeneticsAnalyzer.cs:L127#xml | PopulationGeneticsAnalyzer.CalculateMAF | 1.0.0 | stable |
| 2 | `pop_nucleotide_diversity` | Calculate nucleotide diversity | ✓ | PopulationGeneticsAnalyzer.cs:L168#xml | PopulationGeneticsAnalyzer.CalculateNucleotideDiversity | 1.0.0 | stable |
| 3 | `pop_watterson_theta` | Calculate Watterson's theta | ✓ | PopulationGeneticsAnalyzer.cs:L202#xml | PopulationGeneticsAnalyzer.CalculateWattersonTheta | 1.0.0 | stable |
| 4 | `pop_tajimas_d` | Calculate Tajima's D | ✓ | PopulationGeneticsAnalyzer.cs:L220#xml | PopulationGeneticsAnalyzer.CalculateTajimasD | 1.0.0 | stable |
| 5 | `pop_diversity_stats` | Calculate diversity statistics | ✓ | PopulationGeneticsAnalyzer.cs:L261#xml | PopulationGeneticsAnalyzer.CalculateDiversityStatistics | 1.0.0 | stable |
| 6 | `pop_hardy_weinberg` | Test Hardy-Weinberg equilibrium | ✓ | PopulationGeneticsAnalyzer.cs:L352#xml | PopulationGeneticsAnalyzer.TestHardyWeinberg | 1.0.0 | stable |
| 7 | `pop_fst` | Calculate Fst | ✓ | PopulationGeneticsAnalyzer.cs:L453#xml | PopulationGeneticsAnalyzer.CalculateFst | 1.0.0 | stable |
| 8 | `pop_f_statistics` | Calculate F-statistics | ✓ | PopulationGeneticsAnalyzer.cs:L512#xml | PopulationGeneticsAnalyzer.CalculateFStatistics | 1.0.0 | stable |
| 9 | `pop_ld` | Calculate linkage disequilibrium | ✓ | PopulationGeneticsAnalyzer.cs:L560#xml | PopulationGeneticsAnalyzer.CalculateLD | 1.0.0 | stable |
| 10 | `pop_ihs` | Calculate iHS | ✓ | PopulationGeneticsAnalyzer.cs:L684#xml | PopulationGeneticsAnalyzer.CalculateIHS | 1.0.0 | stable |
| 11 | `pop_inbreeding_roh` | Calculate inbreeding from ROH | ✓ | PopulationGeneticsAnalyzer.cs:L870#xml | PopulationGeneticsAnalyzer.CalculateInbreedingFromROH | 1.0.0 | stable |
| 12 | `meta_taxonomic_profile` | Generate taxonomic profile | ✓ | MetagenomicsAnalyzer.cs:L229#xml | MetagenomicsAnalyzer.GenerateTaxonomicProfile | 1.0.0 | stable |
| 13 | `meta_alpha_diversity` | Calculate alpha diversity | ✓ | MetagenomicsAnalyzer.cs:L288#xml | MetagenomicsAnalyzer.CalculateAlphaDiversity | 1.0.0 | stable |
| 14 | `meta_beta_diversity` | Calculate beta diversity | ✓ | MetagenomicsAnalyzer.cs:L353#xml | MetagenomicsAnalyzer.CalculateBetaDiversity | 1.0.0 | stable |
| 15 | `comparative_reversal_distance` | Calculate reversal distance | ✓ | ComparativeGenomics.cs:L468#xml | ComparativeGenomics.CalculateReversalDistance | 1.0.0 | stable |

---

### Server 10: Epigenetics (14 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `epi_cpg_observed_expected` | CpG observed/expected ratio | ✓ | EpigeneticsAnalyzer.cs:L191#xml | EpigeneticsAnalyzer.CalculateCpGObservedExpected | 1.0.0 | stable |
| 2 | `epi_bisulfite_simulation` | Simulate bisulfite conversion | ✓ | EpigeneticsAnalyzer.cs:L283#xml | EpigeneticsAnalyzer.SimulateBisulfiteConversion | 1.0.0 | stable |
| 3 | `epi_methylation_profile` | Generate methylation profile | ✓ | EpigeneticsAnalyzer.cs:L385#xml | EpigeneticsAnalyzer.GenerateMethylationProfile | 1.0.0 | stable |
| 4 | `epi_chromatin_state` | Predict chromatin state | ✓ | EpigeneticsAnalyzer.cs:L519#xml | EpigeneticsAnalyzer.PredictChromatinState | 1.0.0 | stable |
| 5 | `epi_epigenetic_age` | Calculate epigenetic age | ✓ | EpigeneticsAnalyzer.cs:L721#xml | EpigeneticsAnalyzer.CalculateEpigeneticAge | 1.0.0 | stable |
| 6 | `mirna_seed` | Get miRNA seed sequence | ✓ | MiRnaAnalyzer.cs:L85#xml | MiRnaAnalyzer.GetSeedSequence | 1.0.0 | stable |
| 7 | `mirna_create` | Create miRNA object | ✓ | MiRnaAnalyzer.cs:L96#xml | MiRnaAnalyzer.CreateMiRna | 1.0.0 | stable |
| 8 | `mirna_reverse_complement` | Get miRNA reverse complement | ✓ | MiRnaAnalyzer.cs:L226#xml | MiRnaAnalyzer.GetReverseComplement | 1.0.0 | stable |
| 9 | `mirna_can_pair` | Check base pairing | ✓ | MiRnaAnalyzer.cs:L254#xml | MiRnaAnalyzer.CanPair | 1.0.0 | stable |
| 10 | `mirna_wobble_pair` | Check wobble pairing | ✓ | MiRnaAnalyzer.cs:L267#xml | MiRnaAnalyzer.IsWobblePair | 1.0.0 | stable |
| 11 | `mirna_align_target` | Align miRNA to target | ✓ | MiRnaAnalyzer.cs:L282#xml | MiRnaAnalyzer.AlignMiRnaToTarget | 1.0.0 | stable |
| 12 | `mirna_site_accessibility` | Calculate site accessibility | ✓ | MiRnaAnalyzer.cs:L536#xml | MiRnaAnalyzer.CalculateSiteAccessibility | 1.0.0 | stable |
| 13 | `mirna_gc_content` | Calculate miRNA GC content | ✓ | MiRnaAnalyzer.cs:L615#xml | MiRnaAnalyzer.CalculateGcContent | 1.0.0 | stable |
| 14 | `pangenome_core_alignment` | Create core genome alignment | ✓ | PanGenomeAnalyzer.cs:L487#xml | PanGenomeAnalyzer.CreateCoreGenomeAlignment | 1.0.0 | stable |

---

### Server 11: Structure (14 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `rna_can_pair` | Check RNA base pairing | ✓ | RnaSecondaryStructure.cs:L142#xml | RnaSecondaryStructure.CanPair | 1.0.0 | stable |
| 2 | `rna_complement` | Get RNA complement | ✓ | RnaSecondaryStructure.cs:L174#xml | RnaSecondaryStructure.GetComplement | 1.0.0 | stable |
| 3 | `rna_stem_energy` | Calculate stem energy | ✓ | RnaSecondaryStructure.cs:L286#xml | RnaSecondaryStructure.CalculateStemEnergy | 1.0.0 | stable |
| 4 | `rna_hairpin_energy` | Calculate hairpin loop energy | ✓ | RnaSecondaryStructure.cs:L309#xml | RnaSecondaryStructure.CalculateHairpinLoopEnergy | 1.0.0 | stable |
| 5 | `rna_mfe` | Calculate minimum free energy | ✓ | RnaSecondaryStructure.cs:L365#xml | RnaSecondaryStructure.CalculateMinimumFreeEnergy | 1.0.0 | stable |
| 6 | `rna_predict_structure` | Predict secondary structure | ✓ | RnaSecondaryStructure.cs:L428#xml | RnaSecondaryStructure.PredictStructure | 1.0.0 | stable |
| 7 | `rna_validate_dotbracket` | Validate dot-bracket notation | ✓ | RnaSecondaryStructure.cs:L620#xml | RnaSecondaryStructure.ValidateDotBracket | 1.0.0 | stable |
| 8 | `rna_structure_probability` | Calculate structure probability | ✓ | RnaSecondaryStructure.cs:L684#xml | RnaSecondaryStructure.CalculateStructureProbability | 1.0.0 | stable |
| 9 | `rna_generate_random` | Generate random RNA | ✓ | RnaSecondaryStructure.cs:L698#xml | RnaSecondaryStructure.GenerateRandomRna | 1.0.0 | stable |
| 10 | `protein_motif_find` | Find protein motifs | ✓ | ProteinMotifFinder.cs#xml | ProteinMotifFinder.FindMotifs | 1.0.0 | stable |
| 11 | `protein_prosite_to_regex` | Convert PROSITE to regex | ✓ | ProteinMotifFinder.cs:L201#xml | ProteinMotifFinder.ConvertPrositeToRegex | 1.0.0 | stable |
| 12 | `protein_disorder_predict` | Predict disordered regions | ✓ | DisorderPredictor.cs:L130#xml | DisorderPredictor.PredictDisorder | 1.0.0 | stable |
| 13 | `protein_disorder_propensity` | Get disorder propensity | ✓ | DisorderPredictor.cs:L558#xml | DisorderPredictor.GetDisorderPropensity | 1.0.0 | stable |
| 14 | `protein_disorder_promoting` | Check if disorder-promoting | ✓ | DisorderPredictor.cs:L566#xml | DisorderPredictor.IsDisorderPromoting | 1.0.0 | stable |

---

### Server 12: Annotation (8 tools)

| # | ToolName | Purpose | HasDocs | DocRef | MethodId | Version | Stability |
|---|----------|---------|---------|--------|----------|---------|-----------|
| 1 | `annotation_coding_potential` | Calculate coding potential | ✓ | GenomeAnnotator.cs:L460#xml | GenomeAnnotator.CalculateCodingPotential | 1.0.0 | stable |
| 2 | `splice_predict_structure` | Predict gene structure | ✓ | SpliceSitePredictor.cs:L478#xml | SpliceSitePredictor.PredictGeneStructure | 1.0.0 | stable |
| 3 | `splice_maxent_score` | Calculate MaxEnt score | ✓ | SpliceSitePredictor.cs:L713#xml | SpliceSitePredictor.CalculateMaxEntScore | 1.0.0 | stable |
| 4 | `splice_within_coding` | Check if within coding region | ✓ | SpliceSitePredictor.cs:L756#xml | SpliceSitePredictor.IsWithinCodingRegion | 1.0.0 | stable |
| 5 | `chromosome_karyotype` | Analyze karyotype | ✓ | ChromosomeAnalyzer.cs:L136#xml | ChromosomeAnalyzer.AnalyzeKaryotype | 1.0.0 | stable |
| 6 | `chromosome_telomeres` | Analyze telomeres | ✓ | ChromosomeAnalyzer.cs:L235#xml | ChromosomeAnalyzer.AnalyzeTelomeres | 1.0.0 | stable |
| 7 | `chromosome_centromere` | Analyze centromere | ✓ | ChromosomeAnalyzer.cs:L346#xml | ChromosomeAnalyzer.AnalyzeCentromere | 1.0.0 | stable |
| 8 | `chromosome_arm_ratio` | Calculate arm ratio | ✓ | ChromosomeAnalyzer.cs:L844#xml | ChromosomeAnalyzer.CalculateArmRatio | 1.0.0 | stable |

---

## 5. Implementation Roadmap

### Phase 0: Infrastructure & Traceability
- [ ] Create solution structure in `Seqeron.sln` (12 projects + Shared)
- [ ] Set up build pipeline
- [ ] Create traceability matrix template
- [ ] Define error code catalog
- [ ] Set up documentation templates

### Phase 1: Core Servers (Docs-first)
**Priority**: All tools have HasDocs=true

| Order | Server | Tools | Estimated Effort |
|-------|--------|-------|------------------|
| 1.1 | Core | 12 | S |
| 1.2 | Sequence | 35 | M |
| 1.3 | Parsers | 45 | L |

### Phase 2: Analysis Servers
| Order | Server | Tools | Estimated Effort |
|-------|--------|-------|------------------|
| 2.1 | Alignment | 15 | M |
| 2.2 | Variants | 22 | M |

### Phase 3: Molecular Biology
| Order | Server | Tools | Estimated Effort |
|-------|--------|-------|------------------|
| 3.1 | MolBio | 38 | L |
| 3.2 | Assembly | 13 | S |

### Phase 4: Specialized Analysis
| Order | Server | Tools | Estimated Effort |
|-------|--------|-------|------------------|
| 4.1 | Phylogenetics | 10 | S |
| 4.2 | Population | 15 | M |
| 4.3 | Epigenetics | 14 | M |
| 4.4 | Structure | 14 | M |
| 4.5 | Annotation | 8 | S |

### Phase 5: Stabilization
- [ ] Cross-server integration testing
- [ ] Performance benchmarking
- [ ] Documentation review
- [ ] Version 1.0.0 release

---

## 6. Standards

### Naming Conventions

| Element | Convention | Example |
|---------|------------|---------|
| Server name | `SuffixTree.Mcp.{Domain}` | `SuffixTree.Mcp.Core` |
| Assembly name | `suffixtree-mcp-{domain}` | `suffixtree-mcp-core` |
| Tool name | `{domain}_{action}` or `{action}` | `align_global`, `gc_content` |
| Test class | `{ToolName}Tests` | `AlignGlobalTests` |
| Doc file | `{toolName}.mcp.json`, `{toolName}.md` | `align_global.mcp.json` |

### Schema Rules

| Rule | Specification |
|------|---------------|
| JSON Schema draft | 2020-12 |
| Required fields | Explicitly marked, no implicit |
| Nullability | Use `"type": ["string", "null"]` |
| Defaults | Specified in schema, documented |
| Enums | String enums with descriptions |

### Error Mapping

| Code Range | Category |
|------------|----------|
| 1000-1999 | Input validation errors |
| 2000-2999 | Sequence format errors |
| 3000-3999 | File I/O errors |
| 4000-4999 | Algorithm errors |
| 5000-5999 | Resource limit errors |

### Versioning Policy

- **SemVer 2.0** for all tools
- Breaking changes: Major version bump
- New optional parameters: Minor version bump
- Bug fixes: Patch version bump
- Experimental tools: `0.x.y` version

---

## 7. Testing Policy

### Minimum Tests per Tool: 2

| Test Type | Purpose | Required |
|-----------|---------|----------|
| Schema Test | Validates inputSchema/outputSchema | ✓ |
| Binding Test | Validates method invocation works | ✓ |

### Additional Tests (When Justified)

| Condition | Extra Tests |
|-----------|-------------|
| Multiple overloads | +1 per overload |
| Union input types | +1 per variant |
| Complex validation | +1 for edge cases |

### Test Naming

```
{ToolName}_Schema_ValidatesCorrectly
{ToolName}_Binding_InvokesSuccessfully
{ToolName}_Overload_{Variant}_Binding_InvokesSuccessfully
```

---

## 8. MCP Documentation Policy

### Required Files per Tool

1. **`{toolName}.mcp.json`** — Machine-readable specification
2. **`{toolName}.md`** — Human-readable documentation

### `.mcp.json` Contract (JSON Schema draft 2020-12)

```json
{
  "name": "string (required)",
  "server": "string (required)",
  "description": "string (required)",
  "methodId": "string (required, or 'TBD')",
  "inputSchema": { /* JSON Schema */ },
  "outputSchema": { /* JSON Schema or null */ },
  "references": [
    { "title": "string", "ref": "string (path/URL)" }
  ],
  "examples": [
    {
      "title": "string",
      "userPrompt": "string",
      "expectedToolCall": {
        "name": "string",
        "arguments": { /* object */ }
      }
    }
  ],
  "version": "string (semver)",
  "stability": "stable | experimental"
}
```

### `.md` Structure

```markdown
# {toolName}

## Description
{description}

## Parameters
| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|

## Returns
{output description}

## Errors
| Code | Message | Cause |
|------|---------|-------|

## Examples
### {Example Title}
**User Prompt**: {prompt}
**Tool Call**:
```json
{tool call}
```
**Response**:
```json
{response}
```

## References
- [{title}]({ref})
```

---

## 9. Quality Gates

### G1: Coverage Gate
- [ ] Total tools listed = 241
- [ ] Distribution matches: Core(12), Sequence(35), Parsers(45), Alignment(15), Variants(22), MolBio(38), Assembly(13), Phylogenetics(10), Population(15), Epigenetics(14), Structure(14), Annotation(8)

### G2: Docs-first Gate
- [ ] Each server lists HasDocs=true tools first
- [ ] HasDocs=false tools have remediation plan

### G3: Traceability Gate
- [ ] Every tool has row in traceability matrix
- [ ] MethodId linked or marked TBD with task

### G4: Documentation Gate
- [ ] Every tool has `.mcp.json` file
- [ ] Every tool has `.md` file
- [ ] All references valid
- [ ] All examples complete

### G5: Tests Gate
- [ ] Every tool has 2 minimum tests
- [ ] Additional tests justified in checklist

---

## 10. Definition of Done

### Tool DoD
- [ ] Schema defined and validated
- [ ] Method binding implemented
- [ ] 2 tests passing
- [ ] `.mcp.json` created with references and examples
- [ ] `.md` created with parameters and errors
- [ ] Traceability row added

### Server DoD
- [ ] All tools meet Tool DoD
- [ ] Program.cs configured
- [ ] Server builds successfully
- [ ] Integration test passing
- [ ] README.md for server

### Delivery DoD
- [ ] All 12 servers meet Server DoD
- [ ] Quality gates G1-G5 passing
- [ ] Cross-server tests passing
- [ ] Release notes written
- [ ] Version tags applied

---

## Appendix: Traceability Matrix

| ToolName | Server | MethodId | DocRef | HasDocs | Tests | MCP Docs | Status |
|----------|--------|----------|--------|---------|-------|----------|--------|
| suffix_tree_contains | Core | SuffixTree.Contains | SuffixTree.cs#xml | ✓ | TBD | TBD | Planned |
| suffix_tree_count | Core | SuffixTree.CountOccurrences | SuffixTree.cs#xml | ✓ | TBD | TBD | Planned |
| ... | ... | ... | ... | ... | ... | ... | ... |

*(Full 241-row matrix to be maintained in separate tracking sheet)*
