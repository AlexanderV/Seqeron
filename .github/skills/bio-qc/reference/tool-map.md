# bio-qc tool map (76 tools)

Grouped by sub-task. Each row: **[MCP] tool** · `Method ID` · one-line purpose. Open the linked
per-tool doc for the full I/O schema — do not guess parameters. Servers: **Sequence** (35), **Parsers** (41).

## Parsing by format (Parsers server)

### FASTA
| Tool | Method ID | Purpose |
|---|---|---|
| [`fasta_parse`](../../../../docs/mcp/tools/parsers/fasta_parse.md) | `FastaParser.Parse` | Parse FASTA text into records (id, description, sequence). |
| [`fasta_format`](../../../../docs/mcp/tools/parsers/fasta_format.md) | `FastaParser.ToFasta` | Serialize records back to FASTA text (line-wrapped). |
| [`fasta_write`](../../../../docs/mcp/tools/parsers/fasta_write.md) | `FastaParser.WriteFile` | Write records to a FASTA file. |

### FASTQ
| Tool | Method ID | Purpose |
|---|---|---|
| [`fastq_parse`](../../../../docs/mcp/tools/parsers/fastq_parse.md) | `FastqParser.Parse` | Parse FASTQ into records (id, sequence, quality). |
| [`fastq_statistics`](../../../../docs/mcp/tools/parsers/fastq_statistics.md) | `FastqParser.CalculateStatistics` | Read/base counts, mean Q, Q20/Q30 %, GC content. |
| [`fastq_detect_encoding`](../../../../docs/mcp/tools/parsers/fastq_detect_encoding.md) | `FastqParser.DetectEncoding` | Detect Phred+33 vs +64. **Guarded (PARSE-FASTQ-001)** — Ambiguous → STOP. |
| [`fastq_encode_quality`](../../../../docs/mcp/tools/parsers/fastq_encode_quality.md) | `FastqParser.EncodeQualityScores` | Encode numeric Phred scores to an ASCII quality string. |
| [`fastq_error_to_phred`](../../../../docs/mcp/tools/parsers/fastq_error_to_phred.md) | `FastqParser.ErrorProbabilityToPhred` | Convert error probability → Phred score. |
| [`fastq_phred_to_error`](../../../../docs/mcp/tools/parsers/fastq_phred_to_error.md) | `FastqParser.PhredToErrorProbability` | Convert Phred score → error probability. |
| [`fastq_filter`](../../../../docs/mcp/tools/parsers/fastq_filter.md) | `FastqParser.FilterByQuality` | Keep reads meeting a mean-quality threshold. |
| [`fastq_trim_quality`](../../../../docs/mcp/tools/parsers/fastq_trim_quality.md) | `FastqParser.TrimByQuality` | Trim low-quality ends of reads. |
| [`fastq_trim_adapter`](../../../../docs/mcp/tools/parsers/fastq_trim_adapter.md) | `FastqParser.TrimAdapter` | Remove adapter sequence from reads. |
| [`fastq_format`](../../../../docs/mcp/tools/parsers/fastq_format.md) | `FastqParser.ToFastqString` | Serialize records back to FASTQ text. |
| [`fastq_write`](../../../../docs/mcp/tools/parsers/fastq_write.md) | `FastqParser.WriteToFile` | Write records to a FASTQ file. |

### GenBank
| Tool | Method ID | Purpose |
|---|---|---|
| [`genbank_parse`](../../../../docs/mcp/tools/parsers/genbank_parse.md) | `GenBankParser.Parse` | Parse a GenBank record (metadata + features + sequence). |
| [`genbank_features`](../../../../docs/mcp/tools/parsers/genbank_features.md) | `GenBankParser.GetFeatures` | Extract the feature table (type, location, qualifiers). |
| [`genbank_extract_sequence`](../../../../docs/mcp/tools/parsers/genbank_extract_sequence.md) | `GenBankParser.ExtractSequence` | Extract the sub-sequence for a feature location. |
| [`genbank_parse_location`](../../../../docs/mcp/tools/parsers/genbank_parse_location.md) | `GenBankParser.ParseLocation` | Parse a GenBank location string (join/complement/ranges). |
| [`genbank_statistics`](../../../../docs/mcp/tools/parsers/genbank_statistics.md) | `GenBankParser.Statistics` | Summary counts for a GenBank record. |

### GFF
| Tool | Method ID | Purpose |
|---|---|---|
| [`gff_parse`](../../../../docs/mcp/tools/parsers/gff_parse.md) | `GffParser.Parse` | Parse GFF3 into feature records. |
| [`gff_filter`](../../../../docs/mcp/tools/parsers/gff_filter.md) | `GffParser.Filter*` | Filter features by type / region / attribute. |
| [`gff_statistics`](../../../../docs/mcp/tools/parsers/gff_statistics.md) | `GffParser.CalculateStatistics` | Feature-type counts and length stats. |

### VCF
| Tool | Method ID | Purpose |
|---|---|---|
| [`vcf_parse`](../../../../docs/mcp/tools/parsers/vcf_parse.md) | `VcfParser.Parse` | Parse VCF into variant records. |
| [`vcf_statistics`](../../../../docs/mcp/tools/parsers/vcf_statistics.md) | `VcfParser.CalculateStatistics` | Ts/Tv, SNP/indel counts and other VCF stats. |
| [`vcf_filter`](../../../../docs/mcp/tools/parsers/vcf_filter.md) | `VcfParser.Filter*` | Filter variants (quality, region, type…). |
| [`vcf_classify`](../../../../docs/mcp/tools/parsers/vcf_classify.md) | `VcfParser.ClassifyVariant` | Classify a variant (SNP / insertion / deletion / MNP). |
| [`vcf_variant_length`](../../../../docs/mcp/tools/parsers/vcf_variant_length.md) | `VcfParser.GetVariantLength` | Length of a variant (indel size). |
| [`vcf_is_snp`](../../../../docs/mcp/tools/parsers/vcf_is_snp.md) | `VcfParser.IsSNP` | Is the variant a SNP? |
| [`vcf_is_indel`](../../../../docs/mcp/tools/parsers/vcf_is_indel.md) | `VcfParser.IsIndel` | Is the variant an indel? |
| [`vcf_is_het`](../../../../docs/mcp/tools/parsers/vcf_is_het.md) | `VcfParser.IsHet` | Is the genotype heterozygous? |
| [`vcf_is_hom_ref`](../../../../docs/mcp/tools/parsers/vcf_is_hom_ref.md) | `VcfParser.IsHomRef` | Is the genotype homozygous reference? |
| [`vcf_is_hom_alt`](../../../../docs/mcp/tools/parsers/vcf_is_hom_alt.md) | `VcfParser.IsHomAlt` | Is the genotype homozygous alternate? |
| [`vcf_has_flag`](../../../../docs/mcp/tools/parsers/vcf_has_flag.md) | `VcfParser.HasInfoFlag` | Does the INFO column carry a given flag? |
| [`vcf_write`](../../../../docs/mcp/tools/parsers/vcf_write.md) | `VcfParser.WriteToFile` | Write variant records to a VCF file. |

### BED
| Tool | Method ID | Purpose |
|---|---|---|
| [`bed_parse`](../../../../docs/mcp/tools/parsers/bed_parse.md) | `BedParser.Parse` | Parse BED intervals (0-based, half-open). |
| [`bed_filter`](../../../../docs/mcp/tools/parsers/bed_filter.md) | `BedParser.FilterBy*` | Filter intervals by chrom / score / name. |
| [`bed_intersect`](../../../../docs/mcp/tools/parsers/bed_intersect.md) | `BedParser.Intersect` | Intersect two interval sets. |
| [`bed_merge`](../../../../docs/mcp/tools/parsers/bed_merge.md) | `BedParser.MergeOverlapping` | Merge overlapping/adjacent intervals. |

### EMBL
| Tool | Method ID | Purpose |
|---|---|---|
| [`embl_parse`](../../../../docs/mcp/tools/parsers/embl_parse.md) | `EmblParser.Parse` | Parse an EMBL record. |
| [`embl_features`](../../../../docs/mcp/tools/parsers/embl_features.md) | `EmblParser.GetFeatures` | Extract the EMBL feature table. |
| [`embl_statistics`](../../../../docs/mcp/tools/parsers/embl_statistics.md) | `EmblParser.Statistics` | Summary counts for an EMBL record. |

## Validation (Sequence server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`dna_validate`](../../../../docs/mcp/tools/sequence/dna_validate.md) | `DnaSequence.TryCreate` | Validate DNA; returns `valid`, `length`, `error`. |
| [`rna_validate`](../../../../docs/mcp/tools/sequence/rna_validate.md) | `RnaSequence.TryCreate` | Validate RNA sequence. |
| [`protein_validate`](../../../../docs/mcp/tools/sequence/protein_validate.md) | `ProteinSequence.TryCreate` | Validate protein (amino-acid alphabet). |
| [`is_valid_dna`](../../../../docs/mcp/tools/sequence/is_valid_dna.md) | `SequenceExtensions.IsValidDna` | Quick boolean DNA-alphabet check. |
| [`is_valid_rna`](../../../../docs/mcp/tools/sequence/is_valid_rna.md) | `SequenceExtensions.IsValidRna` | Quick boolean RNA-alphabet check. |

## Composition & statistics (Sequence server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`gc_content`](../../../../docs/mcp/tools/sequence/gc_content.md) | `SequenceExtensions.CalculateGcContentFast` | GC% + gc/total counts. |
| [`nucleotide_composition`](../../../../docs/mcp/tools/sequence/nucleotide_composition.md) | `SequenceStatistics.CalculateNucleotideComposition` | Per-base counts/frequencies. |
| [`amino_acid_composition`](../../../../docs/mcp/tools/sequence/amino_acid_composition.md) | `SequenceStatistics.CalculateAminoAcidComposition` | Per-residue counts/frequencies. |
| [`summarize_sequence`](../../../../docs/mcp/tools/sequence/summarize_sequence.md) | `SequenceStatistics.SummarizeNucleotideSequence` | One-shot summary (length, GC, composition, Tm…). |
| [`melting_temperature`](../../../../docs/mcp/tools/sequence/melting_temperature.md) | `SequenceStatistics.CalculateMeltingTemperature` | Tm (see doc for method/salt/conc params). |
| [`thermodynamics`](../../../../docs/mcp/tools/sequence/thermodynamics.md) | `SequenceStatistics.CalculateThermodynamics` | ΔG/ΔH/ΔS nearest-neighbor thermodynamics. |
| [`molecular_weight_nucleotide`](../../../../docs/mcp/tools/sequence/molecular_weight_nucleotide.md) | `SequenceStatistics.CalculateNucleotideMolecularWeight` | MW of a nucleotide sequence. |
| [`molecular_weight_protein`](../../../../docs/mcp/tools/sequence/molecular_weight_protein.md) | `SequenceStatistics.CalculateMolecularWeight` | MW of a protein sequence. |
| [`hydrophobicity`](../../../../docs/mcp/tools/sequence/hydrophobicity.md) | `SequenceStatistics.CalculateHydrophobicity` | Protein hydrophobicity (GRAVY-style). |
| [`isoelectric_point`](../../../../docs/mcp/tools/sequence/isoelectric_point.md) | `SequenceStatistics.CalculateIsoelectricPoint` | Protein isoelectric point (pI). |
| [`shannon_entropy`](../../../../docs/mcp/tools/sequence/shannon_entropy.md) | `SequenceStatistics.CalculateShannonEntropy` | Shannon entropy of a sequence. |
| [`linguistic_complexity`](../../../../docs/mcp/tools/sequence/linguistic_complexity.md) | `SequenceStatistics.CalculateLinguisticComplexity` | Linguistic complexity (stats variant). |

## Complexity / low-complexity (Sequence server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`complexity_shannon`](../../../../docs/mcp/tools/sequence/complexity_shannon.md) | `SequenceComplexity.CalculateShannonEntropy` | Shannon entropy (complexity module). |
| [`complexity_dust_score`](../../../../docs/mcp/tools/sequence/complexity_dust_score.md) | `SequenceComplexity.CalculateDustScore` | DUST low-complexity score. |
| [`complexity_linguistic`](../../../../docs/mcp/tools/sequence/complexity_linguistic.md) | `SequenceComplexity.CalculateLinguisticComplexity` | Linguistic complexity (complexity module). |
| [`complexity_kmer_entropy`](../../../../docs/mcp/tools/sequence/complexity_kmer_entropy.md) | `SequenceComplexity.CalculateKmerEntropy` | k-mer entropy (complexity module). |
| [`complexity_compression_ratio`](../../../../docs/mcp/tools/sequence/complexity_compression_ratio.md) | `SequenceComplexity.EstimateCompressionRatio` | Lempel-Ziv-style compression ratio. |
| [`complexity_mask_low`](../../../../docs/mcp/tools/sequence/complexity_mask_low.md) | `SequenceComplexity.MaskLowComplexity` | Mask (soft/hard) low-complexity regions. |

## k-mer profiling (Sequence server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`kmer_count`](../../../../docs/mcp/tools/sequence/kmer_count.md) | `KmerAnalyzer.CountKmers` | Count k-mers of a sequence. |
| [`kmer_analyze`](../../../../docs/mcp/tools/sequence/kmer_analyze.md) | `KmerAnalyzer.AnalyzeKmers` | k-mer profile / frequency analysis. |
| [`kmer_entropy`](../../../../docs/mcp/tools/sequence/kmer_entropy.md) | `KmerAnalyzer.CalculateKmerEntropy` | k-mer entropy of a sequence. |
| [`kmer_distance`](../../../../docs/mcp/tools/sequence/kmer_distance.md) | `KmerAnalyzer.KmerDistance` | k-mer-based distance between two sequences. |

> For assembly/de-Bruijn k-mer workflows use **`bio-assembly`**; this skill owns the composition-screen subset.

## Translation & strand (Sequence server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`rna_from_dna`](../../../../docs/mcp/tools/sequence/rna_from_dna.md) | `RnaSequence.FromDna` | Transcribe DNA → RNA (T→U). |
| [`translate_dna`](../../../../docs/mcp/tools/sequence/translate_dna.md) | `Translator.Translate(DnaSequence)` | Translate DNA → protein. |
| [`translate_rna`](../../../../docs/mcp/tools/sequence/translate_rna.md) | `Translator.Translate(RnaSequence)` | Translate RNA → protein. |
| [`dna_reverse_complement`](../../../../docs/mcp/tools/sequence/dna_reverse_complement.md) | `DnaSequence.GetReverseComplementString` | Reverse-complement a DNA sequence. |
| [`complement_base`](../../../../docs/mcp/tools/sequence/complement_base.md) | `SequenceExtensions.GetComplementBase` | Complement of a single base. |

## IUPAC ambiguity codes (Sequence server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`iupac_code`](../../../../docs/mcp/tools/sequence/iupac_code.md) | `IupacDnaSequence.GetIupacCode` | IUPAC code for a set of bases. |
| [`iupac_match`](../../../../docs/mcp/tools/sequence/iupac_match.md) | `IupacDnaSequence.CodesMatch` | Do two IUPAC codes match? |
| [`iupac_matches`](../../../../docs/mcp/tools/sequence/iupac_matches.md) | `IupacHelper.MatchesIupac` | Does a base match an IUPAC code / pattern? |
