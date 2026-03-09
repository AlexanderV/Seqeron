namespace Seqeron.Mcp.Sequence.Tools;

// ================================
// DNA Results
// ================================

/// <summary>Result of dna_validate operation.</summary>
public record DnaValidateResult(bool Valid, int Length, string? Error);

/// <summary>Result of dna_reverse_complement operation.</summary>
public record DnaReverseComplementResult(string ReverseComplement);

// ================================
// RNA Results
// ================================

/// <summary>Result of rna_validate operation.</summary>
public record RnaValidateResult(bool Valid, int Length, string? Error);

/// <summary>Result of rna_from_dna operation.</summary>
public record RnaFromDnaResult(string Rna);

// ================================
// Protein Results
// ================================

/// <summary>Result of protein_validate operation.</summary>
public record ProteinValidateResult(bool Valid, int Length, string? Error);

/// <summary>Result of amino_acid_composition operation.</summary>
public record AminoAcidCompositionResult(
    int Length,
    Dictionary<string, int> Counts,
    double MolecularWeight,
    double IsoelectricPoint,
    double Hydrophobicity,
    double ChargedResidueRatio,
    double AromaticResidueRatio);

/// <summary>Result of molecular_weight_protein operation.</summary>
public record MolecularWeightProteinResult(double MolecularWeight, string Unit);

/// <summary>Result of isoelectric_point operation.</summary>
public record IsoelectricPointResult(double PI);

/// <summary>Result of hydrophobicity operation.</summary>
public record HydrophobicityResult(double Gravy);

// ================================
// Nucleotide / Sequence Results
// ================================

/// <summary>Result of nucleotide_composition operation.</summary>
public record NucleotideCompositionResult(int Length, int A, int T, int G, int C, int U, int Other, double GcContent);

/// <summary>Result of molecular_weight_nucleotide operation.</summary>
public record MolecularWeightNucleotideResult(double MolecularWeight, string Unit, string SequenceType);

/// <summary>Result of thermodynamics operation.</summary>
public record ThermodynamicsResult(double DeltaH, double DeltaS, double DeltaG, double MeltingTemperature);

/// <summary>Result of melting_temperature operation.</summary>
public record MeltingTemperatureResult(double Tm, string Unit);

/// <summary>Result of shannon_entropy operation.</summary>
public record ShannonEntropyResult(double Entropy);

/// <summary>Result of linguistic_complexity operation.</summary>
public record LinguisticComplexityResult(double Complexity);

/// <summary>Result of summarize_sequence operation.</summary>
public record SummarizeSequenceResult(
    int Length,
    double GcContent,
    double Entropy,
    double Complexity,
    double MeltingTemperature,
    Dictionary<string, int> Composition);

/// <summary>Result of gc_content operation.</summary>
public record GcContentResult(double GcContent, int GcCount, int TotalCount);

/// <summary>Result of complement_base operation.</summary>
public record ComplementBaseResult(string Complement, string Original);

/// <summary>Result of is_valid_dna operation.</summary>
public record IsValidDnaResult(bool IsValid, int Length);

/// <summary>Result of is_valid_rna operation.</summary>
public record IsValidRnaResult(bool IsValid, int Length);

// ================================
// K-mer Results
// ================================

/// <summary>Result of kmer_entropy operation.</summary>
public record KmerEntropyResult(double Entropy, int K);

/// <summary>Result of kmer_count operation.</summary>
public record KmerCountResult(Dictionary<string, int> Counts, int K, int UniqueKmers, int TotalKmers);

/// <summary>Result of kmer_distance operation.</summary>
public record KmerDistanceResult(double Distance, int K);

/// <summary>Result of kmer_analyze operation.</summary>
public record KmerAnalyzeResult(
    int TotalKmers,
    int UniqueKmers,
    int MaxCount,
    int MinCount,
    double AverageCount,
    double Entropy,
    int K);

// ================================
// Complexity Results
// ================================

/// <summary>Result of complexity_linguistic operation.</summary>
public record ComplexityLinguisticResult(double Complexity, int MaxWordLength);

/// <summary>Result of complexity_shannon operation.</summary>
public record ComplexityShannonResult(double Entropy);

/// <summary>Result of complexity_kmer_entropy operation.</summary>
public record ComplexityKmerEntropyResult(double Entropy, int K);

/// <summary>Result of complexity_dust_score operation.</summary>
public record ComplexityDustScoreResult(double DustScore, int WordSize);

/// <summary>Result of complexity_mask_low operation.</summary>
public record ComplexityMaskLowResult(string MaskedSequence, int OriginalLength, char MaskChar);

/// <summary>Result of complexity_compression_ratio operation.</summary>
public record ComplexityCompressionRatioResult(double CompressionRatio);

// ================================
// IUPAC Results
// ================================

/// <summary>Result of iupac_code operation.</summary>
public record IupacCodeResult(string Code, string InputBases);

/// <summary>Result of iupac_match operation.</summary>
public record IupacMatchResult(bool Matches, string Code1, string Code2);

/// <summary>Result of iupac_matches operation.</summary>
public record IupacMatchesResult(bool Matches, string Nucleotide, string IupacCode);

// ================================
// Translation Results
// ================================

/// <summary>Result of translate_dna operation.</summary>
public record TranslateDnaResult(string Protein, int Frame, int DnaLength);

/// <summary>Result of translate_rna operation.</summary>
public record TranslateRnaResult(string Protein, int Frame, int RnaLength);
