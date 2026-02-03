using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Sequence.Tools;

/// <summary>
/// MCP tools for DNA/RNA sequence operations.
/// </summary>
[McpServerToolType]
public static class SequenceTools
{
    /// <summary>
    /// Validate a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "dna_validate")]
    [Description("Validate a DNA sequence. Returns whether the sequence contains only valid nucleotides (A, C, G, T).")]
    public static DnaValidateResult DnaValidate(
        [Description("The DNA sequence to validate")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var isValid = global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence, out var dna);

        if (isValid)
        {
            return new DnaValidateResult(true, sequence.Length, null);
        }
        else
        {
            // Find the invalid character for error message
            var upperSeq = sequence.ToUpperInvariant();
            for (int i = 0; i < upperSeq.Length; i++)
            {
                char c = upperSeq[i];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                {
                    return new DnaValidateResult(false, sequence.Length, $"Invalid nucleotide '{sequence[i]}' at position {i}");
                }
            }
            return new DnaValidateResult(false, sequence.Length, "Invalid sequence");
        }
    }

    /// <summary>
    /// Get the reverse complement of a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "dna_reverse_complement")]
    [Description("Get the reverse complement of a DNA sequence. A↔T, C↔G, then reversed.")]
    public static DnaReverseComplementResult DnaReverseComplement(
        [Description("The DNA sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        // Validate DNA first
        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid DNA sequence", nameof(sequence));

        var reverseComplement = global::Seqeron.Genomics.Core.DnaSequence.GetReverseComplementString(sequence);
        return new DnaReverseComplementResult(reverseComplement);
    }

    /// <summary>
    /// Validate an RNA sequence.
    /// </summary>
    [McpServerTool(Name = "rna_validate")]
    [Description("Validate an RNA sequence. Returns whether the sequence contains only valid nucleotides (A, C, G, U).")]
    public static RnaValidateResult RnaValidate(
        [Description("The RNA sequence to validate")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var isValid = global::Seqeron.Genomics.Core.RnaSequence.TryCreate(sequence, out _);

        if (isValid)
        {
            return new RnaValidateResult(true, sequence.Length, null);
        }
        else
        {
            // Find the invalid character for error message
            var upperSeq = sequence.ToUpperInvariant();
            for (int i = 0; i < upperSeq.Length; i++)
            {
                char c = upperSeq[i];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'U')
                {
                    return new RnaValidateResult(false, sequence.Length, $"Invalid nucleotide '{sequence[i]}' at position {i}");
                }
            }
            return new RnaValidateResult(false, sequence.Length, "Invalid sequence");
        }
    }

    /// <summary>
    /// Transcribe DNA to RNA (T→U).
    /// </summary>
    [McpServerTool(Name = "rna_from_dna")]
    [Description("Transcribe DNA to RNA by replacing T (thymine) with U (uracil).")]
    public static RnaFromDnaResult RnaFromDna(
        [Description("The DNA sequence to transcribe")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence, out var dna))
            throw new ArgumentException("Invalid DNA sequence", nameof(sequence));

        var rna = global::Seqeron.Genomics.Core.RnaSequence.FromDna(dna!);
        return new RnaFromDnaResult(rna.ToString());
    }

    /// <summary>
    /// Validate a protein (amino acid) sequence.
    /// </summary>
    [McpServerTool(Name = "protein_validate")]
    [Description("Validate a protein sequence. Returns whether the sequence contains only valid amino acids.")]
    public static ProteinValidateResult ProteinValidate(
        [Description("The protein sequence to validate")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var isValid = global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _);

        if (isValid)
        {
            return new ProteinValidateResult(true, sequence.Length, null);
        }
        else
        {
            // Find the invalid character for error message
            var upperSeq = sequence.ToUpperInvariant();
            var validChars = global::Seqeron.Genomics.Core.ProteinSequence.ValidCharacters;
            for (int i = 0; i < upperSeq.Length; i++)
            {
                if (!validChars.Contains(upperSeq[i]))
                {
                    return new ProteinValidateResult(false, sequence.Length, $"Invalid amino acid '{sequence[i]}' at position {i}");
                }
            }
            return new ProteinValidateResult(false, sequence.Length, "Invalid sequence");
        }
    }

    /// <summary>
    /// Calculate nucleotide composition of a DNA/RNA sequence.
    /// </summary>
    [McpServerTool(Name = "nucleotide_composition")]
    [Description("Calculate nucleotide composition (A, T, G, C, U counts) and GC content of a DNA/RNA sequence.")]
    public static NucleotideCompositionResult NucleotideComposition(
        [Description("The DNA or RNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var comp = global::Seqeron.Genomics.SequenceStatistics.CalculateNucleotideComposition(sequence);
        return new NucleotideCompositionResult(
            comp.Length,
            comp.CountA,
            comp.CountT,
            comp.CountG,
            comp.CountC,
            comp.CountU,
            comp.CountOther,
            comp.GcContent);
    }

    /// <summary>
    /// Calculate amino acid composition of a protein sequence.
    /// </summary>
    [McpServerTool(Name = "amino_acid_composition")]
    [Description("Calculate amino acid composition, molecular weight, and other properties of a protein sequence.")]
    public static AminoAcidCompositionResult AminoAcidComposition(
        [Description("The protein sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var comp = global::Seqeron.Genomics.SequenceStatistics.CalculateAminoAcidComposition(sequence);
        return new AminoAcidCompositionResult(
            comp.Length,
            comp.Counts.ToDictionary(kv => kv.Key.ToString(), kv => kv.Value),
            comp.MolecularWeight,
            comp.IsoelectricPoint,
            comp.Hydrophobicity,
            comp.ChargedResidueRatio,
            comp.AromaticResidueRatio);
    }

    /// <summary>
    /// Calculate molecular weight of a protein sequence.
    /// </summary>
    [McpServerTool(Name = "molecular_weight_protein")]
    [Description("Calculate the molecular weight of a protein sequence in Daltons (Da).")]
    public static MolecularWeightProteinResult MolecularWeightProtein(
        [Description("The protein sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var mw = global::Seqeron.Genomics.SequenceStatistics.CalculateMolecularWeight(sequence);
        return new MolecularWeightProteinResult(mw, "Da");
    }

    /// <summary>
    /// Calculate molecular weight of a DNA or RNA sequence.
    /// </summary>
    [McpServerTool(Name = "molecular_weight_nucleotide")]
    [Description("Calculate the molecular weight of a DNA or RNA sequence in Daltons (Da).")]
    public static MolecularWeightNucleotideResult MolecularWeightNucleotide(
        [Description("The DNA or RNA sequence")] string sequence,
        [Description("True for DNA, false for RNA (default: true)")] bool isDna = true)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var mw = global::Seqeron.Genomics.SequenceStatistics.CalculateNucleotideMolecularWeight(sequence, isDna);
        return new MolecularWeightNucleotideResult(mw, "Da", isDna ? "DNA" : "RNA");
    }

    /// <summary>
    /// Calculate isoelectric point (pI) of a protein sequence.
    /// </summary>
    [McpServerTool(Name = "isoelectric_point")]
    [Description("Calculate the isoelectric point (pI) of a protein sequence. pI is the pH at which the protein has no net charge.")]
    public static IsoelectricPointResult IsoelectricPoint(
        [Description("The protein sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var pI = global::Seqeron.Genomics.SequenceStatistics.CalculateIsoelectricPoint(sequence);
        return new IsoelectricPointResult(pI);
    }

    /// <summary>
    /// Calculate hydrophobicity (GRAVY index) of a protein sequence.
    /// </summary>
    [McpServerTool(Name = "hydrophobicity")]
    [Description("Calculate the grand average of hydropathy (GRAVY) index of a protein sequence. Positive values indicate hydrophobic, negative indicate hydrophilic.")]
    public static HydrophobicityResult Hydrophobicity(
        [Description("The protein sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var gravy = global::Seqeron.Genomics.SequenceStatistics.CalculateHydrophobicity(sequence);
        return new HydrophobicityResult(gravy);
    }

    /// <summary>
    /// Calculate thermodynamic properties of a DNA duplex.
    /// </summary>
    [McpServerTool(Name = "thermodynamics")]
    [Description("Calculate thermodynamic properties (ΔH, ΔS, ΔG, Tm) of a DNA duplex using the nearest-neighbor method.")]
    public static ThermodynamicsResult Thermodynamics(
        [Description("The DNA sequence")] string sequence,
        [Description("Na+ concentration in M (default: 0.05 = 50mM)")] double naConcentration = 0.05,
        [Description("Primer concentration in M (default: 0.00000025 = 250nM)")] double primerConcentration = 0.00000025)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var props = global::Seqeron.Genomics.SequenceStatistics.CalculateThermodynamics(sequence, naConcentration, primerConcentration);
        return new ThermodynamicsResult(props.DeltaH, props.DeltaS, props.DeltaG, props.MeltingTemperature);
    }

    /// <summary>
    /// Calculate simple melting temperature of a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "melting_temperature")]
    [Description("Calculate the melting temperature (Tm) of a DNA sequence using Wallace rule or GC formula.")]
    public static MeltingTemperatureResult MeltingTemperature(
        [Description("The DNA sequence")] string sequence,
        [Description("Use Wallace rule for short oligos (default: true)")] bool useWallaceRule = true)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var tm = global::Seqeron.Genomics.SequenceStatistics.CalculateMeltingTemperature(sequence, useWallaceRule);
        return new MeltingTemperatureResult(tm, "°C");
    }

    /// <summary>
    /// Calculate Shannon entropy of a sequence.
    /// </summary>
    [McpServerTool(Name = "shannon_entropy")]
    [Description("Calculate Shannon entropy of a sequence. Higher values indicate more complexity/randomness.")]
    public static ShannonEntropyResult ShannonEntropy(
        [Description("The sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var entropy = global::Seqeron.Genomics.SequenceStatistics.CalculateShannonEntropy(sequence);
        return new ShannonEntropyResult(entropy);
    }

    /// <summary>
    /// Calculate linguistic complexity of a sequence.
    /// </summary>
    [McpServerTool(Name = "linguistic_complexity")]
    [Description("Calculate linguistic complexity of a sequence based on k-mer diversity. Values range from 0 to 1.")]
    public static LinguisticComplexityResult LinguisticComplexity(
        [Description("The sequence to analyze")] string sequence,
        [Description("Maximum k-mer length to consider (default: 6)")] int maxK = 6)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var complexity = global::Seqeron.Genomics.SequenceStatistics.CalculateLinguisticComplexity(sequence, maxK);
        return new LinguisticComplexityResult(complexity);
    }

    /// <summary>
    /// Generate comprehensive summary statistics for a DNA/RNA sequence.
    /// </summary>
    [McpServerTool(Name = "summarize_sequence")]
    [Description("Generate comprehensive summary statistics for a DNA/RNA sequence including composition, GC content, entropy, complexity, and Tm.")]
    public static SummarizeSequenceResult SummarizeSequence(
        [Description("The DNA or RNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var summary = global::Seqeron.Genomics.SequenceStatistics.SummarizeNucleotideSequence(sequence);
        return new SummarizeSequenceResult(
            summary.Length,
            summary.GcContent,
            summary.Entropy,
            summary.Complexity,
            summary.MeltingTemperature,
            summary.Composition.ToDictionary(kv => kv.Key.ToString(), kv => kv.Value));
    }

    /// <summary>
    /// Calculate GC content of a DNA/RNA sequence.
    /// </summary>
    [McpServerTool(Name = "gc_content")]
    [Description("Calculate the GC content (percentage of G and C nucleotides) of a DNA/RNA sequence.")]
    public static GcContentResult GcContent(
        [Description("The DNA or RNA sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var gcContent = global::Seqeron.Genomics.Core.SequenceExtensions.CalculateGcContentFast(sequence);
        int gcCount = sequence.Count(c => c == 'G' || c == 'C' || c == 'g' || c == 'c');
        return new GcContentResult(gcContent, gcCount, sequence.Length);
    }

    /// <summary>
    /// Get the complement of a single nucleotide base.
    /// </summary>
    [McpServerTool(Name = "complement_base")]
    [Description("Get the Watson-Crick complement of a single nucleotide base (A↔T, C↔G for DNA; A↔U for RNA).")]
    public static ComplementBaseResult ComplementBase(
        [Description("The nucleotide base (A, T, G, C, or U)")] string nucleotide)
    {
        if (string.IsNullOrEmpty(nucleotide) || nucleotide.Length != 1)
            throw new ArgumentException("Must provide exactly one nucleotide character", nameof(nucleotide));

        char input = nucleotide[0];
        char complement = global::Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase(input);
        return new ComplementBaseResult(complement.ToString(), input.ToString());
    }

    /// <summary>
    /// Quick validation if a sequence contains only valid DNA characters.
    /// </summary>
    [McpServerTool(Name = "is_valid_dna")]
    [Description("Quick check if a sequence contains only valid DNA characters (A, T, G, C). Faster than dna_validate but returns less information.")]
    public static IsValidDnaResult IsValidDna(
        [Description("The sequence to validate")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        bool isValid = global::Seqeron.Genomics.Core.SequenceExtensions.IsValidDna(sequence.AsSpan());
        return new IsValidDnaResult(isValid, sequence.Length);
    }

    /// <summary>
    /// Quick validation if a sequence contains only valid RNA characters.
    /// </summary>
    [McpServerTool(Name = "is_valid_rna")]
    [Description("Quick check if a sequence contains only valid RNA characters (A, U, G, C). Faster than rna_validate but returns less information.")]
    public static IsValidRnaResult IsValidRna(
        [Description("The sequence to validate")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        bool isValid = global::Seqeron.Genomics.Core.SequenceExtensions.IsValidRna(sequence.AsSpan());
        return new IsValidRnaResult(isValid, sequence.Length);
    }

    /// <summary>
    /// Calculate k-mer entropy of a sequence.
    /// </summary>
    [McpServerTool(Name = "kmer_entropy")]
    [Description("Calculate Shannon entropy based on k-mer frequencies. Higher values indicate more complexity.")]
    public static KmerEntropyResult KmerEntropy(
        [Description("The sequence to analyze")] string sequence,
        [Description("K-mer length (default: 2)")] int k = 2)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        var entropy = global::Seqeron.Genomics.KmerAnalyzer.CalculateKmerEntropy(sequence, k);
        return new KmerEntropyResult(entropy, k);
    }

    /// <summary>
    /// Calculate linguistic complexity using SequenceComplexity class.
    /// </summary>
    [McpServerTool(Name = "complexity_linguistic")]
    [Description("Calculate DNA linguistic complexity as ratio of observed to possible subwords. LC = 1.0 for maximum complexity.")]
    public static ComplexityLinguisticResult ComplexityLinguistic(
        [Description("The DNA sequence to analyze")] string sequence,
        [Description("Maximum word length to consider (default: 10)")] int maxWordLength = 10)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (maxWordLength < 1)
            throw new ArgumentException("Max word length must be at least 1", nameof(maxWordLength));

        var complexity = global::Seqeron.Genomics.SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength);
        return new ComplexityLinguisticResult(complexity, maxWordLength);
    }

    /// <summary>
    /// Calculate Shannon entropy using SequenceComplexity class.
    /// </summary>
    [McpServerTool(Name = "complexity_shannon")]
    [Description("Calculate DNA Shannon entropy (bits per base). Maximum entropy for DNA is 2 bits.")]
    public static ComplexityShannonResult ComplexityShannon(
        [Description("The DNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var entropy = global::Seqeron.Genomics.SequenceComplexity.CalculateShannonEntropy(sequence);
        return new ComplexityShannonResult(entropy);
    }

    /// <summary>
    /// Calculate k-mer entropy using SequenceComplexity class.
    /// </summary>
    [McpServerTool(Name = "complexity_kmer_entropy")]
    [Description("Calculate k-mer based Shannon entropy for DNA complexity analysis.")]
    public static ComplexityKmerEntropyResult ComplexityKmerEntropy(
        [Description("The DNA sequence to analyze")] string sequence,
        [Description("K-mer size (default: 2 for dinucleotides)")] int k = 2)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence, out var dna))
            throw new ArgumentException("Invalid DNA sequence", nameof(sequence));

        var entropy = global::Seqeron.Genomics.SequenceComplexity.CalculateKmerEntropy(dna!, k);
        return new ComplexityKmerEntropyResult(entropy, k);
    }

    /// <summary>
    /// Calculate DUST score for low-complexity filtering.
    /// </summary>
    [McpServerTool(Name = "complexity_dust_score")]
    [Description("Calculate DUST score for low-complexity filtering (as used in BLAST). Higher scores indicate lower complexity.")]
    public static ComplexityDustScoreResult ComplexityDustScore(
        [Description("The DNA sequence to analyze")] string sequence,
        [Description("Word size for triplet counting (default: 3)")] int wordSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (wordSize < 1)
            throw new ArgumentException("Word size must be at least 1", nameof(wordSize));

        var dustScore = global::Seqeron.Genomics.SequenceComplexity.CalculateDustScore(sequence, wordSize);
        return new ComplexityDustScoreResult(dustScore, wordSize);
    }

    /// <summary>
    /// Mask low-complexity regions using DUST algorithm.
    /// </summary>
    [McpServerTool(Name = "complexity_mask_low")]
    [Description("Mask low-complexity regions in a DNA sequence using the DUST algorithm. Replaces low-complexity bases with mask character.")]
    public static ComplexityMaskLowResult ComplexityMaskLow(
        [Description("The DNA sequence to mask")] string sequence,
        [Description("Window size for analysis (default: 64)")] int windowSize = 64,
        [Description("DUST threshold above which to mask (default: 2.0)")] double threshold = 2.0,
        [Description("Character to use for masking (default: 'N')")] char maskChar = 'N')
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (windowSize < 1)
            throw new ArgumentException("Window size must be at least 1", nameof(windowSize));

        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence, out var dna))
            throw new ArgumentException("Invalid DNA sequence", nameof(sequence));

        var masked = global::Seqeron.Genomics.SequenceComplexity.MaskLowComplexity(dna!, windowSize, threshold, maskChar);
        return new ComplexityMaskLowResult(masked, sequence.Length, maskChar);
    }

    /// <summary>
    /// Estimate sequence complexity using compression ratio.
    /// </summary>
    [McpServerTool(Name = "complexity_compression_ratio")]
    [Description("Estimate sequence complexity using compression ratio. Lower ratios indicate more repetitive/less complex sequences.")]
    public static ComplexityCompressionRatioResult ComplexityCompressionRatio(
        [Description("The sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var ratio = global::Seqeron.Genomics.SequenceComplexity.EstimateCompressionRatio(sequence);
        return new ComplexityCompressionRatioResult(ratio);
    }

    /// <summary>
    /// Count k-mer frequencies in a sequence.
    /// </summary>
    [McpServerTool(Name = "kmer_count")]
    [Description("Count k-mer (substring of length k) frequencies in a sequence. Returns a dictionary of k-mers and their counts.")]
    public static KmerCountResult KmerCount(
        [Description("The sequence to analyze")] string sequence,
        [Description("K-mer length (default: 3)")] int k = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        var counts = global::Seqeron.Genomics.KmerAnalyzer.CountKmers(sequence, k);
        return new KmerCountResult(counts, k, counts.Count, counts.Values.Sum());
    }

    /// <summary>
    /// Calculate k-mer distance between two sequences.
    /// </summary>
    [McpServerTool(Name = "kmer_distance")]
    [Description("Calculate k-mer based distance between two sequences using Euclidean distance of k-mer frequencies. Lower values indicate more similar sequences.")]
    public static KmerDistanceResult KmerDistance(
        [Description("First sequence")] string sequence1,
        [Description("Second sequence")] string sequence2,
        [Description("K-mer length (default: 3)")] int k = 3)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));

        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        var distance = global::Seqeron.Genomics.KmerAnalyzer.KmerDistance(sequence1, sequence2, k);
        return new KmerDistanceResult(distance, k);
    }

    /// <summary>
    /// Analyze k-mer composition of a sequence.
    /// </summary>
    [McpServerTool(Name = "kmer_analyze")]
    [Description("Comprehensive k-mer analysis including statistics about frequency distribution, entropy, and unique k-mers.")]
    public static KmerAnalyzeResult KmerAnalyze(
        [Description("The sequence to analyze")] string sequence,
        [Description("K-mer length (default: 3)")] int k = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        var stats = global::Seqeron.Genomics.KmerAnalyzer.AnalyzeKmers(sequence, k);
        return new KmerAnalyzeResult(
            stats.TotalKmers,
            stats.UniqueKmers,
            stats.MaxCount,
            stats.MinCount,
            stats.AverageCount,
            stats.Entropy,
            k);
    }

    /// <summary>
    /// Get IUPAC ambiguity code for a set of bases.
    /// </summary>
    [McpServerTool(Name = "iupac_code")]
    [Description("Get the IUPAC ambiguity code that represents a set of nucleotide bases.")]
    public static IupacCodeResult IupacCode(
        [Description("Nucleotide bases to encode (e.g., 'AG' for purine R)")] string bases)
    {
        if (string.IsNullOrEmpty(bases))
            throw new ArgumentException("Bases cannot be null or empty", nameof(bases));

        var code = global::Seqeron.Genomics.Core.IupacDnaSequence.GetIupacCode(bases.ToUpperInvariant());
        return new IupacCodeResult(code.ToString(), bases.ToUpperInvariant());
    }

    /// <summary>
    /// Check if two IUPAC codes can match the same base.
    /// </summary>
    [McpServerTool(Name = "iupac_match")]
    [Description("Check if two IUPAC codes can represent the same nucleotide base.")]
    public static IupacMatchResult IupacMatch(
        [Description("First IUPAC code")] string code1,
        [Description("Second IUPAC code")] string code2)
    {
        if (string.IsNullOrEmpty(code1) || code1.Length != 1)
            throw new ArgumentException("Code1 must be a single IUPAC character", nameof(code1));

        if (string.IsNullOrEmpty(code2) || code2.Length != 1)
            throw new ArgumentException("Code2 must be a single IUPAC character", nameof(code2));

        var matches = global::Seqeron.Genomics.Core.IupacDnaSequence.CodesMatch(code1[0], code2[0]);
        return new IupacMatchResult(matches, code1.ToUpperInvariant(), code2.ToUpperInvariant());
    }

    /// <summary>
    /// Check if a nucleotide matches an IUPAC ambiguity code.
    /// </summary>
    [McpServerTool(Name = "iupac_matches")]
    [Description("Check if a specific nucleotide matches an IUPAC ambiguity code.")]
    public static IupacMatchesResult IupacMatches(
        [Description("The nucleotide to check (A, C, G, T)")] string nucleotide,
        [Description("The IUPAC code to match against")] string iupacCode)
    {
        if (string.IsNullOrEmpty(nucleotide) || nucleotide.Length != 1)
            throw new ArgumentException("Nucleotide must be a single character (A, C, G, T)", nameof(nucleotide));

        if (string.IsNullOrEmpty(iupacCode) || iupacCode.Length != 1)
            throw new ArgumentException("IUPAC code must be a single character", nameof(iupacCode));

        var matches = global::Seqeron.Genomics.Core.IupacHelper.MatchesIupac(
            char.ToUpperInvariant(nucleotide[0]),
            char.ToUpperInvariant(iupacCode[0]));
        return new IupacMatchesResult(matches, nucleotide.ToUpperInvariant(), iupacCode.ToUpperInvariant());
    }

    /// <summary>
    /// Translate DNA sequence to protein.
    /// </summary>
    [McpServerTool(Name = "translate_dna")]
    [Description("Translate a DNA sequence to protein using the standard genetic code.")]
    public static TranslateDnaResult TranslateDna(
        [Description("The DNA sequence to translate")] string sequence,
        [Description("Reading frame (0, 1, or 2, default: 0)")] int frame = 0,
        [Description("Stop at first stop codon (default: false)")] bool toFirstStop = false)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (frame < 0 || frame > 2)
            throw new ArgumentException("Frame must be 0, 1, or 2", nameof(frame));

        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence, out var dna))
            throw new ArgumentException("Invalid DNA sequence", nameof(sequence));

        var protein = global::Seqeron.Genomics.Core.Translator.Translate(dna!, null, frame, toFirstStop);
        return new TranslateDnaResult(protein.Sequence, frame, sequence.Length);
    }

    /// <summary>
    /// Translate RNA sequence to protein.
    /// </summary>
    [McpServerTool(Name = "translate_rna")]
    [Description("Translate an RNA sequence to protein using the standard genetic code.")]
    public static TranslateRnaResult TranslateRna(
        [Description("The RNA sequence to translate")] string sequence,
        [Description("Reading frame (0, 1, or 2, default: 0)")] int frame = 0,
        [Description("Stop at first stop codon (default: false)")] bool toFirstStop = false)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (frame < 0 || frame > 2)
            throw new ArgumentException("Frame must be 0, 1, or 2", nameof(frame));

        if (!global::Seqeron.Genomics.Core.RnaSequence.TryCreate(sequence, out var rna))
            throw new ArgumentException("Invalid RNA sequence", nameof(sequence));

        var protein = global::Seqeron.Genomics.Core.Translator.Translate(rna!, null, frame, toFirstStop);
        return new TranslateRnaResult(protein.Sequence, frame, sequence.Length);
    }
}

/// <summary>
/// Result of dna_validate operation.
/// </summary>
public record DnaValidateResult(bool Valid, int Length, string? Error);

/// <summary>
/// Result of dna_reverse_complement operation.
/// </summary>
public record DnaReverseComplementResult(string ReverseComplement);

/// <summary>
/// Result of rna_validate operation.
/// </summary>
public record RnaValidateResult(bool Valid, int Length, string? Error);

/// <summary>
/// Result of rna_from_dna operation.
/// </summary>
public record RnaFromDnaResult(string Rna);

/// <summary>
/// Result of protein_validate operation.
/// </summary>
public record ProteinValidateResult(bool Valid, int Length, string? Error);

/// <summary>
/// Result of nucleotide_composition operation.
/// </summary>
public record NucleotideCompositionResult(int Length, int A, int T, int G, int C, int U, int Other, double GcContent);

/// <summary>
/// Result of amino_acid_composition operation.
/// </summary>
public record AminoAcidCompositionResult(
    int Length,
    Dictionary<string, int> Counts,
    double MolecularWeight,
    double IsoelectricPoint,
    double Hydrophobicity,
    double ChargedResidueRatio,
    double AromaticResidueRatio);

/// <summary>
/// Result of molecular_weight_protein operation.
/// </summary>
public record MolecularWeightProteinResult(double MolecularWeight, string Unit);

/// <summary>
/// Result of molecular_weight_nucleotide operation.
/// </summary>
public record MolecularWeightNucleotideResult(double MolecularWeight, string Unit, string SequenceType);

/// <summary>
/// Result of isoelectric_point operation.
/// </summary>
public record IsoelectricPointResult(double PI);

/// <summary>
/// Result of hydrophobicity operation.
/// </summary>
public record HydrophobicityResult(double Gravy);

/// <summary>
/// Result of thermodynamics operation.
/// </summary>
public record ThermodynamicsResult(double DeltaH, double DeltaS, double DeltaG, double MeltingTemperature);

/// <summary>
/// Result of melting_temperature operation.
/// </summary>
public record MeltingTemperatureResult(double Tm, string Unit);

/// <summary>
/// Result of shannon_entropy operation.
/// </summary>
public record ShannonEntropyResult(double Entropy);

/// <summary>
/// Result of linguistic_complexity operation.
/// </summary>
public record LinguisticComplexityResult(double Complexity);

/// <summary>
/// Result of summarize_sequence operation.
/// </summary>
public record SummarizeSequenceResult(
    int Length,
    double GcContent,
    double Entropy,
    double Complexity,
    double MeltingTemperature,
    Dictionary<string, int> Composition);

/// <summary>
/// Result of gc_content operation.
/// </summary>
public record GcContentResult(double GcContent, int GcCount, int TotalCount);

/// <summary>
/// Result of complement_base operation.
/// </summary>
public record ComplementBaseResult(string Complement, string Original);

/// <summary>
/// Result of is_valid_dna operation.
/// </summary>
public record IsValidDnaResult(bool IsValid, int Length);

/// <summary>
/// Result of is_valid_rna operation.
/// </summary>
public record IsValidRnaResult(bool IsValid, int Length);

/// <summary>
/// Result of kmer_entropy operation.
/// </summary>
public record KmerEntropyResult(double Entropy, int K);

/// <summary>
/// Result of complexity_linguistic operation.
/// </summary>
public record ComplexityLinguisticResult(double Complexity, int MaxWordLength);

/// <summary>
/// Result of complexity_shannon operation.
/// </summary>
public record ComplexityShannonResult(double Entropy);

/// <summary>
/// Result of complexity_kmer_entropy operation.
/// </summary>
public record ComplexityKmerEntropyResult(double Entropy, int K);

/// <summary>
/// Result of complexity_dust_score operation.
/// </summary>
public record ComplexityDustScoreResult(double DustScore, int WordSize);

/// <summary>
/// Result of complexity_mask_low operation.
/// </summary>
public record ComplexityMaskLowResult(string MaskedSequence, int OriginalLength, char MaskChar);

/// <summary>
/// Result of complexity_compression_ratio operation.
/// </summary>
public record ComplexityCompressionRatioResult(double CompressionRatio);

/// <summary>
/// Result of kmer_count operation.
/// </summary>
public record KmerCountResult(Dictionary<string, int> Counts, int K, int UniqueKmers, int TotalKmers);

/// <summary>
/// Result of kmer_distance operation.
/// </summary>
public record KmerDistanceResult(double Distance, int K);

/// <summary>
/// Result of kmer_analyze operation.
/// </summary>
public record KmerAnalyzeResult(
    int TotalKmers,
    int UniqueKmers,
    int MaxCount,
    int MinCount,
    double AverageCount,
    double Entropy,
    int K);

/// <summary>
/// Result of iupac_code operation.
/// </summary>
public record IupacCodeResult(string Code, string InputBases);

/// <summary>
/// Result of iupac_match operation.
/// </summary>
public record IupacMatchResult(bool Matches, string Code1, string Code2);

/// <summary>
/// Result of iupac_matches operation.
/// </summary>
public record IupacMatchesResult(bool Matches, string Nucleotide, string IupacCode);

/// <summary>
/// Result of translate_dna operation.
/// </summary>
public record TranslateDnaResult(string Protein, int Frame, int DnaLength);

/// <summary>
/// Result of translate_rna operation.
/// </summary>
public record TranslateRnaResult(string Protein, int Frame, int RnaLength);
