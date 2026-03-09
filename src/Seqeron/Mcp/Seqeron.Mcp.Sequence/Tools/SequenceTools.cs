using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Sequence.Tools;

/// <summary>
/// MCP tools for DNA/RNA sequence operations.
/// </summary>
[McpServerToolType]
public class SequenceTools
{
    /// <summary>
    /// Validate a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "dna_validate", Title = "DNA — Validate Sequence", ReadOnly = true)]
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
    [McpServerTool(Name = "dna_reverse_complement", Title = "DNA — Reverse Complement", ReadOnly = true)]
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
    [McpServerTool(Name = "rna_validate", Title = "RNA — Validate Sequence", ReadOnly = true)]
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
    [McpServerTool(Name = "rna_from_dna", Title = "RNA — Transcribe from DNA", ReadOnly = true)]
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
    [McpServerTool(Name = "protein_validate", Title = "Protein — Validate Sequence", ReadOnly = true)]
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
    [McpServerTool(Name = "nucleotide_composition", Title = "Sequence — Nucleotide Composition", ReadOnly = true)]
    [Description("Calculate nucleotide composition (A, T, G, C, U counts) and GC content of a DNA/RNA sequence.")]
    public static NucleotideCompositionResult NucleotideComposition(
        [Description("The DNA or RNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var comp = SequenceStatistics.CalculateNucleotideComposition(sequence);
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
    [McpServerTool(Name = "amino_acid_composition", Title = "Protein — Amino Acid Composition", ReadOnly = true)]
    [Description("Calculate amino acid composition, molecular weight, and other properties of a protein sequence.")]
    public static AminoAcidCompositionResult AminoAcidComposition(
        [Description("The protein sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var comp = SequenceStatistics.CalculateAminoAcidComposition(sequence);
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
    [McpServerTool(Name = "molecular_weight_protein", Title = "Protein — Molecular Weight", ReadOnly = true)]
    [Description("Calculate the molecular weight of a protein sequence in Daltons (Da).")]
    public static MolecularWeightProteinResult MolecularWeightProtein(
        [Description("The protein sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var mw = SequenceStatistics.CalculateMolecularWeight(sequence);
        return new MolecularWeightProteinResult(mw, "Da");
    }

    /// <summary>
    /// Calculate molecular weight of a DNA or RNA sequence.
    /// </summary>
    [McpServerTool(Name = "molecular_weight_nucleotide", Title = "Sequence — Nucleotide Molecular Weight", ReadOnly = true)]
    [Description("Calculate the molecular weight of a DNA or RNA sequence in Daltons (Da).")]
    public static MolecularWeightNucleotideResult MolecularWeightNucleotide(
        [Description("The DNA or RNA sequence")] string sequence,
        [Description("True for DNA, false for RNA (default: true)")] bool isDna = true)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var mw = SequenceStatistics.CalculateNucleotideMolecularWeight(sequence, isDna);
        return new MolecularWeightNucleotideResult(mw, "Da", isDna ? "DNA" : "RNA");
    }

    /// <summary>
    /// Calculate isoelectric point (pI) of a protein sequence.
    /// </summary>
    [McpServerTool(Name = "isoelectric_point", Title = "Protein — Isoelectric Point", ReadOnly = true)]
    [Description("Calculate the isoelectric point (pI) of a protein sequence. pI is the pH at which the protein has no net charge.")]
    public static IsoelectricPointResult IsoelectricPoint(
        [Description("The protein sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var pI = SequenceStatistics.CalculateIsoelectricPoint(sequence);
        return new IsoelectricPointResult(pI);
    }

    /// <summary>
    /// Calculate hydrophobicity (GRAVY index) of a protein sequence.
    /// </summary>
    [McpServerTool(Name = "hydrophobicity", Title = "Protein — Hydrophobicity (GRAVY)", ReadOnly = true)]
    [Description("Calculate the grand average of hydropathy (GRAVY) index of a protein sequence. Positive values indicate hydrophobic, negative indicate hydrophilic.")]
    public static HydrophobicityResult Hydrophobicity(
        [Description("The protein sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (!global::Seqeron.Genomics.Core.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var gravy = SequenceStatistics.CalculateHydrophobicity(sequence);
        return new HydrophobicityResult(gravy);
    }

    /// <summary>
    /// Calculate thermodynamic properties of a DNA duplex.
    /// </summary>
    [McpServerTool(Name = "thermodynamics", Title = "DNA — Thermodynamic Properties", ReadOnly = true)]
    [Description("Calculate thermodynamic properties (ΔH, ΔS, ΔG, Tm) of a DNA duplex using the nearest-neighbor method.")]
    public static ThermodynamicsResult Thermodynamics(
        [Description("The DNA sequence")] string sequence,
        [Description("Na+ concentration in M (default: 0.05 = 50mM)")] double naConcentration = 0.05,
        [Description("Primer concentration in M (default: 0.00000025 = 250nM)")] double primerConcentration = 0.00000025)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var props = SequenceStatistics.CalculateThermodynamics(sequence, naConcentration, primerConcentration);
        return new ThermodynamicsResult(props.DeltaH, props.DeltaS, props.DeltaG, props.MeltingTemperature);
    }

    /// <summary>
    /// Calculate simple melting temperature of a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "melting_temperature", Title = "DNA — Melting Temperature", ReadOnly = true)]
    [Description("Calculate the melting temperature (Tm) of a DNA sequence using Wallace rule or GC formula.")]
    public static MeltingTemperatureResult MeltingTemperature(
        [Description("The DNA sequence")] string sequence,
        [Description("Use Wallace rule for short oligos (default: true)")] bool useWallaceRule = true)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var tm = SequenceStatistics.CalculateMeltingTemperature(sequence, useWallaceRule);
        return new MeltingTemperatureResult(tm, "°C");
    }

    /// <summary>
    /// Calculate Shannon entropy of a sequence.
    /// </summary>
    [McpServerTool(Name = "shannon_entropy", Title = "Sequence — Shannon Entropy", ReadOnly = true)]
    [Description("Calculate Shannon entropy of a sequence. Higher values indicate more complexity/randomness.")]
    public static ShannonEntropyResult ShannonEntropy(
        [Description("The sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var entropy = SequenceStatistics.CalculateShannonEntropy(sequence);
        return new ShannonEntropyResult(entropy);
    }

    /// <summary>
    /// Calculate linguistic complexity of a sequence.
    /// </summary>
    [McpServerTool(Name = "linguistic_complexity", Title = "Sequence — Linguistic Complexity", ReadOnly = true)]
    [Description("Calculate linguistic complexity of a sequence based on k-mer diversity. Values range from 0 to 1.")]
    public static LinguisticComplexityResult LinguisticComplexity(
        [Description("The sequence to analyze")] string sequence,
        [Description("Maximum k-mer length to consider (default: 6)")] int maxK = 6)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var complexity = SequenceStatistics.CalculateLinguisticComplexity(sequence, maxK);
        return new LinguisticComplexityResult(complexity);
    }

    /// <summary>
    /// Generate comprehensive summary statistics for a DNA/RNA sequence.
    /// </summary>
    [McpServerTool(Name = "summarize_sequence", Title = "Sequence — Comprehensive Summary", ReadOnly = true)]
    [Description("Generate comprehensive summary statistics for a DNA/RNA sequence including composition, GC content, entropy, complexity, and Tm.")]
    public static SummarizeSequenceResult SummarizeSequence(
        [Description("The DNA or RNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var summary = SequenceStatistics.SummarizeNucleotideSequence(sequence);
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
    [McpServerTool(Name = "gc_content", Title = "Sequence — GC Content", ReadOnly = true)]
    [Description("Calculate the GC content (percentage of G and C nucleotides) of a DNA/RNA sequence.")]
    public static GcContentResult GcContent(
        [Description("The DNA or RNA sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var gcContent = global::Seqeron.Genomics.Core.SequenceExtensions.CalculateGcContentFast(sequence);
        int gcCount = sequence.Count(c => c == 'G' || c == 'C' || c == 'g' || c == 'c');
        int validCount = sequence.Count(c => c is 'A' or 'a' or 'T' or 't' or 'G' or 'g' or 'C' or 'c' or 'U' or 'u');
        return new GcContentResult(gcContent, gcCount, validCount);
    }

    /// <summary>
    /// Get the complement of a single nucleotide base.
    /// </summary>
    [McpServerTool(Name = "complement_base", Title = "Sequence — Complement Base", ReadOnly = true)]
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
    [McpServerTool(Name = "is_valid_dna", Title = "DNA — Quick Validate", ReadOnly = true)]
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
    [McpServerTool(Name = "is_valid_rna", Title = "RNA — Quick Validate", ReadOnly = true)]
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
    [McpServerTool(Name = "kmer_entropy", Title = "K-mer — Entropy", ReadOnly = true)]
    [Description("Calculate Shannon entropy based on k-mer frequencies. Higher values indicate more complexity.")]
    public static KmerEntropyResult KmerEntropy(
        [Description("The sequence to analyze")] string sequence,
        [Description("K-mer length (default: 2)")] int k = 2)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        var entropy = KmerAnalyzer.CalculateKmerEntropy(sequence, k);
        return new KmerEntropyResult(entropy, k);
    }

    /// <summary>
    /// Calculate linguistic complexity using SequenceComplexity class.
    /// </summary>
    [McpServerTool(Name = "complexity_linguistic", Title = "Complexity — Linguistic", ReadOnly = true)]
    [Description("Calculate DNA linguistic complexity as ratio of observed to possible subwords. LC = 1.0 for maximum complexity.")]
    public static ComplexityLinguisticResult ComplexityLinguistic(
        [Description("The DNA sequence to analyze")] string sequence,
        [Description("Maximum word length to consider (default: 10)")] int maxWordLength = 10)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (maxWordLength < 1)
            throw new ArgumentException("Max word length must be at least 1", nameof(maxWordLength));

        var complexity = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength);
        return new ComplexityLinguisticResult(complexity, maxWordLength);
    }

    /// <summary>
    /// Calculate Shannon entropy using SequenceComplexity class.
    /// </summary>
    [McpServerTool(Name = "complexity_shannon", Title = "Complexity — Shannon Entropy", ReadOnly = true)]
    [Description("Calculate DNA Shannon entropy (bits per base). Maximum entropy for DNA is 2 bits.")]
    public static ComplexityShannonResult ComplexityShannon(
        [Description("The DNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var entropy = SequenceComplexity.CalculateShannonEntropy(sequence);
        return new ComplexityShannonResult(entropy);
    }

    /// <summary>
    /// Calculate k-mer entropy using SequenceComplexity class.
    /// </summary>
    [McpServerTool(Name = "complexity_kmer_entropy", Title = "Complexity — K-mer Entropy", ReadOnly = true)]
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

        var entropy = SequenceComplexity.CalculateKmerEntropy(dna!, k);
        return new ComplexityKmerEntropyResult(entropy, k);
    }

    /// <summary>
    /// Calculate DUST score for low-complexity filtering.
    /// </summary>
    [McpServerTool(Name = "complexity_dust_score", Title = "Complexity — DUST Score", ReadOnly = true)]
    [Description("Calculate DUST score for low-complexity filtering (as used in BLAST). Higher scores indicate lower complexity.")]
    public static ComplexityDustScoreResult ComplexityDustScore(
        [Description("The DNA sequence to analyze")] string sequence,
        [Description("Word size for triplet counting (default: 3)")] int wordSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (wordSize < 1)
            throw new ArgumentException("Word size must be at least 1", nameof(wordSize));

        var dustScore = SequenceComplexity.CalculateDustScore(sequence, wordSize);
        return new ComplexityDustScoreResult(dustScore, wordSize);
    }

    /// <summary>
    /// Mask low-complexity regions using DUST algorithm.
    /// </summary>
    [McpServerTool(Name = "complexity_mask_low", Title = "Complexity — Mask Low-Complexity Regions", ReadOnly = true)]
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

        var masked = SequenceComplexity.MaskLowComplexity(dna!, windowSize, threshold, maskChar);
        return new ComplexityMaskLowResult(masked, sequence.Length, maskChar);
    }

    /// <summary>
    /// Estimate sequence complexity using compression ratio.
    /// </summary>
    [McpServerTool(Name = "complexity_compression_ratio", Title = "Complexity — Compression Ratio", ReadOnly = true)]
    [Description("Estimate sequence complexity using compression ratio. Lower ratios indicate more repetitive/less complex sequences.")]
    public static ComplexityCompressionRatioResult ComplexityCompressionRatio(
        [Description("The sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var ratio = SequenceComplexity.EstimateCompressionRatio(sequence);
        return new ComplexityCompressionRatioResult(ratio);
    }

    /// <summary>
    /// Count k-mer frequencies in a sequence.
    /// </summary>
    [McpServerTool(Name = "kmer_count", Title = "K-mer — Count Frequencies", ReadOnly = true)]
    [Description("Count k-mer (substring of length k) frequencies in a sequence. Returns a dictionary of k-mers and their counts.")]
    public static KmerCountResult KmerCount(
        [Description("The sequence to analyze")] string sequence,
        [Description("K-mer length (default: 3)")] int k = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        var counts = KmerAnalyzer.CountKmers(sequence, k);
        return new KmerCountResult(counts, k, counts.Count, counts.Values.Sum());
    }

    /// <summary>
    /// Calculate k-mer distance between two sequences.
    /// </summary>
    [McpServerTool(Name = "kmer_distance", Title = "K-mer — Distance Between Sequences", ReadOnly = true)]
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

        var distance = KmerAnalyzer.KmerDistance(sequence1, sequence2, k);
        return new KmerDistanceResult(distance, k);
    }

    /// <summary>
    /// Analyze k-mer composition of a sequence.
    /// </summary>
    [McpServerTool(Name = "kmer_analyze", Title = "K-mer — Comprehensive Analysis", ReadOnly = true)]
    [Description("Comprehensive k-mer analysis including statistics about frequency distribution, entropy, and unique k-mers.")]
    public static KmerAnalyzeResult KmerAnalyze(
        [Description("The sequence to analyze")] string sequence,
        [Description("K-mer length (default: 3)")] int k = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        if (k < 1)
            throw new ArgumentException("K must be at least 1", nameof(k));

        var stats = KmerAnalyzer.AnalyzeKmers(sequence, k);
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
    [McpServerTool(Name = "iupac_code", Title = "IUPAC — Get Ambiguity Code", ReadOnly = true)]
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
    [McpServerTool(Name = "iupac_match", Title = "IUPAC — Codes Match", ReadOnly = true)]
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
    [McpServerTool(Name = "iupac_matches", Title = "IUPAC — Nucleotide Matches Code", ReadOnly = true)]
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
    [McpServerTool(Name = "translate_dna", Title = "DNA — Translate to Protein", ReadOnly = true)]
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
    [McpServerTool(Name = "translate_rna", Title = "RNA — Translate to Protein", ReadOnly = true)]
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
