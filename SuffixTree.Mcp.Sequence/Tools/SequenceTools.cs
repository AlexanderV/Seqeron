using System.ComponentModel;
using ModelContextProtocol.Server;

namespace SuffixTree.Mcp.Sequence.Tools;

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

        var isValid = global::SuffixTree.Genomics.DnaSequence.TryCreate(sequence, out var dna);

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
        if (!global::SuffixTree.Genomics.DnaSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid DNA sequence", nameof(sequence));

        var reverseComplement = global::SuffixTree.Genomics.DnaSequence.GetReverseComplementString(sequence);
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

        var isValid = global::SuffixTree.Genomics.RnaSequence.TryCreate(sequence, out _);

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

        if (!global::SuffixTree.Genomics.DnaSequence.TryCreate(sequence, out var dna))
            throw new ArgumentException("Invalid DNA sequence", nameof(sequence));

        var rna = global::SuffixTree.Genomics.RnaSequence.FromDna(dna!);
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

        var isValid = global::SuffixTree.Genomics.ProteinSequence.TryCreate(sequence, out _);

        if (isValid)
        {
            return new ProteinValidateResult(true, sequence.Length, null);
        }
        else
        {
            // Find the invalid character for error message
            var upperSeq = sequence.ToUpperInvariant();
            var validChars = global::SuffixTree.Genomics.ProteinSequence.ValidCharacters;
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

        var comp = global::SuffixTree.Genomics.SequenceStatistics.CalculateNucleotideComposition(sequence);
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

        if (!global::SuffixTree.Genomics.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var comp = global::SuffixTree.Genomics.SequenceStatistics.CalculateAminoAcidComposition(sequence);
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

        if (!global::SuffixTree.Genomics.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var mw = global::SuffixTree.Genomics.SequenceStatistics.CalculateMolecularWeight(sequence);
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

        var mw = global::SuffixTree.Genomics.SequenceStatistics.CalculateNucleotideMolecularWeight(sequence, isDna);
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

        if (!global::SuffixTree.Genomics.ProteinSequence.TryCreate(sequence, out _))
            throw new ArgumentException("Invalid protein sequence", nameof(sequence));

        var pI = global::SuffixTree.Genomics.SequenceStatistics.CalculateIsoelectricPoint(sequence);
        return new IsoelectricPointResult(pI);
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
