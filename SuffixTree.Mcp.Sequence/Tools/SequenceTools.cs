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
}

/// <summary>
/// Result of dna_validate operation.
/// </summary>
public record DnaValidateResult(bool Valid, int Length, string? Error);
