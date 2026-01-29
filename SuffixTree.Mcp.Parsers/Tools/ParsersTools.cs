using System.ComponentModel;
using ModelContextProtocol.Server;
using SuffixTree.Genomics;

namespace SuffixTree.Mcp.Parsers.Tools;

[McpServerToolType]
public static class ParsersTools
{
    // ========================
    // FASTA Tools
    // ========================

    [McpServerTool(Name = "fasta_parse")]
    [Description("Parse FASTA format string into sequence entries. Returns list of sequences with their IDs and descriptions.")]
    public static FastaParseResult FastaParse(
        [Description("FASTA format content to parse")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var entries = FastaParser.Parse(content).ToList();
        var results = entries.Select(e => new FastaEntryResult(
            e.Id,
            e.Description,
            e.Sequence.Sequence,
            e.Sequence.Length
        )).ToList();

        return new FastaParseResult(results, results.Count);
    }

    [McpServerTool(Name = "fasta_format")]
    [Description("Format sequence(s) to FASTA string. Accepts ID, optional description, and sequence.")]
    public static FastaFormatResult FastaFormat(
        [Description("Sequence identifier")] string id,
        [Description("DNA sequence")] string sequence,
        [Description("Optional sequence description")] string? description = null,
        [Description("Line width for sequence wrapping (default: 80)")] int lineWidth = 80)
    {
        if (string.IsNullOrEmpty(id))
            throw new ArgumentException("ID cannot be null or empty", nameof(id));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (lineWidth < 1)
            throw new ArgumentException("Line width must be at least 1", nameof(lineWidth));

        var dnaSeq = new DnaSequence(sequence);
        var entry = new FastaEntry(id, description, dnaSeq);
        var fasta = FastaParser.ToFasta(new[] { entry }, lineWidth);

        return new FastaFormatResult(fasta.TrimEnd());
    }
}

// ========================
// Result Records
// ========================

public record FastaEntryResult(string Id, string? Description, string Sequence, int Length);
public record FastaParseResult(List<FastaEntryResult> Entries, int Count);
public record FastaFormatResult(string Fasta);
