using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics;

namespace Seqeron.Mcp.Parsers.Tools;

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

    [McpServerTool(Name = "fasta_write")]
    [Description("Write sequence(s) to a FASTA file. Creates or overwrites the file at specified path.")]
    public static FastaWriteResult FastaWrite(
        [Description("File path to write FASTA output")] string filePath,
        [Description("Sequence identifier")] string id,
        [Description("DNA sequence")] string sequence,
        [Description("Optional sequence description")] string? description = null,
        [Description("Line width for sequence wrapping (default: 80)")] int lineWidth = 80)
    {
        if (string.IsNullOrEmpty(filePath))
            throw new ArgumentException("File path cannot be null or empty", nameof(filePath));
        if (string.IsNullOrEmpty(id))
            throw new ArgumentException("ID cannot be null or empty", nameof(id));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (lineWidth < 1)
            throw new ArgumentException("Line width must be at least 1", nameof(lineWidth));

        var dnaSeq = new DnaSequence(sequence);
        var entry = new FastaEntry(id, description, dnaSeq);
        FastaParser.WriteFile(filePath, new[] { entry }, lineWidth);

        return new FastaWriteResult(filePath, 1, sequence.Length);
    }

    // ========================
    // FASTQ Tools
    // ========================

    [McpServerTool(Name = "fastq_parse")]
    [Description("Parse FASTQ format string into sequence entries with quality scores. Supports Phred+33 and Phred+64 encodings.")]
    public static FastqParseResult FastqParse(
        [Description("FASTQ format content to parse")] string content,
        [Description("Quality encoding: 'phred33' (Sanger/Illumina 1.8+), 'phred64' (Illumina 1.3-1.7), or 'auto' (default)")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var results = records.Select(r => new FastqRecordResult(
            r.Id,
            r.Description,
            r.Sequence,
            r.QualityString,
            r.QualityScores.ToList(),
            r.Sequence.Length
        )).ToList();

        return new FastqParseResult(results, results.Count);
    }

    [McpServerTool(Name = "fastq_statistics")]
    [Description("Calculate quality statistics for FASTQ data. Returns read counts, quality metrics, and GC content.")]
    public static FastqStatisticsResult FastqStatistics(
        [Description("FASTQ format content to analyze")] string content,
        [Description("Quality encoding: 'phred33' (Sanger/Illumina 1.8+), 'phred64' (Illumina 1.3-1.7), or 'auto' (default)")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        return new FastqStatisticsResult(
            stats.TotalReads,
            stats.TotalBases,
            stats.MeanReadLength,
            stats.MeanQuality,
            stats.MinReadLength,
            stats.MaxReadLength,
            stats.Q20Percentage,
            stats.Q30Percentage,
            stats.GcContent);
    }

    [McpServerTool(Name = "fastq_filter")]
    [Description("Filter FASTQ reads by minimum average quality score. Returns reads meeting the quality threshold.")]
    public static FastqFilterResult FastqFilter(
        [Description("FASTQ format content to filter")] string content,
        [Description("Minimum average quality score threshold (e.g., 20 or 30)")] double minQuality,
        [Description("Quality encoding: 'phred33' (Sanger/Illumina 1.8+), 'phred64' (Illumina 1.3-1.7), or 'auto' (default)")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));
        if (minQuality < 0)
            throw new ArgumentException("Minimum quality must be non-negative", nameof(minQuality));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var filtered = FastqParser.FilterByQuality(records, minQuality).ToList();

        var results = filtered.Select(r => new FastqRecordResult(
            r.Id,
            r.Description,
            r.Sequence,
            r.QualityString,
            r.QualityScores.ToList(),
            r.Sequence.Length
        )).ToList();

        return new FastqFilterResult(
            results,
            results.Count,
            records.Count,
            records.Count > 0 ? (double)results.Count / records.Count * 100 : 0);
    }
}

// ========================
// Result Records
// ========================

public record FastaEntryResult(string Id, string? Description, string Sequence, int Length);
public record FastaParseResult(List<FastaEntryResult> Entries, int Count);
public record FastaFormatResult(string Fasta);
public record FastaWriteResult(string FilePath, int EntriesWritten, int TotalBases);
public record FastqRecordResult(string Id, string? Description, string Sequence, string QualityString, List<int> QualityScores, int Length);
public record FastqParseResult(List<FastqRecordResult> Entries, int Count);
public record FastqStatisticsResult(
    int TotalReads,
    long TotalBases,
    double MeanReadLength,
    double MeanQuality,
    int MinReadLength,
    int MaxReadLength,
    double Q20Percentage,
    double Q30Percentage,
    double GcContent);
public record FastqFilterResult(
    List<FastqRecordResult> Entries,
    int PassedCount,
    int TotalCount,
    double PassedPercentage);
