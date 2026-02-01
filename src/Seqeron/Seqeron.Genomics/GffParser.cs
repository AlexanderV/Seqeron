using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics;

/// <summary>
/// Parser for GFF3 (General Feature Format) and GTF (Gene Transfer Format) files.
/// These formats describe gene annotations and other genomic features.
/// </summary>
public static class GffParser
{
    #region Records

    /// <summary>Represents a GFF/GTF feature record</summary>
    public readonly record struct GffRecord(
        string Seqid,
        string Source,
        string Type,
        int Start,
        int End,
        double? Score,
        char Strand,
        int? Phase,
        IReadOnlyDictionary<string, string> Attributes);

    /// <summary>GFF file format version</summary>
    public enum GffFormat
    {
        GFF3,
        GTF,
        GFF2,
        Auto
    }

    /// <summary>Hierarchical gene model</summary>
    public readonly record struct GeneModel(
        GffRecord Gene,
        IReadOnlyList<GffRecord> Transcripts,
        IReadOnlyList<GffRecord> Exons,
        IReadOnlyList<GffRecord> CDS,
        IReadOnlyList<GffRecord> UTRs);

    #endregion

    #region Parsing Methods

    /// <summary>
    /// Parses GFF/GTF records from a file.
    /// </summary>
    public static IEnumerable<GffRecord> ParseFile(string filePath, GffFormat format = GffFormat.Auto)
    {
        if (string.IsNullOrEmpty(filePath) || !File.Exists(filePath))
            yield break;

        using var reader = new StreamReader(filePath);
        foreach (var record in Parse(reader, format))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses GFF/GTF records from text content.
    /// </summary>
    public static IEnumerable<GffRecord> Parse(string content, GffFormat format = GffFormat.Auto)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        using var reader = new StringReader(content);
        foreach (var record in Parse(reader, format))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses GFF/GTF records from a TextReader.
    /// </summary>
    public static IEnumerable<GffRecord> Parse(TextReader reader, GffFormat format = GffFormat.Auto)
    {
        var detectedFormat = format;
        string? line;

        while ((line = reader.ReadLine()) != null)
        {
            // Skip empty lines
            if (string.IsNullOrWhiteSpace(line))
                continue;

            // Handle directives
            if (line.StartsWith("##", StringComparison.Ordinal))
            {
                if (line.StartsWith("##gff-version", StringComparison.OrdinalIgnoreCase))
                {
                    if (line.Contains('3'))
                        detectedFormat = GffFormat.GFF3;
                    else
                        detectedFormat = GffFormat.GFF2;
                }
                continue;
            }

            // Skip comments
            if (line.StartsWith('#'))
                continue;

            var record = ParseLine(line, detectedFormat == GffFormat.Auto ? GffFormat.GFF3 : detectedFormat);
            if (record.HasValue)
                yield return record.Value;
        }
    }

    private static GffRecord? ParseLine(string line, GffFormat format)
    {
        var fields = line.Split('\t');
        if (fields.Length < 8)
            return null;

        var seqid = UnescapeGff(fields[0]);
        var source = UnescapeGff(fields[1]);
        var type = UnescapeGff(fields[2]);

        if (!int.TryParse(fields[3], NumberStyles.Integer, CultureInfo.InvariantCulture, out int start))
            return null;
        if (!int.TryParse(fields[4], NumberStyles.Integer, CultureInfo.InvariantCulture, out int end))
            return null;

        double? score = null;
        if (fields[5] != ".")
        {
            if (double.TryParse(fields[5], NumberStyles.Float, CultureInfo.InvariantCulture, out double s))
                score = s;
        }

        char strand = fields[6].Length > 0 ? fields[6][0] : '.';

        int? phase = null;
        if (fields[7] != "." && int.TryParse(fields[7], out int p))
            phase = p;

        var attributes = fields.Length > 8
            ? ParseAttributes(fields[8], format)
            : new Dictionary<string, string>();

        return new GffRecord(seqid, source, type, start, end, score, strand, phase, attributes);
    }

    private static Dictionary<string, string> ParseAttributes(string attrString, GffFormat format)
    {
        var attributes = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

        if (string.IsNullOrWhiteSpace(attrString))
            return attributes;

        if (format == GffFormat.GTF)
        {
            // GTF format: key "value"; key "value";
            var parts = attrString.Split(';', StringSplitOptions.RemoveEmptyEntries);
            foreach (var part in parts)
            {
                var trimmed = part.Trim();
                var spaceIdx = trimmed.IndexOf(' ');
                if (spaceIdx > 0)
                {
                    var key = trimmed[..spaceIdx].Trim();
                    var value = trimmed[(spaceIdx + 1)..].Trim().Trim('"');
                    attributes[key] = value;
                }
            }
        }
        else
        {
            // GFF3 format: key=value;key=value
            var parts = attrString.Split(';', StringSplitOptions.RemoveEmptyEntries);
            foreach (var part in parts)
            {
                var eqIdx = part.IndexOf('=');
                if (eqIdx > 0)
                {
                    var key = UnescapeGff(part[..eqIdx].Trim());
                    var value = UnescapeGff(part[(eqIdx + 1)..].Trim());
                    attributes[key] = value;
                }
            }
        }

        return attributes;
    }

    private static string UnescapeGff(string value)
    {
        return Uri.UnescapeDataString(value
            .Replace("%09", "\t")
            .Replace("%0A", "\n")
            .Replace("%0D", "\r")
            .Replace("%25", "%")
            .Replace("%3B", ";")
            .Replace("%3D", "=")
            .Replace("%26", "&")
            .Replace("%2C", ","));
    }

    #endregion

    #region Filtering Methods

    /// <summary>
    /// Filters records by feature type.
    /// </summary>
    public static IEnumerable<GffRecord> FilterByType(IEnumerable<GffRecord> records, params string[] types)
    {
        var typeSet = new HashSet<string>(types, StringComparer.OrdinalIgnoreCase);
        return records.Where(r => typeSet.Contains(r.Type));
    }

    /// <summary>
    /// Filters records by sequence ID.
    /// </summary>
    public static IEnumerable<GffRecord> FilterBySeqid(IEnumerable<GffRecord> records, string seqid)
    {
        return records.Where(r => r.Seqid.Equals(seqid, StringComparison.OrdinalIgnoreCase));
    }

    /// <summary>
    /// Filters records by genomic region.
    /// </summary>
    public static IEnumerable<GffRecord> FilterByRegion(
        IEnumerable<GffRecord> records,
        string seqid,
        int start,
        int end)
    {
        return records.Where(r =>
            r.Seqid.Equals(seqid, StringComparison.OrdinalIgnoreCase) &&
            r.Start <= end && r.End >= start);
    }

    /// <summary>
    /// Gets all gene features.
    /// </summary>
    public static IEnumerable<GffRecord> GetGenes(IEnumerable<GffRecord> records)
    {
        return FilterByType(records, "gene");
    }

    /// <summary>
    /// Gets all exon features.
    /// </summary>
    public static IEnumerable<GffRecord> GetExons(IEnumerable<GffRecord> records)
    {
        return FilterByType(records, "exon");
    }

    /// <summary>
    /// Gets all CDS features.
    /// </summary>
    public static IEnumerable<GffRecord> GetCDS(IEnumerable<GffRecord> records)
    {
        return FilterByType(records, "CDS");
    }

    #endregion

    #region Gene Model Building

    /// <summary>
    /// Builds hierarchical gene models from GFF records.
    /// </summary>
    public static IEnumerable<GeneModel> BuildGeneModels(IEnumerable<GffRecord> records)
    {
        var recordsList = records.ToList();

        // Index records by Parent attribute
        var childrenByParent = new Dictionary<string, List<GffRecord>>(StringComparer.OrdinalIgnoreCase);
        var recordsById = new Dictionary<string, GffRecord>(StringComparer.OrdinalIgnoreCase);

        foreach (var record in recordsList)
        {
            if (record.Attributes.TryGetValue("ID", out var id))
            {
                recordsById[id] = record;
            }

            if (record.Attributes.TryGetValue("Parent", out var parent))
            {
                if (!childrenByParent.ContainsKey(parent))
                    childrenByParent[parent] = new List<GffRecord>();
                childrenByParent[parent].Add(record);
            }
        }

        // Find top-level genes
        var genes = recordsList.Where(r =>
            r.Type.Equals("gene", StringComparison.OrdinalIgnoreCase)).ToList();

        foreach (var gene in genes)
        {
            if (!gene.Attributes.TryGetValue("ID", out var geneId))
                continue;

            var transcripts = new List<GffRecord>();
            var exons = new List<GffRecord>();
            var cds = new List<GffRecord>();
            var utrs = new List<GffRecord>();

            // Get direct children of gene
            if (childrenByParent.TryGetValue(geneId, out var geneChildren))
            {
                foreach (var child in geneChildren)
                {
                    var type = child.Type.ToLowerInvariant();

                    if (type == "mrna" || type == "transcript" || type == "ncrna")
                    {
                        transcripts.Add(child);

                        // Get children of transcript
                        if (child.Attributes.TryGetValue("ID", out var transcriptId) &&
                            childrenByParent.TryGetValue(transcriptId, out var transcriptChildren))
                        {
                            foreach (var tChild in transcriptChildren)
                            {
                                var tType = tChild.Type.ToLowerInvariant();
                                if (tType == "exon")
                                    exons.Add(tChild);
                                else if (tType == "cds")
                                    cds.Add(tChild);
                                else if (tType.Contains("utr"))
                                    utrs.Add(tChild);
                            }
                        }
                    }
                    else if (type == "exon")
                        exons.Add(child);
                    else if (type == "cds")
                        cds.Add(child);
                    else if (type.Contains("utr"))
                        utrs.Add(child);
                }
            }

            yield return new GeneModel(gene, transcripts, exons, cds, utrs);
        }
    }

    /// <summary>
    /// Gets attribute value from a record.
    /// </summary>
    public static string? GetAttribute(GffRecord record, string attributeName)
    {
        return record.Attributes.TryGetValue(attributeName, out var value) ? value : null;
    }

    /// <summary>
    /// Gets gene name from a record.
    /// </summary>
    public static string? GetGeneName(GffRecord record)
    {
        return GetAttribute(record, "gene_name") ??
               GetAttribute(record, "Name") ??
               GetAttribute(record, "gene_id") ??
               GetAttribute(record, "ID");
    }

    #endregion

    #region Statistics

    /// <summary>
    /// Calculates statistics for GFF records.
    /// </summary>
    public static GffStatistics CalculateStatistics(IEnumerable<GffRecord> records)
    {
        var recordsList = records.ToList();
        var featureTypeCounts = recordsList
            .GroupBy(r => r.Type, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(g => g.Key, g => g.Count());

        var seqids = recordsList.Select(r => r.Seqid).Distinct().ToList();
        var sources = recordsList.Select(r => r.Source).Distinct().ToList();

        int geneCount = recordsList.Count(r =>
            r.Type.Equals("gene", StringComparison.OrdinalIgnoreCase));
        int exonCount = recordsList.Count(r =>
            r.Type.Equals("exon", StringComparison.OrdinalIgnoreCase));

        return new GffStatistics(
            recordsList.Count,
            featureTypeCounts,
            seqids,
            sources,
            geneCount,
            exonCount);
    }

    /// <summary>GFF file statistics</summary>
    public readonly record struct GffStatistics(
        int TotalFeatures,
        IReadOnlyDictionary<string, int> FeatureTypeCounts,
        IReadOnlyList<string> SequenceIds,
        IReadOnlyList<string> Sources,
        int GeneCount,
        int ExonCount);

    #endregion

    #region Writing

    /// <summary>
    /// Writes GFF records to a file.
    /// </summary>
    public static void WriteToFile(string filePath, IEnumerable<GffRecord> records, GffFormat format = GffFormat.GFF3)
    {
        using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
        WriteToStream(writer, records, format);
    }

    /// <summary>
    /// Writes GFF records to a TextWriter.
    /// </summary>
    public static void WriteToStream(TextWriter writer, IEnumerable<GffRecord> records, GffFormat format = GffFormat.GFF3)
    {
        // Write header
        if (format == GffFormat.GFF3)
            writer.WriteLine("##gff-version 3");

        foreach (var record in records)
        {
            writer.WriteLine(FormatRecord(record, format));
        }
    }

    private static string FormatRecord(GffRecord record, GffFormat format)
    {
        var sb = new StringBuilder();

        sb.Append(EscapeGff(record.Seqid)).Append('\t');
        sb.Append(EscapeGff(record.Source)).Append('\t');
        sb.Append(EscapeGff(record.Type)).Append('\t');
        sb.Append(record.Start.ToString(CultureInfo.InvariantCulture)).Append('\t');
        sb.Append(record.End.ToString(CultureInfo.InvariantCulture)).Append('\t');
        sb.Append(record.Score?.ToString("F2", CultureInfo.InvariantCulture) ?? ".").Append('\t');
        sb.Append(record.Strand).Append('\t');
        sb.Append(record.Phase?.ToString(CultureInfo.InvariantCulture) ?? ".").Append('\t');

        // Format attributes
        if (record.Attributes.Count > 0)
        {
            var attrs = new List<string>();
            foreach (var (key, value) in record.Attributes)
            {
                if (format == GffFormat.GTF)
                    attrs.Add($"{key} \"{value}\"");
                else
                    attrs.Add($"{EscapeGff(key)}={EscapeGff(value)}");
            }
            sb.Append(string.Join(format == GffFormat.GTF ? "; " : ";", attrs));
        }

        return sb.ToString();
    }

    private static string EscapeGff(string value)
    {
        return value
            .Replace("%", "%25")
            .Replace("\t", "%09")
            .Replace("\n", "%0A")
            .Replace("\r", "%0D")
            .Replace(";", "%3B")
            .Replace("=", "%3D")
            .Replace("&", "%26")
            .Replace(",", "%2C");
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Extracts sequence for a feature from reference genome.
    /// </summary>
    public static string ExtractSequence(GffRecord record, string referenceSequence)
    {
        if (string.IsNullOrEmpty(referenceSequence))
            return "";

        int start = Math.Max(0, record.Start - 1); // GFF is 1-based
        int end = Math.Min(referenceSequence.Length, record.End);

        if (start >= end)
            return "";

        var sequence = referenceSequence[start..end];

        // Reverse complement if on minus strand
        if (record.Strand == '-')
        {
            return new DnaSequence(sequence).ReverseComplement().Sequence;
        }

        return sequence;
    }

    /// <summary>
    /// Merges overlapping features.
    /// </summary>
    public static IEnumerable<GffRecord> MergeOverlapping(IEnumerable<GffRecord> records)
    {
        var sorted = records
            .OrderBy(r => r.Seqid)
            .ThenBy(r => r.Start)
            .ToList();

        if (sorted.Count == 0)
            yield break;

        var current = sorted[0];

        for (int i = 1; i < sorted.Count; i++)
        {
            var next = sorted[i];

            if (next.Seqid == current.Seqid &&
                next.Start <= current.End + 1)
            {
                // Merge
                current = new GffRecord(
                    current.Seqid,
                    current.Source,
                    current.Type,
                    current.Start,
                    Math.Max(current.End, next.End),
                    null,
                    current.Strand,
                    null,
                    current.Attributes);
            }
            else
            {
                yield return current;
                current = next;
            }
        }

        yield return current;
    }

    #endregion
}
