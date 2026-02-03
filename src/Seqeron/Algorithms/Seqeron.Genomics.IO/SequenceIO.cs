using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.IO;

/// <summary>
/// Provides parsers and writers for various bioinformatics file formats.
/// Supports GenBank, EMBL, GFF, BED, and other common formats.
/// </summary>
public static class SequenceIO
{
    #region Records and Types

    /// <summary>
    /// Represents a feature annotation on a sequence.
    /// </summary>
    public readonly record struct SequenceFeature(
        string Type,
        int Start,
        int End,
        char Strand,
        IReadOnlyDictionary<string, string> Qualifiers);

    /// <summary>
    /// Represents a reference citation.
    /// </summary>
    public readonly record struct Reference(
        int Number,
        string Authors,
        string Title,
        string Journal,
        string? PubMed = null);

    /// <summary>
    /// Represents a complete annotated sequence record.
    /// </summary>
    public readonly record struct SequenceRecord(
        string Id,
        string? Accession,
        string? Description,
        string Sequence,
        string? Organism,
        string? Taxonomy,
        DateTime? Date,
        IReadOnlyList<SequenceFeature> Features,
        IReadOnlyList<Reference> References,
        IReadOnlyDictionary<string, string> Metadata);

    /// <summary>
    /// Represents a BED format entry.
    /// </summary>
    public readonly record struct BedEntry(
        string Chromosome,
        int Start,
        int End,
        string? Name = null,
        int? Score = null,
        char? Strand = null,
        int? ThickStart = null,
        int? ThickEnd = null,
        string? ItemRgb = null);

    /// <summary>
    /// Represents a GFF/GTF format entry.
    /// </summary>
    public readonly record struct GffEntry(
        string SeqId,
        string Source,
        string Type,
        int Start,
        int End,
        double? Score,
        char Strand,
        int? Phase,
        IReadOnlyDictionary<string, string> Attributes);

    /// <summary>
    /// Represents a SAM/BAM alignment record (simplified).
    /// </summary>
    public readonly record struct SamRecord(
        string ReadName,
        int Flag,
        string ReferenceName,
        int Position,
        int MappingQuality,
        string Cigar,
        string Sequence,
        string Quality);

    #endregion

    #region GenBank Parser

    /// <summary>
    /// Parses a GenBank format file.
    /// </summary>
    public static IEnumerable<SequenceRecord> ParseGenBank(TextReader reader)
    {
        string? line;
        var currentRecord = new StringBuilder();

        while ((line = reader.ReadLine()) != null)
        {
            currentRecord.AppendLine(line);

            if (line.StartsWith("//"))
            {
                var record = ParseGenBankRecord(currentRecord.ToString());
                if (record != null)
                    yield return record.Value;
                currentRecord.Clear();
            }
        }
    }

    /// <summary>
    /// Parses a GenBank format string.
    /// </summary>
    public static IEnumerable<SequenceRecord> ParseGenBankString(string content)
    {
        using var reader = new StringReader(content);
        foreach (var record in ParseGenBank(reader))
            yield return record;
    }

    private static SequenceRecord? ParseGenBankRecord(string content)
    {
        var lines = content.Split('\n');

        string? id = null;
        string? accession = null;
        string? description = null;
        string? organism = null;
        var taxonomy = new StringBuilder();
        DateTime? date = null;
        var features = new List<SequenceFeature>();
        var references = new List<Reference>();
        var metadata = new Dictionary<string, string>();
        var sequence = new StringBuilder();

        bool inFeatures = false;
        bool inSequence = false;
        bool inTaxonomy = false;
        var currentFeature = new StringBuilder();

        foreach (var rawLine in lines)
        {
            string line = rawLine.TrimEnd('\r');

            if (line.StartsWith("LOCUS"))
            {
                var parts = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length > 1) id = parts[1];
                if (parts.Length > 6)
                {
                    // Try parse date (GenBank format: DD-MMM-YYYY)
                    var dateStr = parts[^1];
                    if (DateTime.TryParseExact(dateStr, "dd-MMM-yyyy",
                        CultureInfo.InvariantCulture, DateTimeStyles.None, out var parsedDate))
                    {
                        date = parsedDate;
                    }
                }
            }
            else if (line.StartsWith("ACCESSION"))
            {
                accession = line.Substring(12).Trim().Split(' ')[0];
            }
            else if (line.StartsWith("DEFINITION"))
            {
                description = line.Substring(12).Trim();
            }
            else if (line.StartsWith("  ORGANISM"))
            {
                organism = line.Substring(12).Trim();
                inTaxonomy = true;
            }
            else if (inTaxonomy && line.StartsWith("            "))
            {
                taxonomy.Append(line.Trim());
            }
            else if (line.StartsWith("FEATURES"))
            {
                inFeatures = true;
                inTaxonomy = false;
            }
            else if (line.StartsWith("ORIGIN"))
            {
                inFeatures = false;
                inSequence = true;

                // Process last feature
                if (currentFeature.Length > 0)
                {
                    var feature = ParseGenBankFeature(currentFeature.ToString());
                    if (feature != null) features.Add(feature.Value);
                }
            }
            else if (inFeatures)
            {
                // Feature line - check if it's a new feature type
                // Feature types start at column 5 (0-indexed) and are followed by location
                if (line.Length > 21 && line.StartsWith("     ") && !line.StartsWith("                     "))
                {
                    // New feature type (starts at column 5, not indented qualifier)
                    if (currentFeature.Length > 0)
                    {
                        var feature = ParseGenBankFeature(currentFeature.ToString());
                        if (feature != null) features.Add(feature.Value);
                    }
                    currentFeature.Clear();
                }
                currentFeature.AppendLine(line);
            }
            else if (inSequence && !line.StartsWith("//"))
            {
                // Sequence line: remove numbers and spaces
                var seqLine = Regex.Replace(line, @"[\s\d]", "");
                sequence.Append(seqLine);
            }
        }

        if (id == null)
            return null;

        return new SequenceRecord(
            Id: id,
            Accession: accession,
            Description: description,
            Sequence: sequence.ToString().ToUpperInvariant(),
            Organism: organism,
            Taxonomy: taxonomy.ToString(),
            Date: date,
            Features: features,
            References: references,
            Metadata: metadata);
    }

    private static SequenceFeature? ParseGenBankFeature(string content)
    {
        var lines = content.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        if (lines.Length == 0) return null;

        // First line: type and location
        var firstLine = lines[0].TrimEnd('\r');
        var match = Regex.Match(firstLine, @"^\s{5}(\S+)\s+(.+)$");
        if (!match.Success) return null;

        string type = match.Groups[1].Value;
        string location = match.Groups[2].Value;

        // Parse location
        char strand = '+';
        int start = 0, end = 0;

        if (location.StartsWith("complement("))
        {
            strand = '-';
            location = location.Substring(11, location.Length - 12);
        }

        // Handle join() and other complex locations by taking the overall range
        var numbers = Regex.Matches(location, @"\d+").Cast<Match>().Select(m => int.Parse(m.Value)).ToList();
        if (numbers.Count >= 2)
        {
            start = numbers.Min();
            end = numbers.Max();
        }
        else if (numbers.Count == 1)
        {
            start = end = numbers[0];
        }

        // Parse qualifiers
        var qualifiers = new Dictionary<string, string>();
        var currentQualifier = new StringBuilder();
        string? currentKey = null;

        for (int i = 1; i < lines.Length; i++)
        {
            var line = lines[i].TrimEnd('\r').Trim();

            if (line.StartsWith("/"))
            {
                // Save previous qualifier
                if (currentKey != null)
                {
                    qualifiers[currentKey] = currentQualifier.ToString().Trim('"');
                }

                // Parse new qualifier
                var eqPos = line.IndexOf('=');
                if (eqPos > 0)
                {
                    currentKey = line.Substring(1, eqPos - 1);
                    currentQualifier.Clear();
                    currentQualifier.Append(line.Substring(eqPos + 1));
                }
                else
                {
                    currentKey = line.Substring(1);
                    currentQualifier.Clear();
                    currentQualifier.Append("true");
                }
            }
            else if (currentKey != null)
            {
                currentQualifier.Append(line);
            }
        }

        if (currentKey != null)
        {
            qualifiers[currentKey] = currentQualifier.ToString().Trim('"');
        }

        return new SequenceFeature(type, start, end, strand, qualifiers);
    }

    /// <summary>
    /// Writes a sequence record in GenBank format.
    /// </summary>
    public static string ToGenBank(SequenceRecord record)
    {
        var sb = new StringBuilder();

        // LOCUS line
        string dateStr = record.Date?.ToString("dd-MMM-yyyy", CultureInfo.InvariantCulture).ToUpperInvariant()
                        ?? DateTime.Now.ToString("dd-MMM-yyyy", CultureInfo.InvariantCulture).ToUpperInvariant();
        sb.AppendLine($"LOCUS       {record.Id,-16} {record.Sequence.Length,11} bp    DNA     linear   UNK {dateStr}");

        if (record.Description != null)
            sb.AppendLine($"DEFINITION  {record.Description}");

        if (record.Accession != null)
            sb.AppendLine($"ACCESSION   {record.Accession}");

        if (record.Organism != null)
        {
            sb.AppendLine("SOURCE      " + record.Organism);
            sb.AppendLine("  ORGANISM  " + record.Organism);
            if (record.Taxonomy != null)
                sb.AppendLine("            " + record.Taxonomy);
        }

        // Features
        if (record.Features.Count > 0)
        {
            sb.AppendLine("FEATURES             Location/Qualifiers");
            foreach (var feature in record.Features)
            {
                string location = feature.Strand == '-'
                    ? $"complement({feature.Start}..{feature.End})"
                    : $"{feature.Start}..{feature.End}";
                sb.AppendLine($"     {feature.Type,-15} {location}");

                foreach (var (key, value) in feature.Qualifiers)
                {
                    sb.AppendLine($"                     /{key}=\"{value}\"");
                }
            }
        }

        // Sequence
        sb.AppendLine("ORIGIN");
        for (int i = 0; i < record.Sequence.Length; i += 60)
        {
            sb.Append($"{i + 1,9}");
            for (int j = 0; j < 60 && i + j < record.Sequence.Length; j += 10)
            {
                sb.Append(' ');
                sb.Append(record.Sequence.Substring(i + j, Math.Min(10, record.Sequence.Length - i - j)).ToLowerInvariant());
            }
            sb.AppendLine();
        }
        sb.AppendLine("//");

        return sb.ToString();
    }

    #endregion

    #region EMBL Parser

    /// <summary>
    /// Parses an EMBL format file.
    /// </summary>
    public static IEnumerable<SequenceRecord> ParseEmbl(TextReader reader)
    {
        string? line;
        var currentRecord = new StringBuilder();

        while ((line = reader.ReadLine()) != null)
        {
            currentRecord.AppendLine(line);

            if (line.StartsWith("//"))
            {
                var record = ParseEmblRecord(currentRecord.ToString());
                if (record != null)
                    yield return record.Value;
                currentRecord.Clear();
            }
        }
    }

    /// <summary>
    /// Parses an EMBL format string.
    /// </summary>
    public static IEnumerable<SequenceRecord> ParseEmblString(string content)
    {
        using var reader = new StringReader(content);
        foreach (var record in ParseEmbl(reader))
            yield return record;
    }

    private static SequenceRecord? ParseEmblRecord(string content)
    {
        var lines = content.Split('\n');

        string? id = null;
        string? accession = null;
        string? description = null;
        string? organism = null;
        var features = new List<SequenceFeature>();
        var metadata = new Dictionary<string, string>();
        var sequence = new StringBuilder();

        bool inSequence = false;

        foreach (var rawLine in lines)
        {
            string line = rawLine.TrimEnd('\r');
            if (line.Length < 2) continue;

            string lineType = line.Length >= 2 ? line.Substring(0, 2) : "";
            string lineContent = line.Length > 5 ? line.Substring(5) : "";

            switch (lineType)
            {
                case "ID":
                    var idParts = lineContent.Split(';')[0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    if (idParts.Length > 0) id = idParts[0];
                    break;

                case "AC":
                    accession = lineContent.Split(';')[0].Trim();
                    break;

                case "DE":
                    description = (description ?? "") + lineContent.Trim() + " ";
                    break;

                case "OS":
                    organism = lineContent.Trim();
                    break;

                case "SQ":
                    inSequence = true;
                    break;

                case "  " when inSequence:
                    var seqLine = Regex.Replace(line, @"[\s\d]", "");
                    sequence.Append(seqLine);
                    break;
            }
        }

        if (id == null)
            return null;

        return new SequenceRecord(
            Id: id,
            Accession: accession,
            Description: description?.Trim(),
            Sequence: sequence.ToString().ToUpperInvariant(),
            Organism: organism,
            Taxonomy: null,
            Date: null,
            Features: features,
            References: new List<Reference>(),
            Metadata: metadata);
    }

    #endregion

    #region BED Parser

    /// <summary>
    /// Parses a BED format file.
    /// </summary>
    public static IEnumerable<BedEntry> ParseBed(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            line = line.Trim();
            if (string.IsNullOrEmpty(line) || line.StartsWith("#") || line.StartsWith("track") || line.StartsWith("browser"))
                continue;

            var entry = ParseBedLine(line);
            if (entry != null)
                yield return entry.Value;
        }
    }

    /// <summary>
    /// Parses a BED format string.
    /// </summary>
    public static IEnumerable<BedEntry> ParseBedString(string content)
    {
        using var reader = new StringReader(content);
        foreach (var entry in ParseBed(reader))
            yield return entry;
    }

    private static BedEntry? ParseBedLine(string line)
    {
        var fields = line.Split('\t');
        if (fields.Length < 3) return null;

        if (!int.TryParse(fields[1], out int start)) return null;
        if (!int.TryParse(fields[2], out int end)) return null;

        return new BedEntry(
            Chromosome: fields[0],
            Start: start,
            End: end,
            Name: fields.Length > 3 ? fields[3] : null,
            Score: fields.Length > 4 && int.TryParse(fields[4], out int score) ? score : null,
            Strand: fields.Length > 5 && fields[5].Length > 0 ? fields[5][0] : null,
            ThickStart: fields.Length > 6 && int.TryParse(fields[6], out int ts) ? ts : null,
            ThickEnd: fields.Length > 7 && int.TryParse(fields[7], out int te) ? te : null,
            ItemRgb: fields.Length > 8 ? fields[8] : null);
    }

    /// <summary>
    /// Writes BED entries to a string.
    /// </summary>
    public static string ToBed(IEnumerable<BedEntry> entries)
    {
        var sb = new StringBuilder();

        foreach (var entry in entries)
        {
            sb.Append($"{entry.Chromosome}\t{entry.Start}\t{entry.End}");

            if (entry.Name != null)
            {
                sb.Append($"\t{entry.Name}");
                sb.Append($"\t{entry.Score ?? 0}");
                sb.Append($"\t{entry.Strand ?? '.'}");

                if (entry.ThickStart != null)
                {
                    sb.Append($"\t{entry.ThickStart}\t{entry.ThickEnd ?? entry.End}");
                    if (entry.ItemRgb != null)
                        sb.Append($"\t{entry.ItemRgb}");
                }
            }

            sb.AppendLine();
        }

        return sb.ToString();
    }

    #endregion

    #region GFF Parser

    /// <summary>
    /// Parses a GFF3/GTF format file.
    /// </summary>
    public static IEnumerable<GffEntry> ParseGff(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            line = line.Trim();
            if (string.IsNullOrEmpty(line) || line.StartsWith("#"))
                continue;

            var entry = ParseGffLine(line);
            if (entry != null)
                yield return entry.Value;
        }
    }

    /// <summary>
    /// Parses a GFF3/GTF format string.
    /// </summary>
    public static IEnumerable<GffEntry> ParseGffString(string content)
    {
        using var reader = new StringReader(content);
        foreach (var entry in ParseGff(reader))
            yield return entry;
    }

    private static GffEntry? ParseGffLine(string line)
    {
        var fields = line.Split('\t');
        if (fields.Length < 9) return null;

        if (!int.TryParse(fields[3], out int start)) return null;
        if (!int.TryParse(fields[4], out int end)) return null;

        double? score = fields[5] != "." && double.TryParse(fields[5], out double s) ? s : null;
        char strand = fields[6].Length > 0 ? fields[6][0] : '.';
        int? phase = fields[7] != "." && int.TryParse(fields[7], out int p) ? p : null;

        // Parse attributes
        var attributes = new Dictionary<string, string>();
        if (fields.Length > 8 && !string.IsNullOrEmpty(fields[8]))
        {
            // Handle both GFF3 (key=value;) and GTF (key "value";) formats
            var attrParts = fields[8].Split(new[] { ';' }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var attr in attrParts)
            {
                var trimmed = attr.Trim();

                // GFF3 format: key=value
                var eqPos = trimmed.IndexOf('=');
                if (eqPos > 0)
                {
                    var key = trimmed.Substring(0, eqPos);
                    var value = Uri.UnescapeDataString(trimmed.Substring(eqPos + 1));
                    attributes[key] = value;
                }
                else
                {
                    // GTF format: key "value"
                    var match = Regex.Match(trimmed, @"(\S+)\s+""?([^""]+)""?");
                    if (match.Success)
                    {
                        attributes[match.Groups[1].Value] = match.Groups[2].Value;
                    }
                }
            }
        }

        return new GffEntry(
            SeqId: fields[0],
            Source: fields[1],
            Type: fields[2],
            Start: start,
            End: end,
            Score: score,
            Strand: strand,
            Phase: phase,
            Attributes: attributes);
    }

    /// <summary>
    /// Writes GFF entries to a string.
    /// </summary>
    public static string ToGff(IEnumerable<GffEntry> entries, bool gff3 = true)
    {
        var sb = new StringBuilder();

        if (gff3)
            sb.AppendLine("##gff-version 3");

        foreach (var entry in entries)
        {
            sb.Append($"{entry.SeqId}\t{entry.Source}\t{entry.Type}\t");
            sb.Append($"{entry.Start}\t{entry.End}\t");
            sb.Append(entry.Score?.ToString(CultureInfo.InvariantCulture) ?? ".");
            sb.Append($"\t{entry.Strand}\t");
            sb.Append(entry.Phase?.ToString() ?? ".");
            sb.Append('\t');

            // Attributes
            if (gff3)
            {
                sb.Append(string.Join(";", entry.Attributes.Select(kv =>
                    $"{kv.Key}={Uri.EscapeDataString(kv.Value)}")));
            }
            else
            {
                sb.Append(string.Join("; ", entry.Attributes.Select(kv =>
                    $"{kv.Key} \"{kv.Value}\"")));
            }

            sb.AppendLine();
        }

        return sb.ToString();
    }

    #endregion

    #region SAM Parser

    /// <summary>
    /// Parses a SAM format file (simplified, alignments only).
    /// </summary>
    public static IEnumerable<SamRecord> ParseSam(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (line.StartsWith("@")) continue; // Skip header

            var record = ParseSamLine(line);
            if (record != null)
                yield return record.Value;
        }
    }

    /// <summary>
    /// Parses a SAM format string.
    /// </summary>
    public static IEnumerable<SamRecord> ParseSamString(string content)
    {
        using var reader = new StringReader(content);
        foreach (var record in ParseSam(reader))
            yield return record;
    }

    private static SamRecord? ParseSamLine(string line)
    {
        var fields = line.Split('\t');
        if (fields.Length < 11) return null;

        if (!int.TryParse(fields[1], out int flag)) return null;
        if (!int.TryParse(fields[3], out int pos)) return null;
        if (!int.TryParse(fields[4], out int mapq)) return null;

        return new SamRecord(
            ReadName: fields[0],
            Flag: flag,
            ReferenceName: fields[2],
            Position: pos,
            MappingQuality: mapq,
            Cigar: fields[5],
            Sequence: fields[9],
            Quality: fields[10]);
    }

    /// <summary>
    /// Checks SAM flags for specific properties.
    /// </summary>
    public static bool IsPaired(int flag) => (flag & 0x1) != 0;
    public static bool IsProperPair(int flag) => (flag & 0x2) != 0;
    public static bool IsUnmapped(int flag) => (flag & 0x4) != 0;
    public static bool IsMateUnmapped(int flag) => (flag & 0x8) != 0;
    public static bool IsReverse(int flag) => (flag & 0x10) != 0;
    public static bool IsMateReverse(int flag) => (flag & 0x20) != 0;
    public static bool IsRead1(int flag) => (flag & 0x40) != 0;
    public static bool IsRead2(int flag) => (flag & 0x80) != 0;
    public static bool IsSecondary(int flag) => (flag & 0x100) != 0;
    public static bool IsQcFail(int flag) => (flag & 0x200) != 0;
    public static bool IsDuplicate(int flag) => (flag & 0x400) != 0;
    public static bool IsSupplementary(int flag) => (flag & 0x800) != 0;

    #endregion

    #region VCF Parser

    /// <summary>
    /// Represents a VCF variant record.
    /// </summary>
    public readonly record struct VcfRecord(
        string Chromosome,
        int Position,
        string Id,
        string Reference,
        IReadOnlyList<string> Alternatives,
        double? Quality,
        string Filter,
        IReadOnlyDictionary<string, string> Info,
        IReadOnlyList<string> SampleData);

    /// <summary>
    /// Parses a VCF format file.
    /// </summary>
    public static IEnumerable<VcfRecord> ParseVcf(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (line.StartsWith("#")) continue;

            var record = ParseVcfLine(line);
            if (record != null)
                yield return record.Value;
        }
    }

    /// <summary>
    /// Parses a VCF format string.
    /// </summary>
    public static IEnumerable<VcfRecord> ParseVcfString(string content)
    {
        using var reader = new StringReader(content);
        foreach (var record in ParseVcf(reader))
            yield return record;
    }

    private static VcfRecord? ParseVcfLine(string line)
    {
        var fields = line.Split('\t');
        if (fields.Length < 8) return null;

        if (!int.TryParse(fields[1], out int pos)) return null;

        double? qual = fields[5] != "." && double.TryParse(fields[5], out double q) ? q : null;

        var alts = fields[4].Split(',').ToList();

        // Parse INFO field
        var info = new Dictionary<string, string>();
        if (fields[7] != ".")
        {
            foreach (var item in fields[7].Split(';'))
            {
                var eqPos = item.IndexOf('=');
                if (eqPos > 0)
                    info[item.Substring(0, eqPos)] = item.Substring(eqPos + 1);
                else
                    info[item] = "true";
            }
        }

        // Sample data (if present)
        var sampleData = fields.Length > 9 ? fields.Skip(9).ToList() : new List<string>();

        return new VcfRecord(
            Chromosome: fields[0],
            Position: pos,
            Id: fields[2],
            Reference: fields[3],
            Alternatives: alts,
            Quality: qual,
            Filter: fields[6],
            Info: info,
            SampleData: sampleData);
    }

    #endregion

    #region Phylip Format

    /// <summary>
    /// Parses a Phylip format alignment (sequential or interleaved).
    /// </summary>
    public static IEnumerable<(string Name, string Sequence)> ParsePhylip(TextReader reader)
    {
        string? header = reader.ReadLine()?.Trim();
        if (header == null) yield break;

        var headerParts = header.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
        if (headerParts.Length < 2) yield break;

        int numSeq = int.Parse(headerParts[0]);
        // int seqLen = int.Parse(headerParts[1]); // Not strictly needed

        var names = new List<string>();
        var sequences = new List<StringBuilder>();

        string? line;
        int seqIndex = 0;
        bool firstBlock = true;

        while ((line = reader.ReadLine()) != null)
        {
            line = line.TrimEnd();
            if (string.IsNullOrWhiteSpace(line))
            {
                if (sequences.Count > 0)
                {
                    firstBlock = false;
                    seqIndex = 0;
                }
                continue;
            }

            if (firstBlock)
            {
                // First block: name + sequence
                string name = line.Length >= 10 ? line.Substring(0, 10).Trim() : line.Split(' ')[0];
                string seq = line.Length > 10 ? line.Substring(10).Replace(" ", "") : "";

                names.Add(name);
                sequences.Add(new StringBuilder(seq));
            }
            else
            {
                // Subsequent blocks: sequence continuation
                if (seqIndex < sequences.Count)
                {
                    sequences[seqIndex].Append(line.Replace(" ", ""));
                    seqIndex++;
                }
            }
        }

        for (int i = 0; i < names.Count; i++)
        {
            yield return (names[i], sequences[i].ToString());
        }
    }

    /// <summary>
    /// Writes sequences in Phylip format.
    /// </summary>
    public static string ToPhylip(IEnumerable<(string Name, string Sequence)> sequences, bool interleaved = false)
    {
        var seqList = sequences.ToList();
        if (seqList.Count == 0) return "";

        var sb = new StringBuilder();
        int seqLen = seqList[0].Sequence.Length;

        sb.AppendLine($" {seqList.Count} {seqLen}");

        if (!interleaved)
        {
            // Sequential format
            foreach (var (name, seq) in seqList)
            {
                string paddedName = name.PadRight(10).Substring(0, 10);
                sb.AppendLine($"{paddedName}{seq}");
            }
        }
        else
        {
            // Interleaved format
            int blockSize = 60;
            for (int pos = 0; pos < seqLen; pos += blockSize)
            {
                foreach (var (name, seq) in seqList)
                {
                    if (pos == 0)
                    {
                        string paddedName = name.PadRight(10).Substring(0, 10);
                        sb.Append(paddedName);
                    }
                    else
                    {
                        sb.Append("          "); // 10 spaces
                    }

                    int len = Math.Min(blockSize, seqLen - pos);
                    sb.AppendLine(seq.Substring(pos, len));
                }
                sb.AppendLine();
            }
        }

        return sb.ToString();
    }

    #endregion

    #region Clustal Format

    /// <summary>
    /// Parses a Clustal format alignment.
    /// </summary>
    public static IEnumerable<(string Name, string Sequence)> ParseClustal(TextReader reader)
    {
        var sequences = new Dictionary<string, StringBuilder>();
        var order = new List<string>();

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line) || line.StartsWith("CLUSTAL") ||
                line.Contains("*") || line.Contains(":") || line.Contains("."))
            {
                // Skip header and conservation lines
                continue;
            }

            var parts = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length >= 2)
            {
                string name = parts[0];
                string seq = parts[1];

                if (!sequences.ContainsKey(name))
                {
                    sequences[name] = new StringBuilder();
                    order.Add(name);
                }

                sequences[name].Append(seq);
            }
        }

        foreach (string name in order)
        {
            yield return (name, sequences[name].ToString());
        }
    }

    #endregion
}

