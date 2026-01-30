using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics;

/// <summary>
/// Parser for VCF (Variant Call Format) files.
/// VCF is the standard format for storing genetic variation data.
/// </summary>
public static class VcfParser
{
    #region Records

    /// <summary>Represents a VCF variant record</summary>
    public readonly record struct VcfRecord(
        string Chrom,
        int Pos,
        string Id,
        string Ref,
        string[] Alt,
        double? Qual,
        string[] Filter,
        IReadOnlyDictionary<string, string> Info,
        string[]? Format = null,
        IReadOnlyList<IReadOnlyDictionary<string, string>>? Samples = null);

    /// <summary>Represents VCF header metadata</summary>
    public readonly record struct VcfHeader(
        string FileFormat,
        IReadOnlyList<VcfInfoField> InfoFields,
        IReadOnlyList<VcfFormatField> FormatFields,
        IReadOnlyList<VcfFilterField> FilterFields,
        IReadOnlyList<string> SampleNames,
        IReadOnlyDictionary<string, string> OtherMetadata);

    /// <summary>VCF INFO field definition</summary>
    public readonly record struct VcfInfoField(string Id, string Number, string Type, string Description);

    /// <summary>VCF FORMAT field definition</summary>
    public readonly record struct VcfFormatField(string Id, string Number, string Type, string Description);

    /// <summary>VCF FILTER field definition</summary>
    public readonly record struct VcfFilterField(string Id, string Description);

    /// <summary>Variant type classification</summary>
    public enum VariantType
    {
        SNP,
        MNP,
        Insertion,
        Deletion,
        Complex,
        Symbolic,
        Unknown
    }

    #endregion

    #region Parsing Methods

    /// <summary>
    /// Parses VCF records from a file.
    /// </summary>
    public static IEnumerable<VcfRecord> ParseFile(string filePath)
    {
        if (string.IsNullOrEmpty(filePath) || !File.Exists(filePath))
            return Enumerable.Empty<VcfRecord>();

        var content = File.ReadAllText(filePath);
        return Parse(content);
    }

    /// <summary>
    /// Parses VCF file with header.
    /// </summary>
    public static (VcfHeader Header, IEnumerable<VcfRecord> Records) ParseFileWithHeader(string filePath)
    {
        if (string.IsNullOrEmpty(filePath) || !File.Exists(filePath))
            return (new VcfHeader("VCFv4.3", [], [], [], [], new Dictionary<string, string>()), []);

        var content = File.ReadAllText(filePath);
        return ParseWithHeader(content);
    }

    /// <summary>
    /// Parses VCF content with header.
    /// </summary>
    public static (VcfHeader Header, IEnumerable<VcfRecord> Records) ParseWithHeader(string content)
    {
        using var reader = new StringReader(content);
        var header = ParseHeader(reader);
        var sampleNames = header.SampleNames.ToArray();

        var records = new List<VcfRecord>();
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (!string.IsNullOrWhiteSpace(line) && !line.StartsWith('#'))
            {
                var record = ParseLine(line, sampleNames);
                if (record.HasValue)
                    records.Add(record.Value);
            }
        }

        return (header, records);
    }

    /// <summary>
    /// Parses VCF records from text content.
    /// </summary>
    public static IEnumerable<VcfRecord> Parse(string content)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        using var reader = new StringReader(content);
        string[]? sampleNames = null;

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line))
                continue;

            if (line.StartsWith("##"))
                continue;

            if (line.StartsWith("#CHROM"))
            {
                var fields = line.Split('\t');
                if (fields.Length > 9)
                    sampleNames = fields.Skip(9).ToArray();
                continue;
            }

            if (line.StartsWith('#'))
                continue;

            var record = ParseLine(line, sampleNames);
            if (record.HasValue)
                yield return record.Value;
        }
    }

    private static IEnumerable<VcfRecord> ParseFromReader(TextReader reader, string[]? sampleNames)
    {
        string? line;

        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line))
                continue;

            if (line.StartsWith("##"))
                continue;

            if (line.StartsWith("#CHROM"))
                continue;

            if (line.StartsWith('#'))
                continue;

            var record = ParseLine(line, sampleNames);
            if (record.HasValue)
                yield return record.Value;
        }
    }

    private static VcfHeader ParseHeader(TextReader reader)
    {
        string fileFormat = "VCFv4.3";
        var infoFields = new List<VcfInfoField>();
        var formatFields = new List<VcfFormatField>();
        var filterFields = new List<VcfFilterField>();
        var sampleNames = new List<string>();
        var otherMetadata = new Dictionary<string, string>();

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line))
                continue;

            if (!line.StartsWith('#'))
            {
                // Put line back - but we can't with TextReader, so we'll break
                break;
            }

            if (line.StartsWith("##fileformat="))
            {
                fileFormat = line["##fileformat=".Length..];
            }
            else if (line.StartsWith("##INFO="))
            {
                var info = ParseMetadataLine(line["##INFO=".Length..]);
                if (info != null)
                {
                    infoFields.Add(new VcfInfoField(
                        info.GetValueOrDefault("ID", ""),
                        info.GetValueOrDefault("Number", "."),
                        info.GetValueOrDefault("Type", "String"),
                        info.GetValueOrDefault("Description", "")));
                }
            }
            else if (line.StartsWith("##FORMAT="))
            {
                var format = ParseMetadataLine(line["##FORMAT=".Length..]);
                if (format != null)
                {
                    formatFields.Add(new VcfFormatField(
                        format.GetValueOrDefault("ID", ""),
                        format.GetValueOrDefault("Number", "."),
                        format.GetValueOrDefault("Type", "String"),
                        format.GetValueOrDefault("Description", "")));
                }
            }
            else if (line.StartsWith("##FILTER="))
            {
                var filter = ParseMetadataLine(line["##FILTER=".Length..]);
                if (filter != null)
                {
                    filterFields.Add(new VcfFilterField(
                        filter.GetValueOrDefault("ID", ""),
                        filter.GetValueOrDefault("Description", "")));
                }
            }
            else if (line.StartsWith("#CHROM"))
            {
                var fields = line.Split('\t');
                if (fields.Length > 9)
                    sampleNames.AddRange(fields.Skip(9));
                break;
            }
            else if (line.StartsWith("##"))
            {
                var eqIdx = line.IndexOf('=');
                if (eqIdx > 2)
                {
                    var key = line[2..eqIdx];
                    var value = line[(eqIdx + 1)..];
                    otherMetadata[key] = value;
                }
            }
        }

        return new VcfHeader(fileFormat, infoFields, formatFields, filterFields, sampleNames, otherMetadata);
    }

    private static Dictionary<string, string>? ParseMetadataLine(string content)
    {
        if (!content.StartsWith('<') || !content.EndsWith('>'))
            return null;

        content = content[1..^1];
        var result = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        var parts = new List<string>();
        var current = new StringBuilder();
        bool inQuotes = false;

        foreach (char c in content)
        {
            if (c == '"')
            {
                inQuotes = !inQuotes;
            }
            else if (c == ',' && !inQuotes)
            {
                parts.Add(current.ToString());
                current.Clear();
            }
            else
            {
                current.Append(c);
            }
        }
        if (current.Length > 0)
            parts.Add(current.ToString());

        foreach (var part in parts)
        {
            var eqIdx = part.IndexOf('=');
            if (eqIdx > 0)
            {
                var key = part[..eqIdx].Trim();
                var value = part[(eqIdx + 1)..].Trim().Trim('"');
                result[key] = value;
            }
        }

        return result;
    }

    private static VcfRecord? ParseLine(string line, string[]? sampleNames)
    {
        var fields = line.Split('\t');
        if (fields.Length < 8)
            return null;

        var chrom = fields[0];
        if (!int.TryParse(fields[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int pos))
            return null;

        var id = fields[2];
        var refAllele = fields[3];
        var altAlleles = fields[4].Split(',');

        double? qual = null;
        if (fields[5] != "." && double.TryParse(fields[5], NumberStyles.Float, CultureInfo.InvariantCulture, out double q))
            qual = q;

        var filters = fields[6] == "." ? Array.Empty<string>() : fields[6].Split(';');

        var info = ParseInfo(fields[7]);

        string[]? format = null;
        List<IReadOnlyDictionary<string, string>>? samples = null;

        if (fields.Length > 9 && sampleNames != null)
        {
            format = fields[8].Split(':');
            samples = new List<IReadOnlyDictionary<string, string>>();

            for (int i = 9; i < fields.Length && i - 9 < sampleNames.Length; i++)
            {
                var sampleValues = fields[i].Split(':');
                var sampleDict = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

                for (int j = 0; j < format.Length && j < sampleValues.Length; j++)
                {
                    sampleDict[format[j]] = sampleValues[j];
                }

                samples.Add(sampleDict);
            }
        }

        return new VcfRecord(chrom, pos, id, refAllele, altAlleles, qual, filters, info, format, samples);
    }

    private static Dictionary<string, string> ParseInfo(string infoString)
    {
        var info = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

        if (string.IsNullOrEmpty(infoString) || infoString == ".")
            return info;

        foreach (var part in infoString.Split(';'))
        {
            var eqIdx = part.IndexOf('=');
            if (eqIdx > 0)
            {
                info[part[..eqIdx]] = part[(eqIdx + 1)..];
            }
            else
            {
                info[part] = "true"; // Flag
            }
        }

        return info;
    }

    #endregion

    #region Variant Classification

    /// <summary>
    /// Classifies the variant type.
    /// </summary>
    public static VariantType ClassifyVariant(VcfRecord record, int altIndex = 0)
    {
        if (altIndex >= record.Alt.Length)
            return VariantType.Unknown;

        var alt = record.Alt[altIndex];

        // Symbolic alleles
        if (alt.StartsWith('<') || alt.StartsWith('[') || alt.StartsWith(']'))
            return VariantType.Symbolic;

        var refLen = record.Ref.Length;
        var altLen = alt.Length;

        if (refLen == altLen)
        {
            return refLen == 1 ? VariantType.SNP : VariantType.MNP;
        }

        if (refLen == 1 && altLen > 1)
            return VariantType.Insertion;

        if (refLen > 1 && altLen == 1)
            return VariantType.Deletion;

        return VariantType.Complex;
    }

    /// <summary>
    /// Determines if variant is a SNP.
    /// </summary>
    public static bool IsSNP(VcfRecord record) => ClassifyVariant(record) == VariantType.SNP;

    /// <summary>
    /// Determines if variant is an indel.
    /// </summary>
    public static bool IsIndel(VcfRecord record)
    {
        var type = ClassifyVariant(record);
        return type == VariantType.Insertion || type == VariantType.Deletion;
    }

    /// <summary>
    /// Gets the length of the variant.
    /// </summary>
    public static int GetVariantLength(VcfRecord record, int altIndex = 0)
    {
        if (altIndex >= record.Alt.Length)
            return 0;

        return Math.Abs(record.Alt[altIndex].Length - record.Ref.Length);
    }

    #endregion

    #region Filtering Methods

    /// <summary>
    /// Filters variants by chromosome.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterByChrom(IEnumerable<VcfRecord> records, string chrom)
    {
        return records.Where(r => r.Chrom.Equals(chrom, StringComparison.OrdinalIgnoreCase));
    }

    /// <summary>
    /// Filters variants by region.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterByRegion(
        IEnumerable<VcfRecord> records,
        string chrom,
        int start,
        int end)
    {
        return records.Where(r =>
            r.Chrom.Equals(chrom, StringComparison.OrdinalIgnoreCase) &&
            r.Pos >= start && r.Pos <= end);
    }

    /// <summary>
    /// Filters variants by quality score.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterByQuality(
        IEnumerable<VcfRecord> records,
        double minQuality)
    {
        return records.Where(r => r.Qual.HasValue && r.Qual.Value >= minQuality);
    }

    /// <summary>
    /// Filters to only passing variants.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterPassing(IEnumerable<VcfRecord> records)
    {
        return records.Where(r =>
            r.Filter.Length == 0 ||
            (r.Filter.Length == 1 && r.Filter[0].Equals("PASS", StringComparison.OrdinalIgnoreCase)));
    }

    /// <summary>
    /// Filters by variant type.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterByType(IEnumerable<VcfRecord> records, VariantType type)
    {
        return records.Where(r => ClassifyVariant(r) == type);
    }

    /// <summary>
    /// Filters SNPs only.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterSNPs(IEnumerable<VcfRecord> records)
    {
        return FilterByType(records, VariantType.SNP);
    }

    /// <summary>
    /// Filters indels only.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterIndels(IEnumerable<VcfRecord> records)
    {
        return records.Where(r => IsIndel(r));
    }

    /// <summary>
    /// Filters by INFO field.
    /// </summary>
    public static IEnumerable<VcfRecord> FilterByInfo(
        IEnumerable<VcfRecord> records,
        string infoKey,
        Func<string, bool> predicate)
    {
        return records.Where(r =>
            r.Info.TryGetValue(infoKey, out var value) && predicate(value));
    }

    #endregion

    #region Genotype Analysis

    /// <summary>
    /// Gets genotype for a sample.
    /// </summary>
    public static string? GetGenotype(VcfRecord record, int sampleIndex)
    {
        if (record.Samples == null || sampleIndex >= record.Samples.Count)
            return null;

        return record.Samples[sampleIndex].TryGetValue("GT", out var gt) ? gt : null;
    }

    /// <summary>
    /// Determines if sample is homozygous reference.
    /// </summary>
    public static bool IsHomRef(VcfRecord record, int sampleIndex)
    {
        var gt = GetGenotype(record, sampleIndex);
        if (gt == null) return false;
        return gt == "0/0" || gt == "0|0";
    }

    /// <summary>
    /// Determines if sample is homozygous alternate.
    /// </summary>
    public static bool IsHomAlt(VcfRecord record, int sampleIndex)
    {
        var gt = GetGenotype(record, sampleIndex);
        if (gt == null) return false;

        var alleles = gt.Replace('|', '/').Split('/');
        return alleles.Length == 2 &&
               alleles[0] == alleles[1] &&
               alleles[0] != "0" &&
               alleles[0] != ".";
    }

    /// <summary>
    /// Determines if sample is heterozygous.
    /// </summary>
    public static bool IsHet(VcfRecord record, int sampleIndex)
    {
        var gt = GetGenotype(record, sampleIndex);
        if (gt == null) return false;

        var alleles = gt.Replace('|', '/').Split('/');
        return alleles.Length == 2 && alleles[0] != alleles[1];
    }

    /// <summary>
    /// Gets allele depth for a sample.
    /// </summary>
    public static int[]? GetAlleleDepth(VcfRecord record, int sampleIndex)
    {
        if (record.Samples == null || sampleIndex >= record.Samples.Count)
            return null;

        if (!record.Samples[sampleIndex].TryGetValue("AD", out var ad))
            return null;

        return ad.Split(',')
            .Select(s => int.TryParse(s, out int v) ? v : 0)
            .ToArray();
    }

    /// <summary>
    /// Gets read depth for a sample.
    /// </summary>
    public static int? GetReadDepth(VcfRecord record, int sampleIndex)
    {
        if (record.Samples == null || sampleIndex >= record.Samples.Count)
            return null;

        if (record.Samples[sampleIndex].TryGetValue("DP", out var dp) &&
            int.TryParse(dp, out int depth))
            return depth;

        return null;
    }

    /// <summary>
    /// Calculates allele frequency.
    /// </summary>
    public static double? CalculateAlleleFrequency(IEnumerable<VcfRecord> records, int altIndex = 0)
    {
        int altCount = 0;
        int totalAlleles = 0;

        foreach (var record in records)
        {
            if (record.Samples == null) continue;

            for (int i = 0; i < record.Samples.Count; i++)
            {
                var gt = GetGenotype(record, i);
                if (gt == null || gt.Contains('.')) continue;

                var alleles = gt.Replace('|', '/').Split('/');
                foreach (var allele in alleles)
                {
                    totalAlleles++;
                    if (int.TryParse(allele, out int a) && a == altIndex + 1)
                        altCount++;
                }
            }
        }

        return totalAlleles > 0 ? (double)altCount / totalAlleles : null;
    }

    #endregion

    #region Statistics

    /// <summary>
    /// Calculates statistics for VCF records.
    /// </summary>
    public static VcfStatistics CalculateStatistics(IEnumerable<VcfRecord> records)
    {
        var recordsList = records.ToList();

        var typeCounts = recordsList
            .GroupBy(r => ClassifyVariant(r))
            .ToDictionary(g => g.Key, g => g.Count());

        var chromCounts = recordsList
            .GroupBy(r => r.Chrom, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(g => g.Key, g => g.Count());

        int passingCount = recordsList.Count(r =>
            r.Filter.Length == 0 ||
            (r.Filter.Length == 1 && r.Filter[0].Equals("PASS", StringComparison.OrdinalIgnoreCase)));

        var qualities = recordsList.Where(r => r.Qual.HasValue).Select(r => r.Qual!.Value).ToList();
        double? meanQuality = qualities.Count > 0 ? qualities.Average() : null;

        return new VcfStatistics(
            recordsList.Count,
            typeCounts.GetValueOrDefault(VariantType.SNP),
            typeCounts.GetValueOrDefault(VariantType.Insertion) + typeCounts.GetValueOrDefault(VariantType.Deletion),
            typeCounts.GetValueOrDefault(VariantType.Complex),
            passingCount,
            chromCounts,
            meanQuality);
    }

    /// <summary>VCF file statistics</summary>
    public readonly record struct VcfStatistics(
        int TotalVariants,
        int SnpCount,
        int IndelCount,
        int ComplexCount,
        int PassingCount,
        IReadOnlyDictionary<string, int> ChromosomeCounts,
        double? MeanQuality);

    /// <summary>
    /// Calculates transition/transversion ratio for SNPs.
    /// </summary>
    public static double? CalculateTiTvRatio(IEnumerable<VcfRecord> records)
    {
        int transitions = 0;
        int transversions = 0;

        var transitionPairs = new HashSet<string>(StringComparer.OrdinalIgnoreCase)
        {
            "AG", "GA", "CT", "TC"
        };

        foreach (var record in records.Where(r => IsSNP(r)))
        {
            var pair = $"{record.Ref}{record.Alt[0]}";
            if (transitionPairs.Contains(pair))
                transitions++;
            else
                transversions++;
        }

        return transversions > 0 ? (double)transitions / transversions : null;
    }

    #endregion

    #region Writing

    /// <summary>
    /// Writes VCF records to a file.
    /// </summary>
    public static void WriteToFile(
        string filePath,
        IEnumerable<VcfRecord> records,
        VcfHeader? header = null,
        string[]? sampleNames = null)
    {
        using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
        WriteToStream(writer, records, header, sampleNames);
    }

    /// <summary>
    /// Writes VCF records to a TextWriter.
    /// </summary>
    public static void WriteToStream(
        TextWriter writer,
        IEnumerable<VcfRecord> records,
        VcfHeader? header = null,
        string[]? sampleNames = null)
    {
        // Write header
        writer.WriteLine($"##fileformat={header?.FileFormat ?? "VCFv4.3"}");

        if (header.HasValue)
        {
            foreach (var info in header.Value.InfoFields)
            {
                writer.WriteLine($"##INFO=<ID={info.Id},Number={info.Number},Type={info.Type},Description=\"{info.Description}\">");
            }

            foreach (var format in header.Value.FormatFields)
            {
                writer.WriteLine($"##FORMAT=<ID={format.Id},Number={format.Number},Type={format.Type},Description=\"{format.Description}\">");
            }

            foreach (var filter in header.Value.FilterFields)
            {
                writer.WriteLine($"##FILTER=<ID={filter.Id},Description=\"{filter.Description}\">");
            }
        }

        // Column header
        var columns = new List<string> { "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" };

        var samples = sampleNames ?? header?.SampleNames.ToArray() ?? Array.Empty<string>();
        if (samples.Length > 0)
        {
            columns.Add("FORMAT");
            columns.AddRange(samples);
        }

        writer.WriteLine(string.Join('\t', columns));

        // Write records
        foreach (var record in records)
        {
            writer.WriteLine(FormatRecord(record, samples));
        }
    }

    private static string FormatRecord(VcfRecord record, string[] sampleNames)
    {
        var parts = new List<string>
        {
            record.Chrom,
            record.Pos.ToString(CultureInfo.InvariantCulture),
            string.IsNullOrEmpty(record.Id) ? "." : record.Id,
            record.Ref,
            string.Join(',', record.Alt),
            record.Qual?.ToString("F2", CultureInfo.InvariantCulture) ?? ".",
            record.Filter.Length == 0 ? "." : string.Join(';', record.Filter),
            FormatInfo(record.Info)
        };

        if (sampleNames.Length > 0 && record.Format != null && record.Samples != null)
        {
            parts.Add(string.Join(':', record.Format));

            for (int i = 0; i < sampleNames.Length && i < record.Samples.Count; i++)
            {
                var sampleValues = record.Format
                    .Select(f => record.Samples[i].TryGetValue(f, out var v) ? v : ".")
                    .ToArray();
                parts.Add(string.Join(':', sampleValues));
            }
        }

        return string.Join('\t', parts);
    }

    private static string FormatInfo(IReadOnlyDictionary<string, string> info)
    {
        if (info.Count == 0)
            return ".";

        var parts = info.Select(kvp =>
            kvp.Value == "true" ? kvp.Key : $"{kvp.Key}={kvp.Value}");

        return string.Join(';', parts);
    }

    #endregion

    #region Annotation Helpers

    /// <summary>
    /// Gets INFO field value.
    /// </summary>
    public static string? GetInfoValue(VcfRecord record, string key)
    {
        return record.Info.TryGetValue(key, out var value) ? value : null;
    }

    /// <summary>
    /// Gets INFO field as integer.
    /// </summary>
    public static int? GetInfoInt(VcfRecord record, string key)
    {
        var value = GetInfoValue(record, key);
        return value != null && int.TryParse(value, out int i) ? i : null;
    }

    /// <summary>
    /// Gets INFO field as double.
    /// </summary>
    public static double? GetInfoDouble(VcfRecord record, string key)
    {
        var value = GetInfoValue(record, key);
        return value != null && double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out double d) ? d : null;
    }

    /// <summary>
    /// Checks if INFO flag is set.
    /// </summary>
    public static bool HasInfoFlag(VcfRecord record, string key)
    {
        return record.Info.ContainsKey(key);
    }

    #endregion
}
