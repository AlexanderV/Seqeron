using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics;

/// <summary>
/// Parser for BED (Browser Extensible Data) format files.
/// BED format is used to describe genomic regions.
/// </summary>
public static class BedParser
{
    #region Records

    /// <summary>
    /// Represents a BED record with 3-12 columns.
    /// BED3: chrom, start, end
    /// BED6: + name, score, strand
    /// BED12: + thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
    /// </summary>
    public readonly record struct BedRecord(
        string Chrom,
        int ChromStart,
        int ChromEnd,
        string? Name = null,
        int? Score = null,
        char? Strand = null,
        int? ThickStart = null,
        int? ThickEnd = null,
        string? ItemRgb = null,
        int? BlockCount = null,
        int[]? BlockSizes = null,
        int[]? BlockStarts = null)
    {
        /// <summary>Length of the feature</summary>
        public int Length => ChromEnd - ChromStart;

        /// <summary>Whether this is a BED12 record with blocks</summary>
        public bool HasBlocks => BlockCount.HasValue && BlockSizes != null && BlockStarts != null;
    }

    /// <summary>BED format type based on column count</summary>
    public enum BedFormat
    {
        BED3 = 3,
        BED4 = 4,
        BED5 = 5,
        BED6 = 6,
        BED12 = 12,
        Auto = 0
    }

    /// <summary>Represents a genomic interval</summary>
    public readonly record struct GenomicInterval(string Chrom, int Start, int End)
    {
        public int Length => End - Start;

        public bool Overlaps(GenomicInterval other)
        {
            return Chrom.Equals(other.Chrom, StringComparison.OrdinalIgnoreCase) &&
                   Start < other.End && End > other.Start;
        }

        public GenomicInterval? Intersect(GenomicInterval other)
        {
            if (!Overlaps(other))
                return null;
            return new GenomicInterval(Chrom, Math.Max(Start, other.Start), Math.Min(End, other.End));
        }

        public GenomicInterval Union(GenomicInterval other)
        {
            return new GenomicInterval(Chrom, Math.Min(Start, other.Start), Math.Max(End, other.End));
        }
    }

    #endregion

    #region Parsing Methods

    /// <summary>
    /// Parses BED records from a file.
    /// </summary>
    public static IEnumerable<BedRecord> ParseFile(string filePath, BedFormat format = BedFormat.Auto)
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
    /// Parses BED records from text content.
    /// </summary>
    public static IEnumerable<BedRecord> Parse(string content, BedFormat format = BedFormat.Auto)
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
    /// Parses BED records from a TextReader.
    /// </summary>
    public static IEnumerable<BedRecord> Parse(TextReader reader, BedFormat format = BedFormat.Auto)
    {
        string? line;

        while ((line = reader.ReadLine()) != null)
        {
            // Skip empty lines and header lines
            if (string.IsNullOrWhiteSpace(line))
                continue;
            if (line.StartsWith("track ", StringComparison.OrdinalIgnoreCase))
                continue;
            if (line.StartsWith("browser ", StringComparison.OrdinalIgnoreCase))
                continue;
            if (line.StartsWith('#'))
                continue;

            var record = ParseLine(line, format);
            if (record.HasValue)
                yield return record.Value;
        }
    }

    private static BedRecord? ParseLine(string line, BedFormat format)
    {
        var fields = line.Split('\t');
        if (fields.Length < 3)
        {
            // Try space separation
            fields = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
        }

        if (fields.Length < 3)
            return null;

        // Parse required BED3 fields
        var chrom = fields[0];
        if (!int.TryParse(fields[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int chromStart))
            return null;
        if (!int.TryParse(fields[2], NumberStyles.Integer, CultureInfo.InvariantCulture, out int chromEnd))
            return null;

        // Parse optional fields
        string? name = fields.Length > 3 ? fields[3] : null;

        int? score = null;
        if (fields.Length > 4 && int.TryParse(fields[4], out int s))
            score = Math.Clamp(s, 0, 1000); // BED score is 0-1000

        char? strand = null;
        if (fields.Length > 5 && fields[5].Length > 0)
        {
            var strandChar = fields[5][0];
            if (strandChar == '+' || strandChar == '-' || strandChar == '.')
                strand = strandChar;
        }

        int? thickStart = null;
        if (fields.Length > 6 && int.TryParse(fields[6], out int ts))
            thickStart = ts;

        int? thickEnd = null;
        if (fields.Length > 7 && int.TryParse(fields[7], out int te))
            thickEnd = te;

        string? itemRgb = fields.Length > 8 ? fields[8] : null;

        int? blockCount = null;
        int[]? blockSizes = null;
        int[]? blockStarts = null;

        if (fields.Length >= 12)
        {
            if (int.TryParse(fields[9], out int bc))
            {
                blockCount = bc;
                blockSizes = ParseIntList(fields[10]);
                blockStarts = ParseIntList(fields[11]);
            }
        }

        return new BedRecord(
            chrom, chromStart, chromEnd, name, score, strand,
            thickStart, thickEnd, itemRgb,
            blockCount, blockSizes, blockStarts);
    }

    private static int[] ParseIntList(string value)
    {
        return value.Split(',', StringSplitOptions.RemoveEmptyEntries)
            .Select(s => int.TryParse(s.Trim(), out int v) ? v : 0)
            .ToArray();
    }

    #endregion

    #region Filtering Methods

    /// <summary>
    /// Filters records by chromosome.
    /// </summary>
    public static IEnumerable<BedRecord> FilterByChrom(IEnumerable<BedRecord> records, string chrom)
    {
        return records.Where(r => r.Chrom.Equals(chrom, StringComparison.OrdinalIgnoreCase));
    }

    /// <summary>
    /// Filters records by region.
    /// </summary>
    public static IEnumerable<BedRecord> FilterByRegion(
        IEnumerable<BedRecord> records,
        string chrom,
        int start,
        int end)
    {
        return records.Where(r =>
            r.Chrom.Equals(chrom, StringComparison.OrdinalIgnoreCase) &&
            r.ChromStart < end && r.ChromEnd > start);
    }

    /// <summary>
    /// Filters records by strand.
    /// </summary>
    public static IEnumerable<BedRecord> FilterByStrand(IEnumerable<BedRecord> records, char strand)
    {
        return records.Where(r => r.Strand == strand);
    }

    /// <summary>
    /// Filters records by minimum length.
    /// </summary>
    public static IEnumerable<BedRecord> FilterByLength(
        IEnumerable<BedRecord> records,
        int minLength,
        int? maxLength = null)
    {
        return records.Where(r => r.Length >= minLength && (!maxLength.HasValue || r.Length <= maxLength.Value));
    }

    /// <summary>
    /// Filters records by score.
    /// </summary>
    public static IEnumerable<BedRecord> FilterByScore(
        IEnumerable<BedRecord> records,
        int minScore,
        int? maxScore = null)
    {
        return records.Where(r =>
            r.Score.HasValue &&
            r.Score.Value >= minScore &&
            (!maxScore.HasValue || r.Score.Value <= maxScore.Value));
    }

    #endregion

    #region Interval Operations

    /// <summary>
    /// Converts BED records to genomic intervals.
    /// </summary>
    public static IEnumerable<GenomicInterval> ToIntervals(IEnumerable<BedRecord> records)
    {
        return records.Select(r => new GenomicInterval(r.Chrom, r.ChromStart, r.ChromEnd));
    }

    /// <summary>
    /// Merges overlapping intervals.
    /// </summary>
    public static IEnumerable<BedRecord> MergeOverlapping(IEnumerable<BedRecord> records)
    {
        var sorted = records
            .OrderBy(r => r.Chrom)
            .ThenBy(r => r.ChromStart)
            .ToList();

        if (sorted.Count == 0)
            yield break;

        var current = sorted[0];

        for (int i = 1; i < sorted.Count; i++)
        {
            var next = sorted[i];

            if (next.Chrom.Equals(current.Chrom, StringComparison.OrdinalIgnoreCase) &&
                next.ChromStart <= current.ChromEnd)
            {
                // Merge
                current = current with { ChromEnd = Math.Max(current.ChromEnd, next.ChromEnd) };
            }
            else
            {
                yield return current;
                current = next;
            }
        }

        yield return current;
    }

    /// <summary>
    /// Finds intersections between two sets of intervals.
    /// </summary>
    public static IEnumerable<BedRecord> Intersect(
        IEnumerable<BedRecord> recordsA,
        IEnumerable<BedRecord> recordsB)
    {
        var listA = recordsA.ToList();
        var listB = recordsB
            .GroupBy(r => r.Chrom, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(g => g.Key, g => g.OrderBy(r => r.ChromStart).ToList(), StringComparer.OrdinalIgnoreCase);

        foreach (var a in listA)
        {
            if (!listB.TryGetValue(a.Chrom, out var bList))
                continue;

            foreach (var b in bList)
            {
                if (b.ChromStart >= a.ChromEnd)
                    break;
                if (b.ChromEnd <= a.ChromStart)
                    continue;

                // Found intersection
                int start = Math.Max(a.ChromStart, b.ChromStart);
                int end = Math.Min(a.ChromEnd, b.ChromEnd);

                yield return new BedRecord(
                    a.Chrom,
                    start,
                    end,
                    a.Name ?? b.Name,
                    a.Score ?? b.Score,
                    a.Strand ?? b.Strand);
            }
        }
    }

    /// <summary>
    /// Subtracts intervals in B from intervals in A.
    /// </summary>
    public static IEnumerable<BedRecord> Subtract(
        IEnumerable<BedRecord> recordsA,
        IEnumerable<BedRecord> recordsB)
    {
        var listA = recordsA.ToList();
        var listB = recordsB
            .GroupBy(r => r.Chrom, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(g => g.Key, g => g.OrderBy(r => r.ChromStart).ToList(), StringComparer.OrdinalIgnoreCase);

        foreach (var a in listA)
        {
            if (!listB.TryGetValue(a.Chrom, out var bList))
            {
                yield return a;
                continue;
            }

            var intervals = new List<(int Start, int End)> { (a.ChromStart, a.ChromEnd) };

            foreach (var b in bList)
            {
                if (b.ChromStart >= a.ChromEnd)
                    break;
                if (b.ChromEnd <= a.ChromStart)
                    continue;

                var newIntervals = new List<(int Start, int End)>();
                foreach (var (start, end) in intervals)
                {
                    if (b.ChromEnd <= start || b.ChromStart >= end)
                    {
                        newIntervals.Add((start, end));
                    }
                    else
                    {
                        if (b.ChromStart > start)
                            newIntervals.Add((start, b.ChromStart));
                        if (b.ChromEnd < end)
                            newIntervals.Add((b.ChromEnd, end));
                    }
                }
                intervals = newIntervals;
            }

            foreach (var (start, end) in intervals)
            {
                yield return a with { ChromStart = start, ChromEnd = end };
            }
        }
    }

    /// <summary>
    /// Expands intervals by a specified amount on each side.
    /// </summary>
    public static IEnumerable<BedRecord> ExpandIntervals(
        IEnumerable<BedRecord> records,
        int upstream,
        int downstream)
    {
        foreach (var record in records)
        {
            int newStart, newEnd;

            if (record.Strand == '-')
            {
                newStart = Math.Max(0, record.ChromStart - downstream);
                newEnd = record.ChromEnd + upstream;
            }
            else
            {
                newStart = Math.Max(0, record.ChromStart - upstream);
                newEnd = record.ChromEnd + downstream;
            }

            yield return record with { ChromStart = newStart, ChromEnd = newEnd };
        }
    }

    #endregion

    #region Block Operations (BED12)

    /// <summary>
    /// Expands BED12 blocks to individual exon records.
    /// </summary>
    public static IEnumerable<BedRecord> ExpandBlocks(BedRecord record)
    {
        if (!record.HasBlocks || record.BlockCount == 0)
        {
            yield return record;
            yield break;
        }

        for (int i = 0; i < record.BlockCount && i < record.BlockSizes!.Length && i < record.BlockStarts!.Length; i++)
        {
            int blockStart = record.ChromStart + record.BlockStarts[i];
            int blockEnd = blockStart + record.BlockSizes[i];

            yield return new BedRecord(
                record.Chrom,
                blockStart,
                blockEnd,
                record.Name != null ? $"{record.Name}_block{i + 1}" : null,
                record.Score,
                record.Strand);
        }
    }

    /// <summary>
    /// Calculates total block length for BED12 records.
    /// </summary>
    public static int GetTotalBlockLength(BedRecord record)
    {
        if (!record.HasBlocks || record.BlockSizes == null)
            return record.Length;

        return record.BlockSizes.Sum();
    }

    /// <summary>
    /// Gets intron regions from BED12 blocks.
    /// </summary>
    public static IEnumerable<BedRecord> GetIntrons(BedRecord record)
    {
        if (!record.HasBlocks || record.BlockCount < 2)
            yield break;

        for (int i = 0; i < record.BlockCount - 1 && i < record.BlockSizes!.Length - 1; i++)
        {
            int intronStart = record.ChromStart + record.BlockStarts![i] + record.BlockSizes[i];
            int intronEnd = record.ChromStart + record.BlockStarts[i + 1];

            if (intronEnd > intronStart)
            {
                yield return new BedRecord(
                    record.Chrom,
                    intronStart,
                    intronEnd,
                    record.Name != null ? $"{record.Name}_intron{i + 1}" : null,
                    record.Score,
                    record.Strand);
            }
        }
    }

    #endregion

    #region Statistics

    /// <summary>
    /// Calculates statistics for BED records.
    /// </summary>
    public static BedStatistics CalculateStatistics(IEnumerable<BedRecord> records)
    {
        var recordsList = records.ToList();

        if (recordsList.Count == 0)
        {
            return new BedStatistics(0, 0, 0, 0, 0, new Dictionary<string, int>());
        }

        var lengths = recordsList.Select(r => r.Length).ToList();
        var chromCounts = recordsList
            .GroupBy(r => r.Chrom, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(g => g.Key, g => g.Count());

        return new BedStatistics(
            recordsList.Count,
            lengths.Sum(l => (long)l),
            lengths.Min(),
            lengths.Max(),
            lengths.Average(),
            chromCounts);
    }

    /// <summary>BED file statistics</summary>
    public readonly record struct BedStatistics(
        int RecordCount,
        long TotalBases,
        int MinLength,
        int MaxLength,
        double AverageLength,
        IReadOnlyDictionary<string, int> ChromosomeCounts);

    #endregion

    #region Writing

    /// <summary>
    /// Writes BED records to a file.
    /// </summary>
    public static void WriteToFile(string filePath, IEnumerable<BedRecord> records, BedFormat format = BedFormat.BED6)
    {
        using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
        WriteToStream(writer, records, format);
    }

    /// <summary>
    /// Writes BED records to a TextWriter.
    /// </summary>
    public static void WriteToStream(TextWriter writer, IEnumerable<BedRecord> records, BedFormat format = BedFormat.BED6)
    {
        foreach (var record in records)
        {
            writer.WriteLine(FormatRecord(record, format));
        }
    }

    private static string FormatRecord(BedRecord record, BedFormat format)
    {
        var parts = new List<string>
        {
            record.Chrom,
            record.ChromStart.ToString(CultureInfo.InvariantCulture),
            record.ChromEnd.ToString(CultureInfo.InvariantCulture)
        };

        if (format >= BedFormat.BED4)
            parts.Add(record.Name ?? ".");

        if (format >= BedFormat.BED5)
            parts.Add(record.Score?.ToString(CultureInfo.InvariantCulture) ?? "0");

        if (format >= BedFormat.BED6)
            parts.Add(record.Strand?.ToString() ?? ".");

        if (format == BedFormat.BED12)
        {
            parts.Add(record.ThickStart?.ToString(CultureInfo.InvariantCulture) ?? record.ChromStart.ToString(CultureInfo.InvariantCulture));
            parts.Add(record.ThickEnd?.ToString(CultureInfo.InvariantCulture) ?? record.ChromEnd.ToString(CultureInfo.InvariantCulture));
            parts.Add(record.ItemRgb ?? "0,0,0");
            parts.Add(record.BlockCount?.ToString(CultureInfo.InvariantCulture) ?? "1");
            parts.Add(record.BlockSizes != null ? string.Join(",", record.BlockSizes) : record.Length.ToString(CultureInfo.InvariantCulture));
            parts.Add(record.BlockStarts != null ? string.Join(",", record.BlockStarts) : "0");
        }

        return string.Join('\t', parts);
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Sorts records by genomic position.
    /// </summary>
    public static IEnumerable<BedRecord> Sort(IEnumerable<BedRecord> records)
    {
        return records
            .OrderBy(r => r.Chrom, StringComparer.OrdinalIgnoreCase)
            .ThenBy(r => r.ChromStart)
            .ThenBy(r => r.ChromEnd);
    }

    /// <summary>
    /// Gets coverage depth at each position.
    /// </summary>
    public static IEnumerable<(int Position, int Depth)> CalculateCoverage(
        IEnumerable<BedRecord> records,
        string chrom,
        int start,
        int end)
    {
        var chromRecords = records
            .Where(r => r.Chrom.Equals(chrom, StringComparison.OrdinalIgnoreCase))
            .Where(r => r.ChromStart < end && r.ChromEnd > start)
            .ToList();

        var events = new List<(int Position, int Delta)>();
        foreach (var record in chromRecords)
        {
            events.Add((Math.Max(record.ChromStart, start), 1));
            events.Add((Math.Min(record.ChromEnd, end), -1));
        }

        events = events.OrderBy(e => e.Position).ThenBy(e => e.Delta).ToList();

        int depth = 0;
        int lastPos = start;

        foreach (var (position, delta) in events)
        {
            if (position > lastPos)
            {
                yield return (lastPos, depth);
                lastPos = position;
            }
            depth += delta;
        }

        if (lastPos < end)
            yield return (lastPos, depth);
    }

    /// <summary>
    /// Extracts sequence for a BED record from reference.
    /// </summary>
    public static string ExtractSequence(BedRecord record, string referenceSequence)
    {
        if (string.IsNullOrEmpty(referenceSequence))
            return "";

        int start = Math.Max(0, record.ChromStart);
        int end = Math.Min(referenceSequence.Length, record.ChromEnd);

        if (start >= end)
            return "";

        var sequence = referenceSequence[start..end];

        if (record.Strand == '-')
        {
            return new DnaSequence(sequence).ReverseComplement().Sequence;
        }

        return sequence;
    }

    #endregion
}
