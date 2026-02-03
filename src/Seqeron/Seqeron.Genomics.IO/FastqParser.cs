using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.IO;

/// <summary>
/// Parser for FASTQ format files containing sequences with quality scores.
/// Supports Phred+33 (Sanger/Illumina 1.8+) and Phred+64 (Illumina 1.3-1.7) encodings.
/// </summary>
public static class FastqParser
{
    #region Records

    /// <summary>Represents a FASTQ record with sequence and quality</summary>
    public readonly record struct FastqRecord(
        string Id,
        string Description,
        string Sequence,
        string QualityString,
        IReadOnlyList<int> QualityScores);

    /// <summary>Quality encoding format</summary>
    public enum QualityEncoding
    {
        /// <summary>Phred+33 (Sanger, Illumina 1.8+)</summary>
        Phred33,
        /// <summary>Phred+64 (Illumina 1.3-1.7)</summary>
        Phred64,
        /// <summary>Auto-detect from quality string</summary>
        Auto
    }

    /// <summary>Statistics for a FASTQ file</summary>
    public readonly record struct FastqStatistics(
        int TotalReads,
        long TotalBases,
        double MeanReadLength,
        double MeanQuality,
        int MinReadLength,
        int MaxReadLength,
        double Q20Percentage,
        double Q30Percentage,
        double GcContent);

    #endregion

    #region Parsing Methods

    /// <summary>
    /// Parses FASTQ records from a file.
    /// </summary>
    public static IEnumerable<FastqRecord> ParseFile(string filePath, QualityEncoding encoding = QualityEncoding.Auto)
    {
        if (string.IsNullOrEmpty(filePath) || !File.Exists(filePath))
            yield break;

        using var reader = new StreamReader(filePath);
        foreach (var record in Parse(reader, encoding))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses FASTQ records from text content.
    /// </summary>
    public static IEnumerable<FastqRecord> Parse(string content, QualityEncoding encoding = QualityEncoding.Auto)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        using var reader = new StringReader(content);
        foreach (var record in Parse(reader, encoding))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses FASTQ records from a TextReader.
    /// </summary>
    public static IEnumerable<FastqRecord> Parse(TextReader reader, QualityEncoding encoding = QualityEncoding.Auto)
    {
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            // Skip empty lines
            if (string.IsNullOrWhiteSpace(line))
                continue;

            // Header line (starts with @)
            if (!line.StartsWith('@'))
                continue;

            var header = line[1..];
            var (id, description) = ParseHeader(header);

            // Sequence line(s)
            var sequenceBuilder = new StringBuilder();
            while ((line = reader.ReadLine()) != null && !line.StartsWith('+'))
            {
                sequenceBuilder.Append(line.Trim());
            }
            var sequence = sequenceBuilder.ToString();

            // Quality header line (starts with +)
            // Already consumed by the while loop above

            // Quality line(s) - must be same length as sequence
            var qualityBuilder = new StringBuilder();
            while (qualityBuilder.Length < sequence.Length && (line = reader.ReadLine()) != null)
            {
                qualityBuilder.Append(line.Trim());
            }
            var qualityString = qualityBuilder.ToString();

            // Decode quality scores
            var actualEncoding = encoding == QualityEncoding.Auto
                ? DetectEncoding(qualityString)
                : encoding;
            var qualityScores = DecodeQualityScores(qualityString, actualEncoding);

            yield return new FastqRecord(id, description, sequence, qualityString, qualityScores);
        }
    }

    private static (string Id, string Description) ParseHeader(string header)
    {
        var spaceIndex = header.IndexOf(' ');
        if (spaceIndex > 0)
        {
            return (header[..spaceIndex], header[(spaceIndex + 1)..]);
        }
        return (header, "");
    }

    #endregion

    #region Quality Encoding

    /// <summary>
    /// Detects quality encoding from quality string.
    /// </summary>
    public static QualityEncoding DetectEncoding(string qualityString)
    {
        if (string.IsNullOrEmpty(qualityString))
            return QualityEncoding.Phred33;

        foreach (char c in qualityString)
        {
            // Characters below '@' (64) indicate Phred+33
            if (c < '@')
                return QualityEncoding.Phred33;
            // Characters above 'I' (73) indicate Phred+64
            if (c > 'I')
                return QualityEncoding.Phred64;
        }

        // Default to Phred+33 (most common modern format)
        return QualityEncoding.Phred33;
    }

    /// <summary>
    /// Decodes quality string to Phred scores.
    /// </summary>
    public static IReadOnlyList<int> DecodeQualityScores(string qualityString, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(qualityString))
            return Array.Empty<int>();

        int offset = encoding == QualityEncoding.Phred64 ? 64 : 33;
        var scores = new int[qualityString.Length];

        for (int i = 0; i < qualityString.Length; i++)
        {
            scores[i] = Math.Max(0, qualityString[i] - offset);
        }

        return scores;
    }

    /// <summary>
    /// Encodes Phred scores to quality string.
    /// </summary>
    public static string EncodeQualityScores(IEnumerable<int> scores, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        int offset = encoding == QualityEncoding.Phred64 ? 64 : 33;
        var sb = new StringBuilder();

        foreach (var score in scores)
        {
            var clampedScore = Math.Clamp(score, 0, 41);
            sb.Append((char)(clampedScore + offset));
        }

        return sb.ToString();
    }

    /// <summary>
    /// Converts Phred score to error probability.
    /// </summary>
    public static double PhredToErrorProbability(int phredScore)
    {
        return Math.Pow(10, -phredScore / 10.0);
    }

    /// <summary>
    /// Converts error probability to Phred score.
    /// </summary>
    public static int ErrorProbabilityToPhred(double errorProbability)
    {
        if (errorProbability <= 0)
            return 40; // Max quality
        return (int)Math.Round(-10 * Math.Log10(errorProbability));
    }

    #endregion

    #region Filtering and Trimming

    /// <summary>
    /// Filters records by minimum average quality.
    /// </summary>
    public static IEnumerable<FastqRecord> FilterByQuality(
        IEnumerable<FastqRecord> records,
        double minAverageQuality)
    {
        foreach (var record in records)
        {
            if (record.QualityScores.Count > 0)
            {
                double avgQuality = record.QualityScores.Average();
                if (avgQuality >= minAverageQuality)
                    yield return record;
            }
        }
    }

    /// <summary>
    /// Filters records by minimum length.
    /// </summary>
    public static IEnumerable<FastqRecord> FilterByLength(
        IEnumerable<FastqRecord> records,
        int minLength,
        int? maxLength = null)
    {
        foreach (var record in records)
        {
            if (record.Sequence.Length >= minLength &&
                (!maxLength.HasValue || record.Sequence.Length <= maxLength.Value))
            {
                yield return record;
            }
        }
    }

    /// <summary>
    /// Trims low quality bases from ends.
    /// </summary>
    public static FastqRecord TrimByQuality(FastqRecord record, int minQuality = 20)
    {
        if (record.QualityScores.Count == 0)
            return record;

        // Find trim positions
        int start = 0;
        int end = record.Sequence.Length;

        // Trim from start
        while (start < end && record.QualityScores[start] < minQuality)
            start++;

        // Trim from end
        while (end > start && record.QualityScores[end - 1] < minQuality)
            end--;

        if (start >= end)
        {
            // Entire sequence trimmed
            return new FastqRecord(record.Id, record.Description, "", "", Array.Empty<int>());
        }

        var newSequence = record.Sequence[start..end];
        var newQualityString = record.QualityString[start..end];
        var newScores = record.QualityScores.Skip(start).Take(end - start).ToList();

        return new FastqRecord(record.Id, record.Description, newSequence, newQualityString, newScores);
    }

    /// <summary>
    /// Trims adapters from sequences.
    /// </summary>
    public static FastqRecord TrimAdapter(FastqRecord record, string adapter, int minOverlap = 5)
    {
        if (string.IsNullOrEmpty(adapter) || adapter.Length < minOverlap)
            return record;

        adapter = adapter.ToUpperInvariant();
        var sequence = record.Sequence.ToUpperInvariant();

        // Search for adapter at the end
        for (int overlapLen = adapter.Length; overlapLen >= minOverlap; overlapLen--)
        {
            var adapterPrefix = adapter[..overlapLen];
            var searchStart = record.Sequence.Length - overlapLen;

            if (searchStart >= 0 && sequence.EndsWith(adapterPrefix, StringComparison.Ordinal))
            {
                return new FastqRecord(
                    record.Id,
                    record.Description,
                    record.Sequence[..searchStart],
                    record.QualityString[..searchStart],
                    record.QualityScores.Take(searchStart).ToList());
            }
        }

        // Search for full adapter within sequence
        int adapterPos = sequence.IndexOf(adapter, StringComparison.Ordinal);
        if (adapterPos > 0)
        {
            return new FastqRecord(
                record.Id,
                record.Description,
                record.Sequence[..adapterPos],
                record.QualityString[..adapterPos],
                record.QualityScores.Take(adapterPos).ToList());
        }

        return record;
    }

    #endregion

    #region Statistics

    /// <summary>
    /// Calculates statistics for FASTQ records.
    /// </summary>
    public static FastqStatistics CalculateStatistics(IEnumerable<FastqRecord> records)
    {
        int totalReads = 0;
        long totalBases = 0;
        long totalQuality = 0;
        int minLength = int.MaxValue;
        int maxLength = 0;
        long q20Bases = 0;
        long q30Bases = 0;
        int gcCount = 0;

        foreach (var record in records)
        {
            totalReads++;
            totalBases += record.Sequence.Length;

            if (record.Sequence.Length < minLength)
                minLength = record.Sequence.Length;
            if (record.Sequence.Length > maxLength)
                maxLength = record.Sequence.Length;

            foreach (var score in record.QualityScores)
            {
                totalQuality += score;
                if (score >= 20) q20Bases++;
                if (score >= 30) q30Bases++;
            }

            foreach (var c in record.Sequence.ToUpperInvariant())
            {
                if (c == 'G' || c == 'C')
                    gcCount++;
            }
        }

        if (totalReads == 0)
        {
            return new FastqStatistics(0, 0, 0, 0, 0, 0, 0, 0, 0);
        }

        return new FastqStatistics(
            totalReads,
            totalBases,
            (double)totalBases / totalReads,
            totalBases > 0 ? (double)totalQuality / totalBases : 0,
            minLength == int.MaxValue ? 0 : minLength,
            maxLength,
            totalBases > 0 ? 100.0 * q20Bases / totalBases : 0,
            totalBases > 0 ? 100.0 * q30Bases / totalBases : 0,
            totalBases > 0 ? (double)gcCount / totalBases : 0);
    }

    /// <summary>
    /// Calculates per-position quality statistics.
    /// </summary>
    public static IReadOnlyList<(int Position, double MeanQuality, double StdDev)>
        CalculatePositionQuality(IEnumerable<FastqRecord> records)
    {
        var qualityByPosition = new List<List<int>>();

        foreach (var record in records)
        {
            for (int i = 0; i < record.QualityScores.Count; i++)
            {
                while (qualityByPosition.Count <= i)
                    qualityByPosition.Add(new List<int>());
                qualityByPosition[i].Add(record.QualityScores[i]);
            }
        }

        var result = new List<(int, double, double)>();
        for (int i = 0; i < qualityByPosition.Count; i++)
        {
            var scores = qualityByPosition[i];
            if (scores.Count == 0) continue;

            double mean = scores.Average();
            double variance = scores.Sum(s => (s - mean) * (s - mean)) / scores.Count;
            double stdDev = Math.Sqrt(variance);

            result.Add((i + 1, mean, stdDev));
        }

        return result;
    }

    #endregion

    #region Writing

    /// <summary>
    /// Writes FASTQ records to a file.
    /// </summary>
    public static void WriteToFile(string filePath, IEnumerable<FastqRecord> records)
    {
        using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
        WriteToStream(writer, records);
    }

    /// <summary>
    /// Writes FASTQ records to a TextWriter.
    /// </summary>
    public static void WriteToStream(TextWriter writer, IEnumerable<FastqRecord> records)
    {
        foreach (var record in records)
        {
            // Header
            writer.Write('@');
            writer.Write(record.Id);
            if (!string.IsNullOrEmpty(record.Description))
            {
                writer.Write(' ');
                writer.Write(record.Description);
            }
            writer.WriteLine();

            // Sequence
            writer.WriteLine(record.Sequence);

            // Quality header
            writer.WriteLine('+');

            // Quality string
            writer.WriteLine(record.QualityString);
        }
    }

    /// <summary>
    /// Converts a FastqRecord to string format.
    /// </summary>
    public static string ToFastqString(FastqRecord record)
    {
        var sb = new StringBuilder();
        sb.Append('@').Append(record.Id);
        if (!string.IsNullOrEmpty(record.Description))
            sb.Append(' ').Append(record.Description);
        sb.AppendLine();
        sb.AppendLine(record.Sequence);
        sb.AppendLine("+");
        sb.AppendLine(record.QualityString);
        return sb.ToString();
    }

    #endregion

    #region Paired-End Support

    /// <summary>
    /// Interleaves paired-end reads.
    /// </summary>
    public static IEnumerable<FastqRecord> InterleavePairedReads(
        IEnumerable<FastqRecord> read1,
        IEnumerable<FastqRecord> read2)
    {
        using var enum1 = read1.GetEnumerator();
        using var enum2 = read2.GetEnumerator();

        while (enum1.MoveNext() && enum2.MoveNext())
        {
            yield return enum1.Current;
            yield return enum2.Current;
        }
    }

    /// <summary>
    /// Splits interleaved paired-end reads.
    /// </summary>
    public static (IReadOnlyList<FastqRecord> Read1, IReadOnlyList<FastqRecord> Read2)
        SplitInterleavedReads(IEnumerable<FastqRecord> interleaved)
    {
        var read1 = new List<FastqRecord>();
        var read2 = new List<FastqRecord>();
        bool isRead1 = true;

        foreach (var record in interleaved)
        {
            if (isRead1)
                read1.Add(record);
            else
                read2.Add(record);
            isRead1 = !isRead1;
        }

        return (read1, read2);
    }

    #endregion
}

