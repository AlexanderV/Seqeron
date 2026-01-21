using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree.Genomics;

/// <summary>
/// Provides quality score analysis for sequencing data (FASTQ format).
/// Supports Phred+33 (Sanger/Illumina 1.8+) and Phred+64 (Illumina 1.3-1.7) encodings.
/// </summary>
public static class QualityScoreAnalyzer
{
    /// <summary>
    /// Quality encoding format.
    /// </summary>
    public enum QualityEncoding
    {
        /// <summary>Phred+33 (Sanger, Illumina 1.8+). ASCII 33-126, quality 0-93.</summary>
        Phred33,
        /// <summary>Phred+64 (Illumina 1.3-1.7). ASCII 64-126, quality 0-62.</summary>
        Phred64,
        /// <summary>Auto-detect encoding from data.</summary>
        Auto
    }

    /// <summary>
    /// FASTQ record containing sequence and quality data.
    /// </summary>
    public readonly record struct FastqRecord(
        string Id,
        string Sequence,
        string QualityString,
        string? Description = null);

    /// <summary>
    /// Quality statistics for a sequence or dataset.
    /// </summary>
    public readonly record struct QualityStatistics(
        double MeanQuality,
        double MedianQuality,
        int MinQuality,
        int MaxQuality,
        double StandardDeviation,
        int TotalBases,
        int BasesAboveQ20,
        int BasesAboveQ30,
        double PercentAboveQ20,
        double PercentAboveQ30,
        IReadOnlyList<double> PerPositionMeanQuality);

    /// <summary>
    /// Trimming result with statistics.
    /// </summary>
    public readonly record struct TrimResult(
        string Sequence,
        string QualityString,
        int TrimmedFromStart,
        int TrimmedFromEnd,
        int OriginalLength,
        int FinalLength);

    /// <summary>
    /// Converts a quality character to Phred score.
    /// </summary>
    public static int CharToPhred(char qualityChar, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        int offset = encoding == QualityEncoding.Phred64 ? 64 : 33;
        return qualityChar - offset;
    }

    /// <summary>
    /// Converts a Phred score to quality character.
    /// </summary>
    public static char PhredToChar(int phredScore, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        int offset = encoding == QualityEncoding.Phred64 ? 64 : 33;
        return (char)(phredScore + offset);
    }

    /// <summary>
    /// Converts a quality string to array of Phred scores.
    /// </summary>
    public static int[] QualityStringToPhred(string qualityString, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(qualityString))
            return Array.Empty<int>();

        var actualEncoding = encoding == QualityEncoding.Auto
            ? DetectEncoding(qualityString)
            : encoding;

        return qualityString.Select(c => CharToPhred(c, actualEncoding)).ToArray();
    }

    /// <summary>
    /// Converts Phred scores to quality string.
    /// </summary>
    public static string PhredToQualityString(IEnumerable<int> phredScores, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        return new string(phredScores.Select(p => PhredToChar(p, encoding)).ToArray());
    }

    /// <summary>
    /// Detects the quality encoding from a quality string.
    /// </summary>
    public static QualityEncoding DetectEncoding(string qualityString)
    {
        if (string.IsNullOrEmpty(qualityString))
            return QualityEncoding.Phred33;

        int minChar = qualityString.Min();
        int maxChar = qualityString.Max();

        // Phred+33: ASCII 33-73 typically (! to I)
        // Phred+64: ASCII 64-104 typically (@ to h)
        if (minChar < 59) // Characters below ';' strongly suggest Phred+33
            return QualityEncoding.Phred33;
        if (minChar >= 64 && maxChar > 74) // High range suggests Phred+64
            return QualityEncoding.Phred64;

        return QualityEncoding.Phred33; // Default to modern format
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
        if (errorProbability <= 0) return 60; // Max practical quality
        if (errorProbability >= 1) return 0;
        return (int)Math.Round(-10 * Math.Log10(errorProbability));
    }

    /// <summary>
    /// Calculates comprehensive quality statistics for a quality string.
    /// </summary>
    public static QualityStatistics CalculateStatistics(
        string qualityString,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(qualityString))
        {
            return new QualityStatistics(
                MeanQuality: 0,
                MedianQuality: 0,
                MinQuality: 0,
                MaxQuality: 0,
                StandardDeviation: 0,
                TotalBases: 0,
                BasesAboveQ20: 0,
                BasesAboveQ30: 0,
                PercentAboveQ20: 0,
                PercentAboveQ30: 0,
                PerPositionMeanQuality: Array.Empty<double>());
        }

        var phredScores = QualityStringToPhred(qualityString, encoding);
        return CalculateStatisticsFromPhred(phredScores);
    }

    /// <summary>
    /// Calculates quality statistics for multiple reads.
    /// </summary>
    public static QualityStatistics CalculateStatistics(
        IEnumerable<string> qualityStrings,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        var allScores = new List<int>();
        var positionScores = new Dictionary<int, List<int>>();

        foreach (var qualityString in qualityStrings)
        {
            var phred = QualityStringToPhred(qualityString, encoding);
            allScores.AddRange(phred);

            for (int i = 0; i < phred.Length; i++)
            {
                if (!positionScores.ContainsKey(i))
                    positionScores[i] = new List<int>();
                positionScores[i].Add(phred[i]);
            }
        }

        if (allScores.Count == 0)
        {
            return new QualityStatistics(
                MeanQuality: 0, MedianQuality: 0, MinQuality: 0, MaxQuality: 0,
                StandardDeviation: 0, TotalBases: 0, BasesAboveQ20: 0, BasesAboveQ30: 0,
                PercentAboveQ20: 0, PercentAboveQ30: 0,
                PerPositionMeanQuality: Array.Empty<double>());
        }

        var perPositionMean = positionScores
            .OrderBy(kvp => kvp.Key)
            .Select(kvp => kvp.Value.Average())
            .ToList();

        return CalculateStatisticsFromPhred(allScores.ToArray(), perPositionMean);
    }

    private static QualityStatistics CalculateStatisticsFromPhred(
        int[] phredScores,
        IReadOnlyList<double>? perPositionMean = null)
    {
        if (phredScores.Length == 0)
        {
            return new QualityStatistics(
                MeanQuality: 0, MedianQuality: 0, MinQuality: 0, MaxQuality: 0,
                StandardDeviation: 0, TotalBases: 0, BasesAboveQ20: 0, BasesAboveQ30: 0,
                PercentAboveQ20: 0, PercentAboveQ30: 0,
                PerPositionMeanQuality: Array.Empty<double>());
        }

        double mean = phredScores.Average();
        var sorted = phredScores.OrderBy(x => x).ToArray();
        double median = sorted.Length % 2 == 0
            ? (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]) / 2.0
            : sorted[sorted.Length / 2];

        double variance = phredScores.Select(x => Math.Pow(x - mean, 2)).Average();
        double stdDev = Math.Sqrt(variance);

        int aboveQ20 = phredScores.Count(q => q >= 20);
        int aboveQ30 = phredScores.Count(q => q >= 30);

        var positionMean = perPositionMean ?? phredScores.Select(p => (double)p).ToList();

        return new QualityStatistics(
            MeanQuality: mean,
            MedianQuality: median,
            MinQuality: phredScores.Min(),
            MaxQuality: phredScores.Max(),
            StandardDeviation: stdDev,
            TotalBases: phredScores.Length,
            BasesAboveQ20: aboveQ20,
            BasesAboveQ30: aboveQ30,
            PercentAboveQ20: 100.0 * aboveQ20 / phredScores.Length,
            PercentAboveQ30: 100.0 * aboveQ30 / phredScores.Length,
            PerPositionMeanQuality: positionMean);
    }

    /// <summary>
    /// Trims low-quality bases from both ends of a read.
    /// </summary>
    public static TrimResult QualityTrim(
        string sequence,
        string qualityString,
        int minQuality = 20,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(qualityString))
        {
            return new TrimResult(
                Sequence: "",
                QualityString: "",
                TrimmedFromStart: 0,
                TrimmedFromEnd: 0,
                OriginalLength: sequence?.Length ?? 0,
                FinalLength: 0);
        }

        var phred = QualityStringToPhred(qualityString, encoding);
        int len = Math.Min(sequence.Length, phred.Length);

        // Find first position >= minQuality
        int start = 0;
        while (start < len && phred[start] < minQuality)
            start++;

        // Find last position >= minQuality
        int end = len - 1;
        while (end >= start && phred[end] < minQuality)
            end--;

        if (start > end)
        {
            return new TrimResult(
                Sequence: "",
                QualityString: "",
                TrimmedFromStart: len,
                TrimmedFromEnd: 0,
                OriginalLength: len,
                FinalLength: 0);
        }

        int trimmedLen = end - start + 1;
        return new TrimResult(
            Sequence: sequence.Substring(start, trimmedLen),
            QualityString: qualityString.Substring(start, trimmedLen),
            TrimmedFromStart: start,
            TrimmedFromEnd: len - end - 1,
            OriginalLength: len,
            FinalLength: trimmedLen);
    }

    /// <summary>
    /// Trims using sliding window average quality.
    /// </summary>
    public static TrimResult SlidingWindowTrim(
        string sequence,
        string qualityString,
        int windowSize = 4,
        int minAverageQuality = 20,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(qualityString))
        {
            return new TrimResult("", "", 0, 0, sequence?.Length ?? 0, 0);
        }

        var phred = QualityStringToPhred(qualityString, encoding);
        int len = Math.Min(sequence.Length, phred.Length);

        if (len < windowSize)
        {
            double avg = phred.Take(len).Average();
            if (avg >= minAverageQuality)
                return new TrimResult(sequence, qualityString, 0, 0, len, len);
            return new TrimResult("", "", len, 0, len, 0);
        }

        // Find trim point from end (where window average drops below threshold)
        int cutoff = len;
        for (int i = len - windowSize; i >= 0; i--)
        {
            double windowAvg = 0;
            for (int j = 0; j < windowSize; j++)
                windowAvg += phred[i + j];
            windowAvg /= windowSize;

            if (windowAvg >= minAverageQuality)
            {
                cutoff = i + windowSize;
                break;
            }
            cutoff = i;
        }

        if (cutoff <= 0)
            return new TrimResult("", "", len, 0, len, 0);

        return new TrimResult(
            Sequence: sequence.Substring(0, cutoff),
            QualityString: qualityString.Substring(0, cutoff),
            TrimmedFromStart: 0,
            TrimmedFromEnd: len - cutoff,
            OriginalLength: len,
            FinalLength: cutoff);
    }

    /// <summary>
    /// Filters reads based on quality criteria.
    /// </summary>
    public static IEnumerable<FastqRecord> FilterReads(
        IEnumerable<FastqRecord> reads,
        int minLength = 0,
        int maxLength = int.MaxValue,
        double minMeanQuality = 0,
        double maxExpectedErrors = double.MaxValue,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        foreach (var read in reads)
        {
            // Length filter
            if (read.Sequence.Length < minLength || read.Sequence.Length > maxLength)
                continue;

            var phred = QualityStringToPhred(read.QualityString, encoding);

            // Mean quality filter
            if (minMeanQuality > 0)
            {
                double mean = phred.Average();
                if (mean < minMeanQuality)
                    continue;
            }

            // Expected errors filter
            if (maxExpectedErrors < double.MaxValue)
            {
                double expectedErrors = phred.Sum(p => PhredToErrorProbability(p));
                if (expectedErrors > maxExpectedErrors)
                    continue;
            }

            yield return read;
        }
    }

    /// <summary>
    /// Calculates expected number of errors in a read.
    /// </summary>
    public static double CalculateExpectedErrors(
        string qualityString,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        var phred = QualityStringToPhred(qualityString, encoding);
        return phred.Sum(p => PhredToErrorProbability(p));
    }

    /// <summary>
    /// Masks low-quality bases with 'N'.
    /// </summary>
    public static string MaskLowQualityBases(
        string sequence,
        string qualityString,
        int minQuality = 20,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(qualityString))
            return sequence ?? "";

        var phred = QualityStringToPhred(qualityString, encoding);
        var sb = new StringBuilder(sequence.Length);

        for (int i = 0; i < sequence.Length && i < phred.Length; i++)
        {
            sb.Append(phred[i] >= minQuality ? sequence[i] : 'N');
        }

        return sb.ToString();
    }

    /// <summary>
    /// Parses FASTQ format text into records.
    /// </summary>
    public static IEnumerable<FastqRecord> ParseFastq(IEnumerable<string> lines)
    {
        using var enumerator = lines.GetEnumerator();

        while (enumerator.MoveNext())
        {
            string headerLine = enumerator.Current;
            if (string.IsNullOrWhiteSpace(headerLine) || !headerLine.StartsWith("@"))
                continue;

            // Parse header
            string header = headerLine.Substring(1);
            string id = header.Split(' ', '\t')[0];
            string? description = header.Contains(' ') ? header.Substring(header.IndexOf(' ') + 1) : null;

            // Sequence line
            if (!enumerator.MoveNext()) yield break;
            string sequence = enumerator.Current;

            // Plus line
            if (!enumerator.MoveNext()) yield break;
            // Skip the '+' line

            // Quality line
            if (!enumerator.MoveNext()) yield break;
            string quality = enumerator.Current;

            yield return new FastqRecord(id, sequence, quality, description);
        }
    }

    /// <summary>
    /// Formats a FastqRecord as FASTQ text lines.
    /// </summary>
    public static IEnumerable<string> ToFastq(FastqRecord record)
    {
        string header = string.IsNullOrEmpty(record.Description)
            ? $"@{record.Id}"
            : $"@{record.Id} {record.Description}";

        yield return header;
        yield return record.Sequence;
        yield return "+";
        yield return record.QualityString;
    }

    /// <summary>
    /// Formats multiple records as FASTQ text.
    /// </summary>
    public static IEnumerable<string> ToFastq(IEnumerable<FastqRecord> records)
    {
        foreach (var record in records)
        {
            foreach (var line in ToFastq(record))
                yield return line;
        }
    }

    /// <summary>
    /// Calculates per-base quality distribution.
    /// </summary>
    public static IReadOnlyDictionary<int, int> GetQualityDistribution(
        string qualityString,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        var phred = QualityStringToPhred(qualityString, encoding);
        return phred.GroupBy(q => q).ToDictionary(g => g.Key, g => g.Count());
    }

    /// <summary>
    /// Identifies low-quality regions in a read.
    /// </summary>
    public static IEnumerable<(int start, int end, double meanQuality)> FindLowQualityRegions(
        string qualityString,
        int windowSize = 10,
        int maxQuality = 15,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        var phred = QualityStringToPhred(qualityString, encoding);

        if (phred.Length < windowSize)
            yield break;

        int? regionStart = null;
        double regionSum = 0;

        for (int i = 0; i <= phred.Length - windowSize; i++)
        {
            double windowMean;
            if (i == 0)
            {
                windowMean = phred.Take(windowSize).Average();
            }
            else
            {
                // Sliding window update
                windowMean = 0;
                for (int j = 0; j < windowSize; j++)
                    windowMean += phred[i + j];
                windowMean /= windowSize;
            }

            if (windowMean <= maxQuality)
            {
                if (regionStart == null)
                {
                    regionStart = i;
                    regionSum = windowMean;
                }
                else
                {
                    regionSum += windowMean;
                }
            }
            else if (regionStart != null)
            {
                int regionLen = i - regionStart.Value;
                yield return (regionStart.Value, i + windowSize - 1, regionSum / regionLen);
                regionStart = null;
                regionSum = 0;
            }
        }

        if (regionStart != null)
        {
            int regionLen = phred.Length - windowSize - regionStart.Value + 1;
            yield return (regionStart.Value, phred.Length - 1, regionSum / Math.Max(1, regionLen));
        }
    }

    /// <summary>
    /// Calculates consensus quality from multiple aligned quality strings.
    /// </summary>
    public static string CalculateConsensusQuality(
        IReadOnlyList<string> qualityStrings,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (qualityStrings.Count == 0)
            return "";

        int maxLen = qualityStrings.Max(q => q.Length);
        var consensusPhred = new int[maxLen];

        for (int pos = 0; pos < maxLen; pos++)
        {
            var scoresAtPos = new List<int>();
            foreach (var qs in qualityStrings)
            {
                if (pos < qs.Length)
                {
                    scoresAtPos.Add(CharToPhred(qs[pos], encoding));
                }
            }

            if (scoresAtPos.Count > 0)
            {
                // Use maximum quality (most confident base)
                consensusPhred[pos] = scoresAtPos.Max();
            }
        }

        return PhredToQualityString(consensusPhred, encoding);
    }
}
