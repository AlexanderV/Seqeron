using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.IO;

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
    /// Confidence of a file-level <see cref="DetectEncoding(IEnumerable{string})"/> result.
    /// </summary>
    public enum EncodingConfidence
    {
        /// <summary>
        /// Proven from the data: at least one quality character is below ASCII 64, which is outside the
        /// Phred+64 range (ASCII 64-126), so the encoding can only be Phred+33 (Cock et al., 2010).
        /// </summary>
        Definitive,

        /// <summary>
        /// Inferred heuristically: no character below ASCII 64, and at least one character above ASCII 74
        /// (Phred 41, the Illumina 1.8+ Phred+33 ceiling), which is implausible under Phred+33, so the
        /// encoding is taken to be Phred+64. Phred+64 can never be <i>proven</i>, only inferred, because
        /// every Phred+64 quality string is also a syntactically valid Phred+33 string.
        /// </summary>
        Inferred,

        /// <summary>
        /// Undeterminable: every quality character lies in the Phred+33/Phred+64 overlap range
        /// (ASCII 64-74). The encoding is genuinely ambiguous and the result defaults to Phred+33.
        /// </summary>
        Ambiguous
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
    /// Result of file-level quality-encoding detection: the chosen <see cref="QualityEncoding"/>, how
    /// firmly it is established, and the ASCII span and character count the decision was based on.
    /// </summary>
    public readonly record struct EncodingDetectionResult(
        QualityEncoding Encoding,
        EncodingConfidence Confidence,
        int MinAscii,
        int MaxAscii,
        long CharactersExamined);

    // FASTQ encoding offsets and valid Phred score ranges per Cock et al. (2010),
    // Nucleic Acids Research 38(6):1767-1771, https://doi.org/10.1093/nar/gkp1137
    // Sanger/Phred+33: ASCII 33-126 -> Phred 0-93 (offset 33).
    // Illumina 1.3+/Phred+64: ASCII 64-126 -> Phred 0-62 (offset 64).
    private const int Phred33Offset = 33;
    private const int Phred64Offset = 64;
    private const int Phred33MaxScore = 93;
    private const int Phred64MaxScore = 62;
    private const int PhredMinScore = 0;

    // Highest quality character plausible under Phred+33: ASCII 74 ('J') = Phred 41, the Illumina 1.8+
    // ceiling (Cock et al., 2010). With no character below ASCII 64, a character ABOVE this indicates
    // Phred+64. Mirrors the single-string DetectEncoding threshold so per-read results are unchanged.
    private const int Phred33PlausibleMaxChar = 74;

    // Q20/Q30 quality thresholds (inclusive). Q30 = 99.9% base-call accuracy (1-in-1000 error)
    // and is the NGS run-quality benchmark; Q20 = 99% (1-in-100). Per Illumina, "Sequencing
    // Quality Scores" (Q = -10 log10 e), and Ewing & Green (1998) Genome Research 8(3):186-194.
    // A base whose Phred score is >= the threshold is counted (inclusive comparison).
    private const int Q20Threshold = 20;
    private const int Q30Threshold = 30;
    private const double PercentScale = 100.0;

    private static (int Offset, int MaxScore) GetEncodingParameters(QualityEncoding encoding)
    {
        return encoding == QualityEncoding.Phred64
            ? (Phred64Offset, Phred64MaxScore)
            : (Phred33Offset, Phred33MaxScore);
    }

    /// <summary>
    /// Parses a FASTQ quality string into an array of Phred quality scores using the
    /// specified encoding. Decodes each character as Q = ord(char) - offset, where the
    /// offset is 33 for Phred+33 and 64 for Phred+64.
    /// Per Cock et al. (2010), Phred+33 holds scores 0-93 and Phred+64 holds scores 0-62;
    /// a character that decodes outside the encoding's valid range is rejected as malformed.
    /// </summary>
    /// <param name="qualityString">FASTQ quality line (ASCII-encoded Phred scores).</param>
    /// <param name="encoding">Phred+33 (Sanger/Illumina 1.8+) or Phred+64 (Illumina 1.3-1.7). Auto is resolved by <see cref="DetectEncoding"/>.</param>
    /// <returns>Array of Phred scores; empty for an empty input.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="qualityString"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when a character decodes to a Phred score outside the encoding's valid range.</exception>
    public static int[] ParseQualityString(string qualityString, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (qualityString is null)
            throw new ArgumentNullException(nameof(qualityString));
        if (qualityString.Length == 0)
            return Array.Empty<int>();

        var actualEncoding = encoding == QualityEncoding.Auto
            ? DetectEncoding(qualityString)
            : encoding;
        var (offset, maxScore) = GetEncodingParameters(actualEncoding);

        var scores = new int[qualityString.Length];
        for (int i = 0; i < qualityString.Length; i++)
        {
            int score = qualityString[i] - offset;
            if (score < PhredMinScore || score > maxScore)
                throw new ArgumentOutOfRangeException(
                    nameof(qualityString),
                    $"Character '{qualityString[i]}' (ASCII {(int)qualityString[i]}) at index {i} decodes to Phred score {score}, " +
                    $"outside the valid range [{PhredMinScore}, {maxScore}] for {actualEncoding}.");
            scores[i] = score;
        }

        return scores;
    }

    /// <summary>
    /// Encodes an array of Phred quality scores into a FASTQ quality string using the
    /// specified encoding. Encodes each score as char = chr(Q + offset), where the offset
    /// is 33 for Phred+33 and 64 for Phred+64 (Cock et al., 2010).
    /// </summary>
    /// <param name="scores">Phred quality scores.</param>
    /// <param name="encoding">Phred+33 or Phred+64 (Auto is treated as Phred+33 for encoding, the modern default).</param>
    /// <returns>ASCII quality string; empty for an empty input.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="scores"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when a score is outside the encoding's valid range.</exception>
    public static string ToQualityString(IReadOnlyList<int> scores, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (scores is null)
            throw new ArgumentNullException(nameof(scores));
        if (scores.Count == 0)
            return string.Empty;

        var encodeAs = encoding == QualityEncoding.Phred64 ? QualityEncoding.Phred64 : QualityEncoding.Phred33;
        var (offset, maxScore) = GetEncodingParameters(encodeAs);

        var chars = new char[scores.Count];
        for (int i = 0; i < scores.Count; i++)
        {
            int score = scores[i];
            if (score < PhredMinScore || score > maxScore)
                throw new ArgumentOutOfRangeException(
                    nameof(scores),
                    $"Phred score {score} at index {i} is outside the valid range [{PhredMinScore}, {maxScore}] for {encodeAs}.");
            chars[i] = (char)(score + offset);
        }

        return new string(chars);
    }

    /// <summary>
    /// Converts a FASTQ quality string from one Phred encoding to another by re-offsetting
    /// each character while preserving the underlying Phred score. Because the Phred score is
    /// invariant across the Sanger (Phred+33) and Illumina 1.3+ (Phred+64) variants, conversion
    /// is a pure re-offset (Cock et al., 2010).
    /// </summary>
    /// <param name="qualityString">Source FASTQ quality string.</param>
    /// <param name="fromEncoding">Encoding of the input string.</param>
    /// <param name="toEncoding">Target encoding.</param>
    /// <returns>Quality string re-encoded under <paramref name="toEncoding"/>.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="qualityString"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when a decoded score is invalid for the source encoding, or not representable in the target encoding (e.g. a Phred+33 score &gt; 62 has no Phred+64 representation).</exception>
    public static string ConvertEncoding(string qualityString, QualityEncoding fromEncoding, QualityEncoding toEncoding)
    {
        if (qualityString is null)
            throw new ArgumentNullException(nameof(qualityString));

        var scores = ParseQualityString(qualityString, fromEncoding);
        return ToQualityString(scores, toEncoding);
    }

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
    /// Detects the quality encoding across a whole set of reads (the FASTQ file as a unit), which
    /// disambiguates the case the single-string <see cref="DetectEncoding(string)"/> cannot: a read
    /// confined to the Phred+33/Phred+64 overlap range (ASCII 64-74) is resolved by the rest of the
    /// file. The global minimum and maximum ASCII over every character of every quality string are
    /// scanned in one pass and the Cock et al. (2010) ranges applied:
    /// <list type="bullet">
    /// <item>any character &lt; ASCII 64 ⇒ Phred+33 (<see cref="EncodingConfidence.Definitive"/>; that
    /// character is outside the Phred+64 range);</item>
    /// <item>otherwise a character &gt; ASCII 74 ⇒ Phred+64 (<see cref="EncodingConfidence.Inferred"/>;
    /// above the Phred+33 ceiling with no low character);</item>
    /// <item>otherwise every character lies in the overlap range ⇒ Phred+33 default
    /// (<see cref="EncodingConfidence.Ambiguous"/>).</item>
    /// </list>
    /// Null and empty quality strings are skipped; an input with no characters yields Phred+33 /
    /// <see cref="EncodingConfidence.Ambiguous"/> with <c>CharactersExamined = 0</c>. For a single read
    /// the chosen encoding is identical to <see cref="DetectEncoding(string)"/> — only the file-level
    /// aggregation and the confidence signal are new.
    /// </summary>
    /// <param name="qualityStrings">The quality lines of the reads to consider together.</param>
    /// <returns>The detected encoding with its confidence and the ASCII span examined.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="qualityStrings"/> is null.</exception>
    public static EncodingDetectionResult DetectEncoding(IEnumerable<string> qualityStrings)
    {
        if (qualityStrings is null)
            throw new ArgumentNullException(nameof(qualityStrings));

        int min = int.MaxValue;
        int max = int.MinValue;
        long count = 0;

        foreach (var qualityString in qualityStrings)
        {
            if (string.IsNullOrEmpty(qualityString))
                continue;
            foreach (char c in qualityString)
            {
                if (c < min) min = c;
                if (c > max) max = c;
                count++;
            }
        }

        if (count == 0)
            return new EncodingDetectionResult(QualityEncoding.Phred33, EncodingConfidence.Ambiguous, 0, 0, 0);

        QualityEncoding encoding;
        EncodingConfidence confidence;
        if (min < Phred64Offset)
        {
            // A character below ASCII 64 cannot occur in Phred+64 (range 64-126) -> proven Phred+33.
            encoding = QualityEncoding.Phred33;
            confidence = EncodingConfidence.Definitive;
        }
        else if (max > Phred33PlausibleMaxChar)
        {
            // No sub-64 character and a character above the Phred+33 ceiling -> inferred Phred+64.
            encoding = QualityEncoding.Phred64;
            confidence = EncodingConfidence.Inferred;
        }
        else
        {
            // Every character is in the overlap range (ASCII 64-74) -> genuinely ambiguous.
            encoding = QualityEncoding.Phred33;
            confidence = EncodingConfidence.Ambiguous;
        }

        return new EncodingDetectionResult(encoding, confidence, min, max, count);
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

        int aboveQ20 = phredScores.Count(q => q >= Q20Threshold);
        int aboveQ30 = phredScores.Count(q => q >= Q30Threshold);

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
            PercentAboveQ20: PercentScale * aboveQ20 / phredScores.Length,
            PercentAboveQ30: PercentScale * aboveQ30 / phredScores.Length,
            PerPositionMeanQuality: positionMean);
    }

    /// <summary>
    /// Calculates the Q30 percentage: the percentage of bases in the quality string whose
    /// decoded Phred score is greater than or equal to 30 (the NGS run-quality benchmark; a
    /// Q30 base has an estimated 1-in-1000 error, i.e. 99.9% accuracy). The threshold is
    /// inclusive — a base at exactly Q30 is counted. Equivalent to
    /// <see cref="CalculateStatistics(string, QualityEncoding)"/>.<see cref="QualityStatistics.PercentAboveQ30"/>.
    /// </summary>
    /// <param name="qualityString">FASTQ quality line. Null or empty yields 0.</param>
    /// <param name="encoding">Phred+33 (default), Phred+64, or Auto.</param>
    /// <returns>Percentage of bases with Phred score &gt;= 30, in [0, 100]; 0 for empty/null input.</returns>
    public static double CalculateQ30Percentage(
        string qualityString,
        QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(qualityString))
            return 0.0;

        var phred = QualityStringToPhred(qualityString, encoding);
        if (phred.Length == 0)
            return 0.0;

        int aboveQ30 = phred.Count(q => q >= Q30Threshold);
        return PercentScale * aboveQ30 / phred.Length;
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

