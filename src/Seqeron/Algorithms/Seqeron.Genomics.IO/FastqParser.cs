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
    /// <remarks>
    /// With <see cref="QualityEncoding.Auto"/> the Phred offset is detected once for the whole file
    /// (two streaming passes; see <see cref="Parse(TextReader, QualityEncoding)"/>).
    /// </remarks>
    /// <exception cref="FormatException">Thrown (during enumeration) when a record is malformed.</exception>
    public static IEnumerable<FastqRecord> ParseFile(string filePath, QualityEncoding encoding = QualityEncoding.Auto)
    {
        if (string.IsNullOrEmpty(filePath) || !File.Exists(filePath))
            yield break;

        var actualEncoding = encoding;
        if (encoding == QualityEncoding.Auto)
        {
            using var detectReader = new StreamReader(filePath);
            actualEncoding = DetectFileEncoding(ReadRawRecords(detectReader));
        }

        using var reader = new StreamReader(filePath);
        foreach (var raw in ReadRawRecords(reader))
            yield return CreateRecord(raw, actualEncoding);
    }

    /// <summary>
    /// Parses FASTQ records from text content.
    /// </summary>
    /// <remarks>
    /// With <see cref="QualityEncoding.Auto"/> the Phred offset is detected once for the whole content
    /// (see <see cref="Parse(TextReader, QualityEncoding)"/>).
    /// </remarks>
    /// <exception cref="FormatException">Thrown (during enumeration) when a record is malformed.</exception>
    public static IEnumerable<FastqRecord> Parse(string content, QualityEncoding encoding = QualityEncoding.Auto)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        var actualEncoding = encoding == QualityEncoding.Auto
            ? DetectFileEncoding(ReadRawRecords(new StringReader(content)))
            : encoding;

        foreach (var raw in ReadRawRecords(new StringReader(content)))
            yield return CreateRecord(raw, actualEncoding);
    }

    /// <summary>
    /// Parses FASTQ records from a TextReader, following the Sanger FASTQ definition (Cock et al., 2010,
    /// NAR 38:1767) as realised by Biopython's <c>FastqGeneralIterator</c>:
    /// <list type="bullet">
    /// <item>each record starts with an '@' title line (blank lines between records are skipped);</item>
    /// <item>sequence lines are read up to the '+' line; whitespace inside the sequence is rejected;</item>
    /// <item>the '+' line may repeat the title, and if it does it must be identical;</item>
    /// <item>quality lines are read by length (a quality line may itself begin with '@') and the quality
    /// string must have exactly as many symbols as the sequence has letters.</item>
    /// </list>
    /// With <see cref="QualityEncoding.Auto"/> the Phred offset is a property of the whole file (FastQC
    /// derives it from the lowest quality character over all reads): it is determined once with
    /// <see cref="QualityScoreAnalyzer.DetectEncoding(IEnumerable{string})"/> and applied to every record.
    /// Because a <see cref="TextReader"/> cannot be rewound, Auto mode buffers the raw records of this
    /// reader before yielding; pass an explicit encoding (or use <see cref="ParseFile"/>) to stream.
    /// </summary>
    /// <exception cref="FormatException">Thrown (during enumeration) when a record is malformed or a quality
    /// character is outside the range of the encoding.</exception>
    public static IEnumerable<FastqRecord> Parse(TextReader reader, QualityEncoding encoding = QualityEncoding.Auto)
    {
        ArgumentNullException.ThrowIfNull(reader);
        return ParseReader(reader, encoding);
    }

    private static IEnumerable<FastqRecord> ParseReader(TextReader reader, QualityEncoding encoding)
    {
        if (encoding != QualityEncoding.Auto)
        {
            foreach (var raw in ReadRawRecords(reader))
                yield return CreateRecord(raw, encoding);
            yield break;
        }

        var buffered = ReadRawRecords(reader).ToList();
        var detected = DetectFileEncoding(buffered);
        foreach (var raw in buffered)
            yield return CreateRecord(raw, detected);
    }

    private readonly record struct RawFastqRecord(string Title, string Sequence, string Quality);

    /// <summary>
    /// Line-oriented FASTQ state machine (port of Biopython <c>FastqGeneralIterator</c>, 1.88). Yields
    /// (title, sequence, quality) with line breaks removed; throws <see cref="FormatException"/> on the
    /// same malformations Biopython rejects.
    /// </summary>
    private static IEnumerable<RawFastqRecord> ReadRawRecords(TextReader reader)
    {
        string? line = reader.ReadLine();
        while (true)
        {
            while (line != null && string.IsNullOrWhiteSpace(line))
                line = reader.ReadLine();
            if (line == null)
                yield break;

            if (line[0] != '@')
                throw new FormatException(
                    $"Records in FASTQ files must start with the '@' character (found line: '{Abbreviate(line)}').");
            string title = line.Substring(1).Trim();

            var sequenceBuilder = new StringBuilder();
            while (true)
            {
                line = reader.ReadLine();
                if (line == null)
                    throw new FormatException(sequenceBuilder.Length > 0
                        ? $"End of file without quality information for FASTQ record '{title}'."
                        : $"Unexpected end of file in FASTQ record '{title}'.");
                if (line.Length > 0 && line[0] == '+')
                    break;
                sequenceBuilder.Append(line.TrimEnd());
            }

            string plusTitle = line.Substring(1).Trim();
            if (plusTitle.Length > 0 && plusTitle != title)
                throw new FormatException(
                    $"FASTQ sequence and quality captions differ ('@{title}' vs '+{plusTitle}').");

            string sequence = sequenceBuilder.ToString();
            foreach (char c in sequence)
            {
                if (char.IsWhiteSpace(c))
                    throw new FormatException($"Whitespace is not allowed in the sequence of FASTQ record '{title}'.");
            }

            // At least one line of quality data must follow the '+' line.
            line = reader.ReadLine();
            if (line == null)
                throw new FormatException($"Unexpected end of file in FASTQ record '{title}' (no quality line).");

            var qualityBuilder = new StringBuilder(sequence.Length);
            while (line != null)
            {
                // A line starting with '@' only begins the next record once the quality is complete;
                // before that it is quality data whose first symbol is '@' (Phred+33 Q31 / Phred+64 Q0).
                if (line.Length > 0 && line[0] == '@' && qualityBuilder.Length >= sequence.Length)
                    break;
                qualityBuilder.Append(line.TrimEnd());
                line = reader.ReadLine();
            }

            string quality = qualityBuilder.ToString();
            if (quality.Length != sequence.Length)
                throw new FormatException(
                    $"Lengths of sequence and quality values differ for FASTQ record '{title}' " +
                    $"({sequence.Length} and {quality.Length}).");

            yield return new RawFastqRecord(title, sequence, quality);
        }
    }

    private static string Abbreviate(string line) => line.Length <= 40 ? line : line[..40] + "...";

    private static QualityEncoding DetectFileEncoding(IEnumerable<RawFastqRecord> records)
        => FromCanonical(QualityScoreAnalyzer.DetectEncoding(records.Select(r => r.Quality)).Encoding);

    private static FastqRecord CreateRecord(RawFastqRecord raw, QualityEncoding encoding)
    {
        var (id, description) = SequenceFormatHelper.SplitTitle(raw.Title);
        IReadOnlyList<int> scores;
        try
        {
            scores = DecodeQualityScores(raw.Quality, encoding);
        }
        catch (ArgumentOutOfRangeException ex)
        {
            throw new FormatException($"Invalid quality string in FASTQ record '{id}': {ex.Message}", ex);
        }
        return new FastqRecord(id, description ?? "", raw.Sequence, raw.Quality, scores);
    }

    #endregion

    #region Quality Encoding

    private static QualityScoreAnalyzer.QualityEncoding ToCanonical(QualityEncoding encoding) => encoding switch
    {
        QualityEncoding.Phred64 => QualityScoreAnalyzer.QualityEncoding.Phred64,
        QualityEncoding.Auto => QualityScoreAnalyzer.QualityEncoding.Auto,
        _ => QualityScoreAnalyzer.QualityEncoding.Phred33,
    };

    private static QualityEncoding FromCanonical(QualityScoreAnalyzer.QualityEncoding encoding)
        => encoding == QualityScoreAnalyzer.QualityEncoding.Phred64 ? QualityEncoding.Phred64 : QualityEncoding.Phred33;

    /// <summary>
    /// Detects the quality encoding of a single quality string. Delegates to the canonical
    /// <see cref="QualityScoreAnalyzer.DetectEncoding(string)"/>: a character below ASCII 64 proves
    /// Phred+33 regardless of its position; otherwise a character above ASCII 74 ('J' = Q41, the
    /// Illumina 1.8+ Phred+33 ceiling) infers Phred+64; a string confined to ASCII 64-74 is ambiguous
    /// and defaults to Phred+33 (guarded by <c>LimitationPolicy</c> "PARSE-FASTQ-001").
    /// For a whole file use <see cref="QualityScoreAnalyzer.DetectEncoding(IEnumerable{string})"/>,
    /// which is what <see cref="Parse(string, QualityEncoding)"/> does in Auto mode.
    /// </summary>
    public static QualityEncoding DetectEncoding(string qualityString)
    {
        if (string.IsNullOrEmpty(qualityString))
            return QualityEncoding.Phred33;

        return FromCanonical(QualityScoreAnalyzer.DetectEncoding(qualityString));
    }

    /// <summary>
    /// Decodes a quality string to Phred scores (Q = ASCII - offset; offset 33 or 64). Delegates to the
    /// canonical <see cref="QualityScoreAnalyzer.ParseQualityString"/>. A null or empty string yields no
    /// scores. <see cref="QualityEncoding.Auto"/> is resolved per string.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when a character decodes outside the valid
    /// Phred range of the encoding (0-93 for Phred+33, 0-62 for Phred+64; Cock et al., 2010) — the same
    /// characters Biopython rejects with <c>InvalidCharError</c>.</exception>
    public static IReadOnlyList<int> DecodeQualityScores(string qualityString, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        if (string.IsNullOrEmpty(qualityString))
            return Array.Empty<int>();

        return QualityScoreAnalyzer.ParseQualityString(qualityString, ToCanonical(encoding));
    }

    /// <summary>
    /// Encodes Phred scores to a quality string (char = Q + offset). Scores are first clamped to the
    /// representable range of the encoding (0-93 Phred+33, 0-62 Phred+64) — Biopython likewise truncates
    /// at 93 / 62 when writing — and then encoded by the canonical
    /// <see cref="QualityScoreAnalyzer.ToQualityString"/>. Auto encodes as Phred+33.
    /// </summary>
    public static string EncodeQualityScores(IEnumerable<int> scores, QualityEncoding encoding = QualityEncoding.Phred33)
    {
        ArgumentNullException.ThrowIfNull(scores);
        int maxScore = encoding == QualityEncoding.Phred64 ? MaxPhred64Score : MaxPhred33Score;
        var clamped = scores.Select(score => Math.Clamp(score, 0, maxScore)).ToArray();
        return QualityScoreAnalyzer.ToQualityString(clamped, ToCanonical(encoding));
    }

    // Highest Phred score representable in each encoding (ASCII 126 '~' minus the offset; Cock et al., 2010).
    private const int MaxPhred33Score = 93;
    private const int MaxPhred64Score = 62;

    /// <summary>
    /// Converts a Phred score to an error probability, p = 10^(-Q/10). Delegates to the canonical
    /// <see cref="QualityScoreAnalyzer.PhredToErrorProbability"/>.
    /// </summary>
    public static double PhredToErrorProbability(int phredScore)
        => QualityScoreAnalyzer.PhredToErrorProbability(phredScore);

    /// <summary>
    /// Converts an error probability to a Phred score, Q = round(-10·log10(p)), capped at Q93 — the
    /// highest score representable in Sanger/Phred+33 FASTQ (ASCII 126 − 33). p = 0 maps to Q93, and any
    /// p small enough to exceed Q93 is capped there too, so the mapping stays monotone.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="errorProbability"/> is
    /// NaN or outside [0, 1].</exception>
    public static int ErrorProbabilityToPhred(double errorProbability)
    {
        if (double.IsNaN(errorProbability) || errorProbability < 0 || errorProbability > 1)
            throw new ArgumentOutOfRangeException(nameof(errorProbability), errorProbability,
                "Error probability must be in [0, 1].");
        if (errorProbability == 0)
            return MaxPhred33Score;
        return (int)Math.Min(MaxPhred33Score, Math.Round(-10 * Math.Log10(errorProbability)));
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

        var (start, end) = QualityScoreAnalyzer.FindQualityTrimBounds(
            record.QualityScores, record.Sequence.Length, minQuality);

        if (start >= end)
        {
            // Entire sequence trimmed
            return new FastqRecord(record.Id, record.Description, "", "", Array.Empty<int>());
        }

        return Slice(record, start, end);
    }

    private static FastqRecord Slice(FastqRecord record, int start, int end) => new(
        record.Id,
        record.Description,
        record.Sequence[start..end],
        record.QualityString[start..end],
        record.QualityScores.Skip(start).Take(end - start).ToList());

    /// <summary>
    /// Removes a 3' adapter and everything after it (cutadapt regular 3' adapter, <c>-a ADAPTER -e 0
    /// -O minOverlap</c>): the read is cut at the leftmost position i where the read from i onward matches
    /// the adapter exactly (case-insensitive) — either a full adapter occurrence anywhere, including at
    /// position 0 (the whole read is removed), or an adapter prefix of at least
    /// <paramref name="minOverlap"/> bases running off the 3' end.
    /// </summary>
    /// <remarks>Exact matching only (cutadapt's default error rate 0.1 is not modelled).</remarks>
    public static FastqRecord TrimAdapter(FastqRecord record, string adapter, int minOverlap = 5)
    {
        if (string.IsNullOrEmpty(adapter) || adapter.Length < minOverlap)
            return record;

        var sequence = record.Sequence;
        int minLength = Math.Max(1, minOverlap);
        for (int i = 0; i <= sequence.Length - minLength; i++)
        {
            int overlap = Math.Min(adapter.Length, sequence.Length - i);
            if (string.Compare(sequence, i, adapter, 0, overlap, StringComparison.OrdinalIgnoreCase) == 0)
                return Slice(record, 0, i);
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
        // UTF-8 without a byte-order mark: a BOM before the first '@' makes Biopython and cutadapt reject the file.
        using var writer = new StreamWriter(filePath, false, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
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

