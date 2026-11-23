using System.Globalization;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The limit of the output of the suffix-tree repeat finders (<see cref="GenomicAnalyzer.FindRepeats(DnaSequence, int)"/>, <see cref="GenomicAnalyzer.FindMaximalRepeats(DnaSequence, int)"/>, their string and
/// both-strand variants; audit R14-4b, B09 F292). The output of "every repeated substring" is the sum over the repeats of (label length + occurrences): quadratic for a low-complexity sequence (poly-A of length
/// 20 000: 19 998 repeats, about 200 M occurrences, 30.5 s / 1.92 GB; 50 000 about 12 GB). The size is computed from the suffix tree (leaf counts and label lengths) BEFORE any repeat is materialised.
/// </summary>
/// <param name="MaxBytes">The estimated memory of the materialised repeats (2 bytes per label character, 4 per occurrence, 96 per repeat for the string and the list objects); the default 2^29 = 512 MiB lets poly-A of length 10 000 through (300 MB) and refuses 20 000 (1.2 GB).</param>
/// <param name="MaxOccurrences">The largest number of listed occurrences (positions over all the reported repeats); <see cref="long.MaxValue"/> (default) = no separate limit.</param>
public sealed record RepeatOutputLimits(long MaxBytes = RepeatOutputLimits.DefaultMaxBytes, long MaxOccurrences = long.MaxValue)
{
    /// <summary>The default memory limit of the output: 2^29 bytes (512 MiB).</summary>
    public const long DefaultMaxBytes = 1L << 29;

    /// <summary>The default limits (<see cref="DefaultMaxBytes"/>, no occurrence limit).</summary>
    public static RepeatOutputLimits Default { get; } = new();

    /// <summary>No limit: the caller accepts whatever the output costs.</summary>
    public static RepeatOutputLimits Unlimited { get; } = new(long.MaxValue, long.MaxValue);
}

/// <summary>
/// Thrown (when the enumeration starts) by the repeat finders whose output would exceed <see cref="RepeatOutputLimits"/> (audit R14-4b). Nothing has been materialised. An <see cref="ArgumentException"/>: the arguments
/// (the sequence with this <c>minLength</c> and mode) ask for more than the limits allow; raise <c>minLength</c>, use maximal repeats, or pass larger <see cref="RepeatOutputLimits"/>.
/// </summary>
public sealed class RepeatOutputLimitException : ArgumentException
{
    /// <summary>The number of repeats the call would report.</summary>
    public long Repeats { get; }

    /// <summary>The number of occurrences (positions) over all the repeats.</summary>
    public long Occurrences { get; }

    /// <summary>The number of label characters over all the repeats.</summary>
    public long Characters { get; }

    /// <summary>The estimated memory of the output in bytes.</summary>
    public long EstimatedBytes { get; }

    /// <summary>The limits that were exceeded.</summary>
    public RepeatOutputLimits Limits { get; }

    internal RepeatOutputLimitException(long repeats, long occurrences, long characters, long estimatedBytes, RepeatOutputLimits limits)
        : base(string.Create(CultureInfo.InvariantCulture,
            $"The repeats of this sequence would be {repeats:N0} substrings with {occurrences:N0} occurrences and {characters:N0} characters (about {estimatedBytes / 1048576:N0} MB), which exceeds the limit of {Describe(limits)} (RepeatOutputLimits). The output of every repeated substring of a low-complexity sequence is quadratic in its length (poly-A of length 20 000: 19 998 repeats, 200 M occurrences, 1.9 GB); raise minLength, ask for the maximal repeats only, or pass a larger RepeatOutputLimits (RepeatOutputLimits.Unlimited for none)."))
    {
        Repeats = repeats;
        Occurrences = occurrences;
        Characters = characters;
        EstimatedBytes = estimatedBytes;
        Limits = limits;
    }

    private static string Describe(RepeatOutputLimits limits)
    {
        string occurrences = string.Create(CultureInfo.InvariantCulture, $"{limits.MaxOccurrences:N0} occurrences (MaxOccurrences)");
        if (limits.MaxBytes == long.MaxValue)
            return occurrences;
        string bytes = string.Create(CultureInfo.InvariantCulture, $"{limits.MaxBytes / 1048576:N0} MB (MaxBytes)");
        return limits.MaxOccurrences == long.MaxValue ? bytes : bytes + " and " + occurrences;
    }

    private const long BytesPerCharacter = 2;
    private const long BytesPerOccurrence = 4;
    private const long BytesPerRepeat = 96;

    /// <summary>The check of the engine: the sizes from the loci of the tree (leaf count x number of repeat lengths of the node), thrown before anything is materialised.</summary>
    internal static void Check(IEnumerable<(int EdgeEnd, int ParentDepth, int Depth, int Leaves)> nodes, int effectiveMinLength, bool maximalOnly, RepeatOutputLimits limits)
    {
        long repeats = 0, occurrences = 0, characters = 0;
        foreach (var (_, parentDepth, depth, leaves) in nodes)
        {
            long count;
            long chars;
            if (maximalOnly)
            {
                count = 1;
                chars = depth;
            }
            else
            {
                long first = Math.Max(parentDepth + 1, effectiveMinLength);
                count = depth - first + 1;
                if (count <= 0)
                    continue;
                chars = (first + depth) * count / 2;       // the lengths first .. depth
            }

            repeats += count;
            occurrences = SaturatingAdd(occurrences, leaves * count);
            characters = SaturatingAdd(characters, chars);
        }

        long bytes = SaturatingAdd(SaturatingAdd(BytesPerCharacter * Math.Min(characters, long.MaxValue / 4), BytesPerOccurrence * Math.Min(occurrences, long.MaxValue / 8)), BytesPerRepeat * repeats);
        if (bytes > limits.MaxBytes || occurrences > limits.MaxOccurrences)
            throw new RepeatOutputLimitException(repeats, occurrences, characters, bytes, limits);
    }

    private static long SaturatingAdd(long a, long b) => a > long.MaxValue - b ? long.MaxValue : a + b;
}
