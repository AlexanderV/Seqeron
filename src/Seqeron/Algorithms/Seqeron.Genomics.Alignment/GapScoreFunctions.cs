namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The families of gap-score functions <c>f(position, length)</c> that the interfaces without a callable (MCP, JSON, the oracle of the Waterman-Smith-Beyer aligner) can name for <see cref="PairwiseAligner.SetGapFunction"/> (B09 WP67, F307):
/// <c>lin</c> -(a + b·length), <c>log</c> -(a + b·ln(length + 1)), <c>sq</c> -(a + b·length²), <c>pos</c> -(a + b·position + c·length), <c>step</c> -a for a gap up to c positions and -b beyond,
/// <c>tab</c> -(((position·7 + length·3) mod b)·0.5 + a), <c>zero</c> 0, <c>bonus</c> a·length - b (a gap that scores above zero). The position is the position in the sequence the gap lies in (the target for insertions, the query for deletions).
/// </summary>
public static class GapScoreFunctions
{
    /// <summary>The names of the families.</summary>
    public static IReadOnlyList<string> Kinds { get; } = ["lin", "log", "sq", "pos", "step", "tab", "zero", "bonus"];

    /// <summary>The gap-score function of a family with its parameters (missing parameters are 0).</summary>
    /// <exception cref="ArgumentException">Unknown family.</exception>
    public static Func<int, int, double> Create(string kind, params double[] p)
    {
        ArgumentNullException.ThrowIfNull(kind);
        double P(int i) => p is not null && i < p.Length ? p[i] : 0.0;
        return kind switch
        {
            "lin" => (_, l) => -(P(0) + P(1) * l),
            "log" => (_, l) => -(P(0) + P(1) * Math.Log(l + 1)),
            "sq" => (_, l) => -(P(0) + P(1) * l * l),
            "pos" => (i, l) => -(P(0) + P(1) * i + P(2) * l),
            "step" => (_, l) => l <= (int)P(2) ? -P(0) : -P(1),
            "tab" => (i, l) => -(((i * 7 + l * 3) % Math.Max(1, (int)P(1))) * 0.5 + P(0)),
            "zero" => (_, _) => 0.0,
            "bonus" => (_, l) => P(0) * l - P(1),
            _ => throw new ArgumentException($"Unknown gap-score function '{kind}': {string.Join(", ", Kinds)}.", nameof(kind)),
        };
    }
}
