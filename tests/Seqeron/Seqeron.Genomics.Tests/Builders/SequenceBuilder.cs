using System.Text;

namespace Seqeron.Genomics.Tests.Builders;

/// <summary>
/// Fluent builder for constructing test sequences with domain-meaningful intent.
///
/// Usage:
///   var seq = SequenceBuilder.Create()
///       .Ordered(10)        // 10× W (most order-promoting residue)
///       .Disordered(20)     // 20× P (most disorder-promoting residue)
///       .Ordered(10)
///       .Build();
///
///   var dna = SequenceBuilder.Dna()
///       .GcRich(50)
///       .AtRich(50)
///       .BuildDna();
/// </summary>
public sealed class SequenceBuilder
{
    private readonly StringBuilder _sb = new();
    private int _seed = 42;

    public static SequenceBuilder Create() => new();
    public static SequenceBuilder Dna() => new();
    public static SequenceBuilder Protein() => new();

    // ── DNA helpers ──

    public SequenceBuilder GcRich(int length)
    {
        var rng = new Random(_seed++);
        for (int i = 0; i < length; i++)
            _sb.Append(rng.Next(2) == 0 ? 'G' : 'C');
        return this;
    }

    public SequenceBuilder AtRich(int length)
    {
        var rng = new Random(_seed++);
        for (int i = 0; i < length; i++)
            _sb.Append(rng.Next(2) == 0 ? 'A' : 'T');
        return this;
    }

    public SequenceBuilder CpGIsland(int length)
    {
        for (int i = 0; i < length / 2; i++)
            _sb.Append("CG");
        if (length % 2 != 0)
            _sb.Append('C');
        return this;
    }

    // ── Protein/disorder helpers ──

    /// <summary>W (Trp): lowest TOP-IDP propensity (−0.884) → most order-promoting.</summary>
    public SequenceBuilder Ordered(int length) => Repeat('W', length);

    /// <summary>P (Pro): highest TOP-IDP propensity (0.987) → most disorder-promoting.</summary>
    public SequenceBuilder Disordered(int length) => Repeat('P', length);

    /// <summary>E (Glu): intermediate disorder propensity (0.736).</summary>
    public SequenceBuilder Intermediate(int length) => Repeat('E', length);

    // ── General helpers ──

    public SequenceBuilder Repeat(char residue, int count)
    {
        _sb.Append(new string(residue, count));
        return this;
    }

    public SequenceBuilder Repeat(string motif, int count)
    {
        for (int i = 0; i < count; i++)
            _sb.Append(motif);
        return this;
    }

    public SequenceBuilder Append(string fragment)
    {
        _sb.Append(fragment);
        return this;
    }

    public SequenceBuilder Random(int length, string alphabet = "ACGT")
    {
        var rng = new Random(_seed++);
        for (int i = 0; i < length; i++)
            _sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return this;
    }

    public SequenceBuilder WithSeed(int seed)
    {
        _seed = seed;
        return this;
    }

    // ── Build ──

    public string Build() => _sb.ToString();

    public DnaSequence BuildDna() => new(_sb.ToString());

    public int Length => _sb.Length;
}
