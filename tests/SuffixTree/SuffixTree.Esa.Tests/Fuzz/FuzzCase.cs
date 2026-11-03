using System.Globalization;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>Builds the suffix array (n + 1 layout) and compact LCP of a text: the system under test of the fuzz harness.</summary>
internal delegate (int[] Sa, CompactLcp Lcp) EsaBuilder(string text);

/// <summary>
/// One reproducible fuzz input: the text is a pure function of (<see cref="Family"/>, <see cref="Seed"/>, <see cref="N"/>, <see cref="Sigma"/>,
/// <see cref="Param"/>) (<see cref="FuzzTexts.Generate"/>), so a failing case is fully described by its <see cref="Id"/>
/// <c>family:seed:n:sigma:param</c> (seed in hex), independent of the run, the tier and the time budget that found it.
/// </summary>
internal readonly record struct FuzzCase(string Family, ulong Seed, int N, int Sigma, int Param)
{
    /// <summary>The textual identity <c>family:seed(16 hex digits):n:sigma:param</c>.</summary>
    public string Id => string.Create(CultureInfo.InvariantCulture, $"{Family}:{Seed:x16}:{N}:{Sigma}:{Param}");

    /// <summary>Parses <see cref="Id"/>; throws <see cref="FormatException"/> for anything else.</summary>
    public static FuzzCase Parse(string id)
    {
        ArgumentNullException.ThrowIfNull(id);
        string[] f = id.Split(':');
        if (f.Length != 5
            || f[0].Length == 0
            || !ulong.TryParse(f[1], NumberStyles.AllowHexSpecifier, CultureInfo.InvariantCulture, out ulong seed)
            || !int.TryParse(f[2], NumberStyles.None, CultureInfo.InvariantCulture, out int n)
            || !int.TryParse(f[3], NumberStyles.None, CultureInfo.InvariantCulture, out int sigma)
            || !int.TryParse(f[4], NumberStyles.None, CultureInfo.InvariantCulture, out int param))
            throw new FormatException($"not a fuzz case id: '{id}'");
        return new FuzzCase(f[0], seed, n, sigma, param);
    }

    /// <summary>The text of this case.</summary>
    public string Text() => FuzzTexts.Generate(this);
}
