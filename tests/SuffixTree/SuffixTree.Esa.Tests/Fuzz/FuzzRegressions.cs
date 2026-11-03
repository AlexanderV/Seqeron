using System.Globalization;
using System.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>One line of the committed regression corpus: the original failing case, the failed operation, the minimized text.</summary>
internal sealed record FuzzRegression(FuzzCase Case, string Op, string MinimizedText);

/// <summary>
/// The regression corpus of failing fuzz cases, committed as text (<c>tests/SuffixTree/SuffixTree.Esa.Tests/Fuzz/regressions.txt</c>). Format: one
/// case per line, three TAB-separated fields <c>case-id</c>, <c>op</c>, <c>units</c> where units is the minimized text as space-separated hex code
/// units (<c>0061 0062</c>, <c>-</c> for the empty text); blank lines and lines starting with <c>#</c> are ignored. Every entry is replayed twice by
/// the default test run: the original case regenerated from its id, and the minimized text.
/// </summary>
internal static class FuzzRegressions
{
    /// <summary>Path of the committed corpus.</summary>
    internal static string DefaultPath() =>
        Path.Combine(CorpusManifest.RepositoryRoot(), "tests", "SuffixTree", "SuffixTree.Esa.Tests", "Fuzz", "regressions.txt");

    /// <summary>The line to append to the corpus for a failure.</summary>
    internal static string Format(FuzzFailure failure) => $"{failure.Case.Id}\t{failure.Op}\t{EncodeUnits(failure.MinimizedText)}";

    /// <summary>Hex code units, space separated; "-" for the empty text.</summary>
    internal static string EncodeUnits(string text)
    {
        if (text.Length == 0) return "-";
        var sb = new StringBuilder(text.Length * 5);
        foreach (char c in text)
        {
            if (sb.Length > 0) sb.Append(' ');
            sb.Append(((int)c).ToString("x4", CultureInfo.InvariantCulture));
        }
        return sb.ToString();
    }

    /// <summary>Inverse of <see cref="EncodeUnits"/>.</summary>
    internal static string DecodeUnits(string units)
    {
        if (units == "-") return string.Empty;
        string[] parts = units.Split(' ', StringSplitOptions.RemoveEmptyEntries);
        var c = new char[parts.Length];
        for (int i = 0; i < c.Length; i++)
            c[i] = (char)ushort.Parse(parts[i], NumberStyles.AllowHexSpecifier, CultureInfo.InvariantCulture);
        return new string(c);
    }

    /// <summary>Parses corpus text; throws <see cref="FormatException"/> for a malformed line.</summary>
    internal static List<FuzzRegression> Parse(string text)
    {
        var list = new List<FuzzRegression>();
        int lineNo = 0;
        foreach (string raw in text.Split('\n'))
        {
            lineNo++;
            string line = raw.TrimEnd('\r');
            if (line.Trim().Length == 0 || line.StartsWith('#')) continue;
            string[] f = line.Split('\t');
            if (f.Length != 3) throw new FormatException($"regressions line {lineNo}: expected 3 tab-separated fields, got {f.Length}");
            list.Add(new FuzzRegression(FuzzCase.Parse(f[0]), f[1], DecodeUnits(f[2])));
        }
        return list;
    }

    /// <summary>Loads the committed corpus.</summary>
    internal static List<FuzzRegression> Load() => Parse(File.ReadAllText(DefaultPath()));
}
