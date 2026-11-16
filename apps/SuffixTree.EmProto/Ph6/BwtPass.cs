using System.Globalization;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Ph6;

/// <summary>The side file of a BWT file: the text length, the row of the suffix 0 (whose BWT symbol is the sentinel, stored as a dummy byte), the symbol counts.</summary>
internal sealed record BwtMeta(long N, long PrimaryRow, long[] Counts)
{
    internal int Sigma => Counts.Length;

    internal static string MetaPath(string bwtPath) => bwtPath + ".meta";

    internal string Text() => throw new NotImplementedException(N.ToString(System.Globalization.CultureInfo.InvariantCulture));

    internal static BwtMeta Parse(string text) => throw new NotImplementedException();

    internal static BwtMeta Read(string bwtPath) => throw new NotImplementedException();
}

internal sealed record BwtOptions
{
    /// <summary>Bytes of the suffix array read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>Window of the mapped writer of the BWT file.</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    internal ArenaBudget? Budget { get; init; }
}

/// <summary>
/// The BWT of the text from its suffix array file and the packed text (route (ii) of E9: a separate pass over the SA, the text is read at random through its mapping): <c>BWT[j] = T[SA[j] - 1]</c>, one byte per row
/// (<c>n + 1</c> rows, row 0 is the empty suffix, its symbol is the last symbol of the text); the row of the suffix 0 holds the dummy byte 0 and is recorded in the meta file.
/// </summary>
internal static class BwtPass
{
    internal static BwtMeta Build(string saPath, PackedText text, int sigma, string bwtPath, BwtOptions? options = null) => throw new NotImplementedException();
}
