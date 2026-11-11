using System.Text;

namespace SuffixTree.Esa;

/// <summary>
/// Creates <see cref="EsaSuffixTree"/> instances, as <c>PersistentSuffixTreeFactory</c> does for the v6 tree. A factory that offers the Ukkonen tree, the v6 tree and the ESA together
/// would have to reference all three projects; the cross-backend selection used by the tests lives in the test project (<c>SuffixTreeBackends</c>).
/// </summary>
public static class EsaSuffixTreeFactory
{
    private const int ReadChunk = 1 << 16;

    /// <summary>Creates the tree of a text.</summary>
    /// <param name="text">The text.</param>
    /// <returns>The tree, as <see cref="ISuffixTree"/> (it is <see cref="IDisposable"/>).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    public static ISuffixTree Create(string text) => CreateEsa(text);

    /// <summary>Creates the tree of a text source.</summary>
    /// <param name="text">The text source; a <see cref="StringTextSource"/> is used without a copy, any other source is read once.</param>
    /// <returns>The tree, as <see cref="ISuffixTree"/> (it is <see cref="IDisposable"/>).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="text"/> is longer than <see cref="SuffixTreeLimits.MaxTextLength"/> characters (checked before any content is read).</exception>
    public static ISuffixTree Create(ITextSource text) => CreateEsa(text);

    /// <summary>Creates the tree of a text and returns the concrete disposable type, so callers can write <c>using var tree = ...</c> without a cast.</summary>
    /// <param name="text">The text.</param>
    /// <returns>The tree.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    public static EsaSuffixTree CreateEsa(string text) => EsaSuffixTree.Build(text);

    /// <summary>Creates the tree of a text source and returns the concrete disposable type.</summary>
    /// <param name="text">The text source.</param>
    /// <returns>The tree.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="text"/> is longer than <see cref="SuffixTreeLimits.MaxTextLength"/> characters (checked before any content is read).</exception>
    public static EsaSuffixTree CreateEsa(ITextSource text)
    {
        ArgumentNullException.ThrowIfNull(text);
        SuffixTreeLimits.ThrowIfTextTooLong(text.Length, nameof(text));
        return EsaSuffixTree.Build(text is StringTextSource source ? source.Value : Materialize(text));
    }

    /// <summary>Reads a text source into a string in chunks (nothing is preallocated for the whole length, so a source that cannot be read costs nothing).</summary>
    private static string Materialize(ITextSource text)
    {
        int length = text.Length;
        var builder = new StringBuilder();
        for (int start = 0; start < length; start += ReadChunk)
            builder.Append(text.Slice(start, Math.Min(ReadChunk, length - start)));
        return builder.ToString();
    }
}
