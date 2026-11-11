namespace SuffixTree.Esa;

/// <summary>STUB of work package 1b.6 (tests first).</summary>
public static class EsaSuffixTreeFactory
{
    /// <summary>Creates the tree of a text.</summary>
    /// <param name="text">The text.</param>
    /// <returns>The tree.</returns>
    public static ISuffixTree Create(string text) => throw new NotImplementedException();

    /// <summary>Creates the tree of a text source.</summary>
    /// <param name="text">The text source.</param>
    /// <returns>The tree.</returns>
    public static ISuffixTree Create(ITextSource text) => throw new NotImplementedException();

    /// <summary>Creates the tree and returns the concrete disposable type.</summary>
    /// <param name="text">The text.</param>
    /// <returns>The tree.</returns>
    public static EsaSuffixTree CreateEsa(string text) => throw new NotImplementedException();

    /// <summary>Creates the tree of a text source and returns the concrete disposable type.</summary>
    /// <param name="text">The text source.</param>
    /// <returns>The tree.</returns>
    public static EsaSuffixTree CreateEsa(ITextSource text) => throw new NotImplementedException();
}
