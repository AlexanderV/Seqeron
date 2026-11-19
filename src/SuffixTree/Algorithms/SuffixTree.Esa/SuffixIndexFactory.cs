using SuffixTree.Persistent;

namespace SuffixTree.Esa;

/// <summary>
/// The one entry point that opens a suffix index file of either kind (ADR-003 section 8): a SESA container (this engine) or a v6 persistent suffix tree (<see cref="PersistentSuffixTreeFactory"/>), told apart by the
/// first eight bytes of the file.
/// </summary>
public static class SuffixIndexFactory
{
    /// <summary>Opens the index file at <paramref name="path"/>.</summary>
    /// <param name="path">The file.</param>
    /// <param name="options">Validation level and access mode of a SESA container (default <see cref="EsaOpenOptions"/>); a v6 file ignores them.</param>
    /// <returns>The index as <see cref="ISuffixTree"/>; it is <see cref="IDisposable"/> and holds the file mapped until it is disposed.</returns>
    /// <exception cref="ArgumentException"><paramref name="path"/> is null or empty.</exception>
    /// <exception cref="InvalidDataException">The file is of neither kind (<c>unknown suffix index kind</c>), or a SESA container is corrupt at the requested level.</exception>
    /// <exception cref="NotSupportedException">A SESA container of a larger major version, or one that needs a feature this reader does not have.</exception>
    public static ISuffixTree Open(string path, EsaOpenOptions? options = null) => throw new NotImplementedException();
}
