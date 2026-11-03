namespace SuffixTree.Persistent;

public sealed partial class PersistentSuffixTree
{
    /// <summary>
    /// Strict structural validation of the whole tree: every offset, child array, key, edge, leaf count,
    /// suffix link and header summary is checked once, with bounds checks performed before any dereference.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <see cref="Load(IStorageProvider)"/> validates only the header, which keeps cold loads at a few
    /// milliseconds. Query code reads nodes through raw pointers, so a file that is corrupt behind a valid
    /// header (disk error, truncated copy edited in place, hostile file) must be checked with this method
    /// before it is queried. The check guarantees memory-safe traversal; it cannot prove that the index
    /// matches the text (use the SHA-256 hash of the serializer export for that).
    /// </para>
    /// <para>
    /// Cost: one sequential pass over all nodes (several seconds per 10 million characters on a mapped file)
    /// and a bitmap of one bit per four bytes of storage.
    /// </para>
    /// </remarks>
    /// <exception cref="InvalidOperationException">The structure is inconsistent ("Invalid storage format: ...").</exception>
    /// <exception cref="ObjectDisposedException">The tree is disposed.</exception>
    public void ValidateStructure()
    {
        ThrowIfDisposed();
        PersistentStructureValidator.Validate(
            _storage, _hybrid, _rootOffset, _textSource, _deepestInternalNodeOffset, _lrsDepth);
    }

    /// <summary>
    /// Loads a persistent suffix tree from storage, validating the header and, when requested, the whole structure.
    /// </summary>
    /// <param name="storage">The storage that holds the tree.</param>
    /// <param name="validateStructure">
    /// <see langword="true"/> to run <see cref="ValidateStructure"/> before returning; <see langword="false"/>
    /// is identical to <see cref="Load(IStorageProvider)"/> (header-only validation).
    /// </param>
    /// <returns>The loaded tree.</returns>
    /// <exception cref="InvalidOperationException">
    /// The header or, with <paramref name="validateStructure"/>, the structure is invalid. When the structure
    /// is rejected, the partially opened tree and the <paramref name="storage"/> are disposed.
    /// </exception>
    public static PersistentSuffixTree Load(IStorageProvider storage, bool validateStructure)
    {
        PersistentSuffixTree tree = Load(storage);
        if (!validateStructure)
            return tree;

        try
        {
            tree.ValidateStructure();
        }
        catch
        {
            tree.Dispose();
            throw;
        }
        return tree;
    }
}
