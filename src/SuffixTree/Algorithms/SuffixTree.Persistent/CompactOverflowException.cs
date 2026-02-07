namespace SuffixTree.Persistent;

/// <summary>
/// Thrown internally when the compact (uint32) storage format runs out of
/// addressable space during tree construction. Caught by
/// <see cref="PersistentSuffixTreeFactory"/> and <see cref="SuffixTreeSerializer"/>
/// to trigger an automatic rebuild with the Large (int64) layout.
/// </summary>
internal sealed class CompactOverflowException : Exception
{
    public CompactOverflowException()
        : base("Storage exceeded compact (uint32) address space.") { }
}
