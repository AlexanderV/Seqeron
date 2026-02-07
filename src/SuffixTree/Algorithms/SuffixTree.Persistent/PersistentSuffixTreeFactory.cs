namespace SuffixTree.Persistent;

/// <summary>
/// Provides convenience methods for creating and loading persistent suffix trees.
/// <para>
/// The binary storage format uses hybrid continuation: construction always starts
/// with the Compact (32-bit, 28-byte node) layout for maximum cache efficiency.
/// If the tree outgrows the ~4 GB uint32 address space, the builder transparently
/// switches to the Large (64-bit) layout mid-build, using a jump table for
/// cross-zone references. No rebuild is needed.
/// </para>
/// </summary>
public static class PersistentSuffixTreeFactory
{
    /// <summary>
    /// Creates a new persistent suffix tree from the specified text source.
    /// If a filePath is provided, it uses Memory-Mapped Files; otherwise, it uses heap memory.
    /// <para>
    /// Storage format is chosen automatically: trees start in Compact (32-bit) mode
    /// and seamlessly transition to Large (64-bit) if the storage exceeds the uint32
    /// address space, with a jump table bridging cross-zone references (version 5 Hybrid).
    /// </para>
    /// </summary>
    /// <param name="text">The text source to build the tree from.</param>
    /// <param name="filePath">Optional file path for MMF storage.</param>
    /// <returns>An implementation of ISuffixTree.</returns>
    public static ISuffixTree Create(ITextSource text, string? filePath = null)
    {
        ArgumentNullException.ThrowIfNull(text);
        return CreateCore(text, filePath, NodeLayout.CompactMaxOffset);
    }

    /// <summary>
    /// Core build logic with configurable compact offset limit (for testing).
    /// Starts with Compact layout; if the tree outgrows the limit, the builder
    /// transitions to Large layout mid-build (hybrid continuation).
    /// </summary>
    internal static ISuffixTree CreateCore(ITextSource text, string? filePath, long compactOffsetLimit)
    {
        var storage = CreateStorage(filePath);
        try
        {
            var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
            builder.CompactOffsetLimit = compactOffsetLimit;
            long rootOffset = builder.Build(text);

            if (storage is MappedFileStorageProvider mapped)
                mapped.TrimToSize();

            // All trees start with Compact layout; hybrid v5 uses mixed zones
            // but the root and header are always Compact.
            NodeLayout layout = NodeLayout.Compact;
            return new PersistentSuffixTree(storage, rootOffset, text, layout,
                builder.TransitionOffset, builder.JumpTableStart, builder.JumpTableEnd,
                builder.DeepestInternalNodeOffset);
        }
        catch
        {
            storage.Dispose();
            CleanupFile(filePath);
            throw;
        }
    }

    /// <summary>
    /// Loads an existing persistent suffix tree from a file using Memory-Mapped Files.
    /// The format (Compact v4, Large v3, or Hybrid v5) is detected automatically from the file header.
    /// </summary>
    /// <param name="filePath">The path to the existing suffix tree file.</param>
    /// <returns>An implementation of ISuffixTree.</returns>
    public static ISuffixTree Load(string filePath)
    {
        if (string.IsNullOrEmpty(filePath))
            throw new ArgumentException("File path must be provided.", nameof(filePath));

        var storage = new MappedFileStorageProvider(filePath, readOnly: true);
        try
        {
            return PersistentSuffixTree.Load(storage);
        }
        catch
        {
            storage.Dispose();
            throw;
        }
    }

    private static IStorageProvider CreateStorage(string? filePath)
        => !string.IsNullOrEmpty(filePath)
            ? new MappedFileStorageProvider(filePath)
            : new HeapStorageProvider();

    private static void CleanupFile(string? filePath)
    {
        if (!string.IsNullOrEmpty(filePath))
        {
            try { File.Delete(filePath); }
            catch (IOException) { }
            catch (UnauthorizedAccessException) { }
        }
    }
}
