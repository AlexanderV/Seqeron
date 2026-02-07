using System;

namespace SuffixTree.Persistent;

/// <summary>
/// Provides convenience methods for creating and loading persistent suffix trees.
/// </summary>
public static class PersistentSuffixTreeFactory
{
    /// <summary>
    /// Creates a new persistent suffix tree from the specified text source.
    /// If a filePath is provided, it uses Memory-Mapped Files; otherwise, it uses heap memory.
    /// <para>
    /// The binary format is selected automatically: texts up to ~50 M characters use
    /// a compact 32-bit layout (28-byte nodes, ~30 %% smaller files, better cache locality);
    /// larger texts switch to 64-bit offsets transparently.
    /// </para>
    /// </summary>
    /// <param name="text">The text source to build the tree from.</param>
    /// <param name="filePath">Optional file path for MMF storage.</param>
    /// <returns>An implementation of ISuffixTree.</returns>
    public static ISuffixTree Create(ITextSource text, string? filePath = null)
    {
        ArgumentNullException.ThrowIfNull(text);

        IStorageProvider storage = !string.IsNullOrEmpty(filePath)
            ? new MappedFileStorageProvider(filePath)
            : new HeapStorageProvider();

        try
        {
            var layout = NodeLayout.ForTextLength(text.Length);
            var builder = new PersistentSuffixTreeBuilder(storage, layout);
            long rootOffset = builder.Build(text);

            // Trim MMF file to actual data size (reclaim ~50% unused capacity)
            if (storage is MappedFileStorageProvider mapped)
                mapped.TrimToSize();

            return new PersistentSuffixTree(storage, rootOffset, text, layout);
        }
        catch
        {
            storage.Dispose();
            throw;
        }
    }

    /// <summary>
    /// Loads an existing persistent suffix tree from a file using Memory-Mapped Files.
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
}
