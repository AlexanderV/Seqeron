using System;

namespace SuffixTree.Persistent;

/// <summary>
/// Provides convenience methods for creating and loading persistent suffix trees.
/// </summary>
public static class PersistentSuffixTreeFactory
{
    /// <summary>
    /// Creates a new persistent suffix tree from the specified text.
    /// If a filePath is provided, it uses Memory-Mapped Files; otherwise, it uses heap memory.
    /// </summary>
    /// <param name="text">The text to build the tree from.</param>
    /// <param name="filePath">Optional file path for MMF storage.</param>
    /// <returns>An implementation of ISuffixTree.</returns>
    public static ISuffixTree Create(string text, string? filePath = null)
        => Create(new StringTextSource(text), filePath);

    /// <summary>
    /// Creates a new persistent suffix tree from the specified text source.
    /// If a filePath is provided, it uses Memory-Mapped Files; otherwise, it uses heap memory.
    /// </summary>
    /// <param name="text">The text source to build the tree from.</param>
    /// <param name="filePath">Optional file path for MMF storage.</param>
    /// <returns>An implementation of ISuffixTree.</returns>
    public static ISuffixTree Create(ITextSource text, string? filePath = null)
    {
        IStorageProvider storage = !string.IsNullOrEmpty(filePath)
            ? new MappedFileStorageProvider(filePath)
            : new HeapStorageProvider();

        var builder = new PersistentSuffixTreeBuilder(storage);
        long rootOffset = builder.Build(text);

        return new PersistentSuffixTree(storage, rootOffset, text);
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
        return PersistentSuffixTree.Load(storage);
    }
}
