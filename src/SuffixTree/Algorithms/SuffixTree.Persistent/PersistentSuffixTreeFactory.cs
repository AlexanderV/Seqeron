using System;
using System.IO;

namespace SuffixTree.Persistent;

/// <summary>
/// Provides convenience methods for creating and loading persistent suffix trees.
/// <para>
/// The binary storage format is selected automatically:
/// construction always starts with the Compact (32-bit, 28-byte node) layout
/// for maximum cache efficiency. If the tree outgrows the ~4 GB uint32 address
/// space, construction is transparently restarted with the Large (64-bit) layout.
/// </para>
/// </summary>
public static class PersistentSuffixTreeFactory
{
    /// <summary>
    /// Creates a new persistent suffix tree from the specified text source.
    /// If a filePath is provided, it uses Memory-Mapped Files; otherwise, it uses heap memory.
    /// <para>
    /// Storage format is chosen automatically: trees start in Compact (32-bit) mode
    /// and promote to Large (64-bit) only if the storage exceeds the uint32 address space.
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
    /// Starts with Compact layout; on <see cref="CompactOverflowException"/>,
    /// discards partial storage and rebuilds with Large layout.
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

            return new PersistentSuffixTree(storage, rootOffset, text, NodeLayout.Compact);
        }
        catch (CompactOverflowException)
        {
            // Compact (uint32) address space exhausted — rebuild with Large (int64)
            storage.Dispose();
            CleanupFile(filePath);

            storage = CreateStorage(filePath);
            try
            {
                var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Large);
                long rootOffset = builder.Build(text);

                if (storage is MappedFileStorageProvider mapped)
                    mapped.TrimToSize();

                return new PersistentSuffixTree(storage, rootOffset, text, NodeLayout.Large);
            }
            catch
            {
                storage.Dispose();
                CleanupFile(filePath);
                throw;
            }
        }
        catch
        {
            storage.Dispose();
            throw;
        }
    }

    /// <summary>
    /// Loads an existing persistent suffix tree from a file using Memory-Mapped Files.
    /// The format (Compact v4 or Large v3) is detected automatically from the file header.
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
            try { File.Delete(filePath); } catch { /* best-effort cleanup */ }
    }
}
