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
        var storage = CreateStorage(filePath, text.Length);
        IStorageProvider? childStorage = null;
        IStorageProvider? depthStorage = null;
        string? tempChildPath = null;
        string? tempDepthPath = null;
        try
        {
            // For MMF builds, place child-edge and depth data in separate temp MMFs
            // so that the managed heap stays near-zero even for genome-scale inputs.
            if (!string.IsNullOrEmpty(filePath))
            {
                tempChildPath = filePath + ".children.tmp";
                long childCapacity = (long)text.Length * 40; // ~2 entries/char × 16 bytes + headroom
                childStorage = new MappedFileStorageProvider(tempChildPath, Math.Max(childCapacity, 65536));

                tempDepthPath = filePath + ".depth.tmp";
                long depthCapacity = (long)text.Length * 10; // ~2.5 nodes/char × 4 bytes
                depthStorage = new MappedFileStorageProvider(tempDepthPath, Math.Max(depthCapacity, 65536));
            }

            // v6 Compact: 24-byte nodes (no stored DepthFromRoot)
            var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact, childStorage, depthStorage);
            builder.CompactOffsetLimit = compactOffsetLimit;
            long rootOffset = builder.Build(text);

            if (storage is MappedFileStorageProvider mapped)
                mapped.TrimToSize();

            NodeLayout layout = NodeLayout.Compact;
            return new PersistentSuffixTree(storage, rootOffset, text, layout,
                builder.TransitionOffset, builder.JumpTableStart, builder.JumpTableEnd,
                builder.DeepestInternalNodeOffset, builder.LrsDepth);
        }
        catch
        {
            storage.Dispose();
            CleanupFile(filePath);
            throw;
        }
        finally
        {
            // Dispose and delete temp storage regardless of success/failure.
            // The builder already disposed owned stores, but for externally
            // provided storage (this case) we own the lifecycle.
            (childStorage as IDisposable)?.Dispose();
            CleanupFile(tempChildPath);
            (depthStorage as IDisposable)?.Dispose();
            CleanupFile(tempDepthPath);
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
            storage.Prefetch();
            return PersistentSuffixTree.Load(storage);
        }
        catch
        {
            storage.Dispose();
            throw;
        }
    }

    /// <summary>
    /// Estimates the storage capacity needed for a suffix tree over a text of
    /// the given length. The estimate covers header, nodes, child arrays, and
    /// serialized text. Overshooting is cheap (TrimToSize reclaims it);
    /// undershooting causes expensive remap/resize operations.
    /// </summary>
    internal static long EstimateCapacity(int textLength)
    {
        // SlimCompact layout: ≤2N+1 nodes × 24B + ≤2N+1 child entries × 8B
        //                     + text × 2B + header 88B ≈ 66N + 120.
        // Use 72 bytes/char to leave headroom and avoid any remap in practice.
        const int BytesPerChar = 72;
        const int MinCapacity = 65536;
        long estimate = (long)textLength * BytesPerChar + 256;
        return Math.Max(estimate, MinCapacity);
    }

    private static IStorageProvider CreateStorage(string? filePath, int textLength)
    {
        long capacity = EstimateCapacity(textLength);
        return !string.IsNullOrEmpty(filePath)
            ? new MappedFileStorageProvider(filePath, capacity)
            : new HeapStorageProvider((int)Math.Min(capacity, int.MaxValue));
    }

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
