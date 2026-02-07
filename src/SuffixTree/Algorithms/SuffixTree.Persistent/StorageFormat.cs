namespace SuffixTree.Persistent;

/// <summary>
/// Selects the binary storage format for persistent suffix trees.
/// Internal: format is chosen automatically based on text length.
/// </summary>
internal enum StorageFormat
{
    /// <summary>
    /// 32-bit offsets. Node = 28 bytes, child entry = 8 bytes.
    /// ~30% smaller files, better cache locality.
    /// Maximum file size: 4 GB (~58 M characters).
    /// </summary>
    Compact = 4,

    /// <summary>
    /// 64-bit offsets. Node = 40 bytes, child entry = 12 bytes.
    /// No practical size limit.
    /// </summary>
    Large = 3
}
