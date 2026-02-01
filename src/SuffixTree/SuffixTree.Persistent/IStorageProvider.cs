using System;

namespace SuffixTree.Persistent;

/// <summary>
/// Abstraction for persistent storage access.
/// Allows interchangeable use of Memory-Mapped Files, Byte Arrays, or other low-level buffers.
/// </summary>
public interface IStorageProvider : IDisposable
{
    /// <summary>
    /// Gets the current size of the storage in bytes.
    /// </summary>
    long Size { get; }

    /// <summary>
    /// Ensures that the storage has at least the specified capacity.
    /// </summary>
    void EnsureCapacity(long capacity);

    /// <summary>
    /// Reads an integer from the specified offset.
    /// </summary>
    int ReadInt32(long offset);

    /// <summary>
    /// Writes an integer to the specified offset.
    /// </summary>
    void WriteInt32(long offset, int value);

    /// <summary>
    /// Reads a 64-bit integer from the specified offset.
    /// </summary>
    long ReadInt64(long offset);

    /// <summary>
    /// Writes a 64-bit integer to the specified offset.
    /// </summary>
    void WriteInt64(long offset, long value);

    /// <summary>
    /// Reads a character from the specified offset.
    /// </summary>
    char ReadChar(long offset);

    /// <summary>
    /// Writes a character to the specified offset.
    /// </summary>
    void WriteChar(long offset, char value);

    /// <summary>
    /// Allocates a block of memory of the specified size and returns its offset.
    /// </summary>
    long Allocate(int size);
}
