using System.Buffers.Binary;

namespace SuffixTree.Persistent;

/// <summary>
/// An in-memory implementation of IStorageProvider using a byte array.
/// Useful for testing, small trees, or as a reference implementation.
/// </summary>
public sealed class HeapStorageProvider : IStorageProvider
{
    private byte[] _buffer;
    private long _position;
    private int _disposed;

    /// <summary>Initializes a new <see cref="HeapStorageProvider"/> with the specified initial capacity.</summary>
    public HeapStorageProvider(int initialCapacity = 65536)
    {
        _buffer = new byte[initialCapacity];
        _position = 0;
    }

    /// <inheritdoc />
    public long Size => _position;

    private void ThrowIfDisposed()
    {
        ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);
    }

    private void CheckReadBounds(long offset, int size)
    {
        if (offset < 0 || offset + size > _position)
            throw new ArgumentOutOfRangeException(nameof(offset),
                $"Read at offset {offset} with size {size} exceeds logical size {_position}.");
    }

    private void CheckWriteBounds(long offset, int size)
    {
        if (offset < 0 || offset + size > _position)
            throw new ArgumentOutOfRangeException(nameof(offset),
                $"Write at offset {offset} with size {size} exceeds logical size {_position}.");
    }

    /// <inheritdoc />
    public void EnsureCapacity(long capacity)
    {
        ThrowIfDisposed();

        if (capacity < 0)
            throw new ArgumentOutOfRangeException(nameof(capacity), "Capacity must be non-negative.");

        if (_buffer.Length < capacity)
        {
            if (capacity > int.MaxValue)
                throw new InvalidOperationException("HeapStorageProvider cannot exceed 2GB limit of byte array.");

            long newSize = Math.Max(_buffer.Length * 2L, 64L);
            while (newSize < capacity) newSize *= 2;

            if (newSize > int.MaxValue)
                newSize = int.MaxValue;

            Array.Resize(ref _buffer, (int)newSize);
        }
    }

    /// <inheritdoc />
    public int ReadInt32(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 4);
        return BinaryPrimitives.ReadInt32LittleEndian(_buffer.AsSpan((int)offset, 4));
    }

    /// <inheritdoc />
    public void WriteInt32(long offset, int value)
    {
        ThrowIfDisposed();
        CheckWriteBounds(offset, 4);
        BinaryPrimitives.WriteInt32LittleEndian(_buffer.AsSpan((int)offset, 4), value);
    }

    /// <inheritdoc />
    public uint ReadUInt32(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 4);
        return BinaryPrimitives.ReadUInt32LittleEndian(_buffer.AsSpan((int)offset, 4));
    }

    /// <inheritdoc />
    public void WriteUInt32(long offset, uint value)
    {
        ThrowIfDisposed();
        CheckWriteBounds(offset, 4);
        BinaryPrimitives.WriteUInt32LittleEndian(_buffer.AsSpan((int)offset, 4), value);
    }

    /// <inheritdoc />
    public long ReadInt64(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 8);
        return BinaryPrimitives.ReadInt64LittleEndian(_buffer.AsSpan((int)offset, 8));
    }

    /// <inheritdoc />
    public void WriteInt64(long offset, long value)
    {
        ThrowIfDisposed();
        CheckWriteBounds(offset, 8);
        BinaryPrimitives.WriteInt64LittleEndian(_buffer.AsSpan((int)offset, 8), value);
    }

    /// <inheritdoc />
    public char ReadChar(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 2);
        return (char)BinaryPrimitives.ReadInt16LittleEndian(_buffer.AsSpan((int)offset, 2));
    }

    /// <inheritdoc />
    public void WriteChar(long offset, char value)
    {
        ThrowIfDisposed();
        CheckWriteBounds(offset, 2);
        BinaryPrimitives.WriteInt16LittleEndian(_buffer.AsSpan((int)offset, 2), (short)value);
    }

    /// <inheritdoc />
    public void ReadBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, count);
        Buffer.BlockCopy(_buffer, (int)offset, buffer, start, count);
    }

    /// <inheritdoc />
    public void WriteBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        CheckWriteBounds(offset, count);
        Buffer.BlockCopy(buffer, start, _buffer, (int)offset, count);
    }

    /// <inheritdoc />
    public long Allocate(int size)
    {
        ThrowIfDisposed();
        if (size < 0)
            throw new ArgumentOutOfRangeException(nameof(size), "Allocation size must be non-negative.");

        long newPosition = _position + size;
        EnsureCapacity(newPosition);
        long offset = _position;
        _position = newPosition;
        return offset;
    }

    /// <inheritdoc />
    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;
        GC.SuppressFinalize(this);
        _buffer = null!;
    }
}
