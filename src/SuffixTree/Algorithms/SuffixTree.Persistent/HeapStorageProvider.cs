using System;
using System.Buffers.Binary;

namespace SuffixTree.Persistent;

/// <summary>
/// An in-memory implementation of IStorageProvider using a byte array.
/// Useful for testing, small trees, or as a reference implementation.
/// </summary>
public class HeapStorageProvider : IStorageProvider
{
    private byte[] _buffer;
    private long _position;
    private volatile bool _disposed;

    public HeapStorageProvider(int initialCapacity = 65536)
    {
        _buffer = new byte[initialCapacity];
        _position = 0;
    }

    public long Size => _position;

    private void ThrowIfDisposed()
    {
        if (_disposed)
            throw new ObjectDisposedException(nameof(HeapStorageProvider));
    }

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

    public int ReadInt32(long offset)
    {
        ThrowIfDisposed();
        return BinaryPrimitives.ReadInt32LittleEndian(_buffer.AsSpan((int)offset, 4));
    }

    public void WriteInt32(long offset, int value)
    {
        ThrowIfDisposed();
        BinaryPrimitives.WriteInt32LittleEndian(_buffer.AsSpan((int)offset, 4), value);
    }

    public uint ReadUInt32(long offset)
    {
        ThrowIfDisposed();
        return BinaryPrimitives.ReadUInt32LittleEndian(_buffer.AsSpan((int)offset, 4));
    }

    public void WriteUInt32(long offset, uint value)
    {
        ThrowIfDisposed();
        BinaryPrimitives.WriteUInt32LittleEndian(_buffer.AsSpan((int)offset, 4), value);
    }

    public long ReadInt64(long offset)
    {
        ThrowIfDisposed();
        return BinaryPrimitives.ReadInt64LittleEndian(_buffer.AsSpan((int)offset, 8));
    }

    public void WriteInt64(long offset, long value)
    {
        ThrowIfDisposed();
        BinaryPrimitives.WriteInt64LittleEndian(_buffer.AsSpan((int)offset, 8), value);
    }

    public char ReadChar(long offset)
    {
        ThrowIfDisposed();
        return (char)BinaryPrimitives.ReadInt16LittleEndian(_buffer.AsSpan((int)offset, 2));
    }

    public void WriteChar(long offset, char value)
    {
        ThrowIfDisposed();
        BinaryPrimitives.WriteInt16LittleEndian(_buffer.AsSpan((int)offset, 2), (short)value);
    }

    public void ReadBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        Buffer.BlockCopy(_buffer, (int)offset, buffer, start, count);
    }

    public void WriteBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        Buffer.BlockCopy(buffer, start, _buffer, (int)offset, count);
    }

    public long Allocate(int size)
    {
        ThrowIfDisposed();
        if (size < 0)
            throw new ArgumentOutOfRangeException(nameof(size), "Allocation size must be non-negative.");

        long offset = _position;
        _position += size;
        EnsureCapacity(_position);
        return offset;
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        _buffer = null!;
    }
}
