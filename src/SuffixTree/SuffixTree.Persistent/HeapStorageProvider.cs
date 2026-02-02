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

    public HeapStorageProvider(int initialCapacity = 65536)
    {
        _buffer = new byte[initialCapacity];
        _position = 0;
    }

    public long Size => _position;

    public void EnsureCapacity(long capacity)
    {
        if (_buffer.Length < capacity)
        {
            long newSize = _buffer.Length * 2L;
            while (newSize < capacity) newSize *= 2L;

            if (newSize > int.MaxValue)
                throw new InvalidOperationException("HeapStorageProvider cannot exceed 2GB limit of byte array.");

            Array.Resize(ref _buffer, (int)newSize);
        }
    }

    public int ReadInt32(long offset)
        => BinaryPrimitives.ReadInt32LittleEndian(_buffer.AsSpan((int)offset, 4));

    public void WriteInt32(long offset, int value)
        => BinaryPrimitives.WriteInt32LittleEndian(_buffer.AsSpan((int)offset, 4), value);

    public uint ReadUInt32(long offset)
        => BinaryPrimitives.ReadUInt32LittleEndian(_buffer.AsSpan((int)offset, 4));

    public void WriteUInt32(long offset, uint value)
        => BinaryPrimitives.WriteUInt32LittleEndian(_buffer.AsSpan((int)offset, 4), value);

    public long ReadInt64(long offset)
        => BinaryPrimitives.ReadInt64LittleEndian(_buffer.AsSpan((int)offset, 8));

    public void WriteInt64(long offset, long value)
        => BinaryPrimitives.WriteInt64LittleEndian(_buffer.AsSpan((int)offset, 8), value);

    public char ReadChar(long offset)
        => (char)BinaryPrimitives.ReadInt16LittleEndian(_buffer.AsSpan((int)offset, 2));

    public void WriteChar(long offset, char value)
        => BinaryPrimitives.WriteInt16LittleEndian(_buffer.AsSpan((int)offset, 2), (short)value);

    public long Allocate(int size)
    {
        long offset = _position;
        _position += size;
        EnsureCapacity(_position);
        return offset;
    }

    public void Dispose()
    {
        _buffer = null!;
    }
}
