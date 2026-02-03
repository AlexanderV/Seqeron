using System;
using System.IO;
using System.IO.MemoryMappedFiles;

namespace SuffixTree.Persistent;

/// <summary>
/// A persistent implementation of IStorageProvider using Memory-Mapped Files.
/// This allows the suffix tree to reside on disk and be mapped into the process's address space.
/// </summary>
public class MappedFileStorageProvider : IStorageProvider
{
    private readonly string _filePath;
    private MemoryMappedFile _mmf;
    private MemoryMappedViewAccessor _accessor;
    private readonly bool _readOnly;
    private long _capacity;
    private long _position;

    internal MemoryMappedViewAccessor Accessor => _accessor;

    public MappedFileStorageProvider(string filePath, long initialCapacity = 65536, bool readOnly = false)
    {
        _filePath = filePath;
        _readOnly = readOnly;
        _capacity = initialCapacity;

        // Ensure directory exists
        var dir = Path.GetDirectoryName(filePath);
        if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir))
            Directory.CreateDirectory(dir);

        // Open or create the file with the initial capacity
        var access = _readOnly ? FileAccess.Read : FileAccess.ReadWrite;
        var share = FileShare.ReadWrite;
        var mmfAccess = _readOnly ? MemoryMappedFileAccess.Read : MemoryMappedFileAccess.ReadWrite;

        using (var fs = new FileStream(_filePath, _readOnly ? FileMode.Open : FileMode.OpenOrCreate, access, share))
        {
            if (!_readOnly && fs.Length < _capacity)
                fs.SetLength(_capacity);
            else
                _capacity = fs.Length;
        }

        _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, _capacity, mmfAccess);
        _accessor = _mmf.CreateViewAccessor(0, _capacity, mmfAccess);
        _position = 0; // In a real "load" scenario, we'd need to store the current position/size in the file header.
    }

    public long Size => _position;

    public void EnsureCapacity(long capacity)
    {
        if (_readOnly) throw new InvalidOperationException("Cannot expand capacity in read-only mode.");
        if (capacity > _capacity)
        {
            // Re-mapping required
            _accessor.Dispose();
            _mmf.Dispose();

            _capacity = Math.Max(_capacity * 2, capacity);

            using (var fs = new FileStream(_filePath, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))
            {
                fs.SetLength(_capacity);
            }

            _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, _capacity, MemoryMappedFileAccess.ReadWrite);
            _accessor = _mmf.CreateViewAccessor(0, _capacity, MemoryMappedFileAccess.ReadWrite);
        }
    }

    public int ReadInt32(long offset) => _accessor.ReadInt32(offset);

    public void WriteInt32(long offset, int value)
    {
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public uint ReadUInt32(long offset) => _accessor.ReadUInt32(offset);

    public void WriteUInt32(long offset, uint value)
    {
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public long ReadInt64(long offset) => _accessor.ReadInt64(offset);

    public void WriteInt64(long offset, long value)
    {
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public char ReadChar(long offset) => _accessor.ReadChar(offset);

    public void WriteChar(long offset, char value)
    {
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public long Allocate(int size)
    {
        long offset = _position;
        _position += size;
        EnsureCapacity(_position);
        return offset;
    }

    public void Dispose()
    {
        _accessor?.Dispose();
        _mmf?.Dispose();
    }

    /// <summary>
    /// Sets the current position (used when loading an existing tree).
    /// </summary>
    public void SetSize(long size)
    {
        _position = size;
        EnsureCapacity(_position);
    }
}
