using System;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Threading;

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
    private int _disposed;

    internal MemoryMappedViewAccessor Accessor
    {
        get
        {
            ThrowIfDisposed();
            return _accessor;
        }
    }

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
        _position = _readOnly ? _capacity : 0;
    }

    public long Size => _position;

    public void EnsureCapacity(long capacity)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot expand capacity in read-only mode.");
        if (capacity > _capacity)
        {
            var oldAccessor = _accessor;
            var oldMmf = _mmf;
            var oldCapacity = _capacity;

            _capacity = Math.Max(_capacity * 2, capacity);

            try
            {
                // Must close old mapping before resizing file on Windows
                oldAccessor.Dispose();
                oldMmf.Dispose();

                using (var fs = new FileStream(_filePath, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))
                {
                    fs.SetLength(_capacity);
                }

                _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, _capacity, MemoryMappedFileAccess.ReadWrite);
                _accessor = _mmf.CreateViewAccessor(0, _capacity, MemoryMappedFileAccess.ReadWrite);
            }
            catch
            {
                // Dispose any MMF handle that was successfully created before the error
                // to avoid native handle leak (R3)
                try { _mmf?.Dispose(); } catch { }

                // Try to restore previous mapping at old capacity
                _capacity = oldCapacity;
                try
                {
                    _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, oldCapacity, MemoryMappedFileAccess.ReadWrite);
                    _accessor = _mmf.CreateViewAccessor(0, oldCapacity, MemoryMappedFileAccess.ReadWrite);
                }
                catch
                {
                    // Cannot recover — mark as disposed to prevent zombie state
                    Volatile.Write(ref _disposed, 1);
                }
                throw;
            }
        }
    }

    public int ReadInt32(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 4);
        return _accessor.ReadInt32(offset);
    }

    public void WriteInt32(long offset, int value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public uint ReadUInt32(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 4);
        return _accessor.ReadUInt32(offset);
    }

    public void WriteUInt32(long offset, uint value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public long ReadInt64(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 8);
        return _accessor.ReadInt64(offset);
    }

    public void WriteInt64(long offset, long value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public char ReadChar(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 2);
        return _accessor.ReadChar(offset);
    }

    public void WriteChar(long offset, char value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.Write(offset, value);
    }

    public void ReadBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, count);
        int read = _accessor.ReadArray(offset, buffer, start, count);
        if (read != count)
            throw new InvalidOperationException(
                $"Partial read: requested {count} bytes at offset {offset}, got {read}. Storage may be corrupted.");
    }

    public void WriteBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        _accessor.WriteArray(offset, buffer, start, count);
    }

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

    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;
        _accessor?.Dispose();
        _mmf?.Dispose();
        _accessor = null!;
        _mmf = null!;
    }

    private void ThrowIfDisposed()
    {
        if (Volatile.Read(ref _disposed) != 0)
            throw new ObjectDisposedException(nameof(MappedFileStorageProvider));
    }

    private void CheckReadBounds(long offset, int size)
    {
        if (offset < 0 || offset + size > _position)
            throw new ArgumentOutOfRangeException(nameof(offset),
                $"Read at offset {offset} with size {size} exceeds logical size {_position}.");
    }

    /// <summary>
    /// Sets the current position (used when loading an existing tree).
    /// </summary>
    public void SetSize(long size)
    {
        ThrowIfDisposed();
        if (size < 0)
            throw new ArgumentOutOfRangeException(nameof(size), "Size must be non-negative.");
        _position = size;
        EnsureCapacity(_position);
    }

    /// <summary>
    /// Trims the backing file to the actual data size, reclaiming unused capacity.
    /// Must be called when writing is complete (e.g., after build).
    /// </summary>
    public void TrimToSize()
    {
        ThrowIfDisposed();
        if (_readOnly || _position >= _capacity) return;

        var oldAccessor = _accessor;
        var oldMmf = _mmf;
        var oldCapacity = _capacity;

        _capacity = _position;

        try
        {
            oldAccessor.Dispose();
            oldMmf.Dispose();

            using (var fs = new FileStream(_filePath, FileMode.Open, FileAccess.ReadWrite, FileShare.None))
            {
                fs.SetLength(_capacity);
            }

            _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, _capacity, MemoryMappedFileAccess.ReadWrite);
            _accessor = _mmf.CreateViewAccessor(0, _capacity, MemoryMappedFileAccess.ReadWrite);
        }
        catch
        {
            // Dispose any MMF handle that was successfully created before the error
            // to avoid native handle leak (R4)
            try { _mmf?.Dispose(); } catch { }

            _capacity = oldCapacity;
            try
            {
                _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, oldCapacity, MemoryMappedFileAccess.ReadWrite);
                _accessor = _mmf.CreateViewAccessor(0, oldCapacity, MemoryMappedFileAccess.ReadWrite);
            }
            catch
            {
                Volatile.Write(ref _disposed, 1);
            }
            throw;
        }
    }
}
