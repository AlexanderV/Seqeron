using System.IO.MemoryMappedFiles;
using System.Runtime.InteropServices;

namespace SuffixTree.EmProto.Text;

/// <summary>Owner of a block of memory addressed by a raw pointer: native heap memory or a memory-mapped file.</summary>
internal unsafe interface IMemoryBlock : IDisposable
{
    /// <summary>First byte of the block.</summary>
    byte* Pointer { get; }

    /// <summary>Size of the block in bytes.</summary>
    long Bytes { get; }
}

/// <summary>A block from <see cref="NativeMemory"/> (BCL API, no P/Invoke), zero-initialised, freed on dispose.</summary>
internal sealed unsafe class NativeBlock : IMemoryBlock
{
    private nint _p;

    internal NativeBlock(long bytes)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(bytes);
        Bytes = bytes;
        _p = (nint)NativeMemory.AllocZeroed((nuint)Math.Max(bytes, 1));
    }

    ~NativeBlock() => Free();

    public byte* Pointer => (byte*)_p;

    public long Bytes { get; }

    public void Dispose()
    {
        Free();
        GC.SuppressFinalize(this);
    }

    private void Free()
    {
        nint p = Interlocked.Exchange(ref _p, 0);
        if (p != 0) NativeMemory.Free((void*)p);
    }
}

/// <summary>
/// A whole file mapped read-only (or read-write) into the address space. Pages are file-backed page cache, not anonymous memory: they count in
/// <c>RssFile</c> and never in <c>RssAnon</c> (ADR-004 section 6a). The view covers the full length; the pointer is valid until disposal.
/// </summary>
internal sealed unsafe class MappedFile : IMemoryBlock
{
    private readonly MemoryMappedFile _file;
    private readonly MemoryMappedViewAccessor _view;
    private byte* _base;
    private bool _disposed;

    internal MappedFile(string path, bool writable)
    {
        Path = path;
        Bytes = new FileInfo(path).Length;
        if (Bytes == 0) throw new ArgumentException("cannot map an empty file", nameof(path));
        MemoryMappedFileAccess access = writable ? MemoryMappedFileAccess.ReadWrite : MemoryMappedFileAccess.Read;
        _file = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, access);
        _view = _file.CreateViewAccessor(0, 0, access);
        byte* p = null;
        _view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        _base = p + _view.PointerOffset;
    }

    internal string Path { get; }

    public byte* Pointer => _base;

    public long Bytes { get; }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        _view.SafeMemoryMappedViewHandle.ReleasePointer();
        _base = null;
        _view.Dispose();
        _file.Dispose();
    }
}
