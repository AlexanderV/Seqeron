using System.Diagnostics;
using System.IO.MemoryMappedFiles;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace SuffixTree.Persistent;

/// <summary>
/// A persistent implementation of IStorageProvider using Memory-Mapped Files.
/// This allows the suffix tree to reside on disk and be mapped into the process's address space.
/// Uses unsafe direct pointer access for minimal-overhead reads and writes.
/// </summary>
public sealed unsafe partial class MappedFileStorageProvider : IStorageProvider
{
    private readonly string _filePath;
    private MemoryMappedFile _mmf;
    private MemoryMappedViewAccessor _accessor;
    private byte* _ptr;
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

    /// <summary>Initializes a new <see cref="MappedFileStorageProvider"/> backed by a memory-mapped file.</summary>
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
        AcquirePointer();
        _position = _readOnly ? _capacity : 0;
    }

    /// <inheritdoc />
    public long Size => _position;

    /// <inheritdoc />
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
                // Release pointer and close old mapping before resizing file on Windows
                ReleasePointer();
                oldAccessor.Dispose();
                oldMmf.Dispose();

                using (var fs = new FileStream(_filePath, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))
                {
                    fs.SetLength(_capacity);
                }

                _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, _capacity, MemoryMappedFileAccess.ReadWrite);
                _accessor = _mmf.CreateViewAccessor(0, _capacity, MemoryMappedFileAccess.ReadWrite);
                AcquirePointer();
            }
            catch (IOException) { RecoverMapping(oldCapacity); throw; }
            catch (UnauthorizedAccessException) { RecoverMapping(oldCapacity); throw; }
            catch (OutOfMemoryException) { RecoverMapping(oldCapacity); throw; }
        }
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int ReadInt32(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 4);
        return *(int*)(_ptr + offset);
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void WriteInt32(long offset, int value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        CheckWriteBounds(offset, 4);
        *(int*)(_ptr + offset) = value;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public uint ReadUInt32(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 4);
        return *(uint*)(_ptr + offset);
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void WriteUInt32(long offset, uint value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        CheckWriteBounds(offset, 4);
        *(uint*)(_ptr + offset) = value;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long ReadInt64(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 8);
        return *(long*)(_ptr + offset);
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void WriteInt64(long offset, long value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        CheckWriteBounds(offset, 8);
        *(long*)(_ptr + offset) = value;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public char ReadChar(long offset)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, 2);
        return *(char*)(_ptr + offset);
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void WriteChar(long offset, char value)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        CheckWriteBounds(offset, 2);
        *(char*)(_ptr + offset) = value;
    }

    /// <inheritdoc />
    public void ReadBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        CheckReadBounds(offset, count);
        fixed (byte* dst = &buffer[start])
        {
            Buffer.MemoryCopy(_ptr + offset, dst, count, count);
        }
    }

    /// <inheritdoc />
    public void WriteBytes(long offset, byte[] buffer, int start, int count)
    {
        ThrowIfDisposed();
        if (_readOnly) throw new InvalidOperationException("Cannot write in read-only mode.");
        CheckWriteBounds(offset, count);
        fixed (byte* src = &buffer[start])
        {
            Buffer.MemoryCopy(src, _ptr + offset, count, count);
        }
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
        ReleasePointer();
        _accessor?.Dispose();
        _mmf?.Dispose();
        _accessor = null!;
        _mmf = null!;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ThrowIfDisposed()
    {
        // Plain read: Dispose() uses Interlocked.CompareExchange which
        // already publishes the write. No Volatile fence needed on the
        // read side — the builder is single-threaded.
        if (_disposed != 0)
            ObjectDisposedException.ThrowIf(true, this);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void CheckReadBounds(long offset, int size)
    {
        if (offset < 0 || offset + size > _position)
            throw new ArgumentOutOfRangeException(nameof(offset),
                $"Read at offset {offset} with size {size} exceeds logical size {_position}.");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void CheckWriteBounds(long offset, int size)
    {
        if (offset < 0 || offset + size > _position)
            throw new ArgumentOutOfRangeException(nameof(offset),
                $"Write at offset {offset} with size {size} exceeds logical size {_position}.");
    }

    // ──────────────── Unchecked fast-path for builder hot loop ────────────────
    // These skip ThrowIfDisposed + CheckBounds. The builder guarantees all
    // offsets come from Allocate() and disposal never races with build.

    /// <summary>Raw pointer — use at your own risk.</summary>
    internal byte* RawPointer
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _ptr;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint ReadUInt32Unchecked(long offset)
    {
        Debug.Assert(offset >= 0 && offset + 4 <= _capacity, "Unchecked read out of bounds.");
        return *(uint*)(_ptr + offset);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void WriteUInt32Unchecked(long offset, uint value)
    {
        Debug.Assert(offset >= 0 && offset + 4 <= _capacity, "Unchecked write out of bounds.");
        *(uint*)(_ptr + offset) = value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal int ReadInt32Unchecked(long offset)
    {
        Debug.Assert(offset >= 0 && offset + 4 <= _capacity, "Unchecked read out of bounds.");
        return *(int*)(_ptr + offset);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void WriteInt32Unchecked(long offset, int value)
    {
        Debug.Assert(offset >= 0 && offset + 4 <= _capacity, "Unchecked write out of bounds.");
        *(int*)(_ptr + offset) = value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal long ReadInt64Unchecked(long offset)
    {
        Debug.Assert(offset >= 0 && offset + 8 <= _capacity, "Unchecked read out of bounds.");
        return *(long*)(_ptr + offset);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void WriteInt64Unchecked(long offset, long value)
    {
        Debug.Assert(offset >= 0 && offset + 8 <= _capacity, "Unchecked write out of bounds.");
        *(long*)(_ptr + offset) = value;
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
            ReleasePointer();
            oldAccessor.Dispose();
            oldMmf.Dispose();

            using (var fs = new FileStream(_filePath, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))
            {
                fs.SetLength(_capacity);
            }

            _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, _capacity, MemoryMappedFileAccess.ReadWrite);
            _accessor = _mmf.CreateViewAccessor(0, _capacity, MemoryMappedFileAccess.ReadWrite);
            AcquirePointer();
        }
        catch (IOException) { RecoverMapping(oldCapacity); throw; }
        catch (UnauthorizedAccessException) { RecoverMapping(oldCapacity); throw; }
        catch (OutOfMemoryException) { RecoverMapping(oldCapacity); throw; }
    }

    /// <summary>
    /// Acquires a raw byte pointer to the memory-mapped view for direct access.
    /// </summary>
    private void AcquirePointer()
    {
        byte* basePtr = null;
        _accessor.SafeMemoryMappedViewHandle.AcquirePointer(ref basePtr);
        _ptr = basePtr + _accessor.PointerOffset;
    }

    /// <summary>
    /// Releases the previously acquired pointer. Safe to call when no pointer is held.
    /// </summary>
    private void ReleasePointer()
    {
        if (_ptr != null)
        {
            _ptr = null;
            _accessor?.SafeMemoryMappedViewHandle.ReleasePointer();
        }
    }

    /// <summary>
    /// Attempts to restore the previous memory-mapped file mapping after a failed resize.
    /// Disposes any partially-created handles, then tries to re-map at the old capacity.
    /// If even restoration fails, marks the provider as disposed to prevent zombie state.
    /// </summary>
    private void RecoverMapping(long oldCapacity)
    {
        // Dispose any MMF handle that was partially created before the error
        // to avoid native handle leak. Dispose must not throw — suppress if it does.
#pragma warning disable CA1031 // Safety Dispose: must not throw during recovery
        try { _mmf?.Dispose(); } catch (Exception) { }
#pragma warning restore CA1031

        _capacity = oldCapacity;
        try
        {
            _mmf = MemoryMappedFile.CreateFromFile(_filePath, FileMode.Open, null, oldCapacity, MemoryMappedFileAccess.ReadWrite);
            _accessor = _mmf.CreateViewAccessor(0, oldCapacity, MemoryMappedFileAccess.ReadWrite);
            AcquirePointer();
        }
        catch (IOException)
        {
            // Cannot recover — mark as disposed to prevent zombie state
            Volatile.Write(ref _disposed, 1);
        }
        catch (UnauthorizedAccessException)
        {
            Volatile.Write(ref _disposed, 1);
        }
    }

    // ──────────────── Prefetch ────────────────

    /// <summary>
    /// Hints the OS to prefetch the mapped region into physical memory.
    /// Useful after loading a tree to avoid page faults on first access.
    /// This is advisory and fail-safe — errors are silently ignored.
    /// </summary>
    public void Prefetch()
    {
        ThrowIfDisposed();
        if (_ptr == null || _position <= 0) return;
        PrefetchRegion(_ptr, _position);
    }

    /// <summary>
    /// Hints the OS to prefetch the entire pre-allocated capacity before a build.
    /// Uses <c>PrefetchVirtualMemory</c> (Windows) / <c>posix_madvise</c> (Linux/macOS)
    /// to asynchronously start paging-in / zero-filling the mapped region.
    /// Advisory and fail-safe — returns immediately, errors are silently ignored.
    /// </summary>
    internal void PrefetchForBuild()
    {
        if (_ptr == null || _capacity <= 0) return;
        PrefetchRegion(_ptr, _capacity);
    }

#pragma warning disable CA1031 // Prefetch is advisory — must not throw
    private static void PrefetchRegion(byte* address, long length)
    {
        try
        {
            if (OperatingSystem.IsWindows())
                PrefetchWindows(address, length);
            else if (OperatingSystem.IsLinux() || OperatingSystem.IsMacOS())
                PrefetchPosix(address, length);
        }
        catch (Exception)
        {
            // Prefetch is advisory — ignore all errors
        }
    }
#pragma warning restore CA1031

    // ── Windows: PrefetchVirtualMemory ──

    private static void PrefetchWindows(byte* address, long length)
    {
        var entry = new WIN32_MEMORY_RANGE_ENTRY
        {
            VirtualAddress = address,
            NumberOfBytes = (nuint)length
        };
        PrefetchVirtualMemory(GetCurrentProcess(), 1, &entry, 0);
    }

    [StructLayout(LayoutKind.Sequential)]
    private struct WIN32_MEMORY_RANGE_ENTRY
    {
        public void* VirtualAddress;
        public nuint NumberOfBytes;
    }

    [LibraryImport("kernel32")]
    [DefaultDllImportSearchPaths(DllImportSearchPath.System32)]
    private static partial nint GetCurrentProcess();

    [LibraryImport("kernel32", SetLastError = true)]
    [DefaultDllImportSearchPaths(DllImportSearchPath.System32)]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static partial bool PrefetchVirtualMemory(
        nint hProcess,
        nuint numberOfEntries,
        WIN32_MEMORY_RANGE_ENTRY* virtualAddresses,
        uint flags);

    // ── Linux/macOS: posix_madvise ──

    private static void PrefetchPosix(byte* address, long length)
    {
        const int POSIX_MADV_WILLNEED = 3;
        _ = PosixMadvise(address, (nuint)length, POSIX_MADV_WILLNEED);
    }

    [LibraryImport("libc", EntryPoint = "posix_madvise")]
    [DefaultDllImportSearchPaths(DllImportSearchPath.System32)]
    private static partial int PosixMadvise(void* addr, nuint length, int advice);
}
