using System.IO.MemoryMappedFiles;
using Microsoft.Win32.SafeHandles;

namespace SuffixTree.Esa.External;

/// <summary>How the blocks of a new output file are reserved (ADR-004 section 5a, 1b.7b, has the measurement behind the default).</summary>
internal enum ReservationMode
{
    /// <summary>The file is created at its length without blocks (<c>ftruncate</c>): the layout E9 measured. A full disk then shows as a signal on the page fault of a mapped view, not as an exception.</summary>
    Sparse,

    /// <summary>The blocks are reserved at creation where the platform can (<c>fallocate</c> through <see cref="FileStreamOptions.PreallocationSize"/>): a full disk is an <see cref="IOException"/> at creation. Where it cannot, the free space of the volume is checked explicitly and the file is created sparse.</summary>
    Preallocate,

    /// <summary>The measured default, <see cref="MappedOutputFile.AutoMode"/>.</summary>
    Auto,
}

/// <summary>What the creation of an output file did about the disk space.</summary>
internal enum ReservationOutcome
{
    /// <summary>Created sparse on request: nothing was reserved or checked.</summary>
    Sparse,

    /// <summary>The platform reserved the blocks: a full disk was ruled out at creation.</summary>
    Preallocated,

    /// <summary>The platform cannot reserve; the free space was read and was enough, the file is sparse (advisory: another process can still fill the disk).</summary>
    FreeSpaceChecked,

    /// <summary>The platform cannot reserve and cannot tell the free space; the file is sparse and nothing was checked.</summary>
    Unchecked,

    /// <summary>The range that was added to an existing file (<see cref="MappedOutputFile.Grow"/>) was reserved by writing zeros into it, because the BCL reserves blocks only for a file that it creates: a full disk is an <see cref="IOException"/> from the write, the file is as it was, and a mapped write into the range never needs a new block.</summary>
    Filled,

    /// <summary>The file was not created by this call: the container of an earlier, interrupted run of the same build was opened again (<see cref="MappedOutputFile.OpenExisting"/>, work package 1b.9e). Its blocks were reserved when it was created, and the range that its growth added was reserved by <see cref="Filled"/> when that happened.</summary>
    Reopened,
}

/// <summary>The platform side of the reservation, injectable so that the fallback is tested where the real one does not exist (macOS, Windows).</summary>
internal interface IReservationSource
{
    /// <summary>True where <see cref="CreatePreallocated"/> reserves blocks.</summary>
    bool CanPreallocate { get; }

    /// <summary>Creates (truncates) the file with the blocks for <paramref name="length"/> bytes reserved and returns it open for reading and writing; its length may still be 0. A full disk is an <see cref="IOException"/>, and no file is left.</summary>
    FileStream CreatePreallocated(string path, long length);

    /// <summary>Bytes available to this process on the volume of <paramref name="directory"/>, or -1 when unknown. Never throws.</summary>
    long AvailableFreeBytes(string directory);

    /// <summary>
    /// Reserves the blocks of the range <c>[offset, offset + length)</c> of a file that already exists and has been extended over it (<see cref="MappedOutputFile.Grow"/>): the BCL can ask the file system for blocks only when
    /// it creates a file, so the default writes zeros into the range, which is the one way in managed code that works on every platform and every file system. A full disk is an <see cref="IOException"/>; the caller gives the
    /// length back. A source that can do better (or one that fails on purpose) overrides it.
    /// </summary>
    void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters) => ZeroFill.Write(file, offset, length, counters);
}

/// <summary>What <see cref="ZeroFill.ClearOutside"/> found and did.</summary>
/// <param name="ScannedBytes">Bytes that were read to be compared with zero (everything outside the ranges that were kept).</param>
/// <param name="ClearedBytes">Bytes that were written as zeros: the units of <see cref="ZeroFill.UnitBytes"/> that held a non-zero byte.</param>
internal readonly record struct ClearResult(long ScannedBytes, long ClearedBytes);

/// <summary>The managed reservation of a range of an existing file: zeros written through <see cref="RandomAccess"/> in chunks of <see cref="ChunkBytes"/>; and the normalisation of a reopened file (<see cref="ClearOutside"/>).</summary>
internal static class ZeroFill
{
    /// <summary>Bytes of one write.</summary>
    internal const int ChunkBytes = 1 << 20;

    /// <summary>The unit in which <see cref="ClearOutside"/> decides whether to write: a unit that holds a non-zero byte is written as zeros, a unit that is all zero is not touched (so its page is not dirtied and its block is not read from or written to the device twice).</summary>
    internal const int UnitBytes = 4096;

    private static readonly byte[] s_zeros = new byte[ChunkBytes];

    /// <summary>Writes <paramref name="length"/> zero bytes at <paramref name="offset"/>; a full disk or a quota is an <see cref="IOException"/>, and every byte is counted in <paramref name="counters"/> (<see cref="IoCounters.ReservedWriteBytes"/>) as it is written.</summary>
    internal static void Write(SafeFileHandle file, long offset, long length, IoCounters? counters)
    {
        ArgumentNullException.ThrowIfNull(file);
        ArgumentOutOfRangeException.ThrowIfNegative(offset);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        long end = offset + length;
        for (long at = offset; at < end; at += ChunkBytes)
        {
            int take = (int)Math.Min(ChunkBytes, end - at);
            RandomAccess.Write(file, s_zeros.AsSpan(0, take), at);
            counters?.AddReserved(take);
        }
    }

    /// <summary>
    /// Makes every byte of <c>[0, <paramref name="length"/>)</c> that lies outside the <paramref name="keep"/> ranges zero (work package 1b.9e.1): a new file has nothing but zeros where no writer has written yet, and the writers of the build
    /// rely on it (the child writer sets only the slots that the scan names, nothing writes the padding between the sections, the room of BUILD_INFO behind its text or the header page before PH7); a file that an earlier run left, or that
    /// another process, a failing disk or a power loss changed, has to be given that property again before the build goes on. The bytes are read in chunks of <see cref="ChunkBytes"/> and compared with zero, and only the units of
    /// <see cref="UnitBytes"/> that hold a non-zero byte are written (as zeros, in place, so the blocks that the file has stay the blocks that it has: nothing is allocated or freed, the reservation of the file is what it was). The ranges
    /// may be given in any order, may overlap and may reach beyond <paramref name="length"/>. A file that is shorter than <paramref name="length"/> is an <see cref="IOException"/>.
    /// </summary>
    internal static ClearResult ClearOutside(SafeFileHandle file, long length, IReadOnlyList<(long Offset, long Length)> keep)
    {
        ArgumentNullException.ThrowIfNull(file);
        ArgumentNullException.ThrowIfNull(keep);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        long scanned = 0, cleared = 0;
        byte[]? buffer = null;
        long at = 0;
        foreach ((long offset, long len) in keep.OrderBy(k => k.Offset).Append((length, 0)))
        {
            long gapEnd = Math.Min(offset, length);
            if (gapEnd > at)
            {
                buffer ??= new byte[ChunkBytes];
                ClearSpan(file, at, gapEnd, buffer, ref scanned, ref cleared);
            }
            at = Math.Max(at, offset + len);
        }
        return new ClearResult(scanned, cleared);
    }

    private static void ClearSpan(SafeFileHandle file, long start, long end, byte[] buffer, ref long scanned, ref long cleared)
    {
        for (long chunk = start; chunk < end; chunk += ChunkBytes)
        {
            int take = (int)Math.Min(ChunkBytes, end - chunk);
            for (int got = 0; got < take;)
            {
                int n = RandomAccess.Read(file, buffer.AsSpan(got, take - got), chunk + got);
                if (n <= 0) throw new IOException($"the file ends at {chunk + got}, before the {end} bytes that its length promised");
                got += n;
            }
            scanned += take;
            int dirtyFrom = -1;
            for (int unit = 0; unit < take; unit += UnitBytes)
            {
                if (buffer.AsSpan(unit, Math.Min(UnitBytes, take - unit)).IndexOfAnyExcept((byte)0) >= 0)
                {
                    if (dirtyFrom < 0) dirtyFrom = unit;
                }
                else if (dirtyFrom >= 0)
                {
                    cleared += WriteZeros(file, chunk + dirtyFrom, unit - dirtyFrom);
                    dirtyFrom = -1;
                }
            }
            if (dirtyFrom >= 0) cleared += WriteZeros(file, chunk + dirtyFrom, take - dirtyFrom);
        }
    }

    private static int WriteZeros(SafeFileHandle file, long offset, int length)
    {
        RandomAccess.Write(file, s_zeros.AsSpan(0, length), offset);
        return length;
    }
}

/// <summary>
/// The real source: on Linux the BCL reserves the blocks with <c>fallocate</c> (<see cref="FileStreamOptions.PreallocationSize"/>, no P/Invoke; it raises "the disk was full" as an <see cref="IOException"/> and removes the file);
/// the free space is <see cref="DriveInfo.AvailableFreeSpace"/> of the directory. Elsewhere <see cref="CanPreallocate"/> is false (the BCL has an equivalent on macOS and Windows that was not run here), so the explicit
/// fallback of <see cref="MappedOutputFile"/> applies.
/// </summary>
internal sealed class SystemReservationSource : IReservationSource
{
    internal static SystemReservationSource Instance { get; } = new();

    public bool CanPreallocate => OperatingSystem.IsLinux();

    public FileStream CreatePreallocated(string path, long length)
        => new(path, new FileStreamOptions { Mode = FileMode.Create, Access = FileAccess.ReadWrite, Share = FileShare.ReadWrite, PreallocationSize = length });

    public long AvailableFreeBytes(string directory)
    {
        try
        {
            return new DriveInfo(directory).AvailableFreeSpace;
        }
        catch (Exception e) when (e is ArgumentException or IOException or UnauthorizedAccessException or NotSupportedException)
        {
            return -1;
        }
    }
}

/// <summary>Options of a new output file.</summary>
internal sealed class MappedFileOptions
{
    /// <summary>How the blocks are reserved (default <see cref="ReservationMode.Auto"/>).</summary>
    internal ReservationMode Reservation { get; init; } = ReservationMode.Auto;

    /// <summary>The platform side (default <see cref="SystemReservationSource.Instance"/>).</summary>
    internal IReservationSource? Source { get; init; }

    /// <summary>Test hook: how <see cref="MappedOutputFile.Grow"/> maps the file again once it has its new length (null = the system's read-write mapping of the whole file); a hook that throws is the failed mapping that Grow must undo.</summary>
    internal Func<string, MemoryMappedFile>? Mapper { get; init; }
}

/// <summary>
/// A file of known final length that several <see cref="MappedWindowWriter"/>s fill through bounded windows of one shared mapping (the categorical rule of the owner, 2026-11-04: stream into memory-mapped files):
/// the file is created at its final length, each writer owns a range of it (disjoint ranges, any order, any thread, offsets above 2^32 bytes) and a window of its own, and the file is completed when every writer is
/// disposed (<see cref="Complete"/>: <c>fsync</c>). Creation reserves the disk space according to <see cref="ReservationMode"/>; a failed creation leaves no file at the path.
/// </summary>
internal sealed class MappedOutputFile : IDisposable
{
    /// <summary>What <see cref="ReservationMode.Auto"/> means: the choice measured in work package 1b.7b (ADR-004 section 5a).</summary>
    internal const ReservationMode AutoMode = ReservationMode.Preallocate;

    private MemoryMappedFile? _map;
    private readonly object _gate = new();
    private int _writers;
    private bool _disposed;

    private MappedOutputFile(string path, long length, ReservationOutcome outcome, MemoryMappedFile? map)
    {
        Path = path;
        Length = length;
        Outcome = outcome;
        _map = map;
    }

    /// <summary>The path the file was created at.</summary>
    internal string Path { get; }

    /// <summary>Length of the file in bytes (it grows with <see cref="Grow"/>).</summary>
    internal long Length { get; private set; }

    /// <summary>What the creation did about the disk space.</summary>
    internal ReservationOutcome Outcome { get; }

    /// <summary>What <see cref="OpenExisting"/> did to the bytes outside the ranges it was asked to keep (zero when it was not asked).</summary>
    internal ClearResult Cleared { get; private init; }

    /// <summary>Creates (truncates) the file at <paramref name="path"/> with the length <paramref name="length"/>; on any failure no file is left at the path.</summary>
    internal static MappedOutputFile Create(string path, long length, MappedFileOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(path);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        options ??= new MappedFileOptions();
        ReservationMode mode = options.Reservation == ReservationMode.Auto ? AutoMode : options.Reservation;
        IReservationSource source = options.Source ?? SystemReservationSource.Instance;
        try
        {
            ReservationOutcome outcome = Allocate(path, length, mode, source);
            MemoryMappedFile? map = length > 0 ? MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite) : null;
            return new MappedOutputFile(path, length, outcome, map);
        }
        catch
        {
            DeleteQuietly(path);
            throw;
        }
    }

    /// <summary>
    /// Opens a file that an earlier run created (an interrupted build whose stages are to be resumed, work package 1b.9e) for writing through bounded windows, at the length it has now: nothing is created, truncated or
    /// reserved, so what the file already holds stays. <paramref name="truncateTo"/> cuts a longer file back first (the part behind the reserved length that an interrupted growth added: <see cref="Grow"/> reserves it again, whole,
    /// when PH4 asks for it, so a growth that was cut off half way through its zeros cannot leave a window of blocks that were never reserved). A file shorter than <paramref name="truncateTo"/> is left as it is.
    /// With <paramref name="keep"/> (work package 1b.9e.1) every byte of the file that lies outside those ranges is then made zero (<see cref="ZeroFill.ClearOutside"/>): the file that an earlier run left holds the verified output of the
    /// stages that are kept and, anywhere else, whatever the killed run or the world put there; a build goes on from a file that has nothing but zeros outside the ranges it keeps, as it does from a new one. Null leaves the bytes alone.
    /// </summary>
    /// <exception cref="FileNotFoundException">There is no such file.</exception>
    internal static MappedOutputFile OpenExisting(string path, long? truncateTo = null, IReadOnlyList<(long Offset, long Length)>? keep = null)
    {
        ArgumentNullException.ThrowIfNull(path);
        if (truncateTo is { } negative) ArgumentOutOfRangeException.ThrowIfNegative(negative);
        ClearResult cleared = default;
        if (truncateTo is not null || keep is not null)
        {
            using var handle = File.OpenHandle(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
            if (truncateTo is { } cut && RandomAccess.GetLength(handle) > cut) RandomAccess.SetLength(handle, cut);
            if (keep is not null) cleared = ZeroFill.ClearOutside(handle, RandomAccess.GetLength(handle), keep);
        }
        long length = new FileInfo(path).Length; // a file that is gone is a FileNotFoundException from here (or from the handle above)
        MemoryMappedFile? map = length > 0 ? MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite) : null;
        return new MappedOutputFile(path, length, ReservationOutcome.Reopened, map) { Cleared = cleared };
    }

    /// <summary>
    /// Extends the file to <paramref name="newLength"/> bytes (the container writer places the LCP, whose width is known only after PH4-B2, behind the sections that were known after PH0): every writer must be disposed.
    /// The BCL reserves blocks only for a file that it creates (<see cref="FileStreamOptions.PreallocationSize"/> is refused for an existing file), so unless the mode is <see cref="ReservationMode.Sparse"/> the range that
    /// is added is reserved by <see cref="IReservationSource.ReserveRange"/> (zeros written into it) before the mapping is made again over the whole file: a mapped write into it can then never need a block, whatever
    /// another process does to the volume meanwhile, and a full disk is an <see cref="IOException"/> here. In the mode <see cref="ReservationMode.Sparse"/> the file is extended without blocks, after a check that the
    /// volume has the room. Returns <see cref="ReservationOutcome.Filled"/>, or what the check did in the sparse mode (<see cref="ReservationOutcome.FreeSpaceChecked"/>, or <see cref="ReservationOutcome.Unchecked"/> where
    /// the platform cannot tell), or <see cref="ReservationOutcome.Sparse"/> when there is nothing to add. Every failure leaves the file at its old length (the zeros that were written are given back with it).
    /// </summary>
    /// <exception cref="IOException">The volume has less room than the growth needs (found by the check, or by the write).</exception>
    /// <exception cref="InvalidOperationException">A writer is still open.</exception>
    internal ReservationOutcome Grow(long newLength, MappedFileOptions? options = null, IoCounters? counters = null)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentOutOfRangeException.ThrowIfLessThan(newLength, Length);
        if (newLength == Length) return ReservationOutcome.Sparse;
        int open = Volatile.Read(ref _writers);
        if (open != 0) throw new InvalidOperationException($"{open} writer(s) of {Path} are still open: dispose every writer before growing the file");
        IReservationSource source = options?.Source ?? SystemReservationSource.Instance;
        bool reserve = (options is null || options.Reservation == ReservationMode.Auto ? AutoMode : options.Reservation) != ReservationMode.Sparse;
        long growth = newLength - Length;
        string directory = System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(Path)) ?? ".";
        long free = source.AvailableFreeBytes(directory);
        if (free >= 0 && free < growth) throw new IOException($"not enough free space to grow {Path} by {growth} bytes: {free} bytes available");
        long previous = Length;
        lock (_gate)
        {
            using (SafeFileHandle handle = File.OpenHandle(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))
            {
                RandomAccess.SetLength(handle, newLength);
                if (reserve)
                {
                    try
                    {
                        source.ReserveRange(handle, previous, growth, counters);
                    }
                    catch
                    {
                        RandomAccess.SetLength(handle, previous);
                        throw;
                    }
                }
            }
            MemoryMappedFile map;
            try
            {
                map = options?.Mapper is { } mapper ? mapper(Path) : MemoryMappedFile.CreateFromFile(Path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);
            }
            catch
            {
                using var undo = new FileStream(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
                undo.SetLength(previous);
                throw;
            }
            _map?.Dispose();
            _map = map;
            Length = newLength;
        }
        if (reserve) return ReservationOutcome.Filled;
        return free >= 0 ? ReservationOutcome.FreeSpaceChecked : ReservationOutcome.Unchecked;
    }

    /// <summary>A writer for the bytes <c>[start, start + length)</c> of the file; its offsets are relative to <paramref name="start"/>. Ranges of different writers must not overlap (not checked).</summary>
    internal MappedWindowWriter CreateWriter(long start, long length, int windowBytes = MappedWindowWriter.DefaultWindowBytes, IoCounters? counters = null)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowBytes, 1);
        if (start > Length - length) throw new ArgumentOutOfRangeException(nameof(start), $"range [{start}, {start + length}) outside [0, {Length})");
        Interlocked.Increment(ref _writers);
        return new MappedWindowWriter(this, start, length, windowBytes, counters);
    }

    /// <summary>A writer for the whole file.</summary>
    internal MappedWindowWriter CreateWriter(int windowBytes = MappedWindowWriter.DefaultWindowBytes, IoCounters? counters = null) => CreateWriter(0, Length, windowBytes, counters);

    /// <summary>
    /// Declares the file written: every writer must be disposed (each disposal released its window, which the .NET view accessor flushed); with <paramref name="fsync"/> the file is forced to storage (and counted as a flush in
    /// <paramref name="counters"/>, the one thing a test can see of it). May be repeated.
    /// </summary>
    internal void Complete(bool fsync = true, IoCounters? counters = null)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        int open = Volatile.Read(ref _writers);
        if (open != 0) throw new InvalidOperationException($"{open} writer(s) of {Path} are still open: dispose every writer before completing the file");
        if (!fsync) return;
        using var handle = File.OpenHandle(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
        RandomAccess.FlushToDisk(handle);
        counters?.AddFlush();
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        _map?.Dispose();
    }

    /// <summary>A view of the mapping (serialised: the shared mapping is created from several threads).</summary>
    internal MemoryMappedViewAccessor OpenView(long fileOffset, long size)
    {
        lock (_gate)
        {
            ObjectDisposedException.ThrowIf(_disposed, this);
            return _map!.CreateViewAccessor(fileOffset, size, MemoryMappedFileAccess.Write);
        }
    }

    /// <summary>A writer of this file was disposed.</summary>
    internal void WriterClosed() => Interlocked.Decrement(ref _writers);

    /// <summary>Counts a writer that was not made by <see cref="CreateWriter(long, long, int, IoCounters?)"/> (the owning constructor of <see cref="MappedWindowWriter"/>).</summary>
    internal void WriterOpened() => Interlocked.Increment(ref _writers);

    private static ReservationOutcome Allocate(string path, long length, ReservationMode mode, IReservationSource source)
    {
        if (mode == ReservationMode.Sparse)
        {
            using var sparse = new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);
            sparse.SetLength(length);
            return ReservationOutcome.Sparse;
        }
        if (source.CanPreallocate)
        {
            using FileStream reserved = source.CreatePreallocated(path, length);
            reserved.SetLength(length);
            return ReservationOutcome.Preallocated;
        }
        string directory = System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(path)) ?? ".";
        long free = source.AvailableFreeBytes(directory);
        if (free >= 0 && free < length) throw new IOException($"not enough free space for {path}: {length} bytes needed, {free} bytes available");
        using var checkedFile = new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);
        checkedFile.SetLength(length);
        return free >= 0 ? ReservationOutcome.FreeSpaceChecked : ReservationOutcome.Unchecked;
    }

    private static void DeleteQuietly(string path)
    {
        try
        {
            File.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // the path is not ours to clean up if it cannot be deleted; the original exception is what the caller needs
        }
    }
}
