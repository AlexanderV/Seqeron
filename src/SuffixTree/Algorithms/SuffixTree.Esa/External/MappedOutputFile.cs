using System.IO.MemoryMappedFiles;

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
    /// Extends the file to <paramref name="newLength"/> bytes (the container writer places the LCP, whose width is known only after PH4-B2, behind the sections that were known after PH0): every writer must be disposed; the file
    /// is extended sparse (the BCL reserves blocks only for a file that it creates), after a check that the volume has the room, and the mapping is made again over the whole file. Returns what the check did
    /// (<see cref="ReservationOutcome.FreeSpaceChecked"/>, or <see cref="ReservationOutcome.Unchecked"/> where the platform cannot tell). A failure leaves the file at its old length.
    /// </summary>
    /// <exception cref="IOException">The volume has less room than the growth needs.</exception>
    /// <exception cref="InvalidOperationException">A writer is still open.</exception>
    internal ReservationOutcome Grow(long newLength, MappedFileOptions? options = null)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentOutOfRangeException.ThrowIfLessThan(newLength, Length);
        if (newLength == Length) return ReservationOutcome.Sparse;
        int open = Volatile.Read(ref _writers);
        if (open != 0) throw new InvalidOperationException($"{open} writer(s) of {Path} are still open: dispose every writer before growing the file");
        IReservationSource source = options?.Source ?? SystemReservationSource.Instance;
        long growth = newLength - Length;
        string directory = System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(Path)) ?? ".";
        long free = source.AvailableFreeBytes(directory);
        if (free >= 0 && free < growth) throw new IOException($"not enough free space to grow {Path} by {growth} bytes: {free} bytes available");
        long previous = Length;
        lock (_gate)
        {
            using (var stream = new FileStream(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite)) stream.SetLength(newLength);
            MemoryMappedFile map;
            try
            {
                map = MemoryMappedFile.CreateFromFile(Path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);
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
