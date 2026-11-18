using System.IO.MemoryMappedFiles;

namespace SuffixTree.Esa.External;

/// <summary>
/// Writes a range of a file of known length through bounded sliding windows of a memory mapping (the categorical rule of the owner, 2026-11-04: stream into memory-mapped files, no "build in RAM and then write"): the file
/// is created at its final length, a view of at most <c>windowBytes</c> is mapped around the write position, the bytes are copied into it, and the view is disposed when the writes move on. The pages are file-backed
/// and leave the process as the kernel writes them back, but the .NET view accessor flushes (<c>msync</c>) a window when it is disposed, so the window size is a first-order cost parameter: E9 measured 7.2 B/char of device
/// writes with 64 KB windows against 5.2 with 1 MiB, hence <see cref="DefaultWindowBytes"/>. One writer is not thread-safe; several writers of one <see cref="MappedOutputFile"/> (disjoint ranges) may run on several threads.
/// <see cref="Append"/> continues where the last write ended; <see cref="WriteAt"/> may jump (a jump opens a new window). A file created sparse (the constructor with a path does that; <see cref="MappedOutputFile.Create"/>
/// can reserve the blocks) shows a full disk as a signal on a page fault of a mapped view, not as an exception.
/// </summary>
internal sealed unsafe class MappedWindowWriter : IDisposable
{
    /// <summary>Default window: 1 MiB.</summary>
    internal const int DefaultWindowBytes = 1 << 20;

    private readonly MappedOutputFile _file;
    private readonly bool _ownsFile;
    private readonly long _start;
    private readonly int _window;
    private readonly IoCounters? _counters;
    private MemoryMappedViewAccessor? _view;
    private byte* _base;
    private long _viewStart;
    private long _viewEnd;
    private long _next;
    private bool _closed;

    /// <summary>A writer that owns a new file of <paramref name="length"/> bytes, created sparse (the layout E9 measured; <see cref="MappedOutputFile.Create"/> offers the reservation).</summary>
    internal MappedWindowWriter(string path, long length, int windowBytes = DefaultWindowBytes, IoCounters? counters = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowBytes, 1);
        _file = MappedOutputFile.Create(path, length, new MappedFileOptions { Reservation = ReservationMode.Sparse });
        _ownsFile = true;
        _file.WriterOpened();
        _start = 0;
        Length = length;
        _window = windowBytes;
        _counters = counters;
    }

    /// <summary>A writer for the range <c>[start, start + length)</c> of <paramref name="file"/> (made by <see cref="MappedOutputFile.CreateWriter(long, long, int, IoCounters?)"/>, which has counted it).</summary>
    internal MappedWindowWriter(MappedOutputFile file, long start, long length, int windowBytes, IoCounters? counters)
    {
        _file = file;
        _ownsFile = false;
        _start = start;
        Length = length;
        _window = windowBytes;
        _counters = counters;
    }

    internal string Path => _file.Path;

    /// <summary>Length of the range in bytes.</summary>
    internal long Length { get; }

    /// <summary>Offset of the next <see cref="Append"/>, relative to the start of the range.</summary>
    internal long Position => _next;

    /// <summary>Offset of the first byte of the range in the file.</summary>
    internal long FileOffset => _start;

    /// <summary>Writes <paramref name="data"/> at <paramref name="offset"/> (relative to the range); the range of the write must lie inside the writer's range.</summary>
    internal void WriteAt(long offset, ReadOnlySpan<byte> data)
    {
        if (offset < 0 || offset > Length - data.Length) throw new ArgumentOutOfRangeException(nameof(offset), $"range [{offset}, {offset + data.Length}) outside [0, {Length})");
        while (!data.IsEmpty)
        {
            if (_view is null || offset < _viewStart || offset >= _viewEnd) OpenWindow(offset);
            int take = (int)Math.Min(data.Length, _viewEnd - offset);
            data[..take].CopyTo(new Span<byte>(_base + (offset - _viewStart), take));
            _counters?.AddMappedWrite(take);
            offset += take;
            data = data[take..];
        }
    }

    /// <summary>Writes <paramref name="data"/> right after the previous write.</summary>
    internal void Append(ReadOnlySpan<byte> data)
    {
        WriteAt(_next, data);
        _next += data.Length;
    }

    /// <summary>Releases the current window (the .NET view accessor flushes it when it is disposed); the next write opens another. A writer that is finished with a stretch of the file calls this to bound its dirty pages.</summary>
    internal void Flush() => CloseWindow();

    public void Dispose()
    {
        if (_closed) return;
        _closed = true;
        CloseWindow();
        if (_ownsFile) _file.Dispose();
        else _file.WriterClosed();
    }

    private void OpenWindow(long offset)
    {
        CloseWindow();
        long size = Math.Min(_window, Length - offset);
        _view = _file.OpenView(_start + offset, size);
        byte* p = null;
        _view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        _base = p + _view.PointerOffset;
        _viewStart = offset;
        _viewEnd = offset + size;
        _counters?.AddWindow(size);
    }

    private void CloseWindow()
    {
        if (_view is null) return;
        _view.SafeMemoryMappedViewHandle.ReleasePointer();
        _view.Dispose();
        _view = null;
        _base = null;
        _viewStart = _viewEnd = 0;
    }
}
