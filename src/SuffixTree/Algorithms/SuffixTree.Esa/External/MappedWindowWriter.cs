using System.IO.MemoryMappedFiles;

namespace SuffixTree.Esa.External;

/// <summary>
/// Writes a file of known length through bounded sliding windows of a memory mapping (the categorical rule of the owner, 2026-11-04: stream into memory-mapped files, no "build in RAM and then write"): the file is
/// created at its final length, a view of at most <c>windowBytes</c> is mapped around the write position, the bytes are copied into it, and the view is disposed when the writes move on. There is no explicit flush: the
/// pages are file-backed and leave the process as the kernel writes them back, but the .NET view accessor flushes (<c>msync</c>) a window when it is disposed, so the window size is a first-order cost parameter: E9
/// measured 7.2 B/char of device writes with 64 KB windows against 5.2 with 1 MiB, hence <see cref="DefaultWindowBytes"/>. Not thread-safe. <see cref="Append"/> continues where the last write ended;
/// <see cref="WriteAt"/> may jump (a jump opens a new window). A full disk surfaces as a signal on a page fault of a mapped view, not as an exception: the file is created sparse, nothing here reserves the blocks.
/// </summary>
internal sealed unsafe class MappedWindowWriter : IDisposable
{
    /// <summary>Default window: 1 MiB.</summary>
    internal const int DefaultWindowBytes = 1 << 20;

    private readonly MemoryMappedFile? _file;
    private readonly int _window;
    private readonly IoCounters? _counters;
    private MemoryMappedViewAccessor? _view;
    private byte* _base;
    private long _viewStart;
    private long _viewEnd;
    private long _next;

    internal MappedWindowWriter(string path, long length, int windowBytes = DefaultWindowBytes, IoCounters? counters = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowBytes, 1);
        Path = path;
        Length = length;
        _window = windowBytes;
        _counters = counters;
        using (var fs = new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite)) fs.SetLength(length);
        if (length > 0) _file = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);
    }

    internal string Path { get; }

    /// <summary>Length of the file in bytes.</summary>
    internal long Length { get; }

    /// <summary>Offset of the next <see cref="Append"/>.</summary>
    internal long Position => _next;

#pragma warning disable S2325
    /// <summary>STUB: offset of the first byte of this writer's range in the file.</summary>
    internal long FileOffset => throw new NotImplementedException();

    /// <summary>STUB: releases (and thereby flushes) the current window.</summary>
    internal void Flush() => throw new NotImplementedException();
#pragma warning restore S2325

    /// <summary>Writes <paramref name="data"/> at <paramref name="offset"/>; the range must lie inside the file.</summary>
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

    public void Dispose()
    {
        CloseWindow();
        _file?.Dispose();
    }

    private void OpenWindow(long offset)
    {
        CloseWindow();
        long size = Math.Min(_window, Length - offset);
        _view = _file!.CreateViewAccessor(offset, size, MemoryMappedFileAccess.Write);
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
