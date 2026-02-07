using System.Collections;
using System.IO.MemoryMappedFiles;

namespace SuffixTree.Persistent;

/// <summary>
/// A memory-mapped implementation of <see cref="ITextSource"/>.
/// Provides efficient character access for very large files.
/// </summary>
public sealed unsafe class MemoryMappedTextSource : ITextSource, IDisposable
{
    private readonly MemoryMappedFile? _mmf;
    private readonly MemoryMappedViewAccessor _accessor;
    private char* _ptr;
    private readonly int _length;
    private readonly bool _ownsAccessor;
    private int _disposed; // 0 = active, 1 = disposed (Interlocked for C14)

    private void ThrowIfDisposed()
    {
        ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="MemoryMappedTextSource"/> class from a file.
    /// </summary>
    public MemoryMappedTextSource(string filePath, long offset, int length)
    {
        _mmf = MemoryMappedFile.CreateFromFile(filePath, System.IO.FileMode.Open, null, 0, MemoryMappedFileAccess.Read);
        try
        {
            _accessor = _mmf.CreateViewAccessor(offset, length * sizeof(char), MemoryMappedFileAccess.Read);
            _length = length;
            _ownsAccessor = true;

            byte* basePtr = null;
            _accessor.SafeMemoryMappedViewHandle.AcquirePointer(ref basePtr);
            _ptr = (char*)(basePtr + _accessor.PointerOffset);
        }
        catch
        {
            _accessor?.Dispose();
            _mmf.Dispose();
            throw;
        }
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="MemoryMappedTextSource"/> class from an existing accessor.
    /// </summary>
    public MemoryMappedTextSource(MemoryMappedViewAccessor accessor, long offset, int length)
    {
        _accessor = accessor ?? throw new ArgumentNullException(nameof(accessor));
        if (offset < 0) throw new ArgumentOutOfRangeException(nameof(offset), "Offset must be non-negative.");
        if (length < 0) throw new ArgumentOutOfRangeException(nameof(length), "Length must be non-negative.");
        long requiredBytes = offset + (long)length * sizeof(char);
        if (requiredBytes > accessor.Capacity)
            throw new ArgumentOutOfRangeException(nameof(length),
                $"offset ({offset}) + length ({length}) * 2 = {requiredBytes} exceeds accessor capacity ({accessor.Capacity}).");
        _length = length;
        _ownsAccessor = false;

        byte* basePtr = null;
        try
        {
            _accessor.SafeMemoryMappedViewHandle.AcquirePointer(ref basePtr);
            _ptr = (char*)(basePtr + _accessor.PointerOffset + offset);
        }
        catch
        {
            if (basePtr != null)
                _accessor.SafeMemoryMappedViewHandle.ReleasePointer();
            throw;
        }
    }

    /// <inheritdoc/>
    public int Length => _length;

    /// <inheritdoc/>
    public char this[int index]
    {
        get
        {
            char* p = _ptr;  // C17: snapshot before disposed check
            ThrowIfDisposed();
            ObjectDisposedException.ThrowIf(p == null, this);
            if (index < 0 || index >= _length) throw new ArgumentOutOfRangeException(nameof(index));
            return p[index];
        }
    }

    /// <inheritdoc/>
    public string Substring(int start, int length)
    {
        char* p = _ptr;  // C17: snapshot before disposed check
        ThrowIfDisposed();
        ObjectDisposedException.ThrowIf(p == null, this);
        if ((uint)start > (uint)_length || (uint)length > (uint)(_length - start))
            throw new ArgumentOutOfRangeException(nameof(start));
        return new string(p, start, length);
    }

    /// <inheritdoc/>
    public ReadOnlySpan<char> Slice(int start, int length)
    {
        char* p = _ptr;  // C17: snapshot before disposed check
        ThrowIfDisposed();
        ObjectDisposedException.ThrowIf(p == null, this);
        if ((uint)start > (uint)_length || (uint)length > (uint)(_length - start))
            throw new ArgumentOutOfRangeException(nameof(start));
        return new ReadOnlySpan<char>(p + start, length);
    }

    /// <inheritdoc/>
    public IEnumerator<char> GetEnumerator()
    {
        ThrowIfDisposed();
        return GetEnumeratorImpl();
    }

    private IEnumerator<char> GetEnumeratorImpl()
    {
        for (int i = 0; i < _length; i++)
        {
            yield return this[i];
        }
    }

    /// <inheritdoc/>
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <inheritdoc/>
    public override string ToString()
    {
        char* p = _ptr;  // C17: snapshot before disposed check
        ThrowIfDisposed();
        ObjectDisposedException.ThrowIf(p == null, this);
        return new string(p, 0, _length);
    }

    /// <inheritdoc/>
    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;

        // C17: Null _ptr FIRST so that concurrent readers snapshot null and
        // throw ObjectDisposedException via the p == null guard. Readers that
        // already captured a non-null local pointer are safe: the actual mapped
        // view remains valid until _accessor.Dispose() below.
        _ptr = null;
        _accessor.SafeMemoryMappedViewHandle.ReleasePointer();

        if (_ownsAccessor)
        {
            _accessor.Dispose();
            _mmf?.Dispose();
        }
    }
}
