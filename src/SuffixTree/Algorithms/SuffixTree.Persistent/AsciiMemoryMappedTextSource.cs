using System.Collections;
using System.IO.MemoryMappedFiles;
using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

/// <summary>
/// A memory-mapped <see cref="ITextSource"/> for ASCII-encoded text (1 byte per char).
/// Provides the same API as <see cref="MemoryMappedTextSource"/> but reads bytes
/// and widens them to <see langword="char"/>, halving the on-disk footprint for
/// pure-ASCII inputs such as DNA sequences.
/// </summary>
public sealed unsafe class AsciiMemoryMappedTextSource : ITextSource, ITextPatternMatcher, IDisposable
{
    private readonly MemoryMappedViewAccessor _accessor;
    private byte* _ptr;
    private readonly int _length;
    private readonly bool _ownsAccessor;
    private int _disposed; // 0 = active, 1 = disposed

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ThrowIfDisposed()
    {
        ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);
    }

    /// <summary>
    /// Initializes from an existing accessor (zero-copy within a <see cref="MappedFileStorageProvider"/>).
    /// </summary>
    public AsciiMemoryMappedTextSource(MemoryMappedViewAccessor accessor, long offset, int length)
    {
        _accessor = accessor ?? throw new ArgumentNullException(nameof(accessor));
        if (offset < 0) throw new ArgumentOutOfRangeException(nameof(offset), "Offset must be non-negative.");
        if (length < 0) throw new ArgumentOutOfRangeException(nameof(length), "Length must be non-negative.");
        long requiredBytes = offset + length; // 1 byte per char
        if (requiredBytes > accessor.Capacity)
            throw new ArgumentOutOfRangeException(nameof(length),
                $"offset ({offset}) + length ({length}) = {requiredBytes} exceeds accessor capacity ({accessor.Capacity}).");
        _length = length;
        _ownsAccessor = false;

        byte* basePtr = null;
        try
        {
            _accessor.SafeMemoryMappedViewHandle.AcquirePointer(ref basePtr);
            _ptr = basePtr + _accessor.PointerOffset + offset;
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
            byte* p = _ptr;
            ThrowIfDisposed();
            ObjectDisposedException.ThrowIf(p == null, this);
            if (index < 0 || index >= _length) throw new ArgumentOutOfRangeException(nameof(index));
            return (char)p[index];
        }
    }

    /// <inheritdoc/>
    public string Substring(int start, int length)
    {
        byte* p = _ptr;
        ThrowIfDisposed();
        ObjectDisposedException.ThrowIf(p == null, this);
        if ((uint)start > (uint)_length || (uint)length > (uint)(_length - start))
            throw new ArgumentOutOfRangeException(nameof(start));

        // Widen bytes to chars
        return string.Create(length, (Ptr: (nint)(p + start), Len: length), static (span, state) =>
        {
            byte* src = (byte*)state.Ptr;
            for (int i = 0; i < state.Len; i++)
                span[i] = (char)src[i];
        });
    }

    /// <inheritdoc/>
    public ReadOnlySpan<char> Slice(int start, int length)
    {
        byte* p = _ptr;
        ThrowIfDisposed();
        ObjectDisposedException.ThrowIf(p == null, this);
        if ((uint)start > (uint)_length || (uint)length > (uint)(_length - start))
            throw new ArgumentOutOfRangeException(nameof(start));

        // Must allocate because we can't return a Span<char> over byte* without widening
        char[] result = new char[length];
        for (int i = 0; i < length; i++)
            result[i] = (char)p[start + i];
        return result;
    }

    /// <summary>
    /// Compares a pattern with the ASCII source at the specified start offset
    /// without widening the full slice into a temporary <see cref="char"/> array.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool SequenceEqualAt(int start, ReadOnlySpan<char> pattern)
    {
        byte* p = _ptr;
        ThrowIfDisposed();
        ObjectDisposedException.ThrowIf(p == null, this);

        int length = pattern.Length;
        if ((uint)start > (uint)_length || (uint)length > (uint)(_length - start))
            return false;

        for (int i = 0; i < length; i++)
        {
            if ((char)p[start + i] != pattern[i])
                return false;
        }

        return true;
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
            yield return this[i];
    }

    /// <inheritdoc/>
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <inheritdoc/>
    public override string ToString()
    {
        byte* p = _ptr;
        ThrowIfDisposed();
        ObjectDisposedException.ThrowIf(p == null, this);
        return Substring(0, _length);
    }

    /// <inheritdoc/>
    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;
        _ptr = null;
        _accessor.SafeMemoryMappedViewHandle.ReleasePointer();
        if (_ownsAccessor)
        {
            _accessor.Dispose();
        }
    }
}
