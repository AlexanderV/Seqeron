using System.Collections;
using System.Collections.Generic;
using System.IO.MemoryMappedFiles;
using System.Runtime.InteropServices;

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
    private volatile bool _disposed;

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
        _length = length;
        _ownsAccessor = false;

        byte* basePtr = null;
        _accessor.SafeMemoryMappedViewHandle.AcquirePointer(ref basePtr);
        // Correct pointer calculation: basePtr is the start of the view, 
        // offset is already used in creating the accessor if called from Constructor 1,
        // but for Constructor 2 (shared accessor) we need to add the relative offset.
        _ptr = (char*)(basePtr + _accessor.PointerOffset + offset);
    }

    /// <inheritdoc/>
    public int Length => _length;

    /// <inheritdoc/>
    public char this[int index]
    {
        get
        {
            if (_disposed) throw new ObjectDisposedException(nameof(MemoryMappedTextSource));
            if (index < 0 || index >= _length) throw new IndexOutOfRangeException();
            return _ptr[index];
        }
    }

    /// <inheritdoc/>
    public string Substring(int start, int length)
    {
        if (_disposed) throw new ObjectDisposedException(nameof(MemoryMappedTextSource));
        if (start < 0 || start + length > _length) throw new IndexOutOfRangeException();
        return new string(_ptr, start, length);
    }

    /// <inheritdoc/>
    public ReadOnlySpan<char> Slice(int start, int length)
    {
        if (_disposed) throw new ObjectDisposedException(nameof(MemoryMappedTextSource));
        if (start < 0 || start + length > _length) throw new IndexOutOfRangeException();
        return new ReadOnlySpan<char>(_ptr + start, length);
    }

    /// <inheritdoc/>
    public IEnumerator<char> GetEnumerator()
    {
        if (_disposed) throw new ObjectDisposedException(nameof(MemoryMappedTextSource));
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
    public override string ToString() => new string(_ptr, 0, _length);

    /// <inheritdoc/>
    public void Dispose()
    {
        if (!_disposed)
        {
            _disposed = true; // Mark disposed first to stop readers
            _ptr = null;      // Prevent use-after-free of released pointer

            // Release pointer before disposing accessor
            _accessor.SafeMemoryMappedViewHandle.ReleasePointer();
            if (_ownsAccessor)
            {
                _accessor.Dispose();
                _mmf?.Dispose();
            }
        }
    }
}
