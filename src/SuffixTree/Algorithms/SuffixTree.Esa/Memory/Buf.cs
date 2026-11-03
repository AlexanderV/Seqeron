using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>Build-mode information of the safety infrastructure.</summary>
internal static class EsaBuild
{
    /// <summary>
    /// True when the assembly was compiled with <c>ESA_BOUNDS</c> (<c>dotnet build -p:EsaBounds=true</c>): every <see cref="Buf{T}"/> access is
    /// bounds-checked. False in Release: a <see cref="Buf{T}"/> is a bare pointer and an access is one machine load or store.
    /// </summary>
    internal static bool BoundsChecked
    {
        get
        {
#if ESA_BOUNDS
            return true;
#else
            return false;
#endif
        }
    }
}

/// <summary>
/// The accessor through which every unsafe kernel touches memory (safety infrastructure, work package 1o.0). In a Release build it is a single
/// pointer with an inlined indexer: the generated code is that of a raw pointer access. With the compile-time symbol <c>ESA_BOUNDS</c>
/// it also carries the logical length and a name, and every access, slice and word read is bounds-checked and throws
/// <see cref="EsaMemoryViolationException"/> naming the buffer, the index and the length. Kernels must not extract the pointer and index it
/// themselves; bulk work goes through <see cref="AsSpan"/>, <see cref="Slice(int, int)"/>, <see cref="Fill"/> and <see cref="ReadUnaligned{TWord}"/>.
/// </summary>
/// <typeparam name="T">Element type (unmanaged).</typeparam>
internal readonly unsafe struct Buf<T>
    where T : unmanaged
{
    private readonly T* _p;
#if ESA_BOUNDS
    private readonly int _length;
    private readonly string? _name;

    /// <summary>Wraps <paramref name="length"/> elements at <paramref name="pointer"/> (the caller owns the memory and keeps it alive and pinned).</summary>
    internal Buf(T* pointer, int length, string? name = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        _p = pointer;
        _length = length;
        _name = name;
    }

    /// <summary>The logical length (known in checked builds only; Release code passes lengths explicitly).</summary>
    internal int CheckedLength => _length;

    /// <summary>Element access; bounds-checked in this build.</summary>
    internal ref T this[int index]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if ((uint)index >= (uint)_length) ThrowOutOfBounds(index, 1);
            return ref _p[index];
        }
    }

    /// <summary>A view of <paramref name="length"/> elements starting at <paramref name="start"/>.</summary>
    internal Buf<T> Slice(int start, int length)
    {
        if (start < 0 || length < 0 || start > _length - length) ThrowOutOfBounds(start, length);
        return new Buf<T>(_p + start, length, _name);
    }

    /// <summary>A view from <paramref name="start"/> to the end.</summary>
    internal Buf<T> Slice(int start)
    {
        if ((uint)start > (uint)_length) ThrowOutOfBounds(start, 0);
        return new Buf<T>(_p + start, _length - start, _name);
    }

    /// <summary>Reads a <typeparamref name="TWord"/> from the bytes that start at element <paramref name="index"/> (unaligned); the whole word must lie inside the buffer.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal TWord ReadUnaligned<TWord>(int index)
        where TWord : unmanaged
    {
        if (index < 0 || ((long)index * sizeof(T)) + sizeof(TWord) > (long)_length * sizeof(T)) ThrowOutOfBounds(index, (sizeof(TWord) + sizeof(T) - 1) / sizeof(T));
        return Unsafe.ReadUnaligned<TWord>(_p + index);
    }

    /// <summary>A span over the first <paramref name="length"/> elements (safe-code access to the same memory).</summary>
    internal Span<T> AsSpan(int length)
    {
        if ((uint)length > (uint)_length) ThrowOutOfBounds(0, length);
        return new Span<T>(_p, length);
    }

    private void ThrowOutOfBounds(int index, int count)
        => throw new EsaMemoryViolationException(
            MemoryViolationKind.OutOfBounds,
            $"out-of-bounds access to buffer '{_name ?? typeof(T).Name}': index {index}, count {count}, length {_length}");
#else
    /// <summary>Wraps the memory at <paramref name="pointer"/> (the caller owns it and keeps it alive and pinned); length and name are used by checked builds only.</summary>
    internal Buf(T* pointer, int length, string? name = null)
    {
        _p = pointer;
        _ = length; // Release keeps no length and no name: the struct is one pointer
        _ = name;
    }

    /// <summary>Element access; one load or store.</summary>
    internal ref T this[int index]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => ref _p[index];
    }

    /// <summary>A view of the elements starting at <paramref name="start"/>.</summary>
    internal Buf<T> Slice(int start, int length) => new(_p + start, length);

    /// <summary>A view from <paramref name="start"/> to the end.</summary>
    internal Buf<T> Slice(int start) => new(_p + start, 0);

    /// <summary>Reads a <typeparamref name="TWord"/> from the bytes that start at element <paramref name="index"/> (unaligned).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal TWord ReadUnaligned<TWord>(int index)
        where TWord : unmanaged
        => Unsafe.ReadUnaligned<TWord>(_p + index);

    /// <summary>A span over the first <paramref name="length"/> elements.</summary>
    internal Span<T> AsSpan(int length) => new(_p, length);
#endif

    /// <summary>Sets the first <paramref name="count"/> elements to <paramref name="value"/>.</summary>
    internal void Fill(T value, int count) => AsSpan(count).Fill(value);

    /// <summary>Sets the first <paramref name="count"/> elements to zero.</summary>
    internal void Clear(int count) => AsSpan(count).Clear();

    /// <summary>The wrapped pointer, for pinning bookkeeping and the arena only; kernels do not index it.</summary>
    internal T* Pointer => _p;
}
