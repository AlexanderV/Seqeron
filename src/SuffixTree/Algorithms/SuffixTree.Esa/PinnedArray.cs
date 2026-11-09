using System.Runtime.InteropServices;

namespace SuffixTree.Esa;

/// <summary>
/// An array whose address does not change for the life of the object, viewed as a <see cref="Buf{T}"/> (work package 1b.4: the navigation pack reads the suffix array,
/// the LCP base array, ISA and the pyramid through pointers, and a navigator must stay valid between calls, so a <c>fixed</c> block does not do). Two ways to get one:
/// <see cref="Allocate"/> puts a new array on the pinned object heap (no handle), <see cref="Pin"/> pins an array that already exists (a handle, freed by
/// <see cref="Dispose"/> or by the finalizer). The wrapped array is kept alive by this object; the <see cref="Buf{T}"/> is valid while this object is.
/// </summary>
/// <typeparam name="T">Element type (unmanaged).</typeparam>
internal sealed unsafe class PinnedArray<T> : IDisposable
    where T : unmanaged
{
    private readonly T[] _array;
    private GCHandle _handle;

    private PinnedArray(T[] array, GCHandle handle, string name)
    {
        _array = array;
        _handle = handle;
        Name = name;
        fixed (T* p = array)
        {
            // the address of an element-0 reference of a pinned array is stable; the pointer is taken once
            Buf = new Buf<T>(p, array.Length, name);
        }
    }

    ~PinnedArray() => Release();

    /// <summary>The name given at creation (shown by checked-build violations).</summary>
    internal string Name { get; }

    /// <summary>The elements as a pointer view (bounds-checked in a checked build).</summary>
    internal Buf<T> Buf { get; }

    /// <summary>Number of elements.</summary>
    internal int Length => _array.Length;

    /// <summary>The wrapped array (for managed-code readers and tests).</summary>
    internal T[] Array => _array;

    /// <summary>True when this object holds a handle that <see cref="Dispose"/> frees (an array that was pinned, not allocated pinned).</summary>
    internal bool HoldsHandle => _handle.IsAllocated;

    /// <summary>A new array of <paramref name="length"/> elements on the pinned object heap; zeroed when <paramref name="zero"/>, else uninitialised.</summary>
    internal static PinnedArray<T> Allocate(int length, bool zero, string name)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        T[] array = zero ? GC.AllocateArray<T>(length, pinned: true) : GC.AllocateUninitializedArray<T>(length, pinned: true);
        return new PinnedArray<T>(array, default, name);
    }

    /// <summary>Pins an existing array until <see cref="Dispose"/> (or the finalizer).</summary>
    internal static PinnedArray<T> Pin(T[] array, string name)
    {
        ArgumentNullException.ThrowIfNull(array);
        return new PinnedArray<T>(array, GCHandle.Alloc(array, GCHandleType.Pinned), name);
    }

    /// <summary>Frees the handle of a pinned array; idempotent. The <see cref="Buf"/> must not be used afterwards.</summary>
    public void Dispose()
    {
        Release();
        GC.SuppressFinalize(this);
    }

    private void Release()
    {
        if (_handle.IsAllocated) _handle.Free();
    }
}
