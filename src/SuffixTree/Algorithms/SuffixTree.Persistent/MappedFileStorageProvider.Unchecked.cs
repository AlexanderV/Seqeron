using System.Diagnostics;
using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

public sealed unsafe partial class MappedFileStorageProvider
{
    // ──────────────── Unchecked fast-path for builder hot loop ────────────────
    // These skip ThrowIfDisposed + CheckBounds. The builder guarantees all
    // offsets come from Allocate() and disposal never races with build.

    /// <summary>Raw pointer — use at your own risk.</summary>
    internal byte* RawPointer
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _ptr;
    }

    [Conditional("DEBUG")]
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void AssertUncheckedAccess(long offset, int size, bool isWrite)
    {
        Debug.Assert(_ptr != null, "Unchecked access requires an acquired MMF pointer.");
        Debug.Assert(size > 0, "Unchecked access size must be positive.");
        Debug.Assert(offset >= 0, "Unchecked access offset must be non-negative.");
        Debug.Assert(offset <= long.MaxValue - size, "Unchecked access end offset overflow.");

        long end = offset + size;
        Debug.Assert(end <= _capacity, $"Unchecked access out of physical bounds: end={end}, capacity={_capacity}.");
        Debug.Assert(end <= _position, $"Unchecked access out of logical bounds: end={end}, size={_position}.");

        if (isWrite)
            Debug.Assert(!_readOnly, "Unchecked write is invalid for read-only storage.");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint ReadUInt32Unchecked(long offset)
    {
        AssertUncheckedAccess(offset, 4, isWrite: false);
        return *(uint*)(_ptr + offset);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void WriteUInt32Unchecked(long offset, uint value)
    {
        AssertUncheckedAccess(offset, 4, isWrite: true);
        *(uint*)(_ptr + offset) = value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal int ReadInt32Unchecked(long offset)
    {
        AssertUncheckedAccess(offset, 4, isWrite: false);
        return *(int*)(_ptr + offset);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void WriteInt32Unchecked(long offset, int value)
    {
        AssertUncheckedAccess(offset, 4, isWrite: true);
        *(int*)(_ptr + offset) = value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal long ReadInt64Unchecked(long offset)
    {
        AssertUncheckedAccess(offset, 8, isWrite: false);
        return *(long*)(_ptr + offset);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void WriteInt64Unchecked(long offset, long value)
    {
        AssertUncheckedAccess(offset, 8, isWrite: true);
        *(long*)(_ptr + offset) = value;
    }
}
