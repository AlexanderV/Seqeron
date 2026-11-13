using System.Runtime.CompilerServices;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Ph5;

/// <summary>What a <see cref="SpillStack{T}"/> did with its spill file (ADR-004 section 3.2: "older blocks are appended to a spill file, read back backwards when a pop leaves the window").</summary>
internal readonly record struct SpillStats(long BlocksWritten, long BlocksRead, long BytesWritten, long BytesRead, long PeakFileBytes, long Flushes, long CleanEvictions);

/// <summary>
/// The external stack of the PH5 scan (ADR-004 section 3.2). The stack is cut into blocks of <c>B</c> entries; <b>two block buffers are in RAM</b>: the block of the top entry and either the full block below it or, when the top
/// block has just been entered from below, a spare. When the top block is full and another entry is pushed, the block below is written to the spill file (block <c>j</c> at offset <c>j x B x sizeof(T)</c>, the file is
/// accessed only at its end: LIFO) and its buffer becomes the new top block; when a pop empties the top block the block below becomes the top and, if it is not resident, it is <b>read back</b> from the file (backwards,
/// block by block). The hysteresis is the one of the ADR: after a refill the block above is a spare and the block below is not needed until the top has climbed through a whole block again, so a push / pop alternation at a
/// block boundary costs no I/O; a block that was read back and not modified (its file copy is current) is not written again. Not thread-safe.
/// </summary>
internal sealed unsafe class SpillStack<T> : IDisposable
    where T : unmanaged
{
    private readonly int _block;
    private readonly int _blockBytes;
    private readonly string _path;
    private readonly long _flushBytes;
    private T* _cur;
    private T* _other;
    private bool _curDirty;
    private bool _otherDirty;
    private bool _otherBelow;
    private long _topBlock;
    private int _count;
    private long _height;
    private long _peak;
    private CountedFile? _file;
    private long _fileBytes;
    private long _unflushed;
    private long _blocksWritten;
    private long _blocksRead;
    private long _bytesWritten;
    private long _bytesRead;
    private long _flushes;
    private long _cleanEvictions;
    private bool _disposed;

    /// <param name="arena">Owner of the two block buffers (charged to its budget).</param>
    /// <param name="spillPath">The spill file; created at the first spill, removed on <see cref="Dispose"/>.</param>
    /// <param name="blockEntries">Entries per block (ADR: 256 K); at least 2.</param>
    /// <param name="flushBytes">fsync the file after this many bytes written since the last flush (0 = never): bounds the dirty pages of the spill file.</param>
    internal SpillStack(BudgetedArena arena, string spillPath, int blockEntries, long flushBytes = 0)
    {
        ArgumentNullException.ThrowIfNull(arena);
        ArgumentOutOfRangeException.ThrowIfLessThan(blockEntries, 2);
        _block = blockEntries;
        _blockBytes = checked(blockEntries * sizeof(T));
        _path = spillPath;
        _flushBytes = flushBytes;
        _cur = arena.Allocate<T>(blockEntries, zero: false, "spill.a").Pointer;
        _other = arena.Allocate<T>(blockEntries, zero: false, "spill.b").Pointer;
    }

    /// <summary>Entries on the stack.</summary>
    internal long Height => _height;

    /// <summary>Largest height reached.</summary>
    internal long PeakHeight => _peak;

    internal int BlockEntries => _block;

    internal SpillStats Stats => new(_blocksWritten, _blocksRead, _bytesWritten, _bytesRead, _fileBytes, _flushes, _cleanEvictions);

    /// <summary>The top entry (read access; a write through this reference is not seen by the spill: use <see cref="TopForUpdate"/>).</summary>
    internal ref T Top
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if (_height == 0) ThrowEmpty();
            return ref _cur[_count - 1];
        }
    }

    /// <summary>The top entry for modification: its block is marked dirty, so a later spill writes it.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ref T TopForUpdate()
    {
        if (_height == 0) ThrowEmpty();
        _curDirty = true;
        return ref _cur[_count - 1];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void Push(in T value)
    {
        if (_count == _block) Advance();
        _cur[_count++] = value;
        _curDirty = true;
        if (++_height > _peak) _peak = _height;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void Pop()
    {
        if (_height == 0) ThrowEmpty();
        _height--;
        if (--_count == 0 && _height > 0) Retreat();
    }

    private static void ThrowEmpty() => throw new InvalidOperationException("the stack is empty");

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        if (_file is null) return;
        _file.Dispose();
        _file = null;
        File.Delete(_path);
        ScratchTracker.Remove(_fileBytes);
    }

    /// <summary>The top block is full and one more entry comes: the next block becomes the top block.</summary>
    [MethodImpl(MethodImplOptions.NoInlining)]
    private void Advance()
    {
        if (_otherBelow)
        {
            // the other buffer holds the full block below the top block: it leaves RAM (written unless its file copy is current)
            if (_otherDirty) WriteBlock(_topBlock - 1, _other);
            else _cleanEvictions++;
        }
        // the other buffer (the evicted block, or the spare) becomes the top block; the full top block becomes the block below
        T* below = _cur;
        _cur = _other;
        _other = below;
        _otherDirty = _curDirty;
        _curDirty = false;
        _otherBelow = true;
        _topBlock++;
        _count = 0;
    }

    /// <summary>The top block is empty and entries remain: the full block below becomes the top block (read back when it is not resident).</summary>
    [MethodImpl(MethodImplOptions.NoInlining)]
    private void Retreat()
    {
        if (!_otherBelow)
        {
            ReadBlock(_topBlock - 1, _other);
            _otherDirty = false;
        }
        T* top = _other;
        _other = _cur;
        _cur = top;
        _curDirty = _otherDirty;
        _otherDirty = false;
        _otherBelow = false;
        _topBlock--;
        _count = _block;
    }

    private void WriteBlock(long block, T* buffer)
    {
        _file ??= new CountedFile(_path, FileMode.Create, FileAccess.ReadWrite);
        _file.WriteAt(new ReadOnlySpan<byte>(buffer, _blockBytes), block * _blockBytes);
        _blocksWritten++;
        _bytesWritten += _blockBytes;
        long end = (block + 1) * _blockBytes;
        if (end > _fileBytes)
        {
            ScratchTracker.Add(end - _fileBytes);
            _fileBytes = end;
        }
        if (_flushBytes <= 0) return;
        _unflushed += _blockBytes;
        if (_unflushed < _flushBytes) return;
        _file.Flush();
        _flushes++;
        _unflushed = 0;
    }

    private void ReadBlock(long block, T* buffer)
    {
        _file!.ReadExactAt(new Span<byte>(buffer, _blockBytes), block * _blockBytes);
        _blocksRead++;
        _bytesRead += _blockBytes;
    }
}
