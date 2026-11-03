using System.Runtime.InteropServices;

namespace SuffixTree.Esa;

/// <summary>
/// Options of a <see cref="NativeArena"/>. Production uses <see cref="Production"/> (no guard bands, no poison, no overhead); the tests use
/// <see cref="Guarded"/>: canary guard bands around every block and a poison pattern in blocks that were not zero-initialised, so that a
/// kernel which writes outside its buffer or reads a cell it never wrote is detected.
/// </summary>
/// <param name="GuardBytes">Size of the guard band before and after every block (0 = none); a multiple of 16 keeps the payload aligned.</param>
/// <param name="Poison">Byte value that fills uninitialised blocks, or null to leave them as the allocator returns them.</param>
internal readonly record struct ArenaOptions(int GuardBytes, byte? Poison)
{
    /// <summary>No guard bands, no poison.</summary>
    internal static ArenaOptions Production => new(0, null);

    /// <summary>256-byte guard bands and the given poison byte (default 0xCD).</summary>
    internal static ArenaOptions Guarded(byte poison = 0xCD) => new(256, poison);
}

/// <summary>
/// Owner of native (unmanaged, <see cref="NativeMemory"/>) blocks used by the unsafe kernels for their workspace. Every block is handed out as a
/// <see cref="Buf{T}"/>; <see cref="Dispose"/> frees all blocks (a finalizer does the same if the owner forgot). With guard bands the memory
/// before and after each payload is filled with a position-dependent canary pattern and <see cref="VerifyCanaries"/> proves that nothing wrote
/// outside a payload; tests call it after every kernel call (safety infrastructure, review MAJ-1). Not thread-safe.
/// </summary>
internal sealed unsafe class NativeArena : IDisposable
{
    private readonly ArenaOptions _options;
    private readonly List<Block> _blocks = [];
    private bool _disposed;

    internal NativeArena()
        : this(ArenaOptions.Production)
    {
    }

    internal NativeArena(ArenaOptions options)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(options.GuardBytes);
        _options = options;
    }

    ~NativeArena() => FreeAll();

    /// <summary>Number of blocks allocated so far.</summary>
    internal int BlockCount => _blocks.Count;

    /// <summary>Payload bytes of all blocks allocated so far (guard bands not counted).</summary>
    internal long AllocatedBytes { get; private set; }

    /// <summary>True when blocks carry guard bands.</summary>
    internal bool Guarded => _options.GuardBytes > 0;

    /// <summary>
    /// Allocates <paramref name="length"/> elements. With <paramref name="zero"/> the block is zero-filled; otherwise its content is undefined
    /// (a poison pattern in a guarded arena, so reads before writes show up as wrong results in the tests).
    /// </summary>
    internal Buf<T> Allocate<T>(int length, bool zero, string name)
        where T : unmanaged
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        ArgumentException.ThrowIfNullOrEmpty(name);
        nuint payload = checked((nuint)length * (nuint)sizeof(T));
        nuint guard = (nuint)_options.GuardBytes;
        nuint total = checked(payload + (2 * guard));
        if (total == 0) total = 1; // a valid, non-null pointer even for an empty buffer
        bool zeroedByAllocator = zero && guard == 0; // calloc: fresh pages are zero without touching them
        byte* raw = (byte*)(zeroedByAllocator ? NativeMemory.AllocZeroed(total) : NativeMemory.Alloc(total));
        byte* body = raw + guard;
        if (zero && !zeroedByAllocator) NativeMemory.Clear(body, payload);
        else if (!zero && _options.Poison is byte poison) NativeMemory.Fill(body, payload, poison);
        var block = new Block(raw, body, payload, guard, name, _blocks.Count);
        if (guard > 0)
        {
            block.FillCanary(before: true);
            block.FillCanary(before: false);
        }
        _blocks.Add(block);
        AllocatedBytes += (long)payload;
        return new Buf<T>((T*)body, length, name);
    }

    /// <summary>Throws <see cref="EsaMemoryViolationException"/> (kind <see cref="MemoryViolationKind.CanaryOverwritten"/>) when any guard band was modified. Does nothing in an unguarded arena.</summary>
    internal void VerifyCanaries()
    {
        foreach (Block block in _blocks) block.Verify();
    }

    /// <summary>Frees all blocks; idempotent.</summary>
    public void Dispose()
    {
        FreeAll();
        GC.SuppressFinalize(this);
    }

    private void FreeAll()
    {
        if (_disposed) return;
        _disposed = true;
        foreach (Block block in _blocks) NativeMemory.Free(block.Raw);
        _blocks.Clear();
    }

    private readonly struct Block(byte* raw, byte* body, nuint payload, nuint guard, string name, int id)
    {
        internal byte* Raw { get; } = raw;

        private byte* Body { get; } = body;

        private nuint Payload { get; } = payload;

        private nuint Guard { get; } = guard;

        private string Name { get; } = name;

        private int Id { get; } = id;

        internal void FillCanary(bool before)
        {
            byte* start = before ? Raw : Body + Payload;
            for (nuint i = 0; i < Guard; i++) start[i] = Pattern(before, i);
        }

        internal void Verify()
        {
            if (Guard == 0) return;
            Check(before: true);
            Check(before: false);
        }

        private void Check(bool before)
        {
            byte* start = before ? Raw : Body + Payload;
            for (nuint i = 0; i < Guard; i++)
            {
                if (start[i] == Pattern(before, i)) continue;
                long offset = before ? (long)i - (long)Guard : (long)Payload + (long)i;
                throw new EsaMemoryViolationException(
                    MemoryViolationKind.CanaryOverwritten,
                    $"canary of block '{Name}' (#{Id}, payload {Payload} bytes) overwritten {(before ? "before" : "after")} the payload at byte offset {offset} (value 0x{start[i]:X2}, expected 0x{Pattern(before, i):X2})");
            }
        }

        private byte Pattern(bool before, nuint i)
            => (byte)(((int)(i * 167) + (Id * 31) + (before ? 0x35 : 0x5A)) ^ 0xA5 ^ (int)(i >> 3));
    }
}
