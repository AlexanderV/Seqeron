using System.IO.MemoryMappedFiles;
using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.External;

/// <summary>
/// Division and remainder of a 32-bit value by a fixed 32-bit divisor with two multiplications (Lemire, Kaser, Kurz, "Faster remainder by direct computation", 2019): the comparator of the difference-cover sample does five of
/// them per call and a 64-bit hardware division costs 25-90 cycles on the machines of this class.
/// </summary>
internal readonly struct FastMod32
{
    private readonly ulong _m;
    private readonly uint _d;

    internal FastMod32(uint divisor)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(divisor, 2u);
        _d = divisor;
        _m = (ulong.MaxValue / divisor) + 1;
    }

    internal uint Divisor => _d;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Div(uint a) => (uint)Math.BigMul(_m, a, out _);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Mod(uint a) => (uint)Math.BigMul(_m * a, _d, out _);
}

/// <summary>
/// The difference-cover sample ranks as PH3 reads them: for every position p of <c>[0, n]</c> with <c>p mod X</c> in D the rank of the suffix p among the sample suffixes (a dense <c>uint</c> per slot, slot =
/// <c>(p / X) * |D| + index of (p mod X) in D</c>, see <see cref="SampleGeometry"/>), in the rank file PH1 wrote, through a read-only memory mapping (random reads through the page cache: the cache floor of ADR-004 E8).
/// Opening checks the length of the file in O(1) and reports a file of another length as <see cref="InvalidDataException"/>; <see cref="Verify"/> is the O(S) check of the contents. Reads are thread-safe.
/// </summary>
internal sealed unsafe class RankTable : IDisposable
{
    private readonly FastMod32 _mod;
    private readonly int[] _indexOf;
    private readonly int _classes;
    private readonly uint* _ranks;
    private MemoryMappedFile? _map;
    private MemoryMappedViewAccessor? _view;
#if ESA_BOUNDS
    private string? _registered; // the checked build: the file that this table maps, in the list of the mapped files (MappedFileRegistry)
#endif

    private RankTable(DifferenceCover cover, long n, long slots, MemoryMappedFile map, MemoryMappedViewAccessor view, uint* ranks)
    {
        Cover = cover;
        N = n;
        Slots = slots;
        _mod = new FastMod32((uint)cover.Period);
        _indexOf = cover.IndexOf;
        _classes = cover.Residues.Length;
        _map = map;
        _view = view;
        _ranks = ranks;
    }

    internal DifferenceCover Cover { get; }

    /// <summary>Division and remainder by the period without a hardware division.</summary>
    internal FastMod32 Mod => _mod;

    internal long N { get; }

    /// <summary>Number of slots (including the unused slots of the last block, which lie beyond n).</summary>
    internal long Slots { get; }

    /// <summary>Maps the rank file of a text of <paramref name="n"/> symbols read-only.</summary>
    /// <exception cref="InvalidDataException">The file does not have <c>4 x Slots</c> bytes.</exception>
    internal static RankTable Open(string path, DifferenceCover cover, long n)
    {
        ArgumentNullException.ThrowIfNull(path);
        ArgumentNullException.ThrowIfNull(cover);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        long slots = ((n / cover.Period) + 1) * cover.Residues.Length;
        long bytes = slots * 4;
        using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);
        if (stream.Length != bytes) throw new InvalidDataException($"the rank file has {stream.Length} bytes but a text of {n} symbols with the period {cover.Period} needs {bytes} ({slots} slots)");
        MemoryMappedFile map = MemoryMappedFile.CreateFromFile(stream, null, 0, MemoryMappedFileAccess.Read, HandleInheritability.None, leaveOpen: true);
        MemoryMappedViewAccessor? view = null;
        try
        {
            view = map.CreateViewAccessor(0, bytes, MemoryMappedFileAccess.Read);
            byte* p = null;
            view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            var table = new RankTable(cover, n, slots, map, view, (uint*)(p + view.PointerOffset));
#if ESA_BOUNDS
            table._registered = MappedFileRegistry.Map(path);
#endif
            return table;
        }
        catch
        {
            view?.Dispose();
            map.Dispose();
            throw;
        }
    }

    /// <summary>
    /// The rank of the sample suffix at <paramref name="p"/> (a position of <c>[0, n]</c> whose residue is in D; a position outside the sample reads an arbitrary slot). In the checked build (<c>ESA_BOUNDS</c>, <c>-p:EsaBounds=true</c>) a
    /// position below 0 or above the last block, a residue outside D and a read after <see cref="Dispose"/> throw <see cref="EsaMemoryViolationException"/> where the production build reads through the raw pointer.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Rank(long p)
    {
        uint block = _mod.Div((uint)p);
        uint r = (uint)p - (block * _mod.Divisor);
#if ESA_BOUNDS
        return ReadChecked(p, block, r);
#else
        return _ranks[((long)block * _classes) + _indexOf[r]];
#endif
    }

#if ESA_BOUNDS
    private uint ReadChecked(long p, uint block, uint r)
    {
        int index = _indexOf[r];
        long slot = ((long)block * _classes) + index;
        if ((ulong)p > (ulong)N || index < 0 || (ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(p, r, slot);
        return _ranks[slot];
    }

    private uint ReadSlot(long slot)
    {
        if ((ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(-1, 0, slot);
        return _ranks[slot];
    }

    private void ThrowOutOfBounds(long p, uint r, long slot)
        => throw new EsaMemoryViolationException(
            MemoryViolationKind.OutOfBounds,
            $"out-of-bounds read of the rank table: position {p} (residue {r}) of a text of {N} symbols, slot {slot} of {Slots}{(_view is null ? "; the table is disposed" : string.Empty)}");
#endif

    /// <summary>Reads every slot once: the ranks of the real slots are distinct values in <c>[1, S + |D|]</c> and the unused slots of the last block are zero.</summary>
    /// <exception cref="InvalidDataException">A slot holds something else.</exception>
    internal void Verify()
    {
        ObjectDisposedException.ThrowIf(_view is null, this);
        long samples = Cover.SampleCount(N);
        long top = samples + _classes;
        var seen = new ulong[(top / 64) + 1];
        for (long slot = 0; slot < Slots; slot++)
        {
#if ESA_BOUNDS
            uint v = ReadSlot(slot);
#else
            uint v = _ranks[slot];
#endif
            if (slot >= samples)
            {
                if (v != 0) throw new InvalidDataException($"the unused slot {slot} holds {v}, not 0");
                continue;
            }
            if (v == 0 || v > top) throw new InvalidDataException($"slot {slot} holds the rank {v}, outside [1, {top}]");
            ulong bit = 1UL << (int)(v & 63);
            if ((seen[v >> 6] & bit) != 0) throw new InvalidDataException($"the rank {v} occurs twice (slot {slot})");
            seen[v >> 6] |= bit;
        }
    }

    public void Dispose()
    {
        MemoryMappedViewAccessor? view = Interlocked.Exchange(ref _view, null);
        MemoryMappedFile? map = Interlocked.Exchange(ref _map, null);
#if ESA_BOUNDS
        if (view is not null && _registered is not null) MappedFileRegistry.Unmap(_registered);
#endif
        if (view is not null)
        {
            view.SafeMemoryMappedViewHandle.ReleasePointer();
            view.Dispose();
        }
        map?.Dispose();
    }
}
