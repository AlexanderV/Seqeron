using System.IO.MemoryMappedFiles;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.Esa;

namespace SuffixTree.EmProto.Ph6;

internal sealed record PsiOptions
{
    /// <summary>Bytes of one output window (one per symbol of the pass); rounded down to whole 4-byte entries, at least one entry.</summary>
    internal int WindowBytes { get; init; } = 64 * 1024;

    /// <summary>Symbols whose output cursors are open in one pass (the ADR: 64); a larger alphabet takes ceil(sigma / this) passes over the BWT.</summary>
    internal int SymbolsPerPass { get; init; } = 64;

    /// <summary>Bytes of the BWT read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>PSI[0]: false = 0 (the convention of ADR-003 / ADR-004), true = the row of the suffix 0 (the cyclic Psi, a single cycle through all n + 1 rows).</summary>
    internal bool CyclicEntryZero { get; init; }

    /// <summary>Cap on (symbols per pass) x (window bytes), the mapped working set of the pass, charged to the budget like an arena block (ADR-004: arena cap 160 MiB).</summary>
    internal long WindowCapBytes { get; init; } = 160L << 20;

    internal ArenaBudget? Budget { get; init; }
}

internal sealed record PsiStats(int Passes, int SymbolsPerPass, long Windows, long MappedBoundBytes, long BwtBytesRead, long PsiBytesWritten);

/// <summary>
/// PH6 of ADR-004 section 3.3: Psi from the BWT by sigma-way distribution, without the inverse permutation. With <c>BWT[j] = T[SA[j] - 1]</c> and <c>C[c]</c> = the number of suffixes that start with a symbol below c (the empty suffix
/// included), <c>Psi[C[c] + k]</c> is the row of the k-th (0-based) occurrence of c in the BWT: one scan of the BWT in row order with one output cursor per symbol, every cursor writing its own contiguous range of Psi
/// strictly sequentially through a small mapped window. At most <see cref="PsiOptions.SymbolsPerPass"/> cursors are open at a time; a larger alphabet takes one scan of the BWT per group of symbols. The row of the
/// suffix 0 has no symbol (the sentinel) and is skipped; PSI[0] is 0 by the convention of ADR-003, or that row in the cyclic variant.
/// </summary>
internal static unsafe class PsiStream
{
    internal static PsiStats Build(string bwtPath, BwtMeta meta, string psiPath, PsiOptions? options = null)
    {
        PsiOptions o = options ?? new PsiOptions();
        ArgumentNullException.ThrowIfNull(meta);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.WindowBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.SymbolsPerPass, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.ReadBytes, 1);
        long n = meta.N;
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(meta), n, "negative text length");
        if (n >= uint.MaxValue - 2) throw new NotSupportedException("positions are 32-bit in the prototype");
        long rows = n + 1;
        if (meta.PrimaryRow < 0 || meta.PrimaryRow >= rows) throw new ArgumentOutOfRangeException(nameof(meta), meta.PrimaryRow, $"the row of the suffix 0 lies outside [0, {n}]");
        long[] counts = meta.Counts;
        int sigma = counts.Length;
        if (sigma > 256) throw new ArgumentOutOfRangeException(nameof(meta), sigma, "the BWT is one byte per row: at most 256 symbols");
        if (counts.Any(c => c < 0) || counts.Sum() != n) throw new ArgumentException($"the symbol counts must be non-negative and sum to n = {n}", nameof(meta));

        int entriesPerWindow = Math.Max(1, o.WindowBytes / 4);
        int group = Math.Min(o.SymbolsPerPass, sigma);
        int passes = sigma == 0 ? 0 : (sigma + o.SymbolsPerPass - 1) / o.SymbolsPerPass;
        long bound = group * (long)entriesPerWindow * 4;
        new ArenaBudget(o.WindowCapBytes).Reserve(bound, "psi.windows");
        ArenaBudget budget = o.Budget ?? new ArenaBudget(long.MaxValue);
        budget.Reserve(bound, "psi.windows");
        try
        {
            using var bwt = new CountedFile(bwtPath, FileMode.Open, FileAccess.Read);
            if (bwt.Length != rows) throw new InvalidDataException($"the BWT file has {bwt.Length} bytes, expected n + 1 = {rows}");

            var start = new long[sigma + 1]; // entry index of the first Psi entry of every symbol; row 0, the empty suffix, comes first
            start[0] = 1;
            for (int c = 0; c < sigma; c++) start[c + 1] = start[c] + counts[c];

            using (var fs = new FileStream(psiPath, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite)) fs.SetLength(4 * rows);
            long windows = 0;
            using var map = MemoryMappedFile.CreateFromFile(psiPath, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);
            WriteEntryZero(map, o.CyclicEntryZero ? (uint)meta.PrimaryRow : 0u);

            int readLen = (int)Math.Min(o.ReadBytes, rows);
            using var arena = new BudgetedArena(budget);
            Buf<byte> readBuf = arena.Allocate<byte>(readLen, zero: false, "psi.bwt");
            for (int pass = 0; pass < passes; pass++)
            {
                int lo = pass * o.SymbolsPerPass;
                int hi = Math.Min(sigma, lo + o.SymbolsPerPass);
                using var cursors = new Cursors(map, start, counts, lo, hi, entriesPerWindow);
                for (long first = 0; first < rows; first += readLen)
                {
                    int take = (int)Math.Min(readLen, rows - first);
                    Span<byte> block = readBuf.AsSpan(take);
                    bwt.ReadExactAt(block, first);
                    long skip = meta.PrimaryRow - first; // the sentinel row carries no symbol
                    fixed (byte* b = block)
                    {
                        if (skip >= 0 && skip < take)
                        {
                            cursors.Scan(b, (int)skip, first);
                            cursors.Scan(b + skip + 1, take - (int)skip - 1, first + skip + 1);
                        }
                        else
                        {
                            cursors.Scan(b, take, first);
                        }
                    }
                }
                cursors.Finish();
                windows += cursors.Windows;
            }
            return new PsiStats(passes, group, windows, bound, passes * rows, 4 * rows);
        }
        finally
        {
            budget.Release(bound);
        }
    }

    private static void WriteEntryZero(MemoryMappedFile map, uint value)
    {
        using MemoryMappedViewAccessor view = map.CreateViewAccessor(0, 4, MemoryMappedFileAccess.Write);
        byte* p = null;
        view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        try
        {
            *(uint*)(p + view.PointerOffset) = value;
        }
        finally
        {
            view.SafeMemoryMappedViewHandle.ReleasePointer();
        }
        MappedIo.AddWrite(4);
        MappedIo.AddWindow(4);
    }

    /// <summary>The output cursors of one pass: one per symbol of [lo, hi), each with the mapped window it is writing into.</summary>
    private sealed class Cursors : IDisposable
    {
        private readonly MemoryMappedFile _map;
        private readonly int _lo;
        private readonly uint _span;
        private readonly int _entriesPerWindow;
        private readonly long[] _next;
        private readonly long[] _limit;
        private readonly long[] _first;
        private readonly MemoryMappedViewAccessor?[] _views;
        private readonly nuint* _cur;
        private readonly nuint* _end;

        internal Cursors(MemoryMappedFile map, long[] start, long[] counts, int lo, int hi, int entriesPerWindow)
        {
            _map = map;
            _lo = lo;
            _span = (uint)(hi - lo);
            _entriesPerWindow = entriesPerWindow;
            _next = new long[hi - lo];
            _limit = new long[hi - lo];
            _first = new long[hi - lo];
            _views = new MemoryMappedViewAccessor?[hi - lo];
            for (int s = 0; s < hi - lo; s++)
            {
                _first[s] = _next[s] = start[lo + s];
                _limit[s] = start[lo + s] + counts[lo + s];
            }
            _cur = (nuint*)System.Runtime.InteropServices.NativeMemory.AllocZeroed((nuint)Math.Max(1, hi - lo), (nuint)sizeof(nuint));
            _end = (nuint*)System.Runtime.InteropServices.NativeMemory.AllocZeroed((nuint)Math.Max(1, hi - lo), (nuint)sizeof(nuint));
        }

        /// <summary>Windows opened by the cursors of this pass.</summary>
        internal long Windows { get; private set; }

        /// <summary>Distributes <paramref name="length"/> BWT bytes whose first row is <paramref name="firstRow"/>.</summary>
        internal void Scan(byte* bytes, int length, long firstRow)
        {
            nuint* cur = _cur, end = _end;
            uint span = _span;
            int lo = _lo;
            for (int j = 0; j < length; j++)
            {
                uint s = (uint)(bytes[j] - lo);
                if (s >= span) continue; // a symbol of another pass
                uint* p = (uint*)cur[s];
                if (p == (uint*)end[s]) p = Advance((int)s);
                *p = (uint)(firstRow + j);
                cur[s] = (nuint)(p + 1);
            }
        }

        /// <summary>Every cursor must have used its whole range: the BWT has exactly the counted number of occurrences of every symbol of the pass.</summary>
        internal void Finish()
        {
            for (int s = 0; s < _next.Length; s++)
            {
                if (_next[s] != _limit[s] || _cur[s] != _end[s])
                {
                    long written = _next[s] - ((uint*)_end[s] - (uint*)_cur[s]) - _first[s];
                    throw new InvalidDataException($"symbol {_lo + s} occurs {written} times in the BWT, its count is {_limit[s] - _first[s]}");
                }
            }
        }

        private uint* Advance(int s)
        {
            ReleaseView(s);
            long remaining = _limit[s] - _next[s];
            if (remaining <= 0) throw new InvalidDataException($"symbol {_lo + s} occurs more often in the BWT than its count {_limit[s] - _first[s]}");
            int entries = (int)Math.Min(_entriesPerWindow, remaining);
            long bytes = entries * 4L;
            MemoryMappedViewAccessor view = _map.CreateViewAccessor(_next[s] * 4, bytes, MemoryMappedFileAccess.Write);
            _views[s] = view;
            byte* p = null;
            view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            uint* w = (uint*)(p + view.PointerOffset);
            _cur[s] = (nuint)w;
            _end[s] = (nuint)(w + entries);
            _next[s] += entries;
            Windows++;
            MappedIo.AddWrite(bytes);
            MappedIo.AddWindow(bytes);
            return w;
        }

        private void ReleaseView(int s)
        {
            MemoryMappedViewAccessor? view = _views[s];
            if (view is null) return;
            view.SafeMemoryMappedViewHandle.ReleasePointer();
            view.Dispose();
            _views[s] = null;
        }

        public void Dispose()
        {
            for (int s = 0; s < _views.Length; s++) ReleaseView(s);
            System.Runtime.InteropServices.NativeMemory.Free(_cur);
            System.Runtime.InteropServices.NativeMemory.Free(_end);
        }
    }
}
