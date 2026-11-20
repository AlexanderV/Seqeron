using System.IO.MemoryMappedFiles;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>Options of <see cref="PsiStream"/>.</summary>
internal sealed record PsiOptions
{
    /// <summary>Bytes of one output window (one per symbol of a pass); at least 1 MiB by default (the .NET view accessor flushes every window when it is disposed, ADR-004 section 5a, E9); rounded down to whole 4-byte entries, at least one entry.</summary>
    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    /// <summary>Symbols whose output cursors are open in one pass (ADR-004 section 3.3: 64); a larger alphabet takes <c>ceil(sigma / this)</c> passes over the BWT.</summary>
    internal int SymbolsPerPass { get; init; } = 64;

    /// <summary>Bytes of the BWT read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>Cap on (symbols per pass) x (window bytes), the mapped working set of a pass, charged to the budget like an arena block (ADR-004: arena cap 160 MiB).</summary>
    internal long WindowCapBytes { get; init; } = 160L << 20;

    /// <summary>The shared cap of the native buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Test hook: options of the arena (guard bands, poison pattern); null = production arena.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>The BWT reads and the mapped writes are counted here (null = not counted).</summary>
    internal IoCounters? Counters { get; init; }
}

/// <summary>What one PSI build did: passes over the BWT, symbols per pass, windows mapped, the bound of the mapped windows, BWT bytes read and PSI bytes written.</summary>
internal sealed record PsiStats(int Passes, int SymbolsPerPass, long Windows, long MappedBoundBytes, long BwtBytesRead, long PsiBytesWritten);

/// <summary>
/// PH6 of ADR-004 section 3.3: Psi from the BWT by sigma-way distribution, without the inverse permutation, written into the PSI section of the container. With <c>BWT[j] = T[SA[j] - 1]</c> and <c>C[c]</c> = the number
/// of suffixes that start with a symbol below c (the empty suffix included), <c>Psi[C[c] + k]</c> is the row of the k-th (0-based) occurrence of c in the BWT: one scan of the BWT in row order with one output cursor per
/// symbol, every cursor writing its own contiguous range of the section strictly sequentially through its own mapped window. At most <see cref="PsiOptions.SymbolsPerPass"/> cursors are open at a time. The row of the
/// suffix 0 has no symbol and is skipped; <c>PSI[0]</c> = 0 by convention. The result is <c>PSI[i] = ISA[SA[i] + 1]</c> for <c>SA[i] &lt; n</c> and 0 for <c>SA[i] = n</c>.
/// </summary>
internal static unsafe class PsiStream
{
    private static long s_liveBlocks;

    /// <summary>Native blocks (two per open pass) that the cursors of all builds hold right now: a leak is a count that does not return to where it was (nothing else in the process can see it; the blocks are raw <see cref="NativeMemory"/>).</summary>
    internal static long LiveNativeBlocks => Interlocked.Read(ref s_liveBlocks);

    /// <summary>
    /// Writes the <c>n + 1</c> little-endian <c>uint</c> entries of Psi into the range <c>[psiOffset, psiOffset + 4 (n + 1))</c> of <paramref name="file"/> (every byte of the range is written) from the BWT file
    /// <paramref name="bwtPath"/> (<c>n + 1</c> bytes), the row of the suffix 0 and the symbol counts of the text. The writers are disposed when it returns; the file is not completed.
    /// </summary>
    /// <exception cref="InvalidDataException">The BWT file has the wrong length, or a symbol occurs more or less often than its count.</exception>
    /// <exception cref="ArgumentException">The counts are negative or do not sum to n, or there are more than 256 symbols.</exception>
    internal static PsiStats Build(string bwtPath, long n, long primaryRow, long[] counts, MappedOutputFile file, long psiOffset, PsiOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(bwtPath);
        ArgumentNullException.ThrowIfNull(counts);
        ArgumentNullException.ThrowIfNull(file);
        options ??= new PsiOptions();
        ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.SymbolsPerPass, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ReadBytes, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        if (n > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        long rows = n + 1;
        if (primaryRow < 0 || primaryRow >= rows) throw new ArgumentOutOfRangeException(nameof(primaryRow), primaryRow, $"the row of the suffix 0 lies outside [0, {n}]");
        int sigma = counts.Length;
        if (sigma > 256) throw new ArgumentException($"the BWT is one byte per row: at most 256 symbols, {sigma} counts given", nameof(counts));
        if (counts.Any(c => c < 0) || counts.Sum() != n) throw new ArgumentException($"the symbol counts must be non-negative and sum to n = {n}", nameof(counts));
        ArgumentOutOfRangeException.ThrowIfNegative(psiOffset);
        if (psiOffset > file.Length - (sizeof(uint) * rows)) throw new ArgumentOutOfRangeException(nameof(psiOffset), psiOffset, $"the range of {sizeof(uint) * rows} bytes lies outside the file of {file.Length} bytes");

        int entriesPerWindow = Math.Max(1, options.WindowBytes / sizeof(uint));
        int group = Math.Min(options.SymbolsPerPass, sigma);
        int passes = sigma == 0 ? 0 : (sigma + options.SymbolsPerPass - 1) / options.SymbolsPerPass;
        long bound = group * (long)entriesPerWindow * sizeof(uint);
        new ArenaBudget(options.WindowCapBytes).Reserve(bound, "psi.windows");
        ArenaBudget budget = options.Budget ?? new ArenaBudget(long.MaxValue);
        budget.Reserve(bound, "psi.windows");
        try
        {
            using var bwt = new CountedFile(bwtPath, FileMode.Open, FileAccess.Read, options.Counters);
            if (bwt.Length != rows) throw new InvalidDataException($"the BWT file has {bwt.Length} bytes, expected n + 1 = {rows}");

            var start = new long[sigma + 1]; // entry index of the first Psi entry of every symbol; row 0, the empty suffix, comes first
            start[0] = 1;
            for (int c = 0; c < sigma; c++) start[c + 1] = start[c] + counts[c];

            using (MappedWindowWriter zero = file.CreateWriter(psiOffset, sizeof(uint), sizeof(uint), options.Counters)) zero.Append(new byte[sizeof(uint)]); // PSI[0] = 0 by convention
            long windows = 0;
            int readLen = (int)Math.Min(options.ReadBytes, rows);
            using var arena = new BudgetedArena(budget, options.ArenaGuards);
            Buf<byte> readBuf = arena.Allocate<byte>(readLen, zero: false, "psi.bwt");
            for (int pass = 0; pass < passes; pass++)
            {
                int lo = pass * options.SymbolsPerPass;
                int hi = Math.Min(sigma, lo + options.SymbolsPerPass);
                using var cursors = new Cursors(file, psiOffset, start, counts, lo, hi, entriesPerWindow, options.Counters);
                for (long first = 0; first < rows; first += readLen)
                {
                    int take = (int)Math.Min(readLen, rows - first);
                    Span<byte> block = readBuf.AsSpan(take);
                    bwt.ReadExactAt(block, first);
                    long skip = primaryRow - first; // the sentinel row carries no symbol
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
            return new PsiStats(passes, group, windows, bound, passes * rows, sizeof(uint) * rows);
        }
        finally
        {
            budget.Release(bound);
        }
    }

    /// <summary>The output cursors of one pass: one per symbol of [lo, hi), each with the mapped window it is writing into. Counts itself as a writer of the file until it is disposed.</summary>
    private sealed class Cursors : IDisposable
    {
        private readonly MappedOutputFile _file;
        private readonly long _origin;
        private readonly int _lo;
        private readonly uint _span;
        private readonly int _entriesPerWindow;
        private readonly IoCounters? _counters;
        private readonly long[] _next;
        private readonly long[] _limit;
        private readonly long[] _first;
        private readonly MemoryMappedViewAccessor?[] _views;
        private readonly nuint* _cur;
        private readonly nuint* _end;

        private static nuint* AllocateBlock(int entries)
        {
            nuint* block = (nuint*)NativeMemory.AllocZeroed((nuint)Math.Max(1, entries), (nuint)sizeof(nuint));
            Interlocked.Increment(ref s_liveBlocks);
            return block;
        }

        private static void FreeBlock(nuint* block)
        {
            NativeMemory.Free(block);
            Interlocked.Decrement(ref s_liveBlocks);
        }

        internal Cursors(MappedOutputFile file, long origin, long[] start, long[] counts, int lo, int hi, int entriesPerWindow, IoCounters? counters)
        {
            _file = file;
            _origin = origin;
            _lo = lo;
            _span = (uint)(hi - lo);
            _entriesPerWindow = entriesPerWindow;
            _counters = counters;
            _next = new long[hi - lo];
            _limit = new long[hi - lo];
            _first = new long[hi - lo];
            _views = new MemoryMappedViewAccessor?[hi - lo];
            for (int s = 0; s < hi - lo; s++)
            {
                _first[s] = _next[s] = start[lo + s];
                _limit[s] = start[lo + s] + counts[lo + s];
            }
            _cur = AllocateBlock(hi - lo);
            _end = AllocateBlock(hi - lo);
            file.WriterOpened();
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

        public void Dispose()
        {
            for (int s = 0; s < _views.Length; s++) ReleaseView(s);
            FreeBlock(_cur);
            FreeBlock(_end);
            _file.WriterClosed();
        }

        private uint* Advance(int s)
        {
            ReleaseView(s);
            long remaining = _limit[s] - _next[s];
            if (remaining <= 0) throw new InvalidDataException($"symbol {_lo + s} occurs more often in the BWT than its count {_limit[s] - _first[s]}");
            int entries = (int)Math.Min(_entriesPerWindow, remaining);
            long bytes = entries * (long)sizeof(uint);
            MemoryMappedViewAccessor view = _file.OpenView(_origin + (_next[s] * sizeof(uint)), bytes);
            _views[s] = view;
            byte* p = null;
            view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            uint* w = (uint*)(p + view.PointerOffset);
            _cur[s] = (nuint)w;
            _end[s] = (nuint)(w + entries);
            _next[s] += entries;
            Windows++;
            _counters?.AddMappedWrite(bytes);
            _counters?.AddWindow(bytes);
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
    }
}
