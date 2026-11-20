using System.Diagnostics;
using System.Security.Cryptography;

namespace SuffixTree.Esa.External;

/// <summary>
/// The memory and scratch limits of the structure pass (<see cref="SesaVerifier.CheckStructure"/>): how many text positions one pass over the suffix array serves, how many open intervals the child-table scan keeps in memory
/// before it spills the older ones to a temporary file, and where that file goes. The defaults are two 16 MiB arrays for a chunk, a stack array of 64 MiB of which only as many pages are touched as the tree is deep (500,000 entries, 8 MB, on the 10^8 texts of G2) and the permutation bitmap
/// (n / 8 bytes): the anonymous memory of the pass was 33.6 and 41.6 MB of buffers and a RssAnon peak of about 105 MB at 10^8; the tests set small values to cross every boundary on small texts.
/// </summary>
internal static class SesaVerifierTuning
{
    /// <summary>Text positions per chunk of the Phi pass: two <c>int</c> arrays of this many entries (8 bytes per position).</summary>
    internal static int ChunkSymbols { get; set; } = DefaultChunkSymbols;

    /// <summary>The default chunk: 4 Mi positions, two 16 MiB arrays.</summary>
    internal const int DefaultChunkSymbols = 1 << 22;

    /// <summary>Open intervals of the child-table scan kept in memory (16 bytes each); the older ones go to the spill file.</summary>
    internal static int StackEntries { get; set; } = 1 << 22;

    /// <summary>The directory of the spill file; null means the system's temporary directory.</summary>
    internal static string? SpillDirectory { get; set; }

    /// <summary>Bytes of a section that one call of the CRC takes (a span has an <c>int</c> length): a larger section is chained in pieces of this size; the default is 1 GiB, the tests set a few bytes.</summary>
    internal static long CrcPieceBytes { get; set; } = DefaultCrcPieceBytes;

    /// <summary>The default piece of the checksum pass: 1 GiB.</summary>
    internal const long DefaultCrcPieceBytes = 1L << 30;

    /// <summary>Restores the defaults.</summary>
    internal static void Reset()
    {
        ChunkSymbols = DefaultChunkSymbols;
        StackEntries = 1 << 22;
        SpillDirectory = null;
        CrcPieceBytes = DefaultCrcPieceBytes;
    }
}

/// <summary>
/// The section and structure levels of the validation (ADR-003 section 6, work package 1b.7f). <see cref="CheckSections"/> is the CRC-32C of every payload that has one, in one sequential pass over the mapping, and the
/// rule that every byte outside the payloads is zero. <see cref="CheckStructure"/> is the O(n) pass that proves the arrays agree with each other and with the text, with bounded working memory: the suffix array is a
/// permutation of 0..n with <c>SA[0] = n</c> (a bitmap of n + 1 bits), every adjacent pair is in strictly ascending order (compared from the stored LCP on, which a damaged order cannot satisfy), the LCP array equals what the Phi pass
/// recomputes from the suffix array and the text (Kasai's amortised scan in text order, the predecessor of every position found by sequential scans of the suffix array, one per chunk of the text), the child table equals the one the
/// LCP array defines (the lcp-interval scan with a stack that spills to a temporary file when it is deeper than the memory limit, every cell compared, the cells the scan does not define counted to be zero), the header
/// statistics are recomputed, the Psi array equals <c>ISA[SA[i] + 1]</c> (checked from the ranks of the Phi pass), every packed symbol is below sigma and the symbol counts equal a recount, and BUILD_INFO has the lines it claims.
/// An adversarial container that is consistent in every check cannot be told from a valid one; the threat model of this level is accidental damage and the damage that the tests of the reader inject (a file from an untrusted
/// source stays memory safe in <see cref="EsaAccessMode.Checked"/> mode either way).
/// </summary>
internal static unsafe class SesaVerifier
{
    private const int TextBlock = 1 << 16;
    private const int VerifiedCarry = 64;

    /// <summary>The CRC-32C of every section payload that has one, and the zero bytes between the sections.</summary>
    internal static void CheckSections(SesaFile file, bool? hardwareCrc)
    {
        long t0 = Stopwatch.GetTimestamp();
        bool hardware = hardwareCrc ?? Crc32C.HardwareAccelerated;
        long bytes = 0;
        foreach (SesaSection s in file.Sections)
        {
            if ((s.Flags & SesaFormat.FlagHasCrc) == 0) continue;
            uint crc = Checksum(file.Data + s.Offset, s.Length, hardware);
            bytes += s.Length;
            if (crc != s.Crc) throw new InvalidDataException($"section {Name(s.Kind)} has CRC-32C {crc:X8} but the section table says {s.Crc:X8}: the file is damaged");
        }
        long end = SesaFormat.PageBytes;
        foreach (SesaSection s in file.Sections.OrderBy(x => x.Offset).ThenBy(x => x.Length))
        {
            if (s.Offset > end && new ReadOnlySpan<byte>(file.Data + end, checked((int)(s.Offset - end))).IndexOfAnyExcept((byte)0) >= 0)
                throw new InvalidDataException($"the bytes between offset {end} and the section {Name(s.Kind)} at {s.Offset} are not zero: the file is damaged");
            end = Math.Max(end, s.Offset + s.Length);
        }
        file.SetReport(file.Report with { SectionsSeconds = Seconds(t0), SectionBytesChecksummed = bytes, HardwareCrc = hardware });
    }

    private static uint Checksum(byte* data, long length, bool hardware)
    {
        long piece = Math.Clamp(SesaVerifierTuning.CrcPieceBytes, 1L, SesaVerifierTuning.DefaultCrcPieceBytes);
        uint crc = 0;
        for (long done = 0; done < length; done += piece)
        {
            long take = Math.Min(piece, length - done);
            var span = new ReadOnlySpan<byte>(data + done, (int)take);
            crc = hardware ? Crc32C.AppendHardware(crc, span) : Crc32C.AppendTable(crc, span);
        }
        return crc;
    }

    private static string Name(int kind) => kind switch
    {
        SesaFormat.KindText => "TEXT",
        SesaFormat.KindAlphabet => "ALPHABET",
        SesaFormat.KindSa => "SA",
        SesaFormat.KindLcp => "LCP",
        SesaFormat.KindLcpOverflow => "LCP_OVERFLOW",
        SesaFormat.KindChild => "CHILD",
        SesaFormat.KindBuildInfo => "BUILD_INFO",
        SesaFormat.KindPsi => "PSI",
        SesaFormat.KindCounts => "COUNTS",
        _ => $"kind {kind}",
    };

    /// <summary>The structure pass.</summary>
    /// <exception cref="NotSupportedException">The text has 2^31 - 2 symbols or more (the 32-bit surface; S2).</exception>
    internal static void CheckStructure(SesaFile file)
    {
        long t0 = Stopwatch.GetTimestamp();
        SesaParsed ps = file.Parsed;
        long n64 = ps.Header.TextLength;
        RequireSupportedLength(n64);
        int n = (int)n64;
        long bitmapBytes = ((n + 1L + 63) >> 6) << 3;
        long buffers = 0;

        CheckText(file, n);
        if (n > 0)
        {
            var bits = new ulong[(n + 1 + 63) >> 6];
            CheckSuffixArray(file, n, bits);
            int chunk = Math.Max(1, Math.Min(n, SesaVerifierTuning.ChunkSymbols));
            buffers += 8L * chunk;
            CheckLcpAndPsi(file, n, chunk);
            Array.Clear(bits);
            buffers += CheckChildAndStatistics(file, n, bits);
        }
        else
        {
            CheckEmpty(file);
        }
        CheckBuildInfo(file);
        file.SetReport(file.Report with { StructureSeconds = Seconds(t0), BitmapBytes = bitmapBytes, BufferBytes = buffers + TextBlock });
    }

    /// <summary>The 32-bit structure pass takes at most 2^31 - 3 symbols (the arrays are indexed by <c>int</c>, the bitmap and the chunks by <c>int</c> too); a longer text is S2's.</summary>
    /// <exception cref="NotSupportedException">The text has 2^31 - 2 symbols or more.</exception>
    internal static void RequireSupportedLength(long symbols)
    {
        if (symbols >= int.MaxValue - 1L) throw new NotSupportedException($"the structure pass of a text of {symbols} symbols needs 64-bit indexes (S2); the 32-bit pass takes at most {int.MaxValue - 2} symbols");
    }

    private static void CheckEmpty(SesaFile file)
    {
        SesaParsed ps = file.Parsed;
        int* sa = (int*)(file.Data + ps.Sa.Offset);
        if (sa[0] != 0) throw new InvalidDataException("the suffix array of an empty text must be [0]: the file is damaged");
        if (ps.Psi is { } psi && (*(uint*)(file.Data + psi.Offset) != 0 || psi.Param0 != 0)) throw new InvalidDataException("the PSI array of an empty text must be [0] with the primary row 0: the file is damaged");
        if (file.Lcp[0] != 0) throw new InvalidDataException("LCP[0] must be 0: the file is damaged");
        if (*(int*)(file.Data + ps.Child.Offset) != 0) throw new InvalidDataException("the child table of an empty text must be [0]: the file is damaged");
    }

    // ---------------------------------------------------------------- text

    private static void CheckText(SesaFile file, int n)
    {
        SesaParsed ps = file.Parsed;
        int sigma = ps.Header.Sigma;
        var seen = new long[sigma];
        var block = new byte[TextBlock];
        bool hash = file.Options.VerifyTextSha256 && ps.Header.TextSha256.AsSpan().IndexOfAnyExcept((byte)0) >= 0;
        using IncrementalHash? sha = hash ? IncrementalHash.CreateHash(HashAlgorithmName.SHA256) : null;
        MappedTextReader text = file.TextReader;
        for (long at = 0; at < n; at += block.Length)
        {
            int take = (int)Math.Min(block.Length, n - at);
            text.ReadSymbols(at, block.AsSpan(0, take));
            for (int k = 0; k < take; k++)
            {
                byte s = block[k];
                if (s >= sigma) throw new InvalidDataException($"the text holds the symbol {s} at position {at + k}, which is not below SIGMA = {sigma}: the file is damaged");
                seen[s]++;
                if (sha is not null) block[k] = (byte)ps.Letters[s];
            }
            sha?.AppendData(block.AsSpan(0, take));
        }
        for (int s = 0; s < sigma; s++)
        {
            if (seen[s] == 0) throw new InvalidDataException($"the letter {s} of the alphabet does not occur in the text: the file is damaged");
            if (ps.Counts64 is { } counts && counts[s] != (ulong)seen[s]) throw new InvalidDataException($"COUNTS says the symbol {s} occurs {counts[s]} times, the text has {seen[s]}: the file is damaged");
        }
        if (sha is not null && !sha.GetHashAndReset().AsSpan().SequenceEqual(ps.Header.TextSha256)) throw new InvalidDataException("TEXT_SHA256 differs from the SHA-256 of the text: the file is damaged");
    }

    private static void CheckBuildInfo(SesaFile file)
    {
        if (file.Parsed.Info is not { } info) return;
        ReadOnlySpan<byte> payload = file.Payload(info);
        if (!System.Text.Unicode.Utf8.IsValid(payload)) throw new InvalidDataException("BUILD_INFO is not UTF-8: the file is damaged");
        if (payload.Count((byte)'\n') != info.Count) throw new InvalidDataException($"BUILD_INFO has {payload.Count((byte)'\n')} lines but its entry says {info.Count}: the file is damaged");
    }

    // ---------------------------------------------------------------- suffix array: permutation, order

    private static void CheckSuffixArray(SesaFile file, int n, ulong[] bits)
    {
        SesaParsed ps = file.Parsed;
        int* sa = (int*)(file.Data + ps.Sa.Offset);
        MappedTextReader text = file.TextReader;
        CompactLcp lcp = file.Lcp;
        if (sa[0] != n) throw new InvalidDataException($"SA[0] is {sa[0]}, not n = {n}: the suffix array must start with the empty suffix");
        if (lcp[0] != 0 || lcp[1] != 0) throw new InvalidDataException($"LCP[0] and LCP[1] must be 0, they are {lcp[0]} and {lcp[1]}: the file is damaged");
        bits[n >> 6] |= 1UL << (n & 63);
        CompactLcp.Cursor cursor = lcp.CreateCursor();
        long steps = 0;
        int previous = n;
        for (int r = 1; r <= n; r++)
        {
            int b = sa[r];
            if ((uint)b >= (uint)n) throw new InvalidDataException($"SA[{r}] is {b}: a position of the text is 0..{n - 1}, and {n} is the empty suffix at rank 0 only: the file is damaged");
            ref ulong word = ref bits[b >> 6];
            ulong mask = 1UL << (b & 63);
            if ((word & mask) != 0) throw new InvalidDataException($"the position {b} occurs twice in the suffix array (second time at rank {r}): it is not a permutation, the file is damaged");
            word |= mask;
            int l = cursor.Get(r);
            if (l < 0 || l > n - (previous == n ? 0 : Math.Max(previous, b))) throw new InvalidDataException($"LCP[{r}] is {l}, longer than the shorter of the two suffixes: the file is damaged");
            if (previous == n)
            {
                if (l != 0) throw new InvalidDataException($"LCP[{r}] is {l}: the empty suffix shares nothing with the first suffix");
            }
            else if (text.CompareSuffixes(previous, b, l, ref steps, out long m) >= 0 || m != l)
            {
                throw new InvalidDataException($"the suffixes at ranks {r - 1} and {r} (positions {previous} and {b}) are not in strictly ascending order after their {l} common symbols: the suffix array or the LCP array is damaged");
            }
            previous = b;
        }
    }

    // ---------------------------------------------------------------- LCP by the Phi pass, Psi

    private static void CheckLcpAndPsi(SesaFile file, int n, int chunk)
    {
        SesaParsed ps = file.Parsed;
        int* sa = (int*)(file.Data + ps.Sa.Offset);
        uint* psi = file.PsiPointer;
        MappedTextReader text = file.TextReader;
        CompactLcp lcp = file.Lcp;
        var phi = new int[chunk];
        var rank = new int[chunk];
        long steps = 0;
        int carry = 0;
        int previousRank = 0;
        if (psi is not null && psi[0] != 0) throw new InvalidDataException($"PSI[0] is {psi[0]}: the empty suffix has no successor and the entry must be 0: the file is damaged");
        for (int c0 = 0; c0 < n; c0 += chunk)
        {
            int c1 = (int)Math.Min((long)n, (long)c0 + chunk);
            int span = c1 - c0;
            for (int r = 1; r <= n; r++)
            {
                int p = sa[r];
                if ((uint)(p - c0) < (uint)span)
                {
                    phi[p - c0] = sa[r - 1];
                    rank[p - c0] = r;
                }
            }
            for (int p = c0; p < c1; p++)
            {
                int f = phi[p - c0];
                int r = rank[p - c0];
                int h;
                if (f == n)
                {
                    h = 0;
                }
                else if (carry <= VerifiedCarry)
                {
                    h = (int)text.Lce(p, f, long.MaxValue, ref steps);
                }
                else
                {
                    if (text.Lce(p, f, VerifiedCarry, ref steps) < VerifiedCarry) throw new InvalidDataException($"the suffixes at {f} and {p} (ranks {r - 1} and {r}) do not share the {VerifiedCarry} symbols that the scan in text order carries: the suffix array is damaged");
                    h = carry + (int)text.Lce(p + carry, f + carry, long.MaxValue, ref steps);
                }
                int stored = lcp[r];
                if (stored != h) throw new InvalidDataException($"LCP[{r}] is {stored}, but the suffixes at {f} and {p} that are adjacent in the suffix array have {h} symbols in common: the file is damaged");
                carry = h > 0 ? h - 1 : 0;
                if (psi is not null)
                {
                    if (p == 0)
                    {
                        if (ps.Psi!.Param0 != r) throw new InvalidDataException($"the PSI primary row is {ps.Psi.Param0} but the suffix 0 has rank {r}: the file is damaged");
                    }
                    else if (psi[previousRank] != (uint)r)
                    {
                        throw new InvalidDataException($"PSI[{previousRank}] is {psi[previousRank]}, but the suffix after the one of rank {previousRank} has rank {r}: the file is damaged");
                    }
                }
                previousRank = r;
            }
        }
        if (psi is not null && psi[previousRank] != 0) throw new InvalidDataException($"PSI[{previousRank}] is {psi[previousRank]}: the successor of the last suffix is the empty suffix at rank 0: the file is damaged");
    }

    // ---------------------------------------------------------------- child table and statistics

    /// <summary>Returns the bytes of memory of the scan's stack (peak).</summary>
    private static long CheckChildAndStatistics(SesaFile file, int n, ulong[] written)
    {
        SesaParsed ps = file.Parsed;
        int* cld = (int*)(file.Data + ps.Child.Offset);
        var sink = new ChildCompareSink(cld, n, written, ps.Child.Codec == SesaFormat.ChildCodecTagged);
        long stackBytes;
        using (var stack = new FrameStack(Math.Max(2, SesaVerifierTuning.StackEntries), SesaVerifierTuning.SpillDirectory))
        {
            Scan(file.Lcp, n, ref sink, stack);
            stackBytes = stack.PeakMemoryBytes;
        }
        long nonZero = CountNonZero(cld, n + 1);
        if (nonZero != sink.Writes) throw new InvalidDataException($"the child table has {nonZero} cells in use but the LCP array defines {sink.Writes}: the file is damaged");
        SesaHeader h = ps.Header;
        long nodeCount = 1 + sink.Internal + (n + 1L);
        if (h.NodeCount != nodeCount) throw new InvalidDataException($"NODE_COUNT is {h.NodeCount} but the lcp-interval tree has {nodeCount} nodes: the file is damaged");
        if (h.LrsDepth != sink.LrsDepth) throw new InvalidDataException($"LRS_DEPTH is {h.LrsDepth} but the deepest internal node has depth {sink.LrsDepth}: the file is damaged");
        long expectedRank = sink.LrsDepth == 0 ? -1 : sink.LrsRank;
        if (h.LrsRank != expectedRank) throw new InvalidDataException($"LRS_RANK is {h.LrsRank} but the first deepest interval starts at rank {expectedRank}: the file is damaged");
        return stackBytes;
    }

    private static long CountNonZero(int* cells, long count)
    {
        long total = 0;
        const long piece = 1L << 28;
        for (long done = 0; done < count; done += piece)
        {
            var span = new ReadOnlySpan<int>(cells + done, (int)Math.Min(piece, count - done));
            total += span.Length - span.Count(0);
        }
        return total;
    }

    /// <summary>
    /// The bottom-up scan of <see cref="LcpIntervalScan"/> with a stack that may spill; the same events in the same order (the structure and the sinks are those of the scan the builder uses, this copy exists because the original keeps its stack in
    /// four arrays that grow with the depth of the tree, 16 bytes per level, and the tree of a^n has n levels).
    /// </summary>
    private static void Scan(CompactLcp lcp, int n, ref ChildCompareSink sink, FrameStack stack)
    {
        CompactLcp.Cursor cursor = lcp.CreateCursor();
        stack.Push(new Frame(0, 0, -1, -1));
        int next = 1;
        for (int i = 1; i <= n + 1; i++)
        {
            int cur = i <= n ? cursor.Get(next++) : -1;
            if (i <= n && cur < 0) throw new InvalidDataException($"LCP[{i}] is negative: the file is damaged");
            int left = i - 1;
            while (stack.Count > 0 && stack.Top.Ell > cur)
            {
                Frame f = stack.Pop();
                int below = stack.Count > 0 ? stack.Top.Ell : -1;
                sink.Close(f.Ell, f.Lb, i - 1, f.First, stack.Count == 0 || cur < below);
                left = f.Lb;
            }
            if (i > n) break;
            ref Frame top = ref stack.TopRef;
            if (top.Ell < cur)
            {
                stack.Push(new Frame(cur, left, i, i));
                sink.LIndex(-1, i);
            }
            else if (top.First < 0)
            {
                top.First = i;
                top.Last = i;
                sink.LIndex(-1, i);
            }
            else
            {
                sink.LIndex(top.Last, i);
                top.Last = i;
            }
        }
    }

    private struct Frame(int ell, int lb, int first, int last)
    {
        internal int Ell = ell;
        internal int Lb = lb;
        internal int First = first;
        internal int Last = last;
    }

    /// <summary>The expected writes of the layout of the section's codec (packed: the raw value, a next link only when it is not <c>k + 1</c>; tagged: the kind in bits 30-31, a next link for every l-index that has a successor), compared with the cells of the section.</summary>
    private struct ChildCompareSink : IIntervalSink
    {
        private const int NextKind = 1 << 30;
        private const int UpKind = 2 << 30;
        private const int DownKind = 3 << 30;

        private readonly int* _cells;
        private readonly int _n;
        private readonly ulong[] _seen;
        private readonly bool _tagged;
        internal long Writes;
        internal long Internal;
        internal int LrsDepth;
        internal int LrsRank;

        internal ChildCompareSink(int* cells, int n, ulong[] seen, bool tagged)
        {
            _cells = cells;
            _n = n;
            _seen = seen;
            _tagged = tagged;
        }

        public void LIndex(int previous, int k)
        {
            if (previous < 0) return;
            if (_tagged) Expect(previous, NextKind | k);
            else if (k != previous + 1) Expect(previous, k);
        }

        public void Close(int ell, int lb, int rb, int first, bool lastChild)
        {
            int kind = lastChild ? DownKind : UpKind;
            Expect(lastChild ? lb : rb, _tagged ? kind | first : first);
            if (ell <= 0) return;
            Internal++;
            if (ell > LrsDepth)
            {
                LrsDepth = ell;
                LrsRank = lb;
            }
        }

        private void Expect(int cell, int value)
        {
            if ((uint)cell > (uint)_n) throw new InvalidDataException($"the lcp-interval scan defines the cell {cell} of a child table of {_n + 1} cells: the file is damaged");
            ref ulong word = ref _seen[cell >> 6];
            ulong mask = 1UL << (cell & 63);
            if ((word & mask) != 0) throw new InvalidDataException($"the lcp-interval scan defines the child table cell {cell} twice: the file is damaged");
            word |= mask;
            Writes++;
            if (_cells[cell] != value) throw new InvalidDataException($"the child table cell {cell} holds {_cells[cell]} but the LCP array defines {value}: the file is damaged");
        }
    }

    /// <summary>
    /// A stack of <see cref="Frame"/>s that keeps at most <c>capacity</c> frames in memory: when it is full the older half goes to a temporary file (block after block, last in first out), and a pop that empties the in-memory part
    /// reads the last block back. The file is created at the first spill and deleted when the stack is disposed.
    /// </summary>
    private sealed class FrameStack : IDisposable
    {
        private readonly Frame[] _frames;
        private readonly int _half;
        private readonly string? _directory;
        private int _count;
        private long _spilled;
        private FileStream? _file;
        private long _peak;

        internal FrameStack(int capacity, string? directory)
        {
            _frames = new Frame[capacity];
            _half = capacity / 2;
            _directory = directory;
        }

        /// <summary>Frames on the stack, those on disk included.</summary>
        internal long Count => _count + _spilled;

        internal long PeakMemoryBytes => _peak * 16;

        internal Frame Top => _frames[_count - 1];

        internal ref Frame TopRef => ref _frames[_count - 1];

        internal void Push(Frame frame)
        {
            if (_count == _frames.Length) Spill();
            _frames[_count++] = frame;
            if (_count > _peak) _peak = _count;
        }

        internal Frame Pop()
        {
            Frame top = _frames[--_count];
            if (_count == 0 && _spilled > 0) Refill();
            return top;
        }

        private void Spill()
        {
            if (_file is null)
            {
                string path = Path.Combine(_directory ?? Path.GetTempPath(), "esa-verify-" + Guid.NewGuid().ToString("N") + ".stack");
                _file = new FileStream(path, FileMode.CreateNew, FileAccess.ReadWrite, FileShare.None, 1 << 16, FileOptions.DeleteOnClose);
            }
            _file.Seek(_spilled * 16, SeekOrigin.Begin);
            _file.Write(System.Runtime.InteropServices.MemoryMarshal.AsBytes(_frames.AsSpan(0, _half)));
            _spilled += _half;
            Array.Copy(_frames, _half, _frames, 0, _count - _half);
            _count -= _half;
        }

        private void Refill()
        {
            int take = (int)Math.Min(_half, _spilled);
            _file!.Seek((_spilled - take) * 16, SeekOrigin.Begin);
            _file.ReadExactly(System.Runtime.InteropServices.MemoryMarshal.AsBytes(_frames.AsSpan(0, take)));
            _spilled -= take;
            _count = take;
        }

        public void Dispose() => _file?.Dispose();
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
}
