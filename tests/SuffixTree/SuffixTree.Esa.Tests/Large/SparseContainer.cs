using System.Buffers.Binary;
using System.Globalization;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// What a sparse container is (work package S2.1, deliverable B): the header fields and the shape of the sections of a SESA container of <see cref="N"/> symbols whose big sections are holes. The container is a valid *shell*:
/// its header page and section table are valid for the width (ADR-003, ADR-005 section 6 and 6.1), every section CRC-32C is that of the section as it is (zeros and the patches), so it opens at the <c>Header</c> and <c>Sections</c>
/// levels of a reader; its content is zero except for what a test patches in, so the <c>Structure</c> level (the suffix array is a permutation) must refuse it.
/// </summary>
internal sealed record SparseSpec
{
    /// <summary>The largest text that 4-byte positions hold in the container (what the S1 writer and header reader accept): 2^32 - 4.</summary>
    internal const long MaxContainer4 = uint.MaxValue - 3L;

    /// <summary>The text length n.</summary>
    internal required long N { get; init; }

    /// <summary>4 or 5; 0 chooses by n (<see cref="WidthFor"/>).</summary>
    internal int PosWidth { get; init; }

    /// <summary>Distinct letters (1 for a^n).</summary>
    internal int Sigma { get; init; } = 1;

    /// <summary><c>NODE_COUNT</c>.</summary>
    internal required long NodeCount { get; init; }

    /// <summary><c>LRS_DEPTH</c>.</summary>
    internal long LrsDepth { get; init; }

    /// <summary><c>LRS_RANK</c> (-1 when there is no repeat).</summary>
    internal long LrsRank { get; init; } = -1;

    /// <summary>The largest LCP value (decides the LCP codec when <see cref="LcpCodec"/> is 0).</summary>
    internal long MaxLcp { get; init; }

    /// <summary>1 U8, 2 U16, 3 U32, 4 U64; 0 the narrowest codec that holds <see cref="MaxLcp"/> without an overflow list.</summary>
    internal int LcpCodec { get; init; }

    /// <summary>CHILD codec: 1 (packed). A shell above 2^30 symbols has no other.</summary>
    internal int ChildCodec { get; init; } = SesaFormat.ChildCodecPacked;

    /// <summary>Whether the PSI section is present.</summary>
    internal bool Psi { get; init; }

    /// <summary>Whether the RMQ section (the persisted LCP pyramid, S2.5) is present.</summary>
    internal bool Rmq { get; init; }

    /// <summary><c>PARAM0</c> of PSI: the rank of the suffix at position 0.</summary>
    internal long PsiPrimaryRow { get; init; }

    /// <summary><c>TEXT_SHA256</c> (null = zeros: not requested).</summary>
    internal byte[]? TextSha256 { get; init; }

    /// <summary>The BUILD_INFO text (UTF-8, at most 4096 bytes).</summary>
    internal string? BuildInfo { get; init; }

    /// <summary>The position width: the given one or the one of n.</summary>
    internal int Width => PosWidth != 0 ? PosWidth : WidthFor(N);

    /// <summary>The LCP codec: the given one or the narrowest.</summary>
    internal int Lcp => LcpCodec != 0 ? LcpCodec : LcpCodecFor(MaxLcp);

    /// <summary>Bits per packed symbol.</summary>
    internal int Bits => PackedTextFormat.BitsFor(Sigma);

    /// <summary>4 bytes up to 2^32 - 4 symbols, 5 above (ADR-005 section 6).</summary>
    internal static int WidthFor(long n) => n <= MaxContainer4 ? 4 : 5;

    /// <summary>The narrowest LCP codec that holds the value without an escape (U8 up to 254, U16 up to 65,534, U32 up to 2^32 - 2, else U64: ADR-005 section 6, "a repeat longer than 2^32 - 2").</summary>
    internal static int LcpCodecFor(long maxLcp) => maxLcp <= 254 ? 1 : maxLcp <= 65_534 ? 2 : maxLcp <= uint.MaxValue - 1L ? 3 : 4;

    /// <summary>Bytes of one LCP entry of a codec.</summary>
    internal static int LcpBytes(int codec) => codec switch { 1 => 1, 2 => 2, 3 => 4, 4 => 8, _ => throw new ArgumentOutOfRangeException(nameof(codec), codec, "LCP codec 1..4") };

    /// <summary>The spec of a closed-form text: the header numbers from its statistics.</summary>
    internal static SparseSpec For(ClosedFormText f, int posWidth = 0, bool psi = false, bool rmq = false) => new()
    {
        Rmq = rmq,
        N = f.Length,
        PosWidth = posWidth,
        Sigma = f.Sigma,
        NodeCount = f.NodeCount,
        LrsDepth = f.LrsDepth,
        LrsRank = f.LrsRank,
        MaxLcp = f.LrsDepth,
        Psi = psi,
        PsiPrimaryRow = f.Length == 0 ? 0 : f.RankOf(0),
    };
}

/// <summary>Where the sections of a sparse container are and what they hold.</summary>
internal sealed record SparseLayout(SparseSpec Spec, int PosWidth, int LcpCodec, IReadOnlyList<SesaSection> Sections, long FileSize)
{
    /// <summary>The section of a kind (with a CRC of 0: the CRCs are filled in by the writer).</summary>
    internal SesaSection Section(int kind) => Sections.First(s => s.Kind == kind);

    /// <summary>Whether the layout has a section of the kind.</summary>
    internal bool Has(int kind) => Sections.Any(s => s.Kind == kind);

    /// <summary>Bytes of the section that carry data (the rest is the zero padding that belongs to the payload).</summary>
    internal long DataBytes(int kind)
    {
        SesaSection s = Section(kind);
        return kind switch
        {
            SesaFormat.KindText => s.Length - 64,
            SesaFormat.KindSa or SesaFormat.KindChild or SesaFormat.KindPsi or SesaFormat.KindRmq => s.Count * s.ElemWidth,
            _ => s.Length,
        };
    }
}

/// <summary>A few bytes of real data at <paramref name="Offset"/> of the payload of the section of <paramref name="Kind"/>.</summary>
internal sealed record SparsePatch(int Kind, long Offset, byte[] Data);

/// <summary>
/// Writes a SESA container whose big sections are holes (work package S2.1, deliverable B): <c>SetLength</c> to the final size, the header page and the few patches written, nothing else. The section CRC-32C values are those
/// of the section as it is: the CRC of N zero bytes by the combine of <see cref="Crc32CCombine"/> (never by a pass), the patches folded in. The layout is that of the product's writer (<see cref="SesaLayout"/>: TEXT, ALPHABET,
/// COUNTS, SA, CHILD, PSI, BUILD_INFO, then LCP; 64-byte alignment, 4096 from 1 MiB) at width 4 and that of ADR-005 section 6.1 at width 5 (entries of 5 bytes and 3 bytes of padding in each position section, version 1.1).
/// <para>Safety: no more than <see cref="LargeDisk.MaxDenseBytes"/> of patches, a nominal size of at most <see cref="LargeDisk.MaxNominalBytes"/>, the filesystem probed for sparse files before any file above 1 GiB is made
/// (<see cref="SparseUnavailableException"/> otherwise), and the file removed when anything fails. A sparse file is never mapped and written to as a whole.</para>
/// </summary>
internal static class SparseContainer
{
    private const int BuildInfoCapacity = SesaLayout.BuildInfoCapacity;
    private const int ZeroPadding5 = 3;
    private const int TextPadding = 64;

    /// <summary>The layout of the container of <paramref name="spec"/>, with CRC 0 in every entry.</summary>
    internal static SparseLayout Plan(SparseSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        ArgumentOutOfRangeException.ThrowIfNegative(spec.N);
        int width = spec.Width;
        if (width is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(spec), width, "POS_WIDTH is 4 or 5");
        if (width == 4 && spec.N > SparseSpec.MaxContainer4) throw new ArgumentOutOfRangeException(nameof(spec), spec.N, $"POS_WIDTH 4 holds a text of at most {SparseSpec.MaxContainer4} symbols");
        if (spec.N >= 1L << 40) throw new ArgumentOutOfRangeException(nameof(spec), spec.N, "POS_WIDTH 5 holds a text of fewer than 2^40 symbols");
        int lcpCodec = spec.Lcp;
        int lcpBytes = SparseSpec.LcpBytes(lcpCodec);
        long n = spec.N, rows = n + 1;
        int sigma = spec.Sigma, bits = spec.Bits;
        long textLength = ((n * bits) + 7) / 8 + TextPadding;
        long alphabetLength = 2L * sigma;
        long countsLength = 8L * sigma;
        long arrayLength = (width * rows) + (width == 5 ? ZeroPadding5 : 0);
        byte[] info = InfoBytes(spec);

        long textOffset = SesaFormat.AlignedOffset(SesaFormat.PageBytes, textLength);
        long alphabetOffset = SesaFormat.AlignedOffset(textOffset + textLength, alphabetLength);
        long countsOffset = SesaFormat.AlignedOffset(alphabetOffset + alphabetLength, countsLength);
        long saOffset = SesaFormat.AlignedOffset(countsOffset + countsLength, arrayLength);
        long childOffset = SesaFormat.AlignedOffset(saOffset + arrayLength, arrayLength);
        long end = childOffset + arrayLength;
        long psiOffset = 0;
        if (spec.Psi)
        {
            psiOffset = SesaFormat.AlignedOffset(end, arrayLength);
            end = psiOffset + arrayLength;
        }
        long infoOffset = SesaFormat.AlignedOffset(end, BuildInfoCapacity);
        int widestLcp = width == 5 ? 8 : 4; // the LCP is aligned for the widest codec of the width: its codec is not known when the product places it
        long lcpOffset = SesaFormat.AlignedOffset(infoOffset + BuildInfoCapacity, widestLcp * rows);
        long lcpLength = rows * lcpBytes;
        long fileSize = lcpOffset + lcpLength;
        long rmqCount = spec.Rmq ? RmqLevels(rows).Sum(l => l.Count) : 0;
        long rmqLength = spec.Rmq ? (width * rmqCount) + (width == 5 ? ZeroPadding5 : 0) : 0;
        long rmqOffset = 0;
        if (spec.Rmq)
        {
            rmqOffset = SesaFormat.AlignedOffset(fileSize, rmqLength);
            fileSize = rmqOffset + rmqLength;
        }

        const int required = SesaFormat.FlagRequired | SesaFormat.FlagHasCrc;
        const int optional = SesaFormat.FlagHasCrc;
        var sections = new List<SesaSection>
        {
            new(SesaFormat.KindText, required, textOffset, textLength, n, 0, PackedTextFormat.CodecId, bits, 0, 0),
            new(SesaFormat.KindAlphabet, required, alphabetOffset, alphabetLength, sigma, 2, 1, 0, 0, 0),
            new(SesaFormat.KindCounts, optional, countsOffset, countsLength, sigma, 8, 1, 0, 0, 0),
            new(SesaFormat.KindSa, required, saOffset, arrayLength, rows, width, 1, 0, 0, 0),
            new(SesaFormat.KindChild, required, childOffset, arrayLength, rows, width, spec.ChildCodec, 0, 0, 0),
        };
        if (spec.Psi) sections.Add(new SesaSection(SesaFormat.KindPsi, optional, psiOffset, arrayLength, rows, width, 1, spec.PsiPrimaryRow, 0, 0));
        sections.Add(new SesaSection(SesaFormat.KindBuildInfo, optional, infoOffset, info.Length, info.Count(b => b == (byte)'\n'), 0, 1, 0, 0, 0));
        sections.Add(new SesaSection(SesaFormat.KindLcp, required, lcpOffset, lcpLength, rows, lcpBytes, lcpCodec, 0, 0, 0));
        if (spec.Rmq) sections.Add(new SesaSection(SesaFormat.KindRmq, optional, rmqOffset, rmqLength, rmqCount, width, 1, SesaFormat.RmqFanout, 0, 0));
        return new SparseLayout(spec, width, lcpCodec, sections, fileSize);
    }

    /// <summary>The levels of the pyramid of an LCP array of <paramref name="rows"/> entries: where each starts (in entries of the RMQ section) and how many entries it has (level 0 first).</summary>
    internal static List<(long Start, long Count)> RmqLevels(long rows)
    {
        var levels = new List<(long, long)>();
        long start = 0, size = rows;
        while (size > 1)
        {
            size = (size + SesaFormat.RmqFanout - 1) / SesaFormat.RmqFanout;
            levels.Add((start, size));
            start += size;
        }
        return levels;
    }

    /// <summary>The nominal size of the file (what <c>ls</c> says; the disk takes only the header page and the patches).</summary>
    internal static long NominalBytes(SparseSpec spec) => Plan(spec).FileSize;

    private static byte[] InfoBytes(SparseSpec spec)
    {
        string text = spec.BuildInfo ?? string.Create(CultureInfo.InvariantCulture, $"format=SESA {(spec.Width == 5 ? "1.1" : "1.0")} sparse test shell (S2.1)\nn={spec.N}\npos_width={spec.Width}\n");
        byte[] bytes = System.Text.Encoding.UTF8.GetBytes(text);
        if (bytes.Length > BuildInfoCapacity) throw new ArgumentException($"BUILD_INFO has {bytes.Length} bytes, the capacity is {BuildInfoCapacity}", nameof(spec));
        return bytes;
    }

    /// <summary>
    /// Writes the container to <paramref name="path"/>: the file is made with <c>SetLength</c> (holes), the header page and the <paramref name="patches"/> are written, the BUILD_INFO text is always written.
    /// Everything that can be refused is refused before the file is made.
    /// </summary>
    /// <exception cref="SparseUnavailableException">The filesystem does not make sparse files (a file above 1 GiB is wanted).</exception>
    /// <exception cref="ArgumentException">More than 64 MiB of patches, a patch outside its section or into the padding, overlapping patches, a nominal size above the cap.</exception>
    internal static void Write(string path, SparseSpec spec, IReadOnlyList<SparsePatch> patches)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        ArgumentNullException.ThrowIfNull(patches);
        SparseLayout layout = Plan(spec);
        if (layout.FileSize > LargeDisk.MaxNominalBytes) throw new ArgumentException($"the nominal size {layout.FileSize} is above the cap of {LargeDisk.MaxNominalBytes} bytes", nameof(spec));
        long dense = SesaFormat.PageBytes + patches.Sum(p => (long)p.Data.Length);
        if (dense > LargeDisk.MaxDenseBytes) throw new ArgumentException($"{dense} bytes of real data: no dense region above {LargeDisk.MaxDenseBytes} bytes is ever written", nameof(patches));

        var all = new List<SparsePatch>(patches);
        byte[] info = InfoBytes(spec);
        all.Add(new SparsePatch(SesaFormat.KindBuildInfo, 0, info));
        var perSection = new Dictionary<int, List<SparsePatch>>();
        foreach (SparsePatch patch in all.OrderBy(p => p.Kind).ThenBy(p => p.Offset))
        {
            if (!layout.Has(patch.Kind)) throw new ArgumentException($"the container has no section of kind {patch.Kind}", nameof(patches));
            long data = layout.DataBytes(patch.Kind);
            if (patch.Offset < 0 || patch.Offset > data - patch.Data.Length) throw new ArgumentException($"a patch of {patch.Data.Length} bytes at {patch.Offset} of section {patch.Kind} leaves its {data} bytes of data (the padding stays zero)", nameof(patches));
            if (!perSection.TryGetValue(patch.Kind, out List<SparsePatch>? list)) perSection[patch.Kind] = list = [];
            list.Add(patch);
        }

        var sections = new List<SesaSection>();
        foreach (SesaSection s in layout.Sections)
        {
            List<(long Offset, byte[] Data)> parts = perSection.TryGetValue(s.Kind, out List<SparsePatch>? l) ? [.. l.Select(p => (p.Offset, p.Data))] : [];
            sections.Add(s with { Crc = Crc32CCombine.OfSparse(s.Length, parts) });
        }
        byte[] page = HeaderPage(spec, layout, sections);

        if (layout.FileSize > LargeDisk.ProbeBytes) LargeDisk.RequireSparse(Path.GetDirectoryName(Path.GetFullPath(path))!);
        if (LargeDisk.FreeBytes(Path.GetDirectoryName(Path.GetFullPath(path))!) < LargeDisk.MinimumFreeBytes) throw new SparseUnavailableException("the volume is too full for a sparse test");
        using var stream = new FileStream(path, FileMode.CreateNew, FileAccess.Write, FileShare.None); // an existing file is an IOException here and is never deleted
        try
        {
            stream.SetLength(layout.FileSize);
            stream.Position = 0;
            stream.Write(page);
            foreach (SesaSection s in layout.Sections)
            {
                if (!perSection.TryGetValue(s.Kind, out List<SparsePatch>? l)) continue;
                foreach (SparsePatch p in l)
                {
                    stream.Position = s.Offset + p.Offset;
                    stream.Write(p.Data);
                }
            }
            stream.Flush();
        }
        catch
        {
            stream.Dispose();
            File.Delete(path);
            throw;
        }
    }

    /// <summary>The 4096 bytes of header page and section table of the container.</summary>
    internal static byte[] HeaderPage(SparseSpec spec, SparseLayout layout, IReadOnlyList<SesaSection> sections)
    {
        var header = new SesaHeader
        {
            TextLength = spec.N,
            NodeCount = spec.NodeCount,
            FileSize = layout.FileSize,
            LrsDepth = spec.LrsDepth,
            LrsRank = spec.LrsRank,
            PosWidth = layout.PosWidth,
            Sigma = spec.Sigma,
            TextSha256 = spec.TextSha256 ?? new byte[32],
        };
        byte[] page = SesaFormat.EncodeHeaderPage(header, sections, withMagic: true);
        if (layout.PosWidth == 5)
        {
            // ADR-005 section 6.1: a container with 5-byte positions is a 1.1 container; the header CRC covers the version word, so it is recomputed
            BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);
            ResealHeader(page);
        }
        return page;
    }

    /// <summary>Recomputes the header CRC-32C of a header page (the bytes [0, 124) followed by the section table).</summary>
    internal static void ResealHeader(byte[] page)
    {
        int count = BinaryPrimitives.ReadInt32LittleEndian(page.AsSpan(48));
        uint crc = Crc32C.Append(Crc32C.Append(0, page.AsSpan(0, SesaFormat.CrcOffset)), page.AsSpan(SesaFormat.HeaderBytes, SesaFormat.EntryBytes * count));
        BinaryPrimitives.WriteUInt32LittleEndian(page.AsSpan(SesaFormat.CrcOffset), crc);
    }

    // ------------------------------------------------------------------------------------------- the data of a closed form

    /// <summary>
    /// The real bytes of the container of a closed-form text at the two ends and in the small sections: ALPHABET, COUNTS, the first and the last <paramref name="entriesPerEnd"/> entries of SA, LCP, CHILD (and PSI), and the
    /// first and the last symbols of TEXT. The middle stays a hole, so the file is *not* the index of the text (only its ends are).
    /// </summary>
    internal static IReadOnlyList<SparsePatch> EndPatches(ClosedFormText f, SparseLayout layout, long entriesPerEnd)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(entriesPerEnd);
        long rows = f.Length + 1;
        var patches = new List<SparsePatch>(SmallSections(f, layout));
        long head = Math.Min(entriesPerEnd, rows);
        long tailFrom = Math.Max(head, rows - entriesPerEnd);
        patches.AddRange(Entries(f, layout, 0, head - 1));
        if (tailFrom < rows) patches.AddRange(Entries(f, layout, tailFrom, rows - 1));
        long symbols = Math.Min(f.Length, (entriesPerEnd + 7) / 8 * 8);
        if (f.Length > 0)
        {
            patches.AddRange(TextWindow(f, layout, 0, symbols));
            long tailStart = Math.Max(symbols, (f.Length - entriesPerEnd) / 8 * 8);
            if (tailStart < f.Length) patches.AddRange(TextWindow(f, layout, tailStart, f.Length));
        }
        return patches;
    }

    /// <summary>The entries of SA, CHILD, LCP (and PSI) of the ranks <paramref name="from"/> to <paramref name="to"/>: a window in the middle of a shell (a page of <c>Locate</c> that straddles 2^31 reads real values).</summary>
    internal static IReadOnlyList<SparsePatch> WindowPatches(ClosedFormText f, SparseLayout layout, long from, long to) => [.. Entries(f, layout, from, to)];

    /// <summary>Every entry and every symbol (a small text only: the data are real bytes in the file).</summary>
    internal static IReadOnlyList<SparsePatch> AllPatches(ClosedFormText f, SparseLayout layout)
    {
        var patches = new List<SparsePatch>(SmallSections(f, layout));
        patches.AddRange(Entries(f, layout, 0, f.Length));
        if (f.Length > 0) patches.AddRange(TextWindow(f, layout, 0, f.Length));
        return patches;
    }

    /// <summary>ALPHABET (the letters) and COUNTS (the symbol counts).</summary>
    private static IEnumerable<SparsePatch> SmallSections(ClosedFormText f, SparseLayout layout)
    {
        int sigma = f.Sigma;
        if (sigma == 0) yield break;
        var alphabet = new byte[2 * sigma];
        var counts = new byte[8 * sigma];
        for (int s = 0; s < sigma; s++)
        {
            BinaryPrimitives.WriteUInt16LittleEndian(alphabet.AsSpan(2 * s), (ushort)ClosedFormText.LetterOf(s));
            BinaryPrimitives.WriteInt64LittleEndian(counts.AsSpan(8 * s), f.SymbolCount(s));
        }
        yield return new SparsePatch(SesaFormat.KindAlphabet, 0, alphabet);
        yield return new SparsePatch(SesaFormat.KindCounts, 0, counts);
    }

    /// <summary>The entries of SA, CHILD, LCP (and PSI) of the ranks <paramref name="from"/> to <paramref name="to"/>, in one patch per section.</summary>
    private static IEnumerable<SparsePatch> Entries(ClosedFormText f, SparseLayout layout, long from, long to)
    {
        if (to < from) yield break;
        int w = layout.PosWidth, lw = SparseSpec.LcpBytes(layout.LcpCodec);
        long count = to - from + 1;
        var sa = new byte[checked(count * w)];
        var child = new byte[checked(count * w)];
        var lcp = new byte[checked(count * lw)];
        byte[]? psi = layout.Has(SesaFormat.KindPsi) ? new byte[checked(count * w)] : null;
        for (long r = from; r <= to; r++)
        {
            long at = r - from;
            PutPosition(sa, at, w, f.SuffixAt(r));
            PutPosition(child, at, w, f.PackedCell(r));
            PutLittleEndian(lcp, at * lw, lw, f.Lcp(r));
            if (psi is not null) PutPosition(psi, at, w, PsiOf(f, r));
        }
        yield return new SparsePatch(SesaFormat.KindSa, from * w, sa);
        yield return new SparsePatch(SesaFormat.KindChild, from * w, child);
        yield return new SparsePatch(SesaFormat.KindLcp, from * lw, lcp);
        if (psi is not null) yield return new SparsePatch(SesaFormat.KindPsi, from * w, psi);
    }

    /// <summary>
    /// The patches of the RMQ section for the queries around the given windows of ranks (call it once: its patches are merged and never overlap): on every level the entries of the blocks that contain a window and the 32 on either side
    /// of them (the group of siblings that the climb of a search scans), and every entry of a level of at most 2^16 entries; each entry is the closed-form minimum of its block of the LCP array (<see cref="ClosedFormText.LcpMin"/>).
    /// </summary>
    internal static IReadOnlyList<SparsePatch> RmqPatches(ClosedFormText f, SparseLayout layout, params (long From, long To)[] windows)
    {
        var patches = new List<SparsePatch>();
        if (!layout.Has(SesaFormat.KindRmq)) return patches;
        long rows = f.Length + 1;
        int w = layout.PosWidth;
        long block = SesaFormat.RmqFanout;
        foreach ((long start, long count) in RmqLevels(rows))
        {
            var spans = new List<(long Lo, long Hi)>();
            if (count <= 1 << 16) spans.Add((0, count - 1));
            else
            {
                foreach ((long from, long to) in windows) spans.Add((Math.Max(0, (from / block) - 32), Math.Min(count - 1, (to / block) + 32)));
                spans.Sort();
                for (int i = spans.Count - 1; i > 0; i--)
                {
                    if (spans[i].Lo <= spans[i - 1].Hi + 1)
                    {
                        spans[i - 1] = (spans[i - 1].Lo, Math.Max(spans[i - 1].Hi, spans[i].Hi));
                        spans.RemoveAt(i);
                    }
                }
            }
            foreach ((long lo, long hi) in spans)
            {
                var data = new byte[checked((hi - lo + 1) * w)];
                for (long j = lo; j <= hi; j++)
                {
                    long first = j * block;
                    PutPosition(data, j - lo, w, f.LcpMin(first, Math.Min(rows - 1, first + block - 1)));
                }
                patches.Add(new SparsePatch(SesaFormat.KindRmq, (start + lo) * w, data));
            }
            block *= SesaFormat.RmqFanout;
        }
        return patches;
    }

    /// <summary><c>PSI[r] = ISA[SA[r] + 1]</c>; 0 for the empty suffix and for the suffix of the last symbol (its successor is the empty suffix).</summary>
    internal static long PsiOf(ClosedFormText f, long rank)
    {
        if (rank == 0) return 0;
        long position = f.SuffixAt(rank);
        return position + 1 >= f.Length ? 0 : f.RankOf(position + 1);
    }

    /// <summary>The packed text of the symbols <c>[from, to)</c> (<paramref name="from"/> a multiple of 8, so the window starts on a byte; <paramref name="to"/> a multiple of 8 or the end of the text).</summary>
    private static List<SparsePatch> TextWindow(ClosedFormText f, SparseLayout layout, long from, long to)
    {
        int bits = layout.Spec.Bits;
        if (from % 8 != 0 || (to % 8 != 0 && to != f.Length)) throw new ArgumentException("a window of the text starts and ends on a byte");
        long bytes = ((to - from) * bits + 7) / 8;
        var data = new byte[checked(bytes)];
        for (long s = from; s < to; s++)
        {
            int code = f.CodeAt(s);
            long bit = (s - from) * bits;
            for (int b = 0; b < bits; b++)
                if (((code >> (bits - 1 - b)) & 1) != 0) data[(bit + b) >> 3] |= (byte)(0x80 >> (int)((bit + b) & 7));
        }
        return data.Any(x => x != 0) ? [new SparsePatch(SesaFormat.KindText, from * bits / 8, data)] : [];
    }

    /// <summary>Stores a position as <paramref name="width"/> little-endian bytes (4 or 5) at entry <paramref name="entry"/>.</summary>
    internal static void PutPosition(byte[] target, long entry, int width, long value)
    {
        if (value < 0 || (width == 4 && value > uint.MaxValue) || (width == 5 && value >= 1L << 40)) throw new ArgumentOutOfRangeException(nameof(value), value, $"does not fit {width} bytes");
        PutLittleEndian(target, entry * width, width, value);
    }

    /// <summary>Stores the low <paramref name="width"/> bytes of <paramref name="value"/> little endian.</summary>
    internal static void PutLittleEndian(byte[] target, long offset, int width, long value)
    {
        for (int i = 0; i < width; i++) target[offset + i] = (byte)(value >> (8 * i));
    }
}
