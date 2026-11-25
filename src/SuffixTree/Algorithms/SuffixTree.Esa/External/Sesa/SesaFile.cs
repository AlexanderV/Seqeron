using System.Diagnostics;
using System.IO.MemoryMappedFiles;

namespace SuffixTree.Esa.External;

/// <summary>What the opening of a container did and what it cost (reported by <see cref="SesaFile.Report"/>).</summary>
/// <param name="Level">The validation level that ran.</param>
/// <param name="Access">The access mode.</param>
/// <param name="HeaderSeconds">Time of the header level (mapping included).</param>
/// <param name="SectionsSeconds">Time of the section checksums (0 below that level).</param>
/// <param name="StructureSeconds">Time of the structure pass (0 below that level).</param>
/// <param name="HardwareCrc">Whether the CRC-32C used the hardware instruction.</param>
/// <param name="SectionBytesChecksummed">Payload bytes that the checksum pass read.</param>
/// <param name="BitmapBytes">Size of the permutation bitmap of the structure pass (n + 1 bits, rounded up to bytes).</param>
/// <param name="BufferBytes">Anonymous memory of the structure pass besides the bitmap (buffers, stacks), at its peak.</param>
internal sealed record SesaOpenReport(EsaValidationLevel Level, EsaAccessMode Access, double HeaderSeconds, double SectionsSeconds, double StructureSeconds, bool HardwareCrc, long SectionBytesChecksummed, long BitmapBytes, long BufferBytes);

/// <summary>
/// A SESA container mapped read-only as one view (ADR-003 section 9, work package 1b.7f): <c>MemoryMappedFile.CreateFromFile</c>, one view of the whole file, the pointer acquired once and released in <see cref="Dispose"/>. The header page and
/// the section table are validated at the requested level (<see cref="SesaHeaderReader"/>, <see cref="SesaVerifier"/>); the sections are then exposed as views over the pointer (<see cref="Sa"/>, <see cref="Lcp"/>, <see cref="Child"/>,
/// <see cref="TextView"/>, <see cref="PsiPointer"/>): no section is copied and no array is allocated per section. The mapping is shared by every reader of the file; it is read-only, so the reads need no lock.
/// </summary>
internal sealed unsafe class SesaFile : IDisposable
{
    private readonly MemoryMappedFile? _map;
    private readonly MemoryMappedViewAccessor? _view;
    private readonly byte* _data;
    private readonly long _length;
    private readonly EsaOpenOptions _options;
    private readonly SesaParsed _parsed;
    private readonly MappedTextReader _reader;
    private readonly EsaText _text;
    private readonly CompactLcp _lcp;
    private readonly ChildTable _child;
    private readonly SaView _sa;
    private readonly bool _narrow;
    private int _disposed;
    private SesaOpenReport _report;

    private SesaFile(MemoryMappedFile? map, MemoryMappedViewAccessor? view, byte* data, long length, EsaOpenOptions options, SesaParsed parsed, SesaOpenReport report)
    {
        _map = map;
        _view = view;
        _data = data;
        _length = length;
        _options = options;
        _parsed = parsed;
        _report = report;
        bool isChecked = options.Access == EsaAccessMode.Checked;
        long n = parsed.Header.TextLength;
        int posWidth = parsed.Header.PosWidth;
        _reader = MappedTextReader.FromMemory(data + parsed.Text.Offset, parsed.Text.Length, n, parsed.Bits);
        _text = EsaText.FromPacked(_reader, parsed.Letters, n);
        // the views of every container, at any size and width (work package S2.3; before it a text of 2^31 - 1 symbols or more had none: debt 2 of ADR-005 section 11). The 32-bit forms that the engine of S1 reads (MappedSearch, the pointer
        // readers) are made where they exist: a text below 2^31 - 1 symbols at POS_WIDTH 4 and an LCP array of 8, 16 or 32 bits; every other container gets the wide forms, whose long members work and whose 32-bit members are checked bridges
        _sa = SaView.FromMapped(data + parsed.Sa.Offset, n + 1, posWidth, isChecked);
        byte* pairs = parsed.Overflow is { } o ? data + o.Offset : null;
        bool tagged = parsed.Child.Codec == SesaFormat.ChildCodecTagged;
        bool narrow = n < int.MaxValue && posWidth == 4 && parsed.LcpKind != LcpKind.U64;
        _narrow = narrow;
        if (narrow)
        {
            _lcp = CompactLcp.FromMapped(parsed.LcpKind, (int)(n + 1), data + parsed.Lcp.Offset, (uint*)pairs, (int)parsed.Lcp.Param0, isChecked);
            _child = ChildTable.FromMapped(new MappedChild((int*)(data + parsed.Child.Offset), (int)(n + 1), _lcp, isChecked, tagged));
        }
        else
        {
            _lcp = CompactLcp.FromMappedWide(parsed.LcpKind, n + 1, data + parsed.Lcp.Offset, pairs, parsed.Lcp.Param0, posWidth, isChecked);
            _child = ChildTable.FromMapped(MappedChild.FromMappedWide(data + parsed.Child.Offset, n + 1, posWidth, _lcp, isChecked, tagged));
        }
    }

    /// <summary>Maps <paramref name="path"/> and validates it at the level of <paramref name="options"/>; <paramref name="hardwareCrc"/> forces (true) or forbids (false) the CRC instruction, null lets the CPU decide.</summary>
    /// <exception cref="InvalidDataException">The file is not a SESA container or is damaged at the requested level.</exception>
    /// <exception cref="NotSupportedException">A larger major version, or a feature of a valid container that this reader does not have.</exception>
    internal static SesaFile Open(string path, EsaOpenOptions options, bool? hardwareCrc = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        CheckOptions(options);
        long t0 = Stopwatch.GetTimestamp();
        var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete, 1, FileOptions.None);
        MemoryMappedFile? map = null;
        MemoryMappedViewAccessor? view = null;
        byte* data = null;
        long length;
        try
        {
            if (!stream.CanSeek) throw new InvalidDataException("not a regular file: the stream cannot seek");
            length = stream.Length;
            if (length < SesaFormat.PageBytes) throw new InvalidDataException($"truncated: the file has {length} bytes, shorter than the header page of {SesaFormat.PageBytes}");
            map = MemoryMappedFile.CreateFromFile(stream, null, 0, MemoryMappedFileAccess.Read, HandleInheritability.None, leaveOpen: false);
            view = map.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
            byte* p = null;
            view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            data = p + view.PointerOffset;
        }
        catch
        {
            view?.Dispose();
            map?.Dispose();
            stream.Dispose(); // a no-op once the map owns the handle; the failure of a step before the mapping closes it here
            throw;
        }
        return Create(map, view, data, length, options, hardwareCrc, t0);
    }

    /// <summary>The container as it lies in memory that the caller owns (tests: memory with guard pages): the same checks, nothing is mapped or freed.</summary>
    internal static SesaFile FromMemory(byte* data, long length, EsaOpenOptions options, bool? hardwareCrc = null)
    {
        CheckOptions(options);
        return Create(null, null, data, length, options, hardwareCrc, Stopwatch.GetTimestamp());
    }

    private static void CheckOptions(EsaOpenOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        if (options.Validation is < EsaValidationLevel.Header or > EsaValidationLevel.Structure) throw new ArgumentOutOfRangeException(nameof(options), options.Validation, "not a validation level");
        if (options.Access is < EsaAccessMode.Checked or > EsaAccessMode.Trusted) throw new ArgumentOutOfRangeException(nameof(options), options.Access, "not an access mode");
    }

    private static SesaFile Create(MemoryMappedFile? map, MemoryMappedViewAccessor? view, byte* data, long length, EsaOpenOptions options, bool? hardwareCrc, long t0)
    {
        SesaFile? file = null;
        try
        {
            SesaParsed parsed = SesaHeaderReader.Parse(data, length);
            double header = Seconds(t0);
            var report = new SesaOpenReport(options.Validation, options.Access, header, 0, 0, hardwareCrc ?? Crc32C.HardwareAccelerated, 0, 0, 0);
            file = new SesaFile(map, view, data, length, options, parsed, report);
            if (options.Validation >= EsaValidationLevel.Sections) SesaVerifier.CheckSections(file, hardwareCrc);
            if (options.Validation >= EsaValidationLevel.Structure) SesaVerifier.CheckStructure(file);
            return file;
        }
        catch
        {
            if (file is not null)
            {
                file.Dispose();
            }
            else
            {
                if (view is not null)
                {
                    view.SafeMemoryMappedViewHandle.ReleasePointer();
                    view.Dispose();
                }
                map?.Dispose();
            }
            throw;
        }
    }

    /// <summary>The header fields.</summary>
    internal SesaHeader Header => _parsed.Header;

    /// <summary>The section table, in the order of the file's table.</summary>
    internal IReadOnlyList<SesaSection> Sections => _parsed.Sections;

    /// <summary>What the opening did and cost.</summary>
    internal SesaOpenReport Report => _report;

    /// <summary>The header level's findings: the sections by kind and the derived facts.</summary>
    internal SesaParsed Parsed => _parsed;

    /// <summary>The options the container was opened with.</summary>
    internal EsaOpenOptions Options => _options;

    /// <summary>True in <see cref="EsaAccessMode.Checked"/> mode.</summary>
    internal bool IsChecked => _options.Access == EsaAccessMode.Checked;

    /// <summary>The first byte of the mapped file.</summary>
    internal byte* Data => _data;

    /// <summary>The length of the file.</summary>
    internal long Length => _length;

    /// <summary>The packed text, as code units.</summary>
    internal EsaText TextView => _text;

    /// <summary>The packed text, as dense symbols.</summary>
    internal MappedTextReader TextReader => _reader;

    /// <summary>The suffix array.</summary>
    internal SaView Sa => _sa;

    /// <summary>The LCP array.</summary>
    internal CompactLcp Lcp => _lcp;

    /// <summary>The child table (layout A; codec 1 packed, codec 2 tagged).</summary>
    internal ChildTable Child => _child;

    /// <summary>
    /// True when the 32-bit views that the descent of S1 reads exist (fewer than 2^31 - 1 symbols, 4-byte entries, an LCP array of 8, 16 or 32 bits): <see cref="MappedSearch"/> then runs its 32-bit loop, else the long descent
    /// (<see cref="WideSearcher"/>). The same condition made the 32-bit forms of <see cref="Lcp"/> and <see cref="Child"/>.
    /// </summary>
    internal bool IsNarrow => _narrow;

    /// <summary>The sections that the long descent reads (<see cref="WideSearcher.Create"/>), as pointers into the mapping.</summary>
    internal WideSections WideSections() =>
        new(IsChecked, _parsed.Header.PosWidth, _parsed.LcpKind, _parsed.Child.Codec == SesaFormat.ChildCodecTagged, _data + _parsed.Sa.Offset, _data + _parsed.Child.Offset, _data + _parsed.Lcp.Offset, _lcp, _data + _parsed.Text.Offset, _parsed.Bits, _parsed.Header.TextLength);

    /// <summary>The header statistics in the form of the in-memory index (LRS_RANK is 0 when there is no repeat); every number is 64-bit, as the header has them.</summary>
    internal EsaStatistics Statistics
    {
        get
        {
            SesaHeader h = _parsed.Header;
            long n = h.TextLength;
            if (n == 0) return new EsaStatistics(0, 1, 0, 0, 0, 0, 0, 0);
            return new EsaStatistics(n, h.NodeCount, h.NodeCount - 1 - (n + 1), n, n, h.LrsDepth, h.LrsDepth == 0 ? 0 : h.LrsRank, 0);
        }
    }

    /// <summary>
    /// The sections that the navigation pack reads in place (work package S2.5): PSI (kind 12, checked by the header reader) and the RMQ pyramid (kind 9, which the header reader does not look at: an optional kind it does not implement is
    /// ignored, a container without it is accepted as before). The RMQ entry is checked here, when the pack is built: its <c>COUNT</c> must be the number of entries of the levels of the LCP array at its fanout, <c>ELEM_WIDTH</c> the
    /// <c>POS_WIDTH</c>, the codec 1, <c>PARAM0</c> a power of two from 2 to 1024, <c>PARAM1</c> 0 and <c>LENGTH</c> the entries (and the 3 bytes of padding at width 5); a section that is not so is damage, in both access modes. The CRC is
    /// the <c>Sections</c> level's; the values of the levels are not checked by the verifier (a wrong level is seen by the first walk that it misleads, as a wrong LCP value is).
    /// </summary>
    /// <exception cref="InvalidDataException">The RMQ entry does not describe the pyramid of this container.</exception>
    internal NavSections NavigationSections()
    {
        long n = _parsed.Header.TextLength;
        int posWidth = _parsed.Header.PosWidth;
        SaView? psi = _parsed.Psi is { } p ? SaView.FromMapped(_data + p.Offset, n + 1, posWidth, IsChecked) : null;
        SesaSection? rmq = null;
        foreach (SesaSection s in _parsed.Sections)
        {
            if (s.Kind == SesaFormat.KindRmq) rmq = s;
        }
        if (rmq is null) return new NavSections(psi, null, SesaFormat.RmqFanout, IsChecked);
        long fanout = rmq.Param0;
        if (fanout < 2 || fanout > 1024 || (fanout & (fanout - 1)) != 0) throw new InvalidDataException($"the RMQ section has the fanout {fanout}: a power of two from 2 to 1024 is expected: the file is damaged");
        long expected = n == 0 ? 0 : RmqWriter.CountFor(n + 1, checked((int)fanout));
        long padding = posWidth == 5 ? PosReader.Padding40 : 0;
        if (rmq.Count != expected || rmq.ElemWidth != posWidth || rmq.Codec != 1 || rmq.Param1 != 0 || rmq.Length != (expected * posWidth) + padding)
            throw new InvalidDataException($"the RMQ entry has COUNT {rmq.Count}, ELEM_WIDTH {rmq.ElemWidth}, codec {rmq.Codec}, PARAM1 {rmq.Param1} and LENGTH {rmq.Length}; the pyramid of {n + 1} LCP entries at fanout {fanout} has {expected} entries of {posWidth} bytes: the file is damaged");
        return new NavSections(psi, expected == 0 ? null : SaView.FromMapped(_data + rmq.Offset, expected, posWidth, IsChecked), checked((int)fanout), IsChecked);
    }

    /// <summary>The pointer to the PSI section, or null when the container has none.</summary>
    internal uint* PsiPointer => _parsed.Psi is { } s ? (uint*)(_data + s.Offset) : null;

    /// <summary>The payload of a section.</summary>
    internal ReadOnlySpan<byte> Payload(SesaSection section)
    {
        ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);
        if (section.Length > int.MaxValue) throw new NotSupportedException($"the payload of section {section.Kind} has {section.Length} bytes, more than a span holds ({int.MaxValue}); read it in blocks through Data");
        return new ReadOnlySpan<byte>(_data + section.Offset, (int)section.Length);
    }

    /// <summary>Stores what the verifier measured.</summary>
    internal void SetReport(SesaOpenReport report) => _report = report;

    /// <summary>Releases the mapping (idempotent). Views and trees obtained from this file must not be used afterwards.</summary>
    public void Dispose()
    {
        if (Interlocked.Exchange(ref _disposed, 1) != 0) return;
        _reader.Dispose();
        if (_view is not null)
        {
            _view.SafeMemoryMappedViewHandle.ReleasePointer();
            _view.Dispose();
        }
        _map?.Dispose();
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
}
