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
        _reader = MappedTextReader.FromMemory(data + parsed.Text.Offset, parsed.Text.Length, n, parsed.Bits);
        _text = EsaText.FromPacked(_reader, parsed.Letters, (int)Math.Min(n, int.MaxValue));
        if (n < int.MaxValue)
        {
            // the 32-bit views; a container of 2^31 - 1 symbols or more can be mapped and its header read (the adapter of e8-verify reads such sections itself) but has no tree here
            _sa = new SaView((int*)(data + parsed.Sa.Offset), (int)(n + 1), isChecked);
            uint* pairs = parsed.Overflow is { } o ? (uint*)(data + o.Offset) : null;
            _lcp = CompactLcp.FromMapped(parsed.LcpKind, (int)(n + 1), data + parsed.Lcp.Offset, pairs, (int)parsed.Lcp.Param0, isChecked);
            _child = ChildTable.FromMapped(new MappedChild((int*)(data + parsed.Child.Offset), (int)(n + 1), _lcp, isChecked));
        }
        else
        {
            _lcp = null!;
            _child = null!;
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
        if (!Enum.IsDefined(options.Validation)) throw new ArgumentOutOfRangeException(nameof(options), options.Validation, "not a validation level");
        if (!Enum.IsDefined(options.Access)) throw new ArgumentOutOfRangeException(nameof(options), options.Access, "not an access mode");
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

    /// <summary>The child table (codec 1, layout A).</summary>
    internal ChildTable Child => _child;

    /// <summary>The header statistics in the form of the in-memory index (LRS_RANK is 0 when there is no repeat).</summary>
    internal EsaStatistics Statistics
    {
        get
        {
            SesaHeader h = _parsed.Header;
            int n = (int)h.TextLength;
            if (n == 0) return new EsaStatistics(0, 1, 0, 0, 0, 0, 0, 0);
            return new EsaStatistics(n, h.NodeCount, (int)(h.NodeCount - 1 - (n + 1L)), n, n, (int)h.LrsDepth, h.LrsDepth == 0 ? 0 : (int)h.LrsRank, 0);
        }
    }

    /// <summary>The pointer to the PSI section, or null when the container has none.</summary>
    internal uint* PsiPointer => _parsed.Psi is { } s ? (uint*)(_data + s.Offset) : null;

    /// <summary>The payload of a section.</summary>
    internal ReadOnlySpan<byte> Payload(SesaSection section)
    {
        ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);
        return new ReadOnlySpan<byte>(_data + section.Offset, checked((int)section.Length));
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
