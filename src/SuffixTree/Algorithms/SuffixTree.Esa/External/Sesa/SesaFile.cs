namespace SuffixTree.Esa.External;

/// <summary>What the opening of a container did and what it cost (reported by <see cref="SesaFile.Report"/>).</summary>
/// <param name="Level">The validation level that ran.</param>
/// <param name="Access">The access mode.</param>
/// <param name="HeaderSeconds">Time of the header level (mapping included).</param>
/// <param name="SectionsSeconds">Time of the section checksums (0 below that level).</param>
/// <param name="StructureSeconds">Time of the structure pass (0 below that level).</param>
/// <param name="HardwareCrc">Whether the CRC-32C used the hardware instruction.</param>
/// <param name="SectionBytesChecksummed">Payload bytes that the checksum pass read.</param>
/// <param name="BitmapBytes">Size of the permutation bitmap of the structure pass (n / 8 rounded up).</param>
/// <param name="BufferBytes">Anonymous memory of the structure pass besides the bitmap (buffers, stacks), at its peak.</param>
internal sealed record SesaOpenReport(EsaValidationLevel Level, EsaAccessMode Access, double HeaderSeconds, double SectionsSeconds, double StructureSeconds, bool HardwareCrc, long SectionBytesChecksummed, long BitmapBytes, long BufferBytes);

/// <summary>A SESA container mapped read-only as one view (ADR-003 section 9): the header, the section table and the typed views of the sections. Throwing stub of the tests-first commit.</summary>
internal sealed unsafe class SesaFile : IDisposable
{
    /// <summary>Maps <paramref name="path"/> and validates it at the level of <paramref name="options"/>; <paramref name="hardwareCrc"/> forces (true) or forbids (false) the CRC instruction, null lets the CPU decide.</summary>
    internal static SesaFile Open(string path, EsaOpenOptions options, bool? hardwareCrc = null) => throw new NotImplementedException();

    /// <summary>The container as it lies in memory that the caller owns (tests: memory with guard pages): the same checks, nothing is mapped or freed.</summary>
    internal static SesaFile FromMemory(byte* data, long length, EsaOpenOptions options, bool? hardwareCrc = null) => throw new NotImplementedException();

    private bool _stub;

    internal SesaHeader Header => _stub ? new SesaHeader() : throw new NotImplementedException();

    internal IReadOnlyList<SesaSection> Sections => _stub ? [] : throw new NotImplementedException();

    internal SesaOpenReport Report => _stub ? null! : throw new NotImplementedException();

    public void Dispose() => _stub = !_stub ? throw new NotImplementedException() : _stub;
}
