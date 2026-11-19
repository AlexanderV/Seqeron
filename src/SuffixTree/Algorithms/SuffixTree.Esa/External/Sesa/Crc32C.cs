namespace SuffixTree.Esa.External;

/// <summary>
/// CRC-32C (Castagnoli, polynomial 0x1EDC6F41, reflected 0x82F63B78, initial value and final xor 0xFFFFFFFF) as ADR-003 section 6 wants it for the section and header checksums: the hardware instruction where the CPU
/// has it (SSE 4.2 on x86-64, the CRC extension on arm64), a slice-by-8 managed table path everywhere else; both give the same value. The functions chain: <c>Append(Append(0, a), b)</c> is the CRC of <c>a</c> followed
/// by <c>b</c>, and <c>Append(0, empty)</c> is 0.
/// </summary>
internal static unsafe class Crc32C
{
    /// <summary>True where <see cref="AppendHardware"/> can run.</summary>
    internal static bool HardwareAccelerated => throw new NotImplementedException();

    /// <summary>The CRC of <paramref name="data"/> continued from <paramref name="crc"/> (0 to start), on the hardware path where there is one.</summary>
    internal static uint Append(uint crc, ReadOnlySpan<byte> data) => throw new NotImplementedException();

    /// <summary>The same for memory at <paramref name="data"/> (a mapped view); <paramref name="length"/> may exceed 2^31.</summary>
    internal static uint Append(uint crc, byte* data, long length) => throw new NotImplementedException();

    /// <summary>The managed slice-by-8 path.</summary>
    internal static uint AppendTable(uint crc, ReadOnlySpan<byte> data) => throw new NotImplementedException();

    /// <summary>The hardware path.</summary>
    /// <exception cref="PlatformNotSupportedException">The CPU has no CRC-32C instruction.</exception>
    internal static uint AppendHardware(uint crc, ReadOnlySpan<byte> data) => throw new NotImplementedException();

    /// <summary>The CRC of <paramref name="data"/>.</summary>
    internal static uint Compute(ReadOnlySpan<byte> data) => Append(0, data);
}
