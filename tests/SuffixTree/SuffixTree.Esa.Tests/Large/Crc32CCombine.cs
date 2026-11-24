namespace SuffixTree.Esa.Tests.Large;

/// <summary>STUB of work package S2.1 (tests first): CRC-32C by the algebra of the polynomial is written after its tests.</summary>
internal static class Crc32CCombine
{
    internal static uint Multiply(uint a, uint b) => throw new NotImplementedException("S2.1 stub");

    internal static uint XPower(ulong exponent) => throw new NotImplementedException("S2.1 stub");

    internal static uint Shift(uint crc, long zeroBytes) => throw new NotImplementedException("S2.1 stub");

    internal static uint Combine(uint crcOfFirst, uint crcOfSecond, long lengthOfSecond) => throw new NotImplementedException("S2.1 stub");

    internal static uint OfZeros(long length) => throw new NotImplementedException("S2.1 stub");

    internal static uint OfSparse(long length, IReadOnlyList<(long Offset, byte[] Data)> patches) => throw new NotImplementedException("S2.1 stub");
}
