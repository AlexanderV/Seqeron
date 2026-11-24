using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// CRC-32C by the algebra of its polynomial (work package S2.1, deliverable B): the CRC of a hole of N zero bytes, and the CRC of a section that is a hole with a few patches, without a pass over N bytes. A CRC is a remainder
/// modulo the generator polynomial P (reflected 0x82F63B78): feeding one zero byte multiplies the register by x^8, so N zero bytes multiply it by <c>x^(8N)</c>, computed by square and multiply in at most 64 steps.
/// <c>CRC(A || B) = CRC(A) * x^(8 |B|) xor CRC(B)</c> (the initial and the final inversion cancel), which is zlib's <c>crc32_combine</c>; the exponent is a 64-bit number (no table of 32 powers, no reliance on the order of x).
/// A value is in the reflected representation of the product (<see cref="Crc32C"/>): bit 31 is the coefficient of x^0, bit 0 of x^31.
/// </summary>
internal static class Crc32CCombine
{
    private const uint ReflectedPolynomial = 0x82F63B78;
    private const uint One = 0x80000000;

    /// <summary><c>a * b mod P</c>.</summary>
    internal static uint Multiply(uint a, uint b)
    {
        uint product = 0;
        uint m = One; // the coefficient of x^i of a is tested at the bit 31 - i, while b is multiplied by x each step
        for (int i = 0; i < 32; i++)
        {
            if ((a & m) != 0) product ^= b;
            m >>= 1;
            b = (b & 1) != 0 ? (b >> 1) ^ ReflectedPolynomial : b >> 1;
        }
        return product;
    }

    /// <summary><c>x^exponent mod P</c>.</summary>
    internal static uint XPower(ulong exponent)
    {
        uint result = One;
        uint square = One >> 1; // x
        while (exponent != 0)
        {
            if ((exponent & 1) != 0) result = Multiply(square, result);
            square = Multiply(square, square);
            exponent >>= 1;
        }
        return result;
    }

    /// <summary>The CRC register state after <paramref name="zeroBytes"/> more zero bytes: <c>crc * x^(8 zeroBytes)</c>.</summary>
    internal static uint Shift(uint crc, long zeroBytes)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(zeroBytes);
        return Multiply(XPower(checked(8UL * (ulong)zeroBytes)), crc);
    }

    /// <summary>The CRC of the concatenation of two pieces, from their CRCs and the length of the second.</summary>
    internal static uint Combine(uint crcOfFirst, uint crcOfSecond, long lengthOfSecond) => Shift(crcOfFirst, lengthOfSecond) ^ crcOfSecond;

    /// <summary>The CRC-32C of <paramref name="length"/> zero bytes: the register <c>0xFFFFFFFF * x^(8 length)</c>, inverted.</summary>
    internal static uint OfZeros(long length) => Shift(0xFFFFFFFF, length) ^ 0xFFFFFFFF;

    /// <summary>
    /// The CRC-32C of a section of <paramref name="length"/> bytes that is zero except for the <paramref name="patches"/> (offset, bytes), given in ascending order without overlap; the zero runs between them are combined, never read.
    /// </summary>
    /// <exception cref="ArgumentException">A patch is before the previous one, overlaps it or lies outside the section.</exception>
    internal static uint OfSparse(long length, IReadOnlyList<(long Offset, byte[] Data)> patches)
    {
        uint crc = 0;
        long at = 0;
        foreach ((long offset, byte[] data) in patches)
        {
            if (offset < at || offset > length - data.Length) throw new ArgumentException($"the patch at {offset} (+{data.Length}) overlaps the previous one (end {at}), is out of order or leaves the section of {length} bytes", nameof(patches));
            if (offset > at) crc = Combine(crc, OfZeros(offset - at), offset - at);
            crc = Combine(crc, Crc32C.Compute(data), data.Length);
            at = offset + data.Length;
        }
        if (length > at) crc = Combine(crc, OfZeros(length - at), length - at);
        return crc;
    }
}
