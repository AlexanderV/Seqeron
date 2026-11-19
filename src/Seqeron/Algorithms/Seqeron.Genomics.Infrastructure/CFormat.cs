using System.Globalization;
using System.Numerics;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>Number formatting of the C / C++ standard library, as the reference programs print numbers (shared by the ports that must reproduce their files byte for byte).</summary>
public static class CFormat
{
    /// <summary>
    /// C <c>printf("%.Nf")</c> / <c>ostream &lt;&lt; fixed &lt;&lt; setprecision(N)</c> of a double: the exact binary value rounded to <paramref name="decimals"/> places, ties to even (glibc, default rounding mode).
    /// .NET's "F" format rounds ties away from zero instead (<c>0.015625</c> with five decimals is <c>0.01562</c> in C, <c>0.01563</c> in .NET).
    /// </summary>
    /// <param name="value">The number.</param>
    /// <param name="decimals">Digits after the decimal point (0 or more).</param>
    public static string Fixed(double value, int decimals)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(decimals);
        if (double.IsNaN(value))
            return "nan";
        if (double.IsInfinity(value))
            return value > 0 ? "inf" : "-inf";

        bool negative = double.IsNegative(value);
        long bits = BitConverter.DoubleToInt64Bits(Math.Abs(value));
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = bits & 0xF_FFFF_FFFF_FFFFL;
        if (exponent == 0)
            exponent = 1;
        else
            mantissa |= 1L << 52;
        exponent -= 1075; // |value| = mantissa * 2^exponent

        var scaled = mantissa * BigInteger.Pow(10, decimals);
        BigInteger units;
        if (exponent >= 0)
        {
            units = scaled << exponent;
        }
        else
        {
            var denominator = BigInteger.One << -exponent;
            units = BigInteger.DivRem(scaled, denominator, out var remainder);
            var twice = remainder * 2;
            if (twice > denominator || (twice == denominator && !units.IsEven))
                units += 1;
        }

        string digits = units.ToString(CultureInfo.InvariantCulture).PadLeft(decimals + 1, '0');
        string text = decimals == 0 ? digits : digits[..^decimals] + "." + digits[^decimals..];
        return negative ? "-" + text : text;
    }
}
