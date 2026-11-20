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

    /// <summary>
    /// C <c>printf("%.Ng")</c> / <c>ostream &lt;&lt; setprecision(N)</c> (default float field) of a double: <paramref name="precision"/> significant digits from the exact binary value (ties to even),
    /// scientific notation (<c>e+XX</c>, at least two exponent digits) when the decimal exponent is below -4 or at least the precision, trailing zeros removed.
    /// </summary>
    /// <param name="value">The number.</param>
    /// <param name="precision">Significant digits (0 is taken as 1, as C does).</param>
    public static string General(double value, int precision)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(precision);
        int p = Math.Max(precision, 1);
        if (double.IsNaN(value))
            return double.IsNegative(value) ? "-nan" : "nan"; // glibc prints the sign bit of a NaN (x86 produces the negative default NaN for 0 / 0)
        if (double.IsInfinity(value))
            return value > 0 ? "inf" : "-inf";
        bool negative = double.IsNegative(value);
        if (value == 0)
            return negative ? "-0" : "0";

        long bits = BitConverter.DoubleToInt64Bits(Math.Abs(value));
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = bits & 0xF_FFFF_FFFF_FFFFL;
        if (exponent == 0)
            exponent = 1;
        else
            mantissa |= 1L << 52;
        exponent -= 1075; // |value| = mantissa * 2^exponent

        // the value as a fraction num / den
        BigInteger num = mantissa, den = BigInteger.One;
        if (exponent >= 0)
            num <<= exponent;
        else
            den <<= -exponent;

        // x = floor(log10(value)): start from an estimate and correct it
        int x = (int)Math.Floor(Math.Log10(Math.Abs(value)));
        while (Compare(num, den, x) < 0)
            x--;
        while (Compare(num, den, x + 1) >= 0)
            x++;

        BigInteger digits = RoundScaled(num, den, p - 1 - x);
        if (digits >= BigInteger.Pow(10, p))
        {
            x++;
            digits = RoundScaled(num, den, p - 1 - x);
        }

        string ds = digits.ToString(CultureInfo.InvariantCulture); // exactly p digits
        string text;
        if (x < -4 || x >= p)
        {
            string mant = ds.Length > 1 ? ds[..1] + "." + ds[1..] : ds;
            if (mant.Contains('.'))
                mant = mant.TrimEnd('0').TrimEnd('.');
            text = mant + "e" + (x < 0 ? "-" : "+") + Math.Abs(x).ToString("00", CultureInfo.InvariantCulture);
        }
        else if (x >= 0)
        {
            string intPart = ds[..(x + 1)];
            string frac = ds[(x + 1)..].TrimEnd('0');
            text = frac.Length == 0 ? intPart : intPart + "." + frac;
        }
        else
        {
            string frac = (new string('0', -x - 1) + ds).TrimEnd('0');
            text = "0." + frac;
        }

        return negative ? "-" + text : text;
    }

    // compares num/den with 10^e
    private static int Compare(BigInteger num, BigInteger den, int e)
        => e >= 0 ? num.CompareTo(den * BigInteger.Pow(10, e)) : (num * BigInteger.Pow(10, -e)).CompareTo(den);

    // round-half-even(num / den * 10^k)
    private static BigInteger RoundScaled(BigInteger num, BigInteger den, int k)
    {
        if (k >= 0)
            num *= BigInteger.Pow(10, k);
        else
            den *= BigInteger.Pow(10, -k);
        var q = BigInteger.DivRem(num, den, out var r);
        var twice = r * 2;
        if (twice > den || (twice == den && !q.IsEven))
            q += 1;
        return q;
    }
}
