namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>Enumeration of every string over the first <c>sigma</c> letters of 'a'.. up to a length bound.</summary>
public static class ExhaustiveStrings
{
    /// <summary>All strings of length 0..<paramref name="maxLength"/>, shorter first, lexicographic within a length.</summary>
    public static IEnumerable<string> All(int sigma, int maxLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(sigma, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(maxLength);
        for (int len = 0; len <= maxLength; len++)
        {
            var digits = new int[len];
            var chars = new char[len];
            Array.Fill(chars, 'a');
            while (true)
            {
                yield return new string(chars);
                int pos = len - 1;
                while (pos >= 0 && digits[pos] == sigma - 1)
                {
                    digits[pos] = 0;
                    chars[pos] = 'a';
                    pos--;
                }
                if (pos < 0) break;
                digits[pos]++;
                chars[pos] = (char)('a' + digits[pos]);
            }
        }
    }
}
