using System.Text;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The string formatting of the AJAX library that the sequence writers use to wrap the text of the entries: <c>ajStrFmtWrap</c> and <c>ajStrFmtWrapAt</c> (<c>ajstr.c</c>). A5-8c3h, B09 F267.</summary>
internal static class EmbossStrFmt
{
    private static bool IsSpace(char c) => c is ' ' or '\t' or '\n' or '\v' or '\f' or '\r';

    /// <summary><c>ajStrFmtWrap</c>: newlines are inserted at white space (or at the width when there is none) so that the lines are at most <paramref name="width"/> characters.</summary>
    public static string Wrap(string s, int width) => Wrap(s, width, '\0', false);

    /// <summary><c>ajStrFmtWrapAt</c>: as <see cref="Wrap(string, int)"/>, with a preference for the white space that follows the character <paramref name="ch"/>.</summary>
    public static string WrapAt(string s, int width, char ch) => Wrap(s, width, ch, true);

    private static string Wrap(string s, int width, char ch, bool preferred)
    {
        if (width > s.Length)
            return s;
        var t = new StringBuilder(s);
        long i = 0;
        long imax = t.Length - width;
        char At(long p) => p >= 0 && p < t.Length ? t[(int)p] : '\0';
        while (i < imax)
        {
            long j = i + width + 1;
            if (j > t.Length)
                j = t.Length;
            long k = j, kk = j;
            long cq = i;
            while (i < j)
            {
                if (IsSpace(At(cq)))
                {
                    k = i;
                    if (At(cq) == '\n')
                        break;
                    if (preferred && i != 0 && At(i - 1) == ch)
                        kk = i;
                }

                cq++;
                i++;
            }

            if (At(cq) != '\n')
            {
                if (preferred && kk < j)
                {
                    t[(int)kk] = '\n';
                    k = kk;
                }
                else if (k == j)
                {
                    t.Insert((int)k, '\n');
                    imax++;
                }
                else
                {
                    t[(int)k] = '\n';
                }
            }

            i = k + 1;
        }

        return t.ToString();
    }
}
