namespace Seqeron.Genomics.Tests;

/// <summary>The escape convention of the oracle TSV files (<c>\n</c>, <c>\t</c>, <c>\r</c>, any other escaped character stands for itself); shared by the alignment test projects (linked source file).</summary>
internal static class OracleEscapes
{
    internal static string Unescape(string s)
    {
        var sb = new System.Text.StringBuilder();
        int i = 0;
        while (i < s.Length)
        {
            if (s[i] != '\\' || i + 1 >= s.Length)
            {
                sb.Append(s[i++]);
                continue;
            }

            char n = s[i + 1];
            sb.Append(n switch { 'n' => '\n', 't' => '\t', 'r' => '\r', _ => n });
            i += 2;
        }

        return sb.ToString();
    }
}
