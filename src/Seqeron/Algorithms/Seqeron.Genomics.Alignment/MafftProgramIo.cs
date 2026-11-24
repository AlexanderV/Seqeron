#pragma warning disable IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The input and output text of the MAFFT 7.505 programs, as the C code of <c>io.c</c> reads and writes it: a record starts with a '&gt;' at the start of a line, the name is the first 254 characters of that line
/// (<c>myfgets( name+1, B-2, fp )</c>), the sequence is the rest of the lines up to the next record; the output is FASTA with 60 residues per line (<c>writeData_pointer</c>), the Clustal and PHYLIP formats of <c>f2cl</c>
/// (<c>clustalout_pointer</c>, <c>phylipout_pointer</c>, <c>setmark_clustal</c>). Text is treated as bytes (every <see cref="char"/> is one byte 0 - 255), as the programs do.
/// </summary>
internal static class MafftProgramIo
{
    /// <summary><c>B</c> - 2 of <c>mltaln.h</c>: the number of characters of a title line that the programs keep.</summary>
    internal const int NameLimit = 254;

    /// <summary><c>C</c> of <c>mltaln.h</c>: residues per line of the FASTA output.</summary>
    internal const int LineWidth = 60;

    /// <summary>The text of a file as the programs see it: the .NET string as UTF-8 bytes, one <see cref="char"/> per byte.</summary>
    internal static string ToBytes(string s)
    {
        bool ascii = true;
        foreach (char c in s)
        {
            if (c >= 128)
            {
                ascii = false;
                break;
            }
        }
        if (ascii)
            return s;
        return Encoding.Latin1.GetString(Encoding.UTF8.GetBytes(s));
    }

    /// <summary>The inverse of <see cref="ToBytes"/> (invalid UTF-8 becomes U+FFFD).</summary>
    internal static string FromBytes(string s)
    {
        bool ascii = true;
        foreach (char c in s)
        {
            if (c >= 128)
            {
                ascii = false;
                break;
            }
        }
        if (ascii)
            return s;
        return Encoding.UTF8.GetString(Encoding.Latin1.GetBytes(s));
    }

    /// <summary><c>cat file | tr "\r" "\n"</c>.</summary>
    internal static string CarriageReturnsToNewlines(string s) => s.Replace('\r', '\n');

    /// <summary><c>grep -v "^$"</c> of a text (blank lines are removed; the last line keeps its newline when it had one).</summary>
    internal static string DropEmptyLines(string s)
    {
        var sb = new StringBuilder(s.Length);
        foreach (var line in SplitLines(s))
        {
            if (line.Length == 0)
                continue;
            sb.Append(line).Append('\n');
        }
        return sb.ToString();
    }

    /// <summary>The lines of a text (a final line without a newline counts; a final newline does not start another line).</summary>
    internal static List<string> SplitLines(string s)
    {
        var lines = new List<string>();
        int start = 0;
        for (int i = 0; i < s.Length; i++)
        {
            if (s[i] == '\n')
            {
                lines.Add(s[start..i]);
                start = i + 1;
            }
        }
        if (start < s.Length)
            lines.Add(s[start..]);
        return lines;
    }

    /// <summary>A record of a FASTA file as <c>io.c</c> reads it.</summary>
    internal sealed record Record(string Name, string Raw);

    /// <summary>The records of a FASTA text (<c>searchKUorWA</c>, <c>myfgets</c>, <c>load1SeqWithoutName_realloc</c>): the name is the title line without '&gt;' cut to <see cref="NameLimit"/> characters; the raw sequence is everything up to the next '&gt;' at the start of a line.</summary>
    internal static List<Record> Records(string text)
    {
        var records = new List<Record>();
        int n = text.Length;
        int i = 0;
        // searchKUorWA: the first '>' that follows a newline (or starts the file)
        while (i < n && !(text[i] == '>' && (i == 0 || text[i - 1] == '\n')))
            i++;
        while (i < n)
        {
            int eol = text.IndexOf('\n', i + 1);
            string title;
            int bodyStart;
            if (eol < 0)
            {
                title = text[(i + 1)..];
                bodyStart = n;
            }
            else
            {
                title = text[(i + 1)..eol];
                bodyStart = eol + 1;
            }
            if (title.Length > NameLimit)
                title = title[..NameLimit];
            int j = bodyStart;
            while (j < n && !(text[j] == '>' && text[j - 1] == '\n'))
                j++;
            records.Add(new Record(title, text[bodyStart..j]));
            i = j;
        }
        return records;
    }

    /// <summary>The number of records the script counts (<c>grep -c '^[&gt;|=]'</c>).</summary>
    internal static int CountRecordLines(string text)
    {
        int count = 0;
        foreach (var line in SplitLines(text))
            if (line.Length > 0 && (line[0] is '>' or '|' or '='))
                count++;
        return count;
    }

    /// <summary><c>writeData_pointer</c> / <c>writeData_reorder_pointer</c>.</summary>
    internal static string Fasta(IReadOnlyList<string> names, IReadOnlyList<string> rows)
    {
        var sb = new StringBuilder();
        for (int i = 0; i < rows.Count; i++)
        {
            sb.Append('>').Append(names[i]).Append('\n');
            string r = rows[i];
            for (int j = 0; j < r.Length; j += LineWidth)
                sb.Append(r, j, Math.Min(LineWidth, r.Length - j)).Append('\n');
        }
        return sb.ToString();
    }

    /// <summary><c>extractfirstword</c> of <c>io.c</c> (tabs become blanks, leading blanks are skipped, the word ends at the next blank).</summary>
    internal static string FirstWord(string name)
    {
        name = name.Replace('\t', ' ');
        int i = 0;
        while (i < name.Length && name[i] == ' ')
            i++;
        int j = i;
        while (j < name.Length && name[j] != ' ')
            j++;
        return name[i..j];
    }

    /// <summary><c>printf( "%-*.*s", width, precision, s )</c>: a negative width means left justification, a negative precision is ignored.</summary>
    internal static string PadName(string s, int width, int precision)
    {
        if (precision >= 0 && s.Length > precision)
            s = s[..precision];
        int w = Math.Abs(width);
        return s.Length >= w ? s : s + new string(' ', w - s.Length);
    }

    /// <summary><c>f2cl</c> with <c>-c strategy</c> (Clustal format): the sequences are padded with blanks to the longest, the marks of <c>setmark_clustal</c> are computed from the letters.</summary>
    internal static string Clustal(string? comment, IReadOnlyList<string> names, IReadOnlyList<string> rows, int namelen, bool nucleotide)
    {
        var (seq, maxlen) = Pad(rows);
        string mark = SetMarkClustal(maxlen, seq, nucleotide);
        if (namelen == -1)
            namelen = 15;
        var sb = new StringBuilder();
        if (comment == null)
            sb.Append("CLUSTAL format alignment by MAFFT (v").Append(MafftVersionNumber).Append(")\n\n");
        else
            sb.Append("CLUSTAL format alignment by MAFFT ").Append(comment).Append(" (v").Append(MafftVersionNumber).Append(")\n\n");
        for (int pos = 0; pos < maxlen; pos += 60)
        {
            sb.Append('\n');
            for (int j = 0; j < seq.Length; j++)
            {
                sb.Append(PadName(FirstWord(names[j]), namelen, namelen)).Append(' ');
                sb.Append(Slice(seq[j], pos, 60)).Append('\n');
            }
            sb.Append(PadName("", namelen, namelen)).Append(' ');
            sb.Append(Slice(mark, pos, 60)).Append('\n');
        }
        return sb.ToString();
    }

    /// <summary><c>f2cl -y</c>: the sequential PHYLIP format (blocks of 10 residues, 50 per line).</summary>
    internal static string Phylip(IReadOnlyList<string> names, IReadOnlyList<string> rows, int namelen)
    {
        var (seq, maxlen) = Pad(rows);
        if (namelen == -1)
            namelen = 10;
        var sb = new StringBuilder();
        sb.Append(' ').Append(seq.Length.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(maxlen.ToString(CultureInfo.InvariantCulture)).Append('\n');
        for (int pos = 0; pos < maxlen; pos += 50)
        {
            for (int j = 0; j < seq.Length; j++)
            {
                sb.Append(PadName(pos == 0 ? FirstWord(names[j]) : "", namelen, namelen));
                int pos2 = pos;
                while (pos2 < pos + 41 && pos2 < maxlen)
                {
                    sb.Append(' ').Append(Slice(seq[j], pos2, 10));
                    pos2 += 10;
                }
                sb.Append('\n');
            }
            sb.Append('\n');
        }
        return sb.ToString();
    }

    /// <summary>The version text of the programs (<c>VERSION</c> of <c>mltaln.h</c>).</summary>
    internal const string MafftVersionNumber = "7.505";

    private static string Slice(string s, int pos, int len) => pos >= s.Length ? "" : s.Substring(pos, Math.Min(len, s.Length - pos));

    private static (string[] Seq, int MaxLen) Pad(IReadOnlyList<string> rows)
    {
        int max = 0;
        foreach (var r in rows)
            max = Math.Max(max, r.Length);
        var seq = new string[rows.Count];
        for (int i = 0; i < seq.Length; i++)
            seq[i] = rows[i].PadRight(max, ' ');
        return (seq, max);
    }

    private static readonly string[] StrongProtein = { "STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW" };
    private static readonly string[] WeakerProtein = { "CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY" };
    private static readonly string[] StrongDna = { "TU" };
    private static readonly string[] WeakerDna = { "AG", "CT", "CU" };

    /// <summary><c>setmark_clustal</c>: '*' identical, ':' strong, '.' weaker group; letters outside the first <c>nalpha</c> letters of the program's alphabet give no mark.</summary>
    private static string SetMarkClustal(int nlen, string[] seq, bool nucleotide)
    {
        var g = SequenceAligner.MafftGlobals.Parse(new[] { nucleotide ? "-D" : "-P" });
        g.Dorp = nucleotide ? 'd' : 'p';
        g.Constants(Array.Empty<string>());
        string[] strong = nucleotide ? StrongDna : StrongProtein;
        string[] weaker = nucleotide ? WeakerDna : WeakerProtein;
        int nalpha = nucleotide ? 10 : 20;
        var mark = new char[nlen];
        int nseq = seq.Length;
        for (int i = 0; i < nlen; i++)
        {
            mark[i] = ' ';
            int j;
            for (j = 0; j < nseq; j++)
            {
                char s = seq[j][i];
                if (s == '-' || s == ' ')
                    break;
            }
            if (j != nseq)
                continue;
            char first = char.ToUpperInvariant(seq[0][i]);
            int an = first < 256 ? g.AminoN[first] : -1;
            if (an >= nalpha || an < 0)
                continue;
            for (j = 0; j < nseq; j++)
                if (char.ToUpperInvariant(seq[j][i]) != first)
                    break;
            if (j == nseq)
            {
                mark[i] = '*';
                continue;
            }
            if (AllIn(strong, seq, i))
            {
                mark[i] = ':';
                continue;
            }
            if (AllIn(weaker, seq, i))
                mark[i] = '.';
        }
        return new string(mark);
    }

    private static bool AllIn(string[] groups, string[] seq, int i)
    {
        foreach (var grp in groups)
        {
            int j;
            for (j = 0; j < seq.Length; j++)
                if (grp.IndexOf(char.ToUpperInvariant(seq[j][i])) < 0)
                    break;
            if (j == seq.Length)
                return true;
        }
        return false;
    }
}
