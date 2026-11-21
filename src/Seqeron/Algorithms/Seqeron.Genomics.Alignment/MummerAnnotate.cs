using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Output of <see cref="MummerAnnotate.Run"/>.</summary>
/// <param name="Text">What the program writes to standard output: the gaps file with the alignment of every gap under its line.</param>
/// <param name="WithErrors">The file <c>witherrors.gaps</c> the program writes into the current directory: the lines of the gaps file with the number of errors of each gap appended.</param>
public sealed record MummerAnnotateResult(string Text, string WithErrors);

/// <summary>
/// <c>mummer-annotate</c> (<c>annotate</c>, <c>src/tigr/annotate.cc</c> of MUMmer 3.23, Delcher): reads the output of <see cref="MummerGaps"/> (the first sequence is named on its <c>&gt;</c> lines, the word <c>reverse</c> there
/// reverse-complements the query sequence), and under every gap line prints the edit-distance alignment of the stretch between the two matches that surround it, with 10 bases of the flanking matches and the number of errors;
/// <c>witherrors.gaps</c> receives the gap lines with the number of errors appended.
/// </summary>
/// <remarks>
/// The program reads <c>Data1 [Start1 - 1 + ...]</c> without a bounds test: a window that reaches before the first base or past the end of a sequence reads the memory around the buffer. Beyond the end the buffer holds zero bytes
/// (the port reads 0); a window before the start reads the heap header (the program's result is undefined, valgrind reports an invalid read): the port refuses it with <see cref="InvalidOperationException"/>.
/// </remarks>
public static class MummerAnnotate
{
    private const int MaxLineLength = 200, PrefixLength = 10, Width = 60, MaxAlign = 10000;

    private sealed class Seq
    {
        public byte[] Bytes = [0, 0];   // 1-based, zero padded
        public long Len;

        public byte this[long index]
        {
            get
            {
                if (index < 1)
                    throw new InvalidOperationException("mummer-annotate: a gap window reaches before the first base of a sequence; the program reads memory in front of its buffer (undefined behaviour).");
                return index <= Len ? Bytes[index] : (byte)0;
            }
        }
    }

    private static Seq FirstRecord(string fasta)
    {
        var reader = new MummerFastaReader(fasta);
        var s = new Seq();
        if (reader.Read(out _, out byte[] seq, out long len, new StringBuilder()))
        {
            s.Bytes = seq;
            s.Len = len;
        }

        return s;
    }

    /// <summary>Runs <c>mummer-annotate &lt;gapfile&gt; &lt;datafile&gt;</c>.</summary>
    /// <param name="gaps">The text of the gap file.</param>
    /// <param name="dataFasta">The text of the data file (the query: its first record is used).</param>
    /// <param name="files">The files named on the <c>&gt;</c> lines of the gap file (the first sequence), by name.</param>
    /// <exception cref="FormatException">A file named in the gap file is missing (the program exits with an error).</exception>
    /// <exception cref="InvalidOperationException">A gap window reaches before the start of a sequence (undefined in the program).</exception>
    public static MummerAnnotateResult Run(string gaps, string dataFasta, IReadOnlyDictionary<string, string> files)
    {
        ArgumentNullException.ThrowIfNull(gaps);
        ArgumentNullException.ThrowIfNull(dataFasta);
        ArgumentNullException.ThrowIfNull(files);
        var ci = CultureInfo.InvariantCulture;
        var data2 = FirstRecord(dataFasta);
        var data1 = new Seq();
        var outText = new StringBuilder();
        var errors = new StringBuilder();
        string olap = string.Empty, gap1s = string.Empty, gap2s = string.Empty;
        long start1 = 0, start2 = 0, len = 0;

        foreach (string chunk in Fgets(gaps))
        {
            string line = chunk;
            outText.Append(line);
            int i = 0;
            while (i < line.Length && MummerFastaReader.IsSpace(line[i]) && line[i] != '\n')
                i++;
            char c = i < line.Length ? line[i] : '\0';
            if (c == '>')
            {
                var tok = Words(line, i + 1);
                string name = tok.Count > 0 ? tok[0] : string.Empty;
                string tag = tok.Count > 1 ? tok[1] : string.Empty;
                if (name == "Wrap")
                    continue;
                if (!files.TryGetValue(name, out string? text))
                    throw new FormatException("ERROR: Could not open file " + name);
                data1 = FirstRecord(text);
                if (tag.ToLowerInvariant() == "reverse")
                {
                    var rc = MummerSequence.ReverseComplement(data2.Bytes);
                    data2.Bytes = rc;
                }

                errors.Append(line);
                continue;
            }

            if (!(c >= '0' && c <= '9'))
            {
                errors.Append(line);
                continue;
            }

            // sscanf (Line, "%ld %ld %ld %s %s %s", ...): fields that are missing keep their old values
            var w = Words(line, 0);
            if (w.Count > 0 && long.TryParse(w[0], NumberStyles.AllowLeadingSign, ci, out long v0))
                start1 = v0;
            if (w.Count > 1 && long.TryParse(w[1], NumberStyles.AllowLeadingSign, ci, out long v1))
                start2 = v1;
            if (w.Count > 2 && long.TryParse(w[2], NumberStyles.AllowLeadingSign, ci, out long v2))
                len = v2;
            if (w.Count > 3)
                olap = w[3];
            if (w.Count > 4)
                gap1s = w[4];
            if (w.Count > 5)
                gap2s = w[5];
            _ = len;

            if (gap1s.Length > 0 && gap1s[0] == '-' || gap2s.Length > 0 && gap2s[0] == '-')
            {
                errors.Append(line);
                continue;
            }

            long gap1 = StrToL(gap1s), gap2 = StrToL(gap2s);
            long olapValue = olap.Length > 0 && char.IsAsciiLetter(olap[0]) ? 0 : StrToL(olap);
            start1 -= gap1 + PrefixLength - olapValue;
            start2 -= gap2 + PrefixLength - olapValue;

            string shown = line.Length > 0 ? line[..^1] : line;   // Line [Line_Len - 1] = '\0'
            ShowAlignment(outText, errors, shown, data1, start1 - 1, gap1 - olapValue + 2 * PrefixLength, data2, start2 - 1, gap2 - olapValue + 2 * PrefixLength);
        }

        return new MummerAnnotateResult(outText.ToString(), errors.ToString());

        static long StrToL(string s)
        {
            int p = 0;
            while (p < s.Length && MummerFastaReader.IsSpace(s[p]))
                p++;
            int st = p;
            if (p < s.Length && (s[p] == '+' || s[p] == '-'))
                p++;
            int d = p;
            while (p < s.Length && s[p] >= '0' && s[p] <= '9')
                p++;
            if (p == d)
                return 0;
            return long.TryParse(s.AsSpan(st, p - st), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long r) ? r : 0;
        }
    }

    /// <summary><c>fgets (Line, MAX_LINE_LEN, file)</c>: at most 199 characters, up to and including the newline.</summary>
    private static IEnumerable<string> Fgets(string text)
    {
        int pos = 0;
        while (pos < text.Length)
        {
            int end = Math.Min(text.Length, pos + MaxLineLength - 1);
            int nl = text.IndexOf('\n', pos, end - pos);
            if (nl >= 0)
                end = nl + 1;
            yield return text.Substring(pos, end - pos);
            pos = end;
        }
    }

    private static List<string> Words(string line, int from)
    {
        var words = new List<string>();
        int i = from;
        while (i < line.Length)
        {
            while (i < line.Length && MummerFastaReader.IsSpace(line[i]))
                i++;
            int j = i;
            while (j < line.Length && !MummerFastaReader.IsSpace(line[j]))
                j++;
            if (j > i)
                words.Add(line.Substring(i, j - i));
            i = j;
        }

        return words;
    }

    private static void ShowAlignment(StringBuilder outText, StringBuilder errors, string line, Seq sa, long aOffset, long m, Seq sb2, long bOffset, long n)
    {
        // A [i] = sa [aOffset + i], B [j] = sb2 [bOffset + j]
        if (m >= MaxAlign || n >= MaxAlign)
        {
            outText.Append("\n   *** Too long ***\n\n");
            errors.Append(line).Append(' ').Append("-".PadLeft(7)).Append('\n');
            return;
        }

        if (m < 0 || n < 0)
            throw new InvalidOperationException("mummer-annotate: a negative window length (the program's matrix allocation is undefined).");

        var d = new int[m + 1, n + 1];
        var op = new char[m + 1, n + 1];
        d[0, 0] = 0;
        op[0, 0] = 'a';
        for (long i = 1; i <= m; i++)
        {
            d[i, 0] = (int)i + 1;
            op[i, 0] = 'i';
        }

        for (long j = 1; j <= n; j++)
        {
            d[0, j] = (int)j + 1;
            op[0, j] = 'd';
        }

        for (long i = 1; i <= m; i++)
        {
            for (long j = 1; j <= n; j++)
            {
                d[i, j] = d[i, j - 1] + 1 + (op[i, j - 1] == 'd' ? 0 : 1);
                op[i, j] = 'd';
                int tmp = d[i - 1, j] + 1 + (op[i - 1, j] == 'i' ? 0 : 1);
                if (tmp < d[i, j])
                {
                    d[i, j] = tmp;
                    op[i, j] = 'i';
                }

                tmp = d[i - 1, j - 1] + (sa[aOffset + i] == sb2[bOffset + j] ? 0 : 1);
                if (tmp < d[i, j])
                {
                    d[i, j] = tmp;
                    op[i, j] = 'a';
                }
            }
        }

        var showA = new byte[m + n + 2];
        var showB = new byte[m + n + 2];
        long ct = 0;
        long ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            switch (op[ii, jj])
            {
                case 'a':
                    showA[ct] = sa[aOffset + ii--];
                    showB[ct++] = sb2[bOffset + jj--];
                    break;
                case 'i':
                    showA[ct] = sa[aOffset + ii--];
                    showB[ct++] = (byte)'.';
                    break;
                default:
                    showA[ct] = (byte)'.';
                    showB[ct++] = sb2[bOffset + jj--];
                    break;
            }
        }

        int errs = 0;
        for (long i = 0; i < ct; i++)
            if (showA[i] != showB[i])
                errs++;
        outText.Append("    Errors = ").Append(errs.ToString(CultureInfo.InvariantCulture)).Append('\n');
        errors.Append(line).Append(' ').Append(errs.ToString(CultureInfo.InvariantCulture).PadLeft(7)).Append('\n');
        do
        {
            outText.Append("T:  ");
            for (long i = ct - 1; i >= 0 && i >= ct - Width; i--)
                outText.Append((char)showA[i]);
            outText.Append("\nS:  ");
            for (long i = ct - 1; i >= 0 && i >= ct - Width; i--)
                outText.Append((char)showB[i]);
            outText.Append("\n    ");
            for (long i = ct - 1; i >= 0 && i >= ct - Width; i--)
                outText.Append(showA[i] == showB[i] ? ' ' : '^');
            outText.Append('\n');
            ct -= Width;
        }
        while (ct > 0);
    }
}
