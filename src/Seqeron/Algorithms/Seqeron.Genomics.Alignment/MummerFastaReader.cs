using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary><c>Read_String (fp, T, Size, Name, FALSE)</c> of <c>tigrinc.cc</c>: a sequential reader of FASTA records (the first word of the header, residues lower-cased, a character that is not a letter or <c>*</c> becomes <c>x</c> with a warning).</summary>
internal sealed class MummerFastaReader(string text)
{
    internal static bool IsSpace(int c) => c == ' ' || (c >= '\t' && c <= '\r');

    internal static bool IsAlpha(int c) => (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');

    private int _pos;

    public bool Read(out string name, out byte[] seq, out long len, StringBuilder err)
    {
        name = string.Empty;
        seq = [0, 0];
        len = 0;
        int ch;
        while (_pos < text.Length && text[_pos] != '>')
            _pos++;
        if (_pos >= text.Length)
            return false;
        _pos++; // the '>'

        // fgets (Line, MAX_LINE = 1024): at most 1023 characters up to and including the newline
        int start = _pos;
        int nl = text.IndexOf('\n', _pos);
        if (nl < 0 || nl - start + 1 > 1023)
            throw new FormatException("Assertion `Len > 0 && Line [Len - 1] == '\\n'' failed: the header line of a FASTA record is missing its newline or is longer than 1023 characters.");
        string line = text.Substring(start, nl - start + 1);
        _pos = nl + 1;
        // strtok (Line, " \t\n")
        int i = 0;
        while (i < line.Length && (line[i] == ' ' || line[i] == '\t' || line[i] == '\n'))
            i++;
        int j = i;
        while (j < line.Length && line[j] != ' ' && line[j] != '\t' && line[j] != '\n')
            j++;
        name = line.Substring(i, j - i);

        var bytes = new List<byte>(1024) { 0 };
        while (_pos < text.Length)
        {
            ch = text[_pos];
            if (ch == '>')
                break;
            _pos++;
            if (IsSpace(ch))
                continue;
            ch = ch >= 'A' && ch <= 'Z' ? ch + 32 : ch;
            if (!IsAlpha(ch) && ch != '*')
            {
                err.Append("Unexpected character `").Append((char)ch).Append("' in string ").Append(name).Append('\n');
                ch = 'x';
            }

            bytes.Add((byte)ch);
        }

        bytes.Add(0);
        seq = [.. bytes];
        len = bytes.Count - 2;
        return true;
    }
}
