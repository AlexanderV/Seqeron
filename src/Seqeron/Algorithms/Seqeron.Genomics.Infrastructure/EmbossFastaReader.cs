namespace Seqeron.Genomics.Infrastructure;

/// <summary>One record of a FASTA file as the EMBOSS 6.6 reader (<c>seqReadFasta</c>, <c>ajSeqParseFasta</c>) stores it.</summary>
/// <param name="Name">The sequence name (the first word of the header; blanks of a non-word id become <c>_</c>).</param>
/// <param name="Accession">The accession: the second word of the header when it is an accession number or a sequence version (<c>AB123456.1</c> gives <c>AB123456</c>), else empty.</param>
/// <param name="Description">The rest of the header.</param>
/// <param name="Text">The sequence lines joined (not cleaned: <see cref="EmbossSequenceReader"/> keeps the letters).</param>
public sealed record EmbossFastaRecord(string Name, string Accession, string Description, string Text);

/// <summary>
/// The plain FASTA reader of EMBOSS 6.6 (<c>ajseqread.c</c> <c>seqReadFasta</c> / <c>ajSeqParseFasta</c>, <c>ajseq.c</c> <c>ajSeqtestIsAccession</c> / <c>ajSeqtestIsSeqversion</c>, EMBOSS-6.6.0 via raw.githubusercontent.com):
/// the header is split into name, accession and description. NCBI-style identifiers with <c>|</c> are read by the separate NCBI parser of the program and are not part of this reader.
/// </summary>
public static class EmbossFastaReader
{
    // ajSeqtestIsAccession
    internal static bool IsAccession(string s)
    {
        if (s.Length < 6)
            return false;
        int p = 0;
        if (!char.IsAsciiLetter(s[p++]))
            return false;
        if (char.IsAsciiLetter(s[p]))
        {
            p++;
            if (p < s.Length && s[p] == '_')
                p++;
            for (; p < s.Length; p++)
                if (!char.IsAsciiDigit(s[p]))
                    return false;
            return true;
        }

        if (char.IsAsciiDigit(s[p]))
        {
            p++;
            for (int i = 0; i < 3; i++)
            {
                if (p < s.Length && char.IsAsciiLetterOrDigit(s[p]))
                    p++;
                else
                    return false;
            }

            if (p >= s.Length || !char.IsAsciiDigit(s[p]))
                return false;
            for (; p < s.Length; p++)
                if (!char.IsAsciiDigit(s[p]))
                    return false;
            return true;
        }

        return false;
    }

    // ajSeqtestIsSeqversion: the accession part of "ACCESSION.version", else null
    internal static string? SeqVersionAccession(string s)
    {
        if (s.Length < 8 || !char.IsAsciiLetter(s[0]))
            return null;
        var acc = new System.Text.StringBuilder();
        int p = 0;
        acc.Append(s[p++]);
        if (char.IsAsciiLetter(s[p]))
        {
            acc.Append(s[p++]);
            if (p < s.Length && char.IsAsciiLetter(s[p]))
                acc.Append(s[p++]);
            else if (p < s.Length && s[p] == '_')
                acc.Append(s[p++]);
        }
        else if (char.IsAsciiDigit(s[p]))
        {
            acc.Append(s[p++]);
            for (int i = 0; i < 3; i++)
            {
                if (p < s.Length && char.IsAsciiLetterOrDigit(s[p]))
                    acc.Append(s[p++]);
                else
                    return null;
            }

            if (p >= s.Length || !char.IsAsciiDigit(s[p]))
                return null;
        }
        else
            return null;
        bool dot = false;
        int v = 0;
        for (; p < s.Length; p++)
        {
            char c = s[p];
            if (c == '.')
            {
                if (dot)
                    return null;
                dot = true;
            }
            else if (char.IsAsciiDigit(c))
            {
                if (dot)
                    v++;
                else
                    acc.Append(c);
            }
            else
                return null;
        }

        return v > 0 ? acc.ToString() : null;
    }

    private static readonly char[] HeaderSeparators = [' ', '\t', '\n', '\r'];

    /// <summary>Splits a FASTA header line (starting with <c>&gt;</c>) as <c>ajSeqParseFasta</c> does.</summary>
    public static (string Name, string Accession, string Description) ParseHeader(string header)
    {
        ArgumentNullException.ThrowIfNull(header);
        var ncbi = EmbossNcbiHeader.Parse(header);
        if (ncbi is not null)
            return (SetName(ncbi.Id), ncbi.Accession, ncbi.Description);
        return ParsePlainHeader(header);
    }

    // seqSetName (ajseqread.c): a word is split at ':' (the last non-empty part is kept) and ",/\\" become '_'; a non-word has its blanks / ",;:/\\" turned into '_'
    internal static string SetName(string id)
    {
        if (id.Length == 0)
            return "";
        if (IsWord(id))
        {
            string last = "";
            foreach (var part in id.Split(':'))
                if (part.Length > 0)
                    last = part;
            return last.Replace(',', '_').Replace('/', '_').Replace('\\', '_');
        }

        var words = id.Split([' ', '\t', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
        return string.Join(' ', words).Replace(' ', '_').Replace(',', '_').Replace(';', '_').Replace(':', '_').Replace('/', '_').Replace('\\', '_');
    }

    // the id as ajSeqParseFasta stores it (before seqSetName)
    internal static (string Name, string Accession, string Description) ParsePlainRaw(string header)
    {
        var (_, acc, desc) = ParsePlainHeader(header);
        string line = header.TrimStart('>', ' ');
        int sp = line.IndexOfAny(HeaderSeparators);
        return (sp < 0 ? line : line[..sp], acc, desc);
    }

    private static (string Name, string Accession, string Description) ParsePlainHeader(string header)
    {
        string line = header.TrimStart('>', ' ');
        int sp = line.IndexOfAny(HeaderSeparators);
        string id = sp < 0 ? line : line[..sp];
        string rest = sp < 0 ? "" : line[(sp + 1)..];
        // the next token is delimited by blanks (ajStrTokenNextParse with the default delimiter set); the description of an accession is the remainder of the line
        rest = rest.TrimStart(' ', '\t');
        int sp2 = rest.IndexOfAny(HeaderSeparators);
        string token = sp2 < 0 ? rest : rest[..sp2];
        string afterToken = sp2 < 0 ? "" : rest[(sp2 + 1)..].TrimEnd('\r', '\n');
        string acc = "", desc;
        string token2 = token.Replace("(", "", StringComparison.Ordinal).Replace(")", "", StringComparison.Ordinal);
        string? fromVersion = token.Length > 0 ? SeqVersionAccession(token2) : null;
        if (fromVersion is not null)
        {
            acc = fromVersion;
            desc = afterToken;
        }
        else if (token.Length > 0 && IsAccession(token2))
        {
            acc = token2;
            desc = afterToken;
        }
        else
        {
            desc = Join(token, afterToken);
        }

        return (SetName(id), acc, desc);
    }

    private static string Join(string token, string rest)
    {
        if (token.Length == 0)
            return "";
        return rest.Length > 0 ? token + " " + rest : token;
    }

    private static bool IsWord(string s) => s.All(c => !char.IsWhiteSpace(c));

    /// <summary>Reads every record of a plain FASTA text.</summary>
    public static IReadOnlyList<EmbossFastaRecord> Read(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        var records = new List<EmbossFastaRecord>();
        int unnamed = 0;
        string? header = null;
        var body = new System.Text.StringBuilder();
        foreach (string raw in text.Split('\n'))
        {
            string line = raw.TrimEnd('\r');
            if (line.StartsWith('>'))
            {
                if (header is not null)
                    records.Add(Make(header, body.ToString(), ref unnamed));
                header = line;
                body.Clear();
            }
            else if (header is not null)
                body.Append(line);
        }

        if (header is not null)
            records.Add(Make(header, body.ToString(), ref unnamed));
        return records;
    }

    // a sequence without a name is called EMBOSS_001, EMBOSS_002, ... in the order of the unnamed sequences of the input (measured)
    private static EmbossFastaRecord Make(string header, string body, ref int unnamed)
    {
        var (name, acc, desc) = ParseHeader(header);
        if (name.Length == 0)
            name = "EMBOSS_" + (++unnamed).ToString("000", System.Globalization.CultureInfo.InvariantCulture);
        return new EmbossFastaRecord(name, acc, desc, body);
    }
}
