#pragma warning disable S3358, S3241, S3776 // the library functions keep their shape (nested conditions, results that the callers of other readers use)
using System.Globalization;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>A citation of a sequence entry (<c>AjPSeqRef</c>: the <c>RN</c> / <c>REFERENCE</c> blocks). A5-8c3h, B09 F268.</summary>
internal sealed class EmbossSeqRef
{
    public uint Number;
    public string Position = "", Groupname = "", Authors = "", Title = "", Comment = "", Xref = "", Location = "", Loctype = "";

    public EmbossSeqRef Clone() => (EmbossSeqRef)MemberwiseClone();

    // ajSeqrefAppend*: the text is appended after a blank (a line end after a full stop for the locations and cross references)
    public void AppendAuthors(string s) => Authors = (Authors.Length > 0 ? Authors + " " : "") + s;

    public void AppendComment(string s) => Comment = (Comment.Length > 0 ? Comment + " " : "") + s;

    public void AppendGroupname(string s) => Groupname = (Groupname.Length > 0 ? Groupname + " " : "") + s;

    public void AppendLocation(string s)
    {
        if (Location.Length > 0)
            Location += Location[^1] is '.' or ';' ? "\n" : " ";
        else
            SetLoctype(s);
        Location += s;
    }

    public void AppendPosition(string s) => Position = (Position.Length > 0 ? Position + " " : "") + s;

    public void AppendTitle(string s) => Title = (Title.Length > 0 ? Title + " " : "") + s;

    public void AppendXref(string s)
    {
        if (Xref.Length > 0)
            Xref += Xref[^1] == '.' ? "\n" : " ";
        Xref += s;
    }

    // ajSeqrefSetLoctype
    public void SetLoctype(string s)
    {
        Loctype = s.StartsWith("Submitted ", StringComparison.Ordinal) ? "submission"
            : s.StartsWith("(in) ", StringComparison.Ordinal) ? "book"
            : s.StartsWith("(er) ", StringComparison.Ordinal) ? "electronic"
            : s.StartsWith("Thesis ", StringComparison.Ordinal) ? "thesis"
            : s.StartsWith("Patent ", StringComparison.Ordinal) ? "patent"
            : s.StartsWith("(misc) ", StringComparison.Ordinal) ? "misc" : "journal";
    }

    // ajSeqrefStandard: the final full stops and semicolons that the entry formats add are removed, the title is unquoted
    public void Standard()
    {
        if (Xref.Length > 0 && Xref[^1] == '.')
            Xref = Xref.TrimEnd('.');
        if (Location.Length > 0 && Location[^1] == '.')
            Location = Location.TrimEnd('.');
        if (Title.Length > 0 && Title[^1] == '.')
            Title = Title.TrimEnd('.');
        if (Title.Length > 0 && Title[^1] == ';')
            Title = Title.TrimEnd(';');
        Title = QuoteStrip(Title);
        if (Authors.Length > 0 && Authors[^1] == ';')
            Authors = Authors.TrimEnd(';');
    }

    // ajStrQuoteStrip: a quote at the end and a quote at the start are cut (each on its own), doubled and escaped quotes inside become one
    internal static string QuoteStrip(string s)
    {
        if (s.Length > 0 && s[^1] == '"')
            s = s[..^1];
        if (s.Length > 0 && s[0] == '"')
            s = s[1..];
        s = s.Replace("\"\"", "\"", StringComparison.Ordinal);
        return s.Replace("\\\"", "\"", StringComparison.Ordinal);
    }

    // ajSeqrefFmtAuthorsGb: the last ", " becomes " and", a blank before an initial ("Smith J.") becomes a comma
    public string FmtAuthorsGb()
    {
        string d = Authors;
        int i = d.LastIndexOf(", ", StringComparison.Ordinal);
        if (i != -1)
            d = d[..i] + " and" + d[(i + 1)..];
        char[] cp = d.ToCharArray();
        int imax = cp.Length - 2;
        for (int k = 0; k < imax; k++)
        {
            if (cp[k] == ' ' && cp[k + 2] == '.' && char.IsAsciiLetter(cp[k + 1]))
                cp[k] = ',';
        }

        return new string(cp);
    }

    // ajSeqrefFmtLocationGb
    public string FmtLocationGb()
    {
        string d = Location.Replace("):", "), ", StringComparison.Ordinal);
        d = d.Replace("(", "\u0001", StringComparison.Ordinal);
        d = d.Replace("\u0001", " (", StringComparison.Ordinal);
        return d.Replace("  (", " (", StringComparison.Ordinal);
    }

    // ajSeqrefFmtTitleGb: "Direct Submission" for the submissions that have no title
    public string FmtTitleGb() => Title.Length == 0 ? (Loctype == "submission" ? "Direct Submission" : "") : Title;
}

/// <summary>A cross-reference (<c>AjPSeqXref</c>, the <c>DR</c> lines).</summary>
internal sealed class EmbossSeqXref
{
    /// <summary>The type: 1 = <c>DR</c> line (<c>XREF_DR</c>).</summary>
    public int Type;
    public string Db = "", Id = "", Secid = "", Terid = "", Quatid = "";
    public int Start, End;

    public EmbossSeqXref Clone() => (EmbossSeqXref)MemberwiseClone();
}

/// <summary>A gene of a Swiss-Prot entry (<c>AjPSeqGene</c>, a <c>GN</c> line group).</summary>
internal sealed class EmbossSeqGene
{
    public string Name = "", Synonyms = "", Oln = "", Orf = "";

    public EmbossSeqGene Clone() => (EmbossSeqGene)MemberwiseClone();

    private static string Append(string old, string s) => (old.Length > 0 ? old + " " : "") + s;

    public void AppendName(string s) => Name = Append(Name, s);

    public void AppendSynonyms(string s) => Synonyms = Append(Synonyms, s);

    public void AppendOln(string s) => Oln = Append(Oln, s);

    public void AppendOrf(string s) => Orf = Append(Orf, s);
}

/// <summary>The dates of an entry (<c>AjPSeqDate</c>): the dates are kept as the day, month and year the entry gave (the <c>tm</c> of <c>mktime</c> after <c>ajTimeNewDayFmt</c> with the format <c>dtline</c>: <c>25-MAR-2010</c>).</summary>
internal sealed class EmbossSeqDate
{
    public DateTime? CreDate, ModDate, SeqDate;
    public string CreRel = "", ModRel = "", SeqRel = "", CreVer = "", ModVer = "", SeqVer = "";

    /// <summary>The version strings exist (an object that <c>ajSeqdateClear</c> cleared keeps them, empty): a date set on it gets no default version; a new object has none and the setters make the version <c>1</c>.</summary>
    public bool VersionsExist;

    public EmbossSeqDate Clone() => (EmbossSeqDate)MemberwiseClone();

    // ajSeqdateExists
    public bool Exists => CreDate is not null || ModDate is not null || SeqDate is not null;

    private static readonly string[] Months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"];

    // seqDateSet: "dd-MON-yyyy"; false (and no date) when the month is not one of the twelve or the day is 0. Like mktime, a day beyond the month moves into the next.
    private static bool Parse(string datestr, out DateTime date)
    {
        date = default;
        var tok = new StrTok(datestr, "-");
        tok.Next(out string d);
        CStrToInt(d, out int day);
        tok.Next(out string m);
        int month = Array.IndexOf(Months, m);
        bool ok = true;
        if (month < 0)
        {
            ok = false;
            month = 1;
        }
        else
        {
            month++;
        }

        tok.Next(out string y);
        CStrToInt(y, out int year);
        if (year > 1900)
            year -= 1900;
        if (day == 0)
        {
            ok = false;
            day = 1;
        }

        if (!ok)
            return false;
        if (year > 1899)
            year -= 1900;
        try
        {
            date = new DateTime(1900, 1, 1, 0, 0, 0, DateTimeKind.Unspecified).AddYears(year).AddMonths(month - 1).AddDays(day - 1);
        }
        catch (ArgumentOutOfRangeException)
        {
            return false;
        }

        return true;
    }

    // ajStrToInt: the leading integer, 0 when there is none
    private static bool CStrToInt(string s, out int v)
    {
        v = 0;
        s = s.Trim();
        int i = 0;
        if (i < s.Length && (s[i] is '-' or '+'))
            i++;
        int st = i;
        while (i < s.Length && char.IsAsciiDigit(s[i]))
            i++;
        if (i == st)
            return false;
        return int.TryParse(s[..i], NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out v);
    }

    public bool SetCreate(string datestr)
    {
        if (!Parse(datestr, out var d))
            return false;
        CreDate = d;
        if (!VersionsExist && CreVer.Length == 0)
            CreVer = "1";
        return true;
    }

    public bool SetModify(string datestr)
    {
        if (!Parse(datestr, out var d))
            return false;
        ModDate = d;
        if (!VersionsExist && ModVer.Length == 0)
            ModVer = "1";
        return true;
    }

    public bool SetModseq(string datestr)
    {
        if (!Parse(datestr, out var d))
            return false;
        SeqDate = d;
        if (!VersionsExist && SeqVer.Length == 0)
            SeqVer = "1";
        return true;
    }

    /// <summary><c>%D</c> of a date with the format <c>dtline</c>: <c>25-MAR-2010</c>.</summary>
    public static string Format(DateTime d) => d.Day.ToString("00", CultureInfo.InvariantCulture) + "-" + Months[d.Month - 1] + "-" + d.Year.ToString("0000", CultureInfo.InvariantCulture);
}

/// <summary>The molecule, class and division tables of <c>ajseq.c</c> (<c>seqMolecule</c>, <c>seqClass</c>, <c>seqDivision</c>): the names of the entry formats and the names the library uses inside.</summary>
internal static class EmbossSeqTables
{
    // Name, Embl, Ddbj, Genbank
    private static readonly (string Name, string Embl, string Ddbj, string Gb)[] Mol =
    [
        ("unassigned DNA", "unassigned DNA", "DNA", "DNA"), ("genomic DNA", "genomic DNA", "DNA", "DNA"), ("genomic RNA", "genomic RNA", "DNA", "RNA"), ("mRNA", "mRNA", "mRNA", "mRNA"), ("tRNA", "tRNA", "RNA", "RNA"),
        ("rRNA", "rRNA", "RNA", "RNA"), ("snoRNA", "snoRNA", "RNA", "RNA"), ("snRNA", "snRNA", "RNA", "RNA"), ("snRNA", "snRNA", "RNA", "RNA"), ("scRNA", "scRNA", "RNA", "RNA"), ("pre-RNA", "pre-RNA", "RNA", "RNA"),
        ("viral cRNA", "viral cRNA", "RNA", "RNA"), ("other RNA", "other RNA", "RNA", "RNA"), ("other DNA", "other DNA", "DNA", "DNA"), ("unassigned RNA", "unassigned RNA", "RNA", "RNA"),
    ];

    private static readonly (string Name, string Embl, string Ddbj, string Gb)[] Cls =
    [
        ("STD", "STD", "", ""), ("ANN", "ANN", "", ""), ("CON", "CON", "CON", "CON"), ("EST", "EST", "EST", "EST"), ("GSS", "GSS", "GSS", "GSS"), ("HTC", "HTC", "HTC", ""), ("HTG", "HTG", "HTG", "HTG"), ("MGA", "MGA", "", ""),
        ("PAT", "PAT", "PAT", "PAT"), ("STS", "STS", "STS", "STS"), ("TPA", "TPA", "TPA", ""), ("TSA", "TSA", "TSA", ""), ("WGS", "WGS", "", ""), ("", "", "", ""),
    ];

    private static readonly (string Name, string Embl, string Ddbj, string Gb)[] Div =
    [
        ("default", "UNC", "UNA", "UNA"), ("CON", "CON", "CON", "CON"), ("ENV", "ENV", "", ""), ("FUN", "FUN", "", "PLN"), ("HUM", "HUM", "HUM", "PRI"), ("INV", "INV", "INV", "INV"), ("MAM", "MAM", "MAM", "MAM"), ("MAMPRI", "MAM", "", "PRI"),
        ("MUS", "MUS", "ROD", "ROD"), ("PHG", "PHG", "PHG", "PHG"), ("PLN", "PLN", "PLN", "PLN"), ("PRO", "PRO", "BCT", "BCT"), ("ROD", "ROD", "ROD", "ROD"), ("SYN", "SYN", "SYN", "SYN"), ("TGN", "TGN", "", ""), ("UNC", "UNC", "", ""),
        ("VRL", "VRL", "VRL", "VRL"), ("VRT", "VRT", "VRT", "VRT"), ("UNA", "UNC", "UNA", "UNA"), ("UNCRNA", "UNC", "", "RNA"), ("", "", "", ""),
    ];

    // the tables map the name of a database format to the internal name; with several entries of one name the last one wins (ajTablePut)
    private static string? Lookup((string Name, string Embl, string Ddbj, string Gb)[] t, Func<(string Name, string Embl, string Ddbj, string Gb), string> sel, string key)
    {
        string? found = null;
        foreach (var e in t)
        {
            string k = sel(e);
            if (k.Length > 0 && string.Equals(k, key, StringComparison.OrdinalIgnoreCase))
                found = e.Name;
        }

        return found;
    }

    /// <summary><c>ajSeqmolSetEmbl</c>: the internal name of an EMBL molecule type, null when there is none.</summary>
    public static string? MolFromEmbl(string s) => Lookup(Mol, e => e.Embl, s);

    /// <summary><c>ajSeqmolSetGb</c>.</summary>
    public static string? MolFromGb(string s) => Lookup(Mol, e => e.Gb, s);

    /// <summary><c>ajSeqdivSetGb</c>: GenBank mixes division and class, so a name that is no division is looked up in the classes; the class tables exist only after the first call of a class function of the process (<paramref name="classTables"/>).</summary>
    public static string? DivFromGb(string s, bool classTables) => Lookup(Div, e => e.Gb, s) ?? (classTables ? Lookup(Cls, e => e.Gb, s) : null);

    /// <summary><c>ajSeqclsSetGb</c>.</summary>
    public static string? ClsFromGb(string s) => Lookup(Cls, e => e.Gb, s);

    private static (string Name, string Embl, string Ddbj, string Gb)? Find((string Name, string Embl, string Ddbj, string Gb)[] t, string name)
    {
        (string Name, string Embl, string Ddbj, string Gb)? found = null;
        foreach (var e in t)
        {
            if (string.Equals(e.Name, name, StringComparison.OrdinalIgnoreCase))
                found = e;
        }

        return found;
    }

    /// <summary><c>ajSeqmolGetEmbl</c>: the EMBL name of a molecule type (the default <c>unassigned DNA</c> for an unknown or empty one).</summary>
    public static string MolEmbl(string mol) => (mol.Length > 0 ? Find(Mol, mol) : null)?.Embl ?? Find(Mol, Mol[0].Name)!.Value.Embl;

    /// <summary><c>ajSeqmolGetGb</c>.</summary>
    public static string MolGb(string mol) => (mol.Length > 0 ? Find(Mol, mol) : null)?.Gb ?? Find(Mol, Mol[0].Name)!.Value.Gb;

    /// <summary><c>ajSeqdivGetGb</c>.</summary>
    public static string DivGb(string div) => (div.Length > 0 ? Find(Div, div) : null)?.Gb ?? Find(Div, Div[0].Name)!.Value.Gb;

    /// <summary><c>ajSeqclsGetGb</c>.</summary>
    public static string ClsGb(string cls) => (cls.Length > 0 ? Find(Cls, cls) : null)?.Gb ?? Find(Cls, Cls[0].Name)!.Value.Gb;
}
