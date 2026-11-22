#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481 // the readers keep the control flow of ajseqread.c

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>A sequence as the sequence readers of EMBOSS 6.6 (<c>ajseqread.c</c>) leave it in the AJAX sequence object (<c>AjPSeq</c>): the fields the readers set and the programs of this repository use.</summary>
internal sealed class EmbossSeq
{
    public string Name = "", Acc = "", Sv = "", Gi = "", Desc = "", Setdb = "", Tax = "", Taxid = "", Organelle = "", Taxcommon = "", Db = "", Entryname = "", Filename = "", Formatstr = "", Usa = "";

    /// <summary>The type <c>N</c> / <c>P</c> or empty while it is not known (<c>thys-&gt;Type</c>).</summary>
    public string Type = "";

    public List<string> Acclist = [];
    public List<string> Taxlist = [];

    /// <summary>The <c>Multi</c> list of the full description (<c>thys-&gt;Fulldesc</c>): null when the reader made none (only the swiss and iguspto readers do); the lines of the comments of an iguspto entry.</summary>
    public List<string>? FullDescMulti;
    public StringBuilder Seq = new();
    public EmbossSeqFeatureTable? Fttable;
    public int Format;
    public bool Circular;

    /// <summary>The text of the entry the reader read (<c>thys-&gt;TextPtr</c>; kept when the input stores it).</summary>
    public StringBuilder? TextPtr;

    /// <summary><c>ajSeqSetNuc</c>: the type, the type of the feature table when it has none, and N as the ambiguity code (every x and X of the sequence becomes n and N).</summary>
    public void SetNuc()
    {
        if (Type == "N")
            return;
        Type = "N";
        Fttable?.SetType("N");
        for (int i = 0; i < Seq.Length; i++)
        {
            if (Seq[i] == 'x')
                Seq[i] = 'n';
            else if (Seq[i] == 'X')
                Seq[i] = 'N';
        }
    }

    /// <summary><c>ajSeqSetProt</c>.</summary>
    public void SetProt()
    {
        if (Type == "P")
            return;
        Type = "P";
        Fttable?.SetType("P");
    }

    public string SeqText => Seq.ToString();
}

/// <summary>The sequence input of one program run: the options of the ACD sequence input that the readers consult (<c>AjPSeqin</c>).</summary>
internal sealed class EmbossSeqinOptions
{
    /// <summary>The input format (<c>-sformat</c>) as typed; null = the formats the program tries are tried in turn.</summary>
    public string? Format;

    /// <summary>The <c>format::</c> prefix of the USA (it wins over <see cref="Format"/>); an empty prefix is the format <c>unknown</c>; null = none.</summary>
    public string? PrefixFormat;

    /// <summary>The name of the file (<c>seqin-&gt;Input-&gt;Filename</c>).</summary>
    public string Filename = "in.dat";

    /// <summary>The features of the sequence are wanted (<c>seqin-&gt;Features</c>: the ACD <c>features</c> attribute, e.g. seqret <c>-feature</c>).</summary>
    public bool Features;

    /// <summary>A feature file was named (<c>seqin-&gt;Ufo</c> is not empty): the features of the sequence file are not read.</summary>
    public bool UfoGiven;

    /// <summary><c>-sformat</c> found no such format: the program stops ("Unknown input format").</summary>
    public bool Lower, Upper, IsNuc, IsProt, Circular;

    /// <summary>The ACD type of the sequence input: <c>gapany</c> (seqret; gaps kept) or <c>any</c> (most programs: every character that is no letter and no <c>*</c> is removed).</summary>
    public string SequenceType = "gapany";
    public int Begin, End;
    public bool Rev;
}

/// <summary>
/// The sequence readers of EMBOSS 6.6.0 <c>ajseqread.c</c> (<c>seqRead</c>, <c>seqReadFmt</c>, <c>seqinFormatDef</c> and the <c>seqRead*</c> functions) over the lines of one file: <c>ajSeqRead</c> returns the sequences of the file one after another.
/// The buffer is the one of the feature library (<see cref="EmbossFeatureLines"/>: <c>ajFilebuffClear</c> / <c>Reset</c> / <c>SetUnbuffered</c>). A5-8c3g, B09 F248.
/// </summary>
internal sealed partial class EmbossSeqin
{
    private const int FmtOk = 0, FmtBadtype = 1, FmtFail = 2, FmtEof = 3, FmtEmpty = 5;

    /// <summary>name, EDAM term, alias, tried when no format is named, nucleotide, protein, features, gaps, multiset, binary (the columns of <c>seqinFormatDef</c>).</summary>
    internal static readonly (string Name, string Obo, bool Alias, bool Try, bool Nuc, bool Prot, bool Feat, bool Gap, bool Multi, bool Binary)[] Table =
    [
        ("unknown", "0000", false, false, true, true, false, true, false, false), ("gcg", "1935", false, true, true, true, false, true, false, false),
        ("gcg8", "1935", true, false, true, true, false, true, false, false), ("embl", "1927", false, true, true, false, true, true, false, false),
        ("em", "1927", true, false, true, false, true, true, false, false), ("swiss", "1963", false, true, false, true, true, true, false, false),
        ("sw", "1963", true, false, false, true, true, true, false, false), ("swissprot", "1963", true, false, false, true, true, true, false, false),
        ("uniprot", "2188", true, false, false, true, true, true, false, false), ("nbrf", "1948", false, true, true, true, true, true, false, false),
        ("pir", "1948", true, false, true, true, true, true, false, false), ("pdb", "1950", false, true, false, true, false, false, false, false),
        ("pdbseq", "1953", false, false, false, true, false, false, false, false), ("pdbnuc", "1951", false, false, true, false, false, false, false, false),
        ("pdbnucseq", "1952", false, false, true, false, false, false, false, false), ("fasta", "1929", false, true, true, true, false, true, false, false),
        ("ncbi", "1929", true, false, true, true, false, true, false, false), ("gifasta", "1940", false, false, true, true, false, true, false, false),
        ("pearson", "1954", false, false, true, true, false, true, false, false), ("fastq", "1930", false, true, true, false, false, false, false, false),
        ("fastq-sanger", "1932", false, false, true, false, false, false, false, false), ("fastq-illumina", "1931", false, false, true, false, false, false, false, false),
        ("fastq-solexa", "1933", false, false, true, false, false, false, false, false), ("sam", "2573", false, true, true, false, false, true, false, false),
        ("genbank", "1936", false, true, true, false, true, true, false, false), ("gb", "1936", true, false, true, false, true, true, false, false),
        ("ddbj", "1936", true, false, true, false, true, true, false, false), ("refseq", "1936", true, false, true, false, true, true, false, false),
        ("refseqp", "1958", false, false, false, true, true, true, false, false), ("genpept", "1937", false, false, false, true, false, true, false, false),
        ("codata", "1925", false, true, true, true, true, true, false, false), ("strider", "1962", false, true, true, false, false, true, false, false),
        ("clustal", "1924", false, true, true, true, false, true, false, false), ("aln", "1924", true, false, true, true, false, true, false, false),
        ("phylip", "1955", false, true, true, true, false, true, true, false), ("phylipnon", "1956", false, false, true, true, false, true, true, false),
        ("ace", "3001", false, true, true, false, false, true, false, false), ("consed", "3001", true, true, true, false, false, true, false, false),
        ("acedb", "1923", false, true, true, true, false, true, false, false), ("dbid", "1926", false, false, true, true, false, true, false, false),
        ("msf", "1947", false, true, true, true, false, true, false, false), ("hennig86", "1941", false, true, true, true, false, true, false, false),
        ("jackknifer", "1944", false, true, true, true, false, true, false, false), ("nexus", "1949", false, true, true, true, false, true, false, false),
        ("paup", "1949", true, false, true, true, false, true, false, false), ("treecon", "1965", false, true, true, true, false, true, false, false),
        ("mega", "1946 1971", false, true, true, true, false, true, false, false), ("igstrict", "1943", false, true, true, true, false, true, false, false),
        ("iguspto", "1942", false, true, true, true, false, true, false, false), ("ig", "1942", false, false, true, true, false, true, false, false),
        ("staden", "1960", false, false, true, true, false, true, false, false), ("textonly", "1964", false, false, true, true, false, true, false, false),
        ("plain", "1964", true, false, true, true, false, true, false, false), ("asis", "1964", true, false, true, true, false, true, false, false),
        ("gff2", "1938", false, true, true, true, true, true, false, false), ("gff3", "1939", false, true, true, true, true, true, false, false),
        ("gff", "1939", true, false, true, true, true, true, false, false), ("stockholm", "1961", false, true, false, true, false, true, false, false),
        ("pfam", "1961", true, true, false, true, false, true, false, false), ("selex", "1959", false, false, true, true, false, true, false, false),
        ("fitch", "1934", false, true, true, true, false, true, false, false), ("biomart", "0000", false, false, true, true, false, true, false, false),
        ("mase", "1945", false, false, true, true, false, true, false, false), ("experiment", "1928", false, true, true, false, false, true, false, false),
        ("gde", "0000", false, true, true, true, false, true, false, false), ("raw", "1957", false, true, true, true, false, false, false, true),
        ("nibble", "0000", false, true, true, false, false, false, false, true), ("nib", "0000", true, false, true, false, false, false, false, true),
        ("abi", "1628", false, true, true, false, false, false, false, true), ("bam", "2572", false, true, true, false, false, true, false, true),
        ("ensembl", "0000", false, false, true, true, true, true, false, false), ("das", "1967", false, false, true, true, true, true, false, false),
        ("scf", "2057", false, true, true, false, false, true, false, true)
    ];

    private readonly EmbossFeatureLines _buff;
    private readonly EmbossFeatureReadContext _context;
    private readonly EmbossSeqinOptions _opt;
    private readonly Func<string, string?> _files;
    private int _format;
    private readonly int _specified;
    private int _count;
    private int _records;
    private int _totRecords;
    private string _line = "";

    // seqToken of the library: a static string the readers parse into and that survives from one reader to the next (the reader of strider takes a token the reader of sam left)
    private string _seqToken = "";

    /// <param name="text">The text of the file.</param>
    /// <param name="options">The options of the sequence input.</param>
    /// <param name="context">The state of the feature library of the run (the messages are collected there).</param>
    /// <param name="files">The other files by name (feature files); null when there are none.</param>
    public EmbossSeqin(string text, EmbossSeqinOptions options, EmbossFeatureReadContext context, Func<string, string?>? files = null)
    {
        _buff = new EmbossFeatureLines(SplitLines(text));
        _opt = options;
        _context = context;
        _files = files ?? (_ => null);
        // ajQuerystrParseFormat (the prefix of the USA) then seqinFormatSet (-sformat): a name that is no format is reported and the formats are tried in turn
        bool found = false;
        if (options.PrefixFormat is { } prefix)
        {
            int f = FindFormat(prefix.Length == 0 ? "unknown" : prefix);
            if (f >= 0)
            {
                _format = f;
                found = true;
            }
            else
                context.Error($"Unknown input format '{prefix}'\n");
        }

        if (!found && options.Format is { Length: > 0 } sf)
        {
            int f = FindFormat(sf);
            if (f >= 0)
                _format = f;
            else
                context.Error($"Unknown input format '{sf}'");
        }

        _specified = _format;
    }

    // the lines of a file as ajReadline returns them: every line ends with a line feed, a last line of a file that does not end with one gets it (measured on 6.6.0)
    private static List<string> SplitLines(string text)
    {
        var parts = new List<string>(text.Split('\n'));
        if (parts.Count > 0 && parts[^1].Length == 0)
            parts.RemoveAt(parts.Count - 1);
        for (int i = 0; i < parts.Count; i++)
            parts[i] += "\n";
        return parts;
    }

    /// <summary>The index of a format by its name (any case) in <see cref="Table"/>; -1 when unknown (<c>seqSetInFormat</c> / <c>ajSeqinSetFormat</c>).</summary>
    internal static int FindFormat(string format)
    {
        string lower = format.ToLowerInvariant();
        for (int i = 0; i < Table.Length; i++)
        {
            if (Table[i].Name == lower || Table[i].Obo == format)
                return i;
        }

        return -1;
    }

    /// <summary>The format the input found or was given (<c>textin-&gt;Format</c>; 0 = none yet).</summary>
    public int Format => _format;

    // ---------------------------------------------------------------------------------------------------------------- the buffer (ajTextin*)

    private bool ReadLine()
    {
        string? l = _buff.Next();
        if (l is null)
        {
            _line = "";                                             // a read that hits the end of the file clears the line
            return false;
        }

        _line = l;
        _records++;
        _totRecords++;
        return true;
    }

    private void StoreReset()
    {
        _buff.Reset();
        _records = 0;
    }

    // ajTextinStoreClear: lines < 0 clears all, lines > 0 gives the last lines back
    private void StoreClear(int lines)
    {
        _buff.Clear(lines);
        if (lines < 0)
        {
            _totRecords -= _records;
            _records = 0;
        }
        else if (lines > 0)
        {
            _records = _records >= lines ? _records - lines : 0;
            _totRecords = _totRecords >= lines ? _totRecords - lines : 0;
        }
    }

    private static bool IsWhite(string s) => s.All(EmbossSeqFeatureTable.IsCSpace);

    private static bool IsWord(string s) => s.Length > 0 && !s.Any(EmbossSeqFeatureTable.IsCSpace);

    // ---------------------------------------------------------------------------------------------------------------- seqRead

    /// <summary><c>ajSeqRead</c>: the next sequence of the file, or null when there is none (the input found no more sequence of the format).</summary>
    public EmbossSeq? Next()
    {
        var t = new EmbossSeq();
        _count++;
        bool ok = !_buff.Unbuffered;
        if (!Table[_format].Binary)
        {
            while (ok)
            {
                string? l = _buff.Next();
                ok = l is not null;
                if (ok)
                    _line = l!;
                if (ok && !IsWhite(l!))
                {
                    _buff.Clear(1);
                    break;
                }
            }
        }

        int istat = 0;
        if (_format == 0)
        {
            for (int i = 1; i < Table.Length; i++)
            {
                if (!Table[i].Try)
                    continue;
                istat = ReadFmt(t, i);
                switch (istat)
                {
                    case FmtOk:
                        Define(t);
                        return t;
                    case FmtBadtype:
                    case FmtEof:
                        return null;
                    case FmtEmpty:
                        _context.Warn($"Sequence '{UsaOf(t)}' has zero length, ignored");
                        break;
                }

                if (_format != 0)
                    break;                                          // something was read: the format is known
                t = new EmbossSeq();
            }

            if (_format == 0)
                return null;                                        // all the formats failed
        }
        else
        {
            _buff.SetUnbuffered();
            istat = ReadFmt(t, _format);
            switch (istat)
            {
                case FmtOk:
                    Define(t);
                    return t;
                case FmtBadtype:
                case FmtFail:
                    return null;
                case FmtEof:
                    if (_records > 0)
                        _context.Error($"Error reading file '{_opt.Filename}' with format '{Table[_format].Name}': end-of-file before end of data (read {_records} records, total {_totRecords})");
                    break;
                case FmtEmpty:
                    _context.Warn($"Sequence '{UsaOf(t)}' has zero length, ignored");
                    break;
            }
        }

        // the first try failed or found an empty sequence: the rest of the file is read again while there is any (the search flag of the input)
        while (!_buff.IsEmpty)
        {
            t = new EmbossSeq();
            switch (ReadFmt(t, _format))
            {
                case FmtOk:
                    Define(t);
                    return t;
                case FmtBadtype:
                case FmtFail:
                case FmtEof:
                    return null;
                case FmtEmpty:
                    if (istat != FmtEmpty)
                        _context.Warn($"Sequence '{UsaOf(t)}' has zero length, ignored");
                    break;
            }
        }

        return null;
    }

    // seqMakeUsa
    private string UsaOf(EmbossSeq t)
    {
        if (t.Db.Length > 0)
            return $"{t.Db}-id:{t.Entryname}";
        string usa = t.Formatstr == "asis" ? "asis::" + t.SeqText : t.Entryname.Length > 0 ? $"{t.Formatstr}::{t.Filename}:{t.Entryname}" : $"{t.Formatstr}::{t.Filename}";
        return usa;
    }

    private int ReadFmt(EmbossSeq t, int format)
    {
        _records = 0;
        if (Dispatch(t, format))
        {
            _format = format;
            t.Formatstr = Table[format].Name;
            t.Filename = _opt.Filename;
            if (t.Entryname.Length == 0)
                t.Entryname = t.Name;
            if (_opt.Features && t.Fttable is null)
            {
                var table = EmbossFeatureUfo.Read("", null, null, t.Entryname, _files, _context);
                if (table is not null)
                {
                    table.Length = t.Seq.Length;
                    t.Fttable = table;
                }
            }

            if (t.Seq.Length == 0)
                return FmtEmpty;
            if (TypeCheckIn(t))
            {
                if (t.Seq.Length == 0)
                    return FmtEmpty;
                if (_opt.Upper)
                    ToCase(t, true);
                if (_opt.Lower)
                    ToCase(t, false);
                return FmtOk;
            }

            return FmtBadtype;
        }

        if (_buff.Unbuffered && _buff.End)
            return FmtEof;
        _buff.Reset();
        _records = 0;
        return FmtFail;
    }

    // ajSeqFmtUpper / ajSeqFmtLower
    private static void ToCase(EmbossSeq t, bool upper)
    {
        for (int i = 0; i < t.Seq.Length; i++)
        {
            char c = t.Seq[i];
            if (upper && char.IsAsciiLetterLower(c))
                t.Seq[i] = (char)(c - 32);
            else if (!upper && char.IsAsciiLetterUpper(c))
                t.Seq[i] = (char)(c + 32);
        }
    }

    // seqDefine: the type from the format when the reader set none (ajSeqSetNuc / ajSeqSetProt), else from the letters (ajSeqType), the circular flag of the table
    private void Define(EmbossSeq t)
    {
        if (t.Type.Length == 0 && _specified != 0)
        {
            var f = Table[_specified];                     // thys->Format is set from the format NAMED for the input (seqinFormatSet), not from the one a search found
            if (f.Nuc && !f.Prot)
                t.SetNuc();
            if (!f.Nuc && f.Prot)
                t.SetProt();
        }

        if (t.Type.Length == 0)
        {
            // the library looks at the residues as a C string: the text from a NUL character on is not seen
            string text = t.SeqText;
            int nul = text.IndexOf('\0');
            if (nul >= 0)
                text = text[..nul];
            if (IsNucleotideText(text))
                t.SetNuc();
            else if (!AnyProteinBad(text))
                t.SetProt();
        }

        if (_opt.Circular)
            t.Circular = true;
        if (t.Fttable is { } table)
        {
            if (t.Circular)
                table.Circular = true;
            else if (table.Circular)
                t.Circular = true;
        }
    }

    // ajSeqTypeAnyprotS: a character that is not a protein code or a gap
    private static bool AnyProteinBad(string s) => s.Any(c => !(char.IsAsciiLetter(c) || c is '*' or '-' or '.' or '~' or '?'));

    // ajSeqIsNuc of a sequence under reading
    private static bool IsNucleotideText(string s) => EmbossSequenceReader.IsNucleotide(s);

    // ajSeqTypeCheckIn of the type gapany (the type of the sequence input of seqret): -snucleotide / -sprotein, the gap characters '.' and '~' become '-', a sequence that holds only letters, '*', '?' and gaps keeps them
    // ('?' becomes X), any other character (# + ...) becomes N (a nucleotide) or X
    private bool TypeCheckIn(EmbossSeq t)
    {
        if (_opt.IsNuc)
            t.SetNuc();
        if (_opt.IsProt)
            t.SetProt();
        // the library checks the residues as a C string: the text behind a NUL character (the mega reader can leave one) is not looked at
        string whole = t.SeqText;
        int nul = whole.IndexOf('\0');
        string tail = nul >= 0 ? whole[nul..] : "";
        string text = nul >= 0 ? whole[..nul] : whole;
        if (_opt.SequenceType == "any")
        {
            // ungapped type: ajStrRemoveGap keeps the letters and '*'
            var kept = new StringBuilder(t.Seq.Length);
            foreach (char c in text)
            {
                if (char.IsAsciiLetter(c) || c == '*')
                    kept.Append(c);
            }

            t.Seq.Clear().Append(kept).Append(tail);
            return true;
        }

        var sb = new StringBuilder(t.Seq.Length);
        bool allGood = true;
        foreach (char ch in text)
        {
            char c = ch is '~' or '.' ? '-' : ch;
            if (!(char.IsAsciiLetter(c) || c is '*' or '-' or '?'))
                allGood = false;
            sb.Append(c);
        }

        string s = sb.ToString();
        if (!allGood)
        {
            // seqTypeFix: the characters that are not good ones become N (a nucleotide) or X
            bool nuc = t.Type == "N" || (t.Type.Length == 0 && IsNucleotideText(s));
            var fixedSeq = new StringBuilder(s.Length);
            foreach (char c in s)
                fixedSeq.Append(char.IsAsciiLetter(c) || c is '*' or '-' or '?' ? c : nuc ? 'N' : 'X');
            s = fixedSeq.ToString();
        }

        t.Seq.Clear().Append(s.Replace('?', 'X')).Append(tail);
        return true;
    }

    // ---------------------------------------------------------------------------------------------------------------- the readers

    private bool Dispatch(EmbossSeq t, int format)
    {
        switch (Table[format].Name)
        {
            case "unknown":
            case "textonly":
            case "plain":
            case "asis":
                return ReadText(t);
            case "gcg":
            case "gcg8":
                return ReadGcg(t);
            case "fasta":
            case "ncbi":
                return ReadNcbi(t);
            case "gifasta":
                return ReadGifasta(t);
            case "pearson":
                return ReadFasta(t);
            case "dbid":
                return ReadDbId(t);
            case "staden":
                return ReadStaden(t);
            case "strider":
                return ReadStrider(t);
            case "codata":
                return ReadCodata(t);
            case "nbrf":
            case "pir":
                return ReadNbrf(t);
            case "embl":
            case "em":
                return ReadEmbl(t);
            case "swiss":
            case "sw":
            case "swissprot":
            case "uniprot":
                return ReadSwiss(t);
            case "genbank":
            case "gb":
            case "ddbj":
            case "refseq":
                return ReadGenbank(t);
            case "genpept":
                return ReadProteinEntry(t, false);
            case "refseqp":
                return ReadProteinEntry(t, true);
            case "pdb":
                return ReadPdbAtoms(t, false);
            case "pdbnuc":
                return ReadPdbAtoms(t, true);
            case "pdbseq":
                return ReadPdbSeqres(t, false);
            case "pdbnucseq":
                return ReadPdbSeqres(t, true);
            case "ig":
                return ReadIg(t);
            case "igstrict":
                return ReadIgstrict(t);
            case "iguspto":
                return ReadIguspto(t);
            case "raw":
                return ReadRaw(t);
            case "gde":
                return ReadGde(t);
            case "clustal":
            case "aln":
                return ReadClustal(t);
            case "msf":
                return ReadMsf(t);
            case "fitch":
                return ReadFitch(t);
            case "mase":
                return ReadMase(t);
            case "biomart":
                return ReadBiomart(t);
            case "treecon":
                return ReadTreecon(t);
            case "acedb":
                return ReadAcedb(t);
            case "phylip":
                return ReadPhylip(t);
            case "phylipnon":
                return ReadPhylipnon(t);
            case "jackknifer":
                return ReadJackknifer(t);
            case "hennig86":
                return ReadHennig86(t);
            case "nexus":
            case "paup":
                return ReadNexus(t);
            case "mega":
                return ReadMega(t);
            case "experiment":
                return ReadExperiment(t);
            case "ace":
            case "consed":
                return ReadAce(t);
            case "selex":
                return ReadSelex(t);
            case "stockholm":
            case "pfam":
                return ReadStockholm(t);
            case "fastq":
                return ReadFastq(t);
            case "fastq-sanger":
                return ReadFastqQuality(t, "SANGER", 33);
            case "fastq-illumina":
                return ReadFastqQuality(t, "ILLUMINA", 64);
            case "fastq-solexa":
                return ReadFastqQuality(t, "SOLEXA", 59);
            case "sam":
                return ReadSam(t);
            case "gff2":
                return ReadGff2(t);
            case "gff3":
            case "gff":
                return ReadGff3(t);
            default:
                return ReadUnported(t, Table[format].Name);
        }
    }

    // seqAppend: the letters and the characters * . ~ ? # + - of the line
    private static void Append(EmbossSeq t, string line)
    {
        foreach (char c in line)
        {
            if (char.IsAsciiLetter(c) || c is '*' or '.' or '~' or '?' or '#' or '+' or '-')
                t.Seq.Append(c);
        }
    }

    // seqAppendWarn: with a format given (textin->Format; EMBOSS_SEQWARN is not set) the letters and the set are kept and the other characters that are no blanks are reported; with none they are dropped silently
    private void AppendWarn(EmbossSeq t, string line)
    {
        if (_format == 0)
        {
            Append(t, line);
            return;
        }

        var rest = new StringBuilder();
        foreach (char c in line)
        {
            if (char.IsAsciiLetter(c) || c is '*' or '.' or '~' or '?' or '#' or '+' or '-')
                t.Seq.Append(c);
            else if (!EmbossSeqFeatureTable.IsCSpace(c))
                rest.Append(c);
        }

        if (rest.Length > 0)
            _context.Warn($"Sequence '{t.Name}' has bad character(s) '{rest}'");
    }

    // seqAppendCommented: <...> (Staden) and <...< / >...> (GCG) comments are cut out of the sequence
    private static void AppendCommented(EmbossSeq t, ref bool incomment, string line)
    {
        var sb = new StringBuilder();
        foreach (char c in line)
        {
            if (char.IsAsciiLetter(c) || c is '*' or '.' or '~' or '?' or '#' or '+' or '-' or '<' or '>')
                sb.Append(c);
        }

        string tmp = sb.ToString();
        while (tmp.Length > 0)
        {
            int i = tmp.IndexOfAny(['<', '>']);
            if (incomment)
            {
                if (i >= 0)
                {
                    tmp = tmp[(i + 1)..];
                    incomment = false;
                }
                else
                    tmp = "";
            }
            else
            {
                if (i >= 0)
                {
                    if (i > 0)
                        t.Seq.Append(tmp, 0, i);
                    tmp = tmp[(i + 1)..];
                    incomment = true;
                }
                else
                {
                    t.Seq.Append(tmp);
                    tmp = "";
                }
            }
        }
    }

    // seqSetName
    private void SetName(EmbossSeq t, string str)
    {
        if (str.Length == 0)
        {
            // ajSeqSetNameMulti(thys, NULL): a sequence with no name is called EMBOSS_001, EMBOSS_002, ... in the order of the unnamed sequences of the run
            if (t.Name.Length == 0)
                t.Name = "EMBOSS_" + (++_context.SeqNameCount).ToString("000", System.Globalization.CultureInfo.InvariantCulture);
            return;
        }

        t.Name = EmbossFastaReader.SetName(str);
    }

    // seqSetNameNospace
    private static string NameNospace(string str)
    {
        if (IsWord(str))
            return str;
        var words = str.Replace('\t', ' ').Replace('\n', ' ').Replace('\r', ' ').Split(' ', StringSplitOptions.RemoveEmptyEntries);
        return string.Join(' ', words).Replace(' ', '_');
    }

    private static void AccSave(EmbossSeq t, string acc)
    {
        t.Acclist.Add(acc);
        if (t.Acc.Length == 0)
            t.Acc = acc;
    }

    private static void SvSave(EmbossSeq t, string sv)
    {
        if (t.Sv.Length == 0)
            t.Sv = sv;
    }

    private static void TaxSave(EmbossSeq t, string tax, int level)
    {
        switch (level)
        {
            case 1:
                if (t.Tax.Length == 0)
                    t.Tax = tax;
                break;
            case 2:
                if (t.Organelle.Length == 0)
                    t.Organelle = tax;
                break;
            case 3:
                if (t.Taxcommon.Length == 0)
                    t.Taxcommon = tax;
                break;
            default:
                t.Taxlist.Add(tax);
                break;
        }
    }

    // ajSeqParseFasta: (id, acc, sv, desc) of a header line that starts with '>'
    private static bool ParseFasta(string instr, out string id, out string acc, out string sv, out string desc)
    {
        id = acc = sv = desc = "";
        if (!instr.StartsWith('>'))
            return false;
        var h = new StrTok(instr, "> ");
        h.Next(out id, " \t\n\r");
        bool ok = h.Next(out string token);
        string token2 = token.Replace("(", "", StringComparison.Ordinal).Replace(")", "", StringComparison.Ordinal);
        if (ok && EmbossFastaReader.SeqVersionAccession(token2) is { } a)
        {
            acc = a;
            sv = token2;
            h.Next(out desc, "\n\r");
        }
        else if (ok && EmbossFastaReader.IsAccession(token2))
        {
            acc = token2;
            h.Next(out desc, "\n\r");
        }
        else if (ok)
        {
            desc = token;
            if (h.Next(out string more, "\n\r"))
                desc = desc + " " + more;
        }

        return true;
    }

    // the sequence lines of a FASTA-like entry up to the next '>' line (kept in the buffer)
    private void ReadFastaBody(EmbossSeq t, bool warn)
    {
        _buff.SetUnbuffered();
        bool ok = ReadLine();
        while (ok && !_line.StartsWith('>'))
        {
            if (warn)
                AppendWarn(t, _line);
            else
                Append(t, _line);
            ok = ReadLine();
        }

        if (ok)
            StoreClear(1);
        else
            _buff.Clear(0);
    }

    // seqReadFasta (pearson): the id is the first word, not parsed further
    private bool ReadFasta(EmbossSeq t)
    {
        if (!ReadLine())
            return false;
        if (_line.Length > 3 && _line[3] == ';')
            _line = ">" + _line[4..];
        if (!_line.StartsWith('>'))
        {
            StoreReset();
            return false;
        }

        if (!ParseFasta(_line, out string id, out string acc, out string sv, out string desc))
        {
            StoreReset();
            return false;
        }

        _buff.SetUnbuffered();
        t.Name = NameNospace(id);
        if (sv.Length > 0)
            SvSave(t, sv);
        if (acc.Length > 0)
            AccSave(t, acc);
        t.Desc = desc;
        ReadFastaBody(t, true);
        return true;
    }

    // ajSeqParseNcbi: an NCBI-style header, else the plain FASTA header; false for a header with ';' at position 3 (PIR) or none that starts with '>'
    private static EmbossNcbiHeader? ParseNcbi(string line)
    {
        if (line.Length > 3 && line[3] == ';')
            return null;
        if (line.Length == 0 || line[0] != '>')
            return null;
        var h = EmbossNcbiHeader.Parse(line);
        if (h is not null)
            return h;
        return ParseFasta(line, out string id, out string acc, out string sv, out string desc) ? new EmbossNcbiHeader(id, acc, sv, "", "", desc) : null;
    }

    // seqReadNcbi (fasta, ncbi)
    private bool ReadNcbi(EmbossSeq t)
    {
        if (!ReadLine())
            return false;
        var h = ParseNcbi(_line);
        if (h is null)
        {
            StoreReset();
            return false;
        }

        t.Setdb = h.Db;
        if (h.Gi.Length > 0)
            t.Gi = h.Gi;
        if (h.SequenceVersion.Length > 0)
            SvSave(t, h.SequenceVersion);
        if (h.Accession.Length > 0)
            AccSave(t, h.Accession);
        SetName(t, h.Id);
        t.Desc = h.Description;
        ReadFastaBody(t, true);
        return true;
    }

    // seqReadGifasta
    private bool ReadGifasta(EmbossSeq t)
    {
        if (!ReadLine())
            return false;
        var h = ParseNcbi(_line);
        if (h is null || h.Gi.Length == 0)
        {
            StoreReset();
            return false;
        }

        t.Gi = h.Gi;
        t.Setdb = h.Db;
        if (h.SequenceVersion.Length > 0)
            SvSave(t, h.SequenceVersion);
        if (h.Accession.Length > 0)
            AccSave(t, h.Accession);
        SetName(t, h.Gi);
        t.Desc = h.Description;
        ReadFastaBody(t, false);
        return true;
    }

    // seqReadDbId
    private bool ReadDbId(EmbossSeq t)
    {
        if (!ReadLine())
            return false;
        if (_line.Length > 3 && _line[3] == ';')
            return false;
        if (!_line.StartsWith('>'))
        {
            StoreReset();
            return false;
        }

        var h = new StrTok(_line, "> ");
        h.Next(out _, " \t\n\r");                         // the database name
        h.Next(out string token, " \t\n\r");
        SetName(t, token);
        h.Next(out token);
        if (EmbossFastaReader.SeqVersionAccession(token) is { } vacc)
        {
            SvSave(t, token);
            AccSave(t, vacc);
            h.Next(out string d, "\n\r");
            t.Desc = d;
        }
        else if (EmbossFastaReader.IsAccession(token))
        {
            AccSave(t, token);
            h.Next(out string d, "\n\r");
            t.Desc = d;
        }
        else
        {
            t.Desc = token;
            if (h.Next(out string more, "\n\r"))
                t.Desc = t.Desc + " " + more;
        }

        ReadFastaBody(t, false);
        return true;
    }

    // seqReadText (textonly, plain, asis, unknown)
    private bool ReadText(EmbossSeq t)
    {
        while (ReadLine())
            Append(t, _line);
        _buff.Clear(0);
        if (_records == 0)
            return false;
        // seqSetNameFile
        t.Name = EmbossFastaReader.SetName(_opt.Filename);
        if (_count > 1)
            t.Name += "_" + _count.ToString(System.Globalization.CultureInfo.InvariantCulture).PadLeft(3);
        return true;
    }

    // seqReadGcg
    private bool ReadGcg(EmbossSeq t)
    {
        bool ok = ReadLine();
        while (ok && IsWhite(_line))
        {
            StoreReset();
            ok = ReadLine();
        }

        if (!ok)
            return false;
        if (_line.StartsWith("!!NA_SEQUENCE", StringComparison.Ordinal))
            t.SetNuc();
        else if (_line.StartsWith("!!AA_SEQUENCE", StringComparison.Ordinal))
            t.SetProt();
        uint len = 0;
        if (!GcgDots(t, 5000, ref len))
        {
            StoreReset();
            return false;
        }

        bool incomment = false;
        while (ok && (uint)t.Seq.Length < len)
        {
            ok = ReadLine();
            if (ok)
                AppendCommented(t, ref incomment, _line);
        }

        _buff.Clear(0);
        return ok;
    }

    private static readonly Regex RegGcgDot = new(@"[.][.]", RegexOptions.CultureInvariant);
    private static readonly Regex RegGcgChk = new(@"[Cc][Hh][Ee][Cc][Kk]:[ \t]*([0-9]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegGcgLen = new(@"[Ll][Ee][Nn][Gg][Tt][Hh]:[ \t]*([0-9]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegGcgTyp = new(@"[Tt][Yy][Pp][Ee]:[ \t]*([NP])", RegexOptions.CultureInvariant);
    private static readonly Regex RegGcgNam = new(@"[^ \t>]+", RegexOptions.CultureInvariant);
    private static readonly Regex RegGcgMsf = new(@"[Mm][Ss][Ff]:[ \t]*([0-9]+)", RegexOptions.CultureInvariant);

    // seqGcgDots: the line with ".." and "Check:" (searched over at most maxlines lines)
    private bool GcgDots(EmbossSeq t, int maxlines, ref uint len)
    {
        int nlines = 0;
        while (nlines < maxlines)
        {
            if (nlines++ > 0)
            {
                if (!ReadLine())
                    return false;
            }

            if (nlines > maxlines)
                return false;
            if (!RegGcgDot.IsMatch(_line))
                continue;
            var chk = RegGcgChk.Match(_line);
            if (!chk.Success)
                return false;
            if (RegGcgMsf.IsMatch(_line))
                return false;
            CStr.ToUint(chk.Groups[1].Value, out _);
            var l = RegGcgLen.Match(_line);
            if (l.Success)
            {
                CStr.ToUint(l.Groups[1].Value, out uint v);
                len = v;
            }

            var nam = RegGcgNam.Match(_line);
            if (nam.Success)
                t.Name = nam.Value;
            var typ = RegGcgTyp.Match(_line);
            if (typ.Success)
                t.Type = typ.Groups[1].Value;
            return true;
        }

        return false;
    }

    private static readonly Regex RegStadenId = new(@"^[<]([^>-]+)[-]*[>]", RegexOptions.CultureInvariant);

    // seqReadStaden
    private bool ReadStaden(EmbossSeq t)
    {
        bool incomment = false;
        if (!ReadLine())
            return false;
        var m = RegStadenId.Match(_line);
        if (m.Success)
        {
            SetName(t, m.Groups[1].Value);
            AppendCommented(t, ref incomment, _line[(m.Index + m.Length)..]);
        }
        else
        {
            SetName(t, FilenameTrimAll(_opt.Filename));
            AppendCommented(t, ref incomment, _line);
        }

        while (ReadLine())
            AppendCommented(t, ref incomment, _line);
        if (_records == 0)
            return false;
        _buff.Clear(0);
        return true;
    }

    private static readonly Regex RegEntryExp = new(@":([A-Za-z0-9_-]+)$", RegexOptions.CultureInvariant);
    private static readonly Regex RegFileExp = new(@"([A-Za-z0-9_-]+)[.][A-Za-z0-9_-]+$", RegexOptions.CultureInvariant);
    private static readonly Regex RegRestExp = new(@"([A-Za-z0-9_-]+)[^A-Za-z0-9_-]*$", RegexOptions.CultureInvariant);

    // ajFilenameTrimAll
    internal static string FilenameTrimAll(string name)
    {
        var m = RegEntryExp.Match(name);
        if (m.Success)
            return m.Groups[1].Value;
        m = RegFileExp.Match(name);
        if (m.Success)
            return m.Groups[1].Value;
        m = RegRestExp.Match(name);
        return m.Success ? m.Groups[1].Value : "";
    }

    // seqReadStrider
    private bool ReadStrider(EmbossSeq t)
    {
        bool ok;
        do
        {
            ok = ReadLine();
            if (ok && _line.StartsWith("; DNA sequence", StringComparison.Ordinal))
            {
                var h = new StrTok(_line, " ;\t,\n");
                h.Next(out _);
                h.Next(out _);
                h.Next(out string token);
                _seqToken = token;
            }
        }
        while (ok && _line.StartsWith(';'));

        if (!ok || _seqToken.Length == 0)
        {
            StoreReset();
            _seqToken = "";
            return false;
        }

        _buff.SetUnbuffered();
        SetName(t, _seqToken);
        while (ok && !_line.StartsWith("//", StringComparison.Ordinal))
        {
            Append(t, _line);
            ok = ReadLine();
        }

        _buff.Clear(0);
        _seqToken = "";
        return true;
    }

    // seqReadCodata
    private bool ReadCodata(EmbossSeq t)
    {
        if (!ReadLine())
            return false;
        if (!_line.StartsWith("ENTRY ", StringComparison.Ordinal))
        {
            StoreReset();
            return false;
        }

        var h = new StrTok(_line, " \n\r");
        h.Next(out _);
        h.Next(out string token);
        SetName(t, token);
        bool ok = ReadLine();
        while (ok && !_line.StartsWith("SEQUENCE", StringComparison.Ordinal))
        {
            bool done = false;
            if (_line.StartsWith("ACCESSION ", StringComparison.Ordinal))
            {
                var a = new StrTok(_line, " ;\n\r");
                a.Next(out _);
                a.Next(out token);
                AccSave(t, token);
            }

            if (_line.StartsWith("TITLE ", StringComparison.Ordinal))
            {
                var d = new StrTok(_line, " ");
                d.Next(out _);
                d.Next(out string desc, "\n\r");
                var descText = new StringBuilder(desc);
                ok = ReadLine();
                done = true;
                while (ok && _line.StartsWith(' '))
                {
                    var c = new StrTok(_line, " ");
                    c.Next(out token, "\n\r");
                    descText.Append(' ').Append(token);
                    ok = ReadLine();
                }

                t.Desc = descText.ToString();
            }

            if (!done)
                ok = ReadLine();
        }

        ok = ReadLine();
        while (ok && !_line.StartsWith("///", StringComparison.Ordinal))
        {
            Append(t, _line);
            ok = ReadLine();
        }

        _buff.Clear(0);
        return true;
    }

    private static readonly Regex RegNbrfId = new(@"^>(..)[>;]([^ \t\n]+)", RegexOptions.CultureInvariant);

    // seqReadNbrf (pir)
    private bool ReadNbrf(EmbossSeq t)
    {
        bool skipheader = true;
        while (skipheader)
        {
            if (!ReadLine())
                return false;
            if (!_line.StartsWith("C;", StringComparison.Ordinal) && !IsWhite(_line))
                skipheader = false;
        }

        var m = RegNbrfId.Match(_line);
        if (!m.Success)
        {
            StoreReset();
            return false;
        }

        string token = m.Groups[1].Value;
        t.Name = m.Groups[2].Value;
        string idline = _line.Length > 4 ? _line[4..] : "";
        switch (char.ToUpperInvariant(token[0]))
        {
            case 'P':
            case 'F':
                t.SetProt();
                break;
            case 'B':
            case 'D':
            case 'R':
                t.SetNuc();
                break;
            default:
                _context.Warn($"Unknown NBRF sequence type '{token}'");
                break;
        }

        if (!ReadLine())
        {
            StoreReset();
            return false;
        }

        t.Desc = _line.EndsWith('\n') ? _line[..^1] : _line;
        bool ok = ReadLine();
        bool dofeat = false;
        bool tryfeat = _opt.Features && !_opt.UfoGiven;
        List<string>? featureLines = null;
        while (ok && !_line.StartsWith('>'))
        {
            if (_line.Length < 2 || _line[1] != ';')
                Append(t, _line);
            else
            {
                if (_line.StartsWith("C;Accession:", StringComparison.Ordinal))
                {
                    var tk = new StrTok(_line.Length > 13 ? _line[13..] : "", " ;\n\r");
                    while (tk.Next(out string a))
                        AccSave(t, a);
                }

                if (_line.StartsWith("C;Species:", StringComparison.Ordinal))
                {
                    var tk = new StrTok(_line.Length > 11 ? _line[11..] : "", ";.\n\r");
                    while (tk.Next(out string a))
                        TaxSave(t, a, 1);
                }

                if (_line.Length > 0 && _line[0] == 'R')
                {
                    while ((ok = ReadLine()))
                    {
                        if ((_line.Length > 1 && _line[1] == ';') || (_line.Length > 0 && _line[0] == '>'))
                            break;
                    }

                    if (ok)
                        continue;
                }
                else if (tryfeat && _line.Length > 0 && _line[0] == 'F')
                {
                    dofeat = true;
                    (featureLines ??= []).Add(Chomp(_line));
                }
            }

            if (ok)
                ok = ReadLine();
            if (ok && t.Seq.Length == 0 && _line.Contains('>'))
            {
                string tmp = _line.TrimStart('0', '1', '2', '3', '4', '5', '6', '7', '8', '9');
                tmp = tmp.Length > 4 ? tmp[4..] : "";
                if (tmp == idline)
                {
                    ok = ReadLine();
                    if (ok && !IsWhite(_line))
                        ok = ReadLine();
                }
            }
        }

        if (t.Seq.Length > 0 && t.Seq[^1] == '*')
            t.Seq.Length--;
        if (ok)
            StoreClear(1);
        else
            _buff.Clear(0);
        if (dofeat && featureLines is not null)
        {
            var table = new EmbossSeqFeatureTable(_context) { SeqId = t.Name };
            table.Found = table.ReadPir(new EmbossFeatureLines(featureLines));
            if (table.Found)
            {
                table.Length = t.Seq.Length;
                table.SeqId = t.Name;
                t.Fttable = table;
            }
        }

        return true;
    }
}

/// <summary>The program died while reading its sequence input (the library message, <c>Died: ...</c>).</summary>
public sealed class EmbossSeqinException(IReadOnlyList<string> messages, bool fatal) : Exception(string.Join("\n", messages))
{
    public IReadOnlyList<string> Messages { get; } = messages;

    public bool Fatal { get; } = fatal;
}
