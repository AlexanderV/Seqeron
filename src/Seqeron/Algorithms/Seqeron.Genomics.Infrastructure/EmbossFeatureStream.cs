#pragma warning disable S1121, S3358, S2325 // the readers keep the control flow of ajfeatread.c

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The input formats of feature tables (<c>featinformatDef</c> of <c>ajfeatread.c</c>).</summary>
internal static class EmbossFeatureFormats
{
    /// <summary>name, EDAM term, alias, tried when the format is unknown, nucleotide, protein.</summary>
    public static readonly (string Name, string Obo, bool Alias, bool Try, bool Nuc, bool Prot)[] Table =
    [
        ("unknown", "0000", false, false, false, false), ("embl", "1927", false, true, true, false), ("em", "1927", true, false, true, false), ("genbank", "1936", true, false, true, false),
        ("gb", "1936", true, false, true, false), ("ddbj", "1936", true, false, true, false), ("refseq", "1958", false, true, true, false), ("refseqp", "0000", false, true, false, true),
        ("genpept", "0000", false, true, false, true), ("gff3", "1939", false, true, true, true), ("gff3emboss", "1939", false, false, true, true), ("gff2", "1938", false, true, true, true),
        ("gff", "1939", true, false, true, true), ("swiss", "1963", false, true, false, true), ("sw", "1963", true, false, false, true), ("uniprot", "2188", true, false, false, true),
        ("swissprot", "0000", true, true, false, true), ("pir", "1948", false, true, false, true), ("nbrf", "1948", true, false, false, true), ("dasgff", "1978", false, true, true, false),
        ("chado", "0000", false, true, true, false),
    ];

    /// <summary><c>feattabinformatFind</c>: the index of a format by its name (any case) or its EDAM term; -1 when unknown ("Unknown input format").</summary>
    public static int Find(string format)
    {
        if (format.Length == 0)
            return -1;
        string lower = format.ToLowerInvariant();
        for (int i = 0; i < Table.Length; i++)
        {
            if (lower == Table[i].Name || format == Table[i].Obo)
                return i;
        }

        return -1;
    }

    /// <summary>The reader of the format on a table (<c>featinformatDef[i].Read</c>).</summary>
    public static bool Read(EmbossSeqFeatureTable table, int format, EmbossFeatureLines file)
    {
        switch (Table[format].Name)
        {
            case "unknown" or "gff3" or "gff":
                return table.ReadGff3(file);
            case "embl" or "em" or "genbank" or "gb" or "ddbj":
                return table.ReadEmbl(file);
            case "refseq":
                return table.ReadFiveBlanks(file, false, false);
            case "refseqp":
                return table.ReadFiveBlanks(file, true, true);
            case "genpept":
                return table.ReadFiveBlanks(file, true, false);
            case "gff3emboss":
                return table.ReadGff3Old(file);
            case "gff2":
                return table.ReadGff2(file);
            case "swiss" or "sw" or "uniprot" or "swissprot":
                return table.ReadSwiss(file);
            case "pir" or "nbrf":
                return table.ReadPir(file);
            case "dasgff":
                return table.ReadDasgff(file);
            default:
                return false;                                  // chado: the data comes from a database only
        }
    }
}

/// <summary>
/// A feature input stream (<c>AjPFeattaball</c> / <c>ajFeattabinRead</c> of <c>ajfeatread.c</c>): the tables of ONE file, read one after the other; with no format given every format that is tried for an unknown input is tried in turn on the start of the file
/// (<c>feattabinRead</c>) and the one that finds features is kept for the next tables. A5-8c3e, B09 F241.
/// </summary>
internal sealed class EmbossFeatureStream
{
    private readonly EmbossFeatureLines _file;
    private readonly EmbossFeatureReadContext _context;
    private int _format;
    private bool _buffered = true;
    private int _count;

    /// <param name="text">The text of the file.</param>
    /// <param name="format">The index of the input format in <see cref="EmbossFeatureFormats.Table"/>; 0 = find it.</param>
    /// <param name="context">The state of the feature library of the run.</param>
    public EmbossFeatureStream(string text, int format, EmbossFeatureReadContext context)
    {
        _file = EmbossFeatureLines.FromText(text);
        _format = format;
        _context = context;
    }

    /// <summary>The index of the format the stream reads (0 until it is found).</summary>
    public int Format => _format;

    /// <summary>The number of tables read.</summary>
    public int Count => _count;

    /// <summary><c>ajFeattabinRead</c>: the next table, or null when there is none (or the text is of no format).</summary>
    public EmbossSeqFeatureTable? Next()
    {
        _count++;
        if (_buffered)
        {
            // skip blank lines: the first line that is not white is given back
            string? line;
            while ((line = _file.Next()) is not null)
            {
                if (!IsWhite(line))
                {
                    _file.Unread();
                    break;
                }
            }
        }

        if (_format == 0)
        {
            for (int i = 1; i < EmbossFeatureFormats.Table.Length; i++)
            {
                if (!EmbossFeatureFormats.Table[i].Try)
                    continue;
                var table = TryRead(i);
                if (table is not null)
                {
                    _format = i;
                    return table;
                }

                _file.Reset();
            }

            return null;
        }

        _buffered = false;
        bool atEnd = _file.AtEnd;
        var next = TryRead(_format);
        if (next is not null && atEnd)
            throw new InvalidOperationException("The reference program never stops here: a reader returns a table at the end of the file again and again (the last feature of a static pointer of the library).");
        return next;
    }

    private EmbossSeqFeatureTable? TryRead(int format)
    {
        var table = new EmbossSeqFeatureTable(_context);
        if (!EmbossFeatureFormats.Read(table, format, _file))
            return null;
        table.FormatName = EmbossFeatureFormats.Table[format].Name;
        return table;
    }

    private static bool IsWhite(string s) => s.All(EmbossSeqFeatureTable.IsCSpace);
}
