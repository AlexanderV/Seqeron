using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The program stopped as the library does ("Died: ..."): <see cref="Exception.Message"/> is the text of the last line, <see cref="Errors"/> the lines printed before it.</summary>
public sealed class EmbossFeatureDiedException(string message, IReadOnlyList<string> errors) : Exception(message)
{
    public IReadOnlyList<string> Errors { get; } = errors;
}

/// <summary>
/// The <c>features</c> / <c>featuresall</c> input of an ACD program (<c>acdSetFeatures</c> of <c>ajacd.c</c>, <c>feattabinQryProcess</c> / <c>ajQuerystrParse*</c> of <c>ajfeatread.c</c> / <c>ajquery.c</c>): the value is a file name with an optional <c>format::</c>
/// prefix and <c>[begin:end:r]</c> range suffix (the <c>-fformat</c> qualifier is set on the wrong field of the input object and has no effect), the tables are read from it one after the other, the first with the range
/// <c>-fbegin</c> / <c>-fend</c> applied. A5-8c3e, B09 F241.
/// </summary>
internal sealed class EmbossFeatureInput
{
    private static readonly Regex RegFmt = new("^([A-Za-z0-9-]*)::(.*)$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
    private static readonly Regex RegRange = new(@"^(.*)[[](-?[0-9]*):(-?[0-9]*)(:(r)?)?[]]$", RegexOptions.CultureInvariant | RegexOptions.IgnoreCase);

    private readonly EmbossFeatureStream _stream = null!;
    private readonly EmbossFeatureReadContext _context;
    private readonly EmbossSeqFeatureTable? _first;
    private bool _firstReturned;

    /// <summary>The messages the library printed (<c>Warning: ...</c>, <c>Error: ...</c>) so far, in order.</summary>
    public IReadOnlyList<string> Errors => _context.Messages;

    /// <param name="value">The value of the <c>-features</c> parameter.</param>
    /// <param name="files">The files by name: null when the file does not exist.</param>
    /// <param name="fbegin">The value of <c>-fbegin</c> (default 1).</param>
    /// <param name="fend">The value of <c>-fend</c> (default: the length of the first table).</param>
    /// <param name="context">The state of the feature library of the run.</param>
    /// <exception cref="EmbossFeatureDiedException">The value is no feature file ("Bad value for '-features' and no more retries"): the ACD asks twice, so every message of the reading is printed twice.</exception>
    public EmbossFeatureInput(string value, Func<string, string?> files, int fbegin, int? fend, EmbossFeatureReadContext context)
    {
        _context = context;
        for (int attempt = 0; attempt < 2 && _first is null; attempt++)
        {
            var stream = Open(value, files);
            if (stream is null)
                continue;
            _stream = stream;
            _first = stream.Next();
        }

        if (_first is null)
            throw new EmbossFeatureDiedException("Died: featcopy terminated: Bad value for '-features' and no more retries", _context.Messages);
        int len = _first.Length;
        int end = fend ?? len;
        _first.Start = unchecked((uint)PositionOf(len, 1, fbegin));
        _first.End = unchecked((uint)PositionOf(len, (int)_first.Start, end));
    }

    // one attempt of the ACD: the query string is parsed (ajQuerystrParseFormat: a format name that is none is an error, twice) and the file opened
    private EmbossFeatureStream? Open(string value, Func<string, string?> files)
    {
        string qry = value.Trim(' ', '\t', '\n');
        int format = 0;
        var m = RegFmt.Match(qry);
        if (m.Success)
        {
            string name = m.Groups[1].Value.Length > 0 ? m.Groups[1].Value : "unknown";
            qry = m.Groups[2].Value;
            int f = EmbossFeatureFormats.Find(name);
            if (f >= 0)
            {
                format = f;
            }
            else
            {
                _context.Error($"Unknown input format '{name}'");
                _context.Error($"Unknown input format '{name}'\n");
            }
        }

        var mr = RegRange.Match(qry);
        if (mr.Success)
            qry = mr.Groups[1].Value;
        string? text = files(qry);
        if (text is null)
        {
            _context.Error($"Failed to open filename '{qry}'");
            return null;
        }

        return new EmbossFeatureStream(text, format, _context);
    }

    // ajFeattablePosII (unsigned: a negative result is a huge number, which is then cut to the length)
    private static int PositionOf(int ilen, int imin, int ipos)
    {
        uint jpos;
        if (ipos < 0)
            jpos = unchecked((uint)(ilen + ipos + 1));
        else
            jpos = ipos != 0 ? (uint)ipos : 1u;
        if (jpos > unchecked((uint)ilen))
            jpos = unchecked((uint)ilen);
        if (jpos < unchecked((uint)imin))
            jpos = unchecked((uint)imin);
        return unchecked((int)jpos);
    }

    /// <summary><c>ajFeattaballNext</c>: the next table (the first is the one the ACD read), or null.</summary>
    public EmbossSeqFeatureTable? Next()
    {
        if (!_firstReturned)
        {
            _firstReturned = true;
            return _first;
        }

        return _stream.Next();
    }

    /// <summary>The state of the feature library the tables are read with.</summary>
    public EmbossFeatureReadContext Context => _context;
}

/// <summary>
/// The features of a sequence read from a separate file (<c>ajFeattableNewReadUfo</c> of <c>ajfeatread.c</c>, called by <c>seqReadFmt</c> when the sequence input has features, none came with the sequence and <c>-ufo</c>, <c>-fformat</c>
/// or <c>-fopenfile</c> name them): the UFO is <c>[format:]file</c>, the format of <c>-fformat</c> is used when the UFO has no prefix, the file is <c>-fopenfile</c> or <c>NAME.FORMAT</c> when the UFO names none; with no format every
/// format of the table is tried in turn. A5-8c3e, B09 F241.
/// </summary>
internal static class EmbossFeatureUfo
{
    private static readonly Regex RegFmt = new("^([A-Za-z0-9][A-Za-z0-9]+):+(.*)$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
    private static readonly Regex RegFile = new("^(([A-Za-z]:)?[^:]+)$", RegexOptions.CultureInvariant | RegexOptions.Singleline);

    /// <summary>The table of the features, or null (no file, no features, a format that reads none).</summary>
    /// <param name="ufo">The UFO (<c>-ufo</c>), may be empty.</param>
    /// <param name="formatOption">The value of <c>-fformat</c>, null / empty = none.</param>
    /// <param name="fileOption">The value of <c>-fopenfile</c>, null / empty = none.</param>
    /// <param name="seqname">The name of the sequence (the name the file is taken from when nothing names it).</param>
    /// <param name="files">The files by name: null when the file does not exist.</param>
    /// <param name="context">The state of the feature library of the run.</param>
    public static EmbossSeqFeatureTable? Read(string ufo, string? formatOption, string? fileOption, string seqname, Func<string, string?> files, EmbossFeatureReadContext context)
    {
        string formatstr = formatOption ?? "";
        string filename = fileOption ?? "";
        string test = ufo;
        int format = 0;
        if (ufo.Length > 0)
        {
            var m = RegFmt.Match(ufo);
            if (m.Success)
            {
                formatstr = m.Groups[1].Value.Length > 0 ? m.Groups[1].Value : "unknown";
                test = m.Groups[2].Value;
                int f = EmbossFeatureFormats.Find(formatstr);
                if (f >= 0)
                {
                    format = f;
                }
                else
                {
                    context.Error($"Unknown input feat format '{formatstr}'");
                    context.Error($"unknown input feature table format '{formatstr}'\nNO Features will be read in");
                }
            }
        }

        // featFormatSet
        if (formatstr.Length > 0)
        {
            int f = EmbossFeatureFormats.Find(formatstr);
            if (f >= 0)
                format = f;
            else
                context.Error($"Unknown input feat format '{formatstr}'");
        }

        var mf = RegFile.Match(test);
        if (mf.Success)
        {
            filename = mf.Groups[1].Value;
        }
        else if (seqname.Length > 0 && formatstr.Length > 0)
        {
            if (filename.Length == 0)
                filename = seqname + "." + formatstr;
        }
        else
        {
            return null;
        }

        string? text = files(filename);
        if (text is null)
            return null;
        var file = EmbossFeatureLines.FromText(text);
        if (format != 0)
            return ReadTable(file, format, seqname, context);
        for (int i = 1; i < EmbossFeatureFormats.Table.Length; i++)
        {
            var ret = ReadTable(file, i, seqname, context);
            if (ret is not null)
                return ret;
            file.Reset();
        }

        return null;
    }

    // ajFeattableNewRead: blank lines are skipped, then the reader of the format reads the table
    private static EmbossSeqFeatureTable? ReadTable(EmbossFeatureLines file, int format, string seqname, EmbossFeatureReadContext context)
    {
        string? line;
        while ((line = file.Next()) is not null)
        {
            if (!line.All(EmbossSeqFeatureTable.IsCSpace))
            {
                file.Clear(1);
                break;
            }
        }

        var table = new EmbossSeqFeatureTable(context) { SeqId = seqname };
        return EmbossFeatureFormats.Read(table, format, file) ? table : null;
    }
}
