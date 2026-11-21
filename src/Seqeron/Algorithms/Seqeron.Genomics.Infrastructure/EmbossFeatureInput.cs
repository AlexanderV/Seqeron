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

    // ajFeattablePosII
    private static int PositionOf(int ilen, int imin, int ipos)
    {
        int jpos;
        if (ipos < 0)
            jpos = ilen + ipos + 1;
        else
            jpos = ipos != 0 ? ipos : 1;
        if (jpos > ilen)
            jpos = ilen;
        if (jpos < imin)
            jpos = imin;
        return jpos;
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
