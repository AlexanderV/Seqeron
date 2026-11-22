#pragma warning disable S1121, S3776, S107, S1854

using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The <c>seqall</c> input of an EMBOSS program: the ACD item (<c>acdSetSeqall</c>: the USA, the associated qualifiers <c>-sformat -sdbname -sid -iquery -ioffset -squick -scircular -supper -slower -sbegin -send -sreverse -snucleotide -sprotein</c>),
/// <c>ajSeqAllRead</c> / <c>ajSeqRead</c> (the USA of <c>seqinUsaProcess</c>: <c>format::</c> prefix, <c>[begin:end:r]</c> range, <c>list::</c> / <c>@</c> list files of <c>seqinListProcess</c>, queries on files, the text of <c>asis::</c>), the readers of
/// <see cref="EmbossSeqin"/> for every sequence of the file, <c>ajSeqallNext</c> (<c>ajSeqSetRange</c> / <c>ajSeqSetRangeRev</c>) and the messages of the library. A5-8c3g, B09 F248; the USA forms, ranges and queries A5-8c3h, B09 F255.
/// </summary>
internal sealed class EmbossSeqall
{
    private sealed class UsaNode
    {
        public string Usa = "", Formatstr = "";
        public int Begin, End;
        public bool Rev, Features;
        public long Fpos;
    }

    private readonly EmbossFeatureReadContext _context;
    private readonly EmbossSeqinOptions _opt = new();
    private readonly EmbossTextin _textin = new();
    private readonly EmbossQueryEnv _env;
    private readonly IReadOnlyDictionary<string, string> _files;
    private List<UsaNode>? _usalist;
    private EmbossSeqin? _input;
    private int _depth;
    private readonly Dictionary<string, int> _expansions = [];
    private int _count;

    // the seqall object
    private int _begin, _end;
    private bool _rev;
    private EmbossSeq? _first;

    private EmbossSeqall(EmbossFeatureReadContext context, IReadOnlyDictionary<string, string> files, EmbossNames? names)
    {
        _context = context;
        _files = files;
        _env = new EmbossQueryEnv
        {
            Files = name => files.TryGetValue(name, out string? t) ? t : null,
            Names = names,
            Directory = dir => ListDirectory(files, dir),
        };
    }

    // the names of the files of a directory (the keys of the files that have the directory as their path); null when no file lies in it
    private static IReadOnlyList<string>? ListDirectory(IReadOnlyDictionary<string, string> files, string dir)
    {
        string prefix = dir == "." ? "" : dir + "/";
        List<string>? names = null;
        foreach (string key in files.Keys)
        {
            if (!key.StartsWith(prefix, StringComparison.Ordinal))
                continue;
            string rest = key[prefix.Length..];
            names ??= [];
            if (!rest.Contains('/'))
                names.Add(rest);
        }

        return names;
    }

    /// <summary>
    /// Opens the input and reads its first sequence (<c>acdSetSeqall</c>). Null: the program dies ("Unable to read sequence", the messages of the library are in <paramref name="context"/>; the caller adds the <c>Died:</c> line).
    /// </summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    /// <param name="single">The input is an ACD <c>sequence</c> (<c>acdSetSeq</c>: <c>ajSeqRead</c>, one sequence, reversed at once) and not a <c>seqall</c> (<c>acdSetSeqall</c>: <c>ajSeqAllRead</c>).</param>
    /// <param name="storeText">The readers store the text of every entry (the ACD attribute <c>entry: Y</c> of <c>entret</c>).</param>
    /// <param name="definitions">The texts of the configuration files (<c>.embossrc</c>) that define databases for the USAs <c>db:id</c>; null = none.</param>
    public static EmbossSeqall? Open(string program, string spec, IReadOnlyDictionary<string, string> files, EmbossSequenceQualifiers q, bool features, EmbossFeatureReadContext context, string sequenceType = "gapany", bool single = false, IReadOnlyList<string>? definitions = null, bool storeText = false)
    {
        context.Program = program;
        EmbossNames? names = definitions is { Count: > 0 } ? EmbossNames.Parse([.. definitions]) : null;
        if (names is not null)
            context.Messages.AddRange(names.Messages);
        var all = new EmbossSeqall(context, files, names);
        var o = all._opt;
        // ajSeqinUsa, then the associated qualifiers
        all._textin.Qry = spec;
        o.Features = features;
        all._textin.Formatstr = q.Format ?? "";
        o.Db = q.DbName ?? "";
        o.Entryname = q.EntryId ?? "";
        all._textin.QryFields = q.QueryFields ?? "";
        all._textin.Fpos = q.Offset;
        o.Circular = q.Circular;
        o.Upper = q.Upper;
        o.Lower = q.Lower;
        o.Begin = q.Begin;
        o.End = q.End;
        o.UfoGiven = !string.IsNullOrEmpty(q.FeaturesUfo);
        o.SequenceType = sequenceType;
        o.IsNuc = q.Nucleotide;
        o.IsProt = q.Protein;
        o.Query = all._textin.Query;
        o.StoreText = storeText;
        if (o.UfoGiven)
            o.Features = true;

        EmbossSeq? first;
        try
        {
            first = single ? all.ReadNext() : all.ReadFirst();
        }
        catch (EmbossFatalException e)
        {
            context.Messages.Add(e.Message);
            context.Fatal = true;
            return null;
        }

        if (first is null)
        {
            context.Error($"Unable to read sequence '{spec}'");
            return null;
        }

        if (single)
        {
            // acdSetSeq: the range comes with the sequence (the begin and end of the input); a nucleotide sequence that is to be reversed is reversed at once
            if (first.Type == "N" && q.Reverse)
            {
                o.Rev = true;
                first.Rev = true;
            }

            if (first.Rev)
                EmbossSeqRange.ReverseDo(first);
            all._first = first;
            return all;
        }

        // now the begin, end and reverse options
        if (o.Begin != 0)
        {
            all._begin = first.Begin = o.Begin;
        }

        if (o.End != 0)
        {
            all._end = first.End = o.End;
        }

        if (first.Type == "N" && q.Reverse)
        {
            o.Rev = true;
            first.Rev = true;
            all._rev = true;
        }

        all._first = first;
        return all;
    }

    /// <summary>The first sequence (the one the ACD read when it opened the input) without taking it.</summary>
    public EmbossSeq Peek() => _first ?? throw new InvalidOperationException("The first sequence was taken.");

    /// <summary><c>ajSeqallNext</c>: the next sequence with its range set (and reversed when it is to be), or null (the library reads it, and prints its messages, when it is asked).</summary>
    public EmbossSeq? Next()
    {
        EmbossSeq? seq;
        if (_count == 0)
        {
            _count = 1;
            seq = _first;
            _first = null;
        }
        else
        {
            seq = ReadNext();
            if (seq is null)
                return null;
            _count++;
        }

        if (_rev)
            EmbossSeqRange.SetRangeRev(seq!, _begin, _end);
        else
            EmbossSeqRange.SetRange(seq!, _begin, _end);
        return seq;
    }

    // ---------------------------------------------------------------------------------------------------------------- ajSeqAllRead / ajSeqRead

    private bool UsalistAny => _usalist is { Count: > 0 };

    private UsaNode Pop()
    {
        var node = _usalist![0];
        _usalist.RemoveAt(0);
        return node;
    }

    private void ClearPos()
    {
        _opt.Rev = false;
        _opt.Begin = 0;
        _opt.End = 0;
    }

    // ajTextinClearNewfile
    private void ClearNewfile() => _textin.Fpos = 0;

    // ajSeqinUsa: the input is cleared (a new USA) and the query set
    private void SeqinUsa(string usa)
    {
        _textin.Clear();
        _input = null;
        _opt.Entryname = "";
        _opt.Rev = false;
        _opt.Begin = 0;
        _opt.End = 0;
        _textin.Fpos = 0;
        _opt.Db = "";
        _textin.Qry = usa;
    }

    private void UsaRestore(UsaNode node)
    {
        _opt.Begin = node.Begin;
        _opt.End = node.End;
        _opt.Rev = node.Rev;
        _textin.Fpos = node.Fpos;
        _opt.Features = node.Features;
        _textin.Formatstr = node.Formatstr;
    }

    private UsaNode UsaSave() => new() { Begin = _opt.Begin, End = _opt.End, Rev = _opt.Rev, Fpos = _textin.Fpos, Features = _opt.Features, Formatstr = _textin.Formatstr };

    // ajSeqAllRead: the first sequence
    private EmbossSeq? ReadFirst()
    {
        bool listdata = false;
        EmbossSeq? t = null;
        if (_input is null)
        {
            if (!UsaProcess() && !UsalistAny)
                return null;
            if (UsalistAny)
                listdata = true;
            ClearNewfile();
        }

        bool ret = SeqRead(out t);
        if (ret)
        {
            if (_opt.Db.Length > 0)
                t!.Db = _opt.Db;
            if (_opt.Entryname.Length > 0)
                t!.Entryname = _opt.Entryname;
        }

        ret = RetryList(ref t, ret, ref listdata, true);
        if (!ret)
        {
            if (listdata)
                _context.Error($"Failed to read sequence '{_textin.Qry}'");
            return null;
        }

        if (_usalist is not null)
            ClearPos();
        return t;
    }

    // ajSeqRead: the next sequence
    private EmbossSeq? ReadNext()
    {
        bool ret = false, listdata = false;
        EmbossSeq? t = null;
        if (_input is not null)
        {
            ret = SeqRead(out t);
        }
        else if (UsalistAny)
        {
            listdata = true;
            var node = Pop();
            SeqinUsa(node.Usa);
            UsaRestore(node);
            if (!UsaProcess() && !UsalistAny)
                return null;
            ClearNewfile();
            ret = SeqRead(out t);
        }
        else
        {
            if (!UsaProcess() && !UsalistAny)
                return null;
            if (UsalistAny)
                listdata = true;
            ClearNewfile();
            ret = SeqRead(out t);
        }

        ret = RetryList(ref t, ret, ref listdata, false);
        if (!ret)
        {
            if (listdata)
                _context.Error($"Failed to read sequence '{_textin.Qry}'");
            return null;
        }

        return t;
    }

    // the loop of ajSeqAllRead / ajSeqRead: a read failed but the list has more USAs
    private bool RetryList(ref EmbossSeq? t, bool ret, ref bool listdata, bool first)
    {
        while (!ret && UsalistAny)
        {
            if (first || listdata)
                _context.Error($"Failed to read sequence '{_textin.Qry}'");
            listdata = true;
            var node = Pop();
            SeqinUsa(node.Usa);
            UsaRestore(node);
            if (!UsaProcess())
                continue;
            ClearNewfile();
            ret = SeqRead(out t);
        }

        return ret;
    }

    // seqRead: the next sequence of the file that is open
    private bool SeqRead(out EmbossSeq? t)
    {
        t = null;
        if (_input is null)
            return false;
        t = _input.Next();
        return t is not null;
    }

    // ---------------------------------------------------------------------------------------------------------------- seqinUsaProcess

    private static bool FormatKnown(string format) => EmbossSeqin.FindFormat(format) >= 0;

    // seqinUsaProcess
    private bool UsaProcess()
    {
        string qrystr = _textin.Qry;
        EmbossQueryParse.ParseFormat(ref qrystr, _textin, FormatKnown, _context);
        EmbossQueryParse.ParseRange(ref qrystr, ref _opt.Begin, ref _opt.End, ref _opt.Rev);
        if (EmbossQueryParse.ParseListfile(ref qrystr))
            return ListProcess(qrystr);
        var r = EmbossQueryParse.ParseRead(ref qrystr, _textin, FormatKnown, _env, _context);
        if (!r.Ok)
            return false;
        // seqinFormatSet: the format named for the input (an unknown one is reported by the reader, which then tries the formats in turn)
        _opt.Format = _textin.Formatstr.Length > 0 ? _textin.Formatstr : null;
        _opt.PrefixFormat = null;
        _opt.Filename = r.Filename;
        string text = r.Text!;
        if (r.DbName is not null)
            _opt.Db = r.DbName;                                     // textAccessDirect: the sequences of a database carry its name
        long fpos = _textin.Query.Fpos;
        _opt.FullText = null;
        _opt.AsisText = r.Asis;
        if (!r.Asis && fpos > 0)
        {
            _opt.FullText = text;
            text = fpos < text.Length ? text[(int)fpos..] : "";
        }
        Func<string, string?> lookup = name => _files.TryGetValue(name, out string? s) ? s : null;
        _input = r.Texts is { Count: > 1 } texts ? new EmbossSeqin(texts, _opt, _context, lookup) : new EmbossSeqin(text, _opt, _context, lookup);
        return true;
    }

    // seqinListProcess
    private bool ListProcess(string listfile)
    {
        _depth++;
        if (_depth > 16)
            throw new EmbossFatalException("\n   EMBOSS An error in ajseqread.c at line 17817:\nUSA List too deep");
        _usalist ??= [];
        // a list that names itself (directly or through another list) is expanded again and again: the reference program never stops and writes without end
        _expansions[listfile] = _expansions.GetValueOrDefault(listfile) + 1;
        if (_expansions[listfile] > 2000)
            throw new InvalidOperationException("The reference program never stops on this input: a list file names itself.");
        var list = new List<UsaNode>();
        if (!_files.TryGetValue(listfile, out string? text))
        {
            _context.Error($"Failed to open list file '{listfile}'");
            _depth--;
            return false;
        }

        var lines = text.Split('\n').ToList();
        if (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);
        int recnum = 0;
        foreach (string raw in lines)
        {
            string line = raw;
            if (line.EndsWith('\r'))
                line = line[..^1];
            ++recnum;
            int hash = line.IndexOf('#');
            if (hash >= 0)
                line = line[..hash];
            int p = 0;
            while (p < line.Length && IsSpace(line[p]))
                p++;
            if (p >= line.Length)
                continue;
            int e = p;
            while (e < line.Length && !IsSpace(line[e]))
                e++;
            string token = line[p..e];
            int j = e;
            while (j < line.Length && IsSpace(line[j]))
                j++;
            string rest = j < line.Length ? line[j..] : "";
            if (rest.Length > 0)
            {
                _context.Error($"Bad record {recnum} in list file '{listfile}'\n'{line}'");
            }
            else if (token.Length > 0)
            {
                var node = UsaSave();
                node.Usa = token;
                list.Add(node);
            }
        }

        // ajListPushlist: the new entries go to the START of the list
        _usalist.InsertRange(0, list);

        // now try the first item on the list (this can descend recursively if it is also a list)
        bool ret = false;
        if (_usalist.Count > 0)
        {
            var node = Pop();
            SeqinUsa(node.Usa);
            UsaRestore(node);
            ret = UsaProcess();
        }

        _depth--;
        return ret;
    }

    private static bool IsSpace(char c) => c is ' ' or '\t' or '\n' or '\v' or '\f' or '\r';
}
