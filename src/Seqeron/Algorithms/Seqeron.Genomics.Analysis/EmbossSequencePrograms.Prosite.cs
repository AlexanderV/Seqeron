#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905, S2681
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The PROSITE programs of EMBOSS 6.6.0 (B09 WP67, F306): <c>embPatPrositeToRegExp</c> (a PROSITE pattern as the regular expression of the library), <c>embPatMatchFind</c> (the matches of an expression in a string),
/// <c>prosextract</c>, <c>patmatmotifs</c> and <c>patmatdb</c>.
/// </summary>
internal static partial class EmbossSequencePrograms
{
    private const string PrositeAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    private static string FatalText(string file, int line, string message) => $"\n   EMBOSS An error in {file} at line {line.ToString(CultureInfo.InvariantCulture)}:\n{message}";

    // ---------------------------------------------------------------------------------------------------------------- ajStrRemoveWhiteExcess on the buffer of the library

    /// <summary>
    /// <c>ajStrRemoveWhiteExcess</c> on the buffer of an <c>AjPStr</c>: the text is cleaned in place, so the bytes after the end of the new text keep what the old text had there. <paramref name="buf"/> holds the text, a NUL and
    /// room for one more; the length of the new text is the first NUL.
    /// </summary>
    private static void RemoveWhiteExcessInPlace(char[] p, int len)
    {
        if (len == 0)
            return;
        for (int i = 0; i < len; ++i)
            if (p[i] is '\t' or '\n' or '\r')
                p[i] = ' ';
        int k = 0;
        while (p[k] != '\0')
        {
            if (p[k] != ' ')
                break;
            ++k;
            --len;
        }

        if (k > 0)
        {
            len++;
            Array.Copy(p, k, p, 0, len);                                       // memmove(p, &p[i], len): the text and its NUL
            len = CStrLen(p);
            if (len == 0)
                return;
        }

        if (p[len - 1] == '\n')
        {
            --len;
            p[len] = '\0';
        }

        if (len == 0)
            return;
        int e;
        for (e = len; e > 0; --e)
            if (p[e - 1] != ' ')
                break;
        p[e] = '\0';
        len = CStrLen(p);
        int j = 0;
        for (int i = 0; i < len; ++i)
        {
            if (p[i] != ' ')
            {
                p[j++] = p[i];
                continue;
            }

            p[j++] = ' ';
            ++i;
            while (p[i] == ' ')
                ++i;
            p[j++] = p[i];
        }

        p[j] = '\0';
    }

    private static int CStrLen(char[] p)
    {
        int n = 0;
        while (p[n] != '\0')
            n++;
        return n;
    }

    // ---------------------------------------------------------------------------------------------------------------- embPatPrositeToRegExp

    /// <summary>
    /// <c>embPatPrositeToRegExp</c> (<c>embpat.c</c>): a PROSITE pattern as the regular expression of the library: upper case, <c>&lt;</c> / <c>&gt;</c> anchor the expression, <c>x</c> is any residue but B J O U X Z, <c>?</c> any but <c>#</c>,
    /// <c>(n,m)</c> a count, <c>[..]</c> a class, <c>{..}</c> an excluded class; the text between the brackets is copied as it stands (the characters <c>&lt;</c> and <c>&gt;</c> of it dropped).
    /// </summary>
    /// <exception cref="EmbossProgramFatalException">An unmatched bracket or a character that is no part of the PROSITE syntax (<c>ajFatal</c>).</exception>
    internal static string PrositeToRegExp(string s)
    {
        int len = s.Length;
        if (len == 0)
            return "";
        var buf = new char[len + 3];
        for (int i = 0; i < len; i++)
            buf[i] = s[i] is >= 'a' and <= 'z' ? (char)(s[i] - 32) : s[i];
        RemoveWhiteExcessInPlace(buf, len);                                    // the cleaned string; the scan for the anchors below reads as many bytes as the ORIGINAL string had (the stale bytes after the new end)
        bool isnt = false, isct = false;
        for (int i = 0; i < s.Length; ++i)
        {
            if (buf[i] == '>')
                isct = true;
            if (buf[i] == '<')
                isnt = true;
        }

        var t = new StringBuilder();
        if (isnt)
            t.Append('^');
        int p = 0;
        while (buf[p] != '\0')
        {
            char ch = buf[p];
            if (ch == '?')
            {
                t.Append("[^#]");
                ++p;
                continue;
            }

            if (ch == 'X')
            {
                t.Append("[^BJOUXZ]");
                ++p;
                continue;
            }

            if (ch == '(')
            {
                ++p;
                t.Append('{');
                while (buf[p] != ')')
                {
                    if (buf[p] == '\0')
                        throw new EmbossProgramFatalException(FatalText("embpat.c", 906, "Unmatched '(' in " + s + "\n"));
                    if (buf[p] is '<' or '>')
                    {
                        ++p;
                        continue;
                    }

                    t.Append(buf[p]);
                    ++p;
                }

                t.Append('}');
                ++p;
                continue;
            }

            if (ch == '[')
            {
                while (buf[p] != ']')
                {
                    if (buf[p] == '\0')
                        throw new EmbossProgramFatalException(FatalText("embpat.c", 929, "Unmatched '[' in " + s + "\n"));
                    if (buf[p] is '<' or '>')
                    {
                        ++p;
                        continue;
                    }

                    t.Append(buf[p]);
                    ++p;
                }

                t.Append(']');
                ++p;
                continue;
            }

            if (ch == '{')
            {
                ++p;
                t.Append("[^");
                while (buf[p] != '}')
                {
                    if (buf[p] == '\0')
                        throw new EmbossProgramFatalException(FatalText("embpat.c", 955, "Unmatched '{' in " + s + "\n"));
                    if (buf[p] is '<' or '>')
                    {
                        ++p;
                        continue;
                    }

                    t.Append(buf[p]);
                    ++p;
                }

                t.Append(']');
                ++p;
                continue;
            }

            if (PrositeAminoAcids.Contains(ch, StringComparison.Ordinal))
            {
                t.Append(ch);
                ++p;
                continue;
            }

            if (!(ch is ' ' or '-' or '>' or '<'))
                throw new EmbossProgramFatalException(FatalText("embpat.c", 982, "Unrecognised character in " + s + "\n"));
            ++p;
        }

        if (isct)
            t.Append('$');
        return t.ToString();
    }

    // ---------------------------------------------------------------------------------------------------------------- embPatMatchFind

    /// <summary>The matches of <c>embPatMatchFind</c>: the start (0-based, an unsigned value) and the length of each, in the order the library finds them.</summary>
    internal static List<(uint Start, uint Length)> PatMatchFind(string regexp, string text)
    {
        EmbossPcre rx;
        try
        {
            rx = EmbossPcre.Compile(regexp);
        }
        catch (EmbossPcreException e)
        {
            throw new InvalidOperationException($"ajRegComp refuses the expression '{regexp}' at position {e.Position.ToString(CultureInfo.InvariantCulture)}: {e.Message}; the reference program then calls pcre_exec with a NULL expression and dies by a segmentation fault: it has no defined result.");
        }

        bool nterm = regexp.Length > 0 && regexp[0] == '^';
        var found = new List<(uint, uint)>();
        int off = 0;
        while (off < text.Length)
        {
            var m = rx.Exec(off == 0 ? text : text[off..]);
            if (m is null)
                break;
            found.Add(((uint)(m.Value.Offset + off), (uint)m.Value.Length));
            off += m.Value.Offset + 1;
            if (nterm)
                break;
        }

        found.Reverse();                                                          // the library pushes every match on the front of a list: the results run from the last match to the first
        return found;
    }

    // the feature a match makes: start = 1 + start, end = 1 + (len + start - 1) (unsigned arithmetic), score = the length
    private static (int Start, int End, int Length) PatMatchCoordinates((uint Start, uint Length) m)
    {
        uint start = 1u + m.Start;
        uint end = unchecked(1u + unchecked(m.Length + m.Start - 1u));
        return (unchecked((int)start), unchecked((int)end), unchecked((int)m.Length));
    }

    // ---------------------------------------------------------------------------------------------------------------- the data files

    // ajDatafileNewInNameS: the work directory, ./.embossdata, the data directory (the home directory of the run has no files)
    private static (string Text, string Path)? DatafileOpen(IReadOnlyDictionary<string, string> files, string name)
    {
        if (files.TryGetValue(name, out string? text))
            return (text, name);
        string dot = ".embossdata/" + name;
        if (files.TryGetValue(dot, out text))
            return (text, dot);
        string data = EmbossRunEnvironment.Current.DataDir + name;
        if (files.TryGetValue(data, out text))
            return (text, data);
        return null;
    }

    // a path of the work directory as the files of a run name it: a leading ./ and doubled slashes dropped
    private static string NormalizePath(string path)
    {
        while (path.StartsWith("./", StringComparison.Ordinal))
            path = path[2..];
        return path.Replace("//", "/", StringComparison.Ordinal);
    }

    // ajFilenameExistsRead of a directory: the path of the directory as the library opens it (the work directory, ./.embossdata or the data directory), null when there is none
    private static string? DatafileDirectory(IReadOnlyDictionary<string, string> files, string name)
    {
        foreach (string path in new[] { name, ".embossdata/" + name, EmbossRunEnvironment.Current.DataDir + name })
        {
            if (files.Keys.Any(k => k.StartsWith(path, StringComparison.Ordinal)))
                return path;
        }

        return null;
    }

    private static List<string> TrimmedLines(string text)
    {
        var lines = ReadLines(text);
        for (int i = 0; i < lines.Count; i++)
        {
            string l = lines[i];
            if (l.EndsWith('\n'))
                l = l[..^1];
            if (l.EndsWith('\r'))
                l = l[..^1];
            lines[i] = l;
        }

        return lines;
    }

    // ajSysFuncStrtok: a token of the text, the delimiters in front of it skipped; the state is kept between the calls
    private sealed class SysStrtok
    {
        private string _s = "";
        private int _p;

        public string? Next(string? src, string delims)
        {
            if (src is not null)
            {
                _s = src;
                _p = 0;
            }

            if (_p >= _s.Length)
                return null;
            while (_p < _s.Length && delims.Contains(_s[_p], StringComparison.Ordinal))
                _p++;
            if (_p >= _s.Length)
                return null;
            int start = _p;
            while (_p < _s.Length && !delims.Contains(_s[_p], StringComparison.Ordinal))
                _p++;
            return _s[start.._p];
        }
    }

    // ---------------------------------------------------------------------------------------------------------------- patmatmotifs

    private static readonly string[] PrunedMotifs = ["myristyl", "asn_glycosylation", "camp_phospho_site", "pkc_phospho_site", "ck2_phospho_site", "tyr_phospho_site"];

    /// <summary>
    /// <c>patmatmotifs</c>: the protein sequence is searched with every pattern of the PROSITE database made by <c>prosextract</c> (<c>PROSITE/prosite.lines</c> in the work directory or in the data directory; the documentation files
    /// <c>PROSITE/PSnnnnn</c> of <c>-full</c> next to it): the hits are reported (<c>-rformat dbmotif</c> by default), the simple patterns of post-translational modifications are skipped with <c>-prune</c> (the default).
    /// </summary>
    internal static EmbossRunResult Patmatmotifs(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("patmatmotifs");
        var all = EmbossProgramBase.Open("patmatmotifs", a, files, ctx, "protein", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var seq = all.Peek();
        bool full, prune;
        EmbossReportArgs.Opened? report;
        try
        {
            full = a.Bool("full", false);
            prune = a.Bool("prune", true);
            report = EmbossReportArgs.Open(a, ctx, "outfile", "patmatmotifs", "dbmotif", [new("str", "motif", "Motif")], 0, "", a.Raw("outfile"),
                EmbossOutPath.DefaultBaseName(seq.Name, a.Raw("sequence") ?? ""), EmbossRunEnvironment.Current.RunDate ?? DateTime.Now, files);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (report is null)
            return EmbossProgramBase.Died(ctx, "outfile");
        string upper = seq.SeqText.ToUpperInvariant();                          // ajSeqFmtUpper: the sequence object itself is in upper case
        AssignUpper(seq, upper);
        var data = DatafileOpen(files, "PROSITE/prosite.lines");
        if (data is null)
            return Fatal(ctx, report, FatalText("patmatmotifs.c", 101, "Either EMBOSS_DATA undefined or PROSEXTRACT needs running"));
        var features = new List<EmbossFeature>();
        var tail = new StringBuilder();
        var tok = new SysStrtok();
        var lines = TrimmedLines(data.Value.Text);
        int pos = 0;
        string name = "", accession = "";
        bool Read(out string line)
        {
            if (pos < lines.Count)
            {
                line = lines[pos++];
                return true;
            }

            line = "";
            return false;
        }

        while (Read(out string regexp))
        {
            if (regexp.Length > 0 && regexp[0] != ' ' && regexp[0] != '^')
            {
                name = tok.Next(regexp, " ") ?? "";
                if (prune && PrunedMotifs.Any(m => string.Equals(m, name, StringComparison.OrdinalIgnoreCase)))
                {
                    for (int i = 0; i < 4; ++i)
                        Read(out _);
                    continue;
                }

                accession = tok.Next(null, " ") ?? "";
            }

            if (regexp.StartsWith('^'))
            {
                var matches = PatMatchFind(regexp[1..], upper);
                foreach (var m in matches)
                {
                    var (start, end, length) = PatMatchCoordinates(m);
                    features.Add(new EmbossFeature("SO:0001067", "patmatmotifs", start, end, length, '\0', 0, [new EmbossFeatureTag("note", "*motif " + name)]));
                }

                if (full && matches.Count > 0)
                {
                    string docName = "PROSITE/" + accession;
                    var doc = DatafileOpen(files, docName);
                    if (doc is null)
                    {
                        // an empty accession names the directory PROSITE/: it opens, and reading it fails
                        if (accession.Length == 0 && DatafileDirectory(files, docName) is { } directory)
                            return Fatal(ctx, report, FatalText("ajfileio.c", 289, $"Error reading from file '{directory}' ferror:1 21: Is a directory\n"));
                        continue;
                    }

                    tail.Append("Motif: ").Append(name).Append('\n');
                    tail.Append("Count: ").Append(matches.Count.ToString(CultureInfo.InvariantCulture)).Append("\n\n");
                    foreach (string text in TrimmedLines(doc.Value.Text))
                        tail.Append(text).Append('\n');
                    tail.Append("\n***************\n\n");
                }
            }
        }

        string header = $"Full: {(full ? "Yes" : "No")}\nPrune: {(prune ? "Yes" : "No")}\nData_file: {data.Value.Path}\n";
        var options = report.Options with { Header = header, Tail = tail.ToString(), WriteFooter = false };
        var sequence = new EmbossReportSequence(seq.Name, EmbossFeatureReports.MakeUsa(seq), seq.Acc, seq.Desc, upper, EmbossSeqRange.GetBegin(seq), EmbossSeqRange.GetEnd(seq), seq.Offset, seq.Reversed, false);
        var table = new EmbossReport(sequence, features, true, "P", TableLength: sequence.Length, Tags: options.Tags);
        string text2 = EmbossFeatureReports.Write(options, [table], ctx, null, out _);
        return EmbossReportArgs.Result(ctx, report.File, text2);
    }

    // the program stops with a message of ajFatal after its report file was opened: an empty file stays when it is not the standard output
    private static EmbossRunResult Fatal(EmbossFeatureReadContext ctx, EmbossReportArgs.Opened report, string message)
    {
        var died = EmbossProgramBase.Fatal(ctx, message);
        return report.File == "stdout" ? died : died with { Files = new Dictionary<string, string> { [report.File] = "" } };
    }

    // ---------------------------------------------------------------------------------------------------------------- patmatdb

    /// <summary><c>patmatdb</c>: the protein sequences are searched with one PROSITE pattern (<c>-motif</c>, case-independent); the hits are reported (<c>-rformat dbmotif</c> by default).</summary>
    internal static EmbossRunResult Patmatdb(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("patmatdb");
        var all = EmbossProgramBase.Open("patmatdb", a, files, ctx, "protein");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        string motif;
        EmbossReportArgs.Opened? report;
        try
        {
            motif = a.Raw("motif") ?? "";
            if (motif.Length < 2)
                throw a.Die("motif", "Error: Too short (motif) - minimum length is 2 characters");
            report = EmbossReportArgs.Open(a, ctx, "outfile", "patmatdb", "dbmotif", [], 0, "", a.Raw("outfile"),
                EmbossOutPath.DefaultBaseName(all.Peek().Name, a.Raw("sequence") ?? ""), EmbossRunEnvironment.Current.RunDate ?? DateTime.Now, files);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (report is null)
            return EmbossProgramBase.Died(ctx, "outfile");
        motif = motif.ToUpperInvariant();
        string regexp;
        try
        {
            regexp = PrositeToRegExp(motif);
        }
        catch (EmbossProgramFatalException e)
        {
            return Fatal(ctx, report, e.Message);
        }

        var options = report.Options with { Header = "Motif: " + motif + "\n" };
        var tables = new List<EmbossReport>();
        long count = 0, total = 0;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            count++;
            total += EmbossSeqRange.GetEnd(seq) - EmbossSeqRange.GetBegin(seq) + 1;
            string str = seq.SeqText.ToUpperInvariant();
            var matches = PatMatchFind(regexp, str);
            if (matches.Count == 0)
                continue;
            var features = new List<EmbossFeature>();
            foreach (var m in matches)
            {
                var (start, end, length) = PatMatchCoordinates(m);
                features.Add(new EmbossFeature("SO:0001067", "patmatdb", start, end, length, '\0', 0, [], Protein: true));
            }

            var sequence = new EmbossReportSequence(seq.Name, EmbossFeatureReports.MakeUsa(seq), seq.Acc, seq.Desc, seq.Seq.ToString(), EmbossSeqRange.GetBegin(seq), EmbossSeqRange.GetEnd(seq), seq.Offset, seq.Reversed, false);
            tables.Add(new EmbossReport(sequence, features, true, "P", TableLength: sequence.Length, Tags: options.Tags));
        }

        string text = EmbossReportArgs.Write(report with { Options = options }, tables, ctx, count, total);
        return EmbossReportArgs.Result(ctx, report.File, text);
    }

    // ---------------------------------------------------------------------------------------------------------------- prosextract

    /// <summary>
    /// <c>prosextract</c>: <c>prosite.dat</c> and <c>prosite.doc</c> of the directory <c>-prositedir</c> become the files <c>PROSITE/prosite.lines</c> (for every PATTERN entry: identifier, accession, the first word of the description, the pattern,
    /// its regular expression) and <c>PROSITE/PSnnnnn</c> (the documentation of every accession) of the data directory.
    /// </summary>
    internal static EmbossRunResult Prosextract(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = new EmbossFeatureReadContext { Program = "prosextract" };
        string dirname;
        try
        {
            dirname = a.Raw("prositedir") ?? "./";
            if (dirname.Length == 0)
            {
                throw a.Die("prositedir", "Error: Directory path is required");
            }

            if (!dirname.EndsWith('/'))
                dirname += "/";
            string bare = dirname.TrimEnd('/');
            bool exists = bare.Length == 0 || bare == "." || EmbossReportArgs.ExistingDirectories(files).Contains(bare);
            if (!exists)
                throw a.Die("prositedir", $"Error: Unable to open directory '{dirname}'");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var written = new Dictionary<string, string>();
        string dataDir = EmbossRunEnvironment.Current.DataDir;
        string linesName = dataDir + "PROSITE/prosite.lines";
        if (!files.TryGetValue(NormalizePath(dirname + "prosite.dat"), out string? dat))
            return EmbossProgramBase.Fatal(ctx, FatalText("prosextract.c", 84, "Cannot open file " + dirname + "prosite.dat"));
        var outf = new StringBuilder();
        var tok = new SysStrtok();
        var datLines = TrimmedLines(dat);
        bool haspattern = false;
        int pos = 0;
        string line = "";
        bool Read()
        {
            if (pos < datLines.Count)
            {
                line = datLines[pos++];
                return true;
            }

            line = "";
            return false;
        }

        // the second token of the line; a line without one gives NULL: ajStrAssignC accepts it, ajStrAppendC (the pattern lines) takes strlen(NULL)
        string Second(string delims, bool append = false)
        {
            tok.Next(line, delims);
            string? t = tok.Next(null, delims);
            if (t is null && append)
                throw new InvalidOperationException("prosextract appends the second token of a PA line to the pattern with ajStrAppendC; a PA line without a second token gives NULL and strlen(NULL) is a segmentation fault: the reference program has no defined result.");
            return t ?? "";
        }

        while (Read())
        {
            if (line.StartsWith("ID", StringComparison.Ordinal))
            {
                if (line.EndsWith("PATTERN.", StringComparison.Ordinal))
                {
                    haspattern = true;
                    outf.Append(Second(" \t;")).Append(' ');
                    continue;
                }

                haspattern = false;
                continue;
            }

            if (!haspattern)
                continue;
            if (line.StartsWith("AC", StringComparison.Ordinal))
            {
                outf.Append(Second(" \t;")).Append("\n ");
                continue;
            }

            if (line.StartsWith("DE", StringComparison.Ordinal))
            {
                outf.Append(Second(" \t.")).Append("\n ");
                continue;
            }

            if (line.StartsWith("PA", StringComparison.Ordinal))
            {
                var pa = new StringBuilder();
                while (line.StartsWith("PA", StringComparison.Ordinal))
                {
                    pa.Append(Second(" \t.", true));
                    Read();
                }

                outf.Append(pa).Append('\n');
                string re;
                try
                {
                    re = PrositeToRegExp(pa.ToString());
                }
                catch (EmbossProgramFatalException e)
                {
                    written[linesName] = outf.ToString();
                    return EmbossProgramBase.Fatal(ctx, e.Message) with { Files = written };
                }

                outf.Append('^').Append(re).Append("\n\n");
                continue;
            }
        }

        written[linesName] = outf.ToString();
        if (!files.TryGetValue(NormalizePath(dirname + "prosite.doc"), out string? doc))
            return EmbossProgramBase.Fatal(ctx, FatalText("prosextract.c", 173, "Cannot open file " + dirname + "prosite.doc")) with { Files = written };
        var docLines = TrimmedLines(doc);
        int dp = 0;
        bool flag = false, isopen = false, goback = false;
        string? current = null;
        int storepos = 0;
        string text = "";
        bool ReadDoc()
        {
            if (dp < docLines.Count)
            {
                text = docLines[dp++];
                return true;
            }

            text = "";
            return false;
        }

        while (ReadDoc())
        {
            if (text.StartsWith("{PS", StringComparison.Ordinal) && isopen && !goback)
                goback = true;
            if (text.StartsWith("{PS", StringComparison.Ordinal) && !isopen)
            {
                storepos = dp;
                var t2 = new SysStrtok();
                string p = t2.Next(text[1..], ";") ?? "";
                current = dataDir + "PROSITE/" + p;
                written[current] = "";
                flag = true;
                isopen = true;
                continue;
            }

            if (text.StartsWith("{BEGIN}", StringComparison.Ordinal) && flag)
            {
                var body = new StringBuilder();
                while (ReadDoc())
                {
                    if (text.StartsWith("{END}", StringComparison.Ordinal))
                        break;
                    body.Append(text).Append('\n');
                }

                if (current is not null)
                    written[current] += body.ToString();                         // the file of the entry is closed after the first {BEGIN}: a later one writes to a closed file (nothing)
                current = null;
                isopen = false;
                if (goback)
                {
                    goback = false;
                    dp = storepos;
                }
            }
        }

        return new EmbossRunResult("", EmbossRunResult.Join(ctx.Messages), 0, written);
    }
}
