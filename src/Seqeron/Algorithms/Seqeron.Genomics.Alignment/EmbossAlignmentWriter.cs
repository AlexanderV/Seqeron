using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The EMBOSS 6.6.0 alignment output formats (<c>-aformat</c>), byte for byte (<c>ajax/core/ajalign.c</c> <c>alignWrite*</c>, <c>ajseqwrite.c</c> <c>seqWrite*</c>,
/// <c>nucleus/embaln.c</c> <c>embAlignReportGlobal</c> / <c>embAlignReportLocal</c>): <c>pair</c> / <c>simple</c> / <c>srspair</c> / <c>srs</c> / <c>markx0..3, 10</c> / <c>match</c> /
/// <c>score</c> / <c>sam</c> / <c>tcoffee</c> for pairwise alignments (the output of needle, water, stretcher, matcher, supermatcher) and <c>fasta</c> / <c>a2m</c> / <c>msf</c> /
/// <c>clustal</c> / <c>mega</c> / <c>meganon</c> / <c>nexus</c> / <c>nexusnon</c> / <c>phylip</c> / <c>phylipnon</c> / <c>selex</c> / <c>treecon</c> for any number of rows. The header (program, run date, command
/// line, matrix, gap penalties, Length / Identity / Similarity / Gaps / Score) and the tail are the ones of <c>ajAlignWriteHeader</c> / <c>ajAlignWriteTail</c>.
/// <para>
/// Use one instance per output file: <see cref="Write(EmbossPairAlignment)"/> once per <c>ajAlignWrite</c> (needle writes one alignment per <c>-bsequence</c> sequence), then <see cref="Close"/>.
/// The pair-markup formats (<c>pair</c>, <c>markx*</c>, ...) of three or more sequences need the EMBOSS <c>cons</c> consensus and are not implemented: they take exactly two sequences,
/// as every EMBOSS program that writes them (needle, water, stretcher, matcher, supermatcher) produces pairs.
/// </para>
/// </summary>
public sealed class EmbossAlignmentWriter
{
    private const string GapChars = "-~.? ";
    private const string GapTest = " .~-";

    private readonly StringBuilder _out = new();
    private readonly EmbossAlignFormat _format;
    private readonly EmbossAlignOptions _o;
    private int _count;
    private int _markCount;
    private bool _seqOnly;
    private bool _closed;

    /// <summary>A writer for one output file.</summary>
    public EmbossAlignmentWriter(EmbossAlignFormat format, EmbossAlignOptions? options = null)
    {
        _format = format;
        _o = options ?? new EmbossAlignOptions();
    }

    /// <summary>Parses an <c>-aformat</c> name (case-insensitive; the EMBOSS names and aliases).</summary>
    public static bool TryParseFormat(string? name, out EmbossAlignFormat format)
    {
        format = default;
        if (name is null)
            return false;
        switch (name.Trim().ToLowerInvariant())
        {
            case "unknown": format = EmbossAlignFormat.Unknown; return true;
            case "fasta": format = EmbossAlignFormat.Fasta; return true;
            case "a2m": format = EmbossAlignFormat.A2m; return true;
            case "msf": format = EmbossAlignFormat.Msf; return true;
            case "clustal": format = EmbossAlignFormat.Clustal; return true;
            case "mega": format = EmbossAlignFormat.Mega; return true;
            case "meganon": format = EmbossAlignFormat.Meganon; return true;
            case "nexus": format = EmbossAlignFormat.Nexus; return true;
            case "paup": format = EmbossAlignFormat.Paup; return true;
            case "nexusnon": format = EmbossAlignFormat.Nexusnon; return true;
            case "paupnon": format = EmbossAlignFormat.Paupnon; return true;
            case "phylip": format = EmbossAlignFormat.Phylip; return true;
            case "phylipnon": format = EmbossAlignFormat.Phylipnon; return true;
            case "selex": format = EmbossAlignFormat.Selex; return true;
            case "treecon": format = EmbossAlignFormat.Treecon; return true;
            case "markx0": format = EmbossAlignFormat.Markx0; return true;
            case "markx1": format = EmbossAlignFormat.Markx1; return true;
            case "markx2": format = EmbossAlignFormat.Markx2; return true;
            case "markx3": format = EmbossAlignFormat.Markx3; return true;
            case "markx10": format = EmbossAlignFormat.Markx10; return true;
            case "match": format = EmbossAlignFormat.Match; return true;
            case "multiple": format = EmbossAlignFormat.Multiple; return true;
            case "pair": format = EmbossAlignFormat.Pair; return true;
            case "simple": format = EmbossAlignFormat.Simple; return true;
            case "sam": format = EmbossAlignFormat.Sam; return true;
            case "score": format = EmbossAlignFormat.Score; return true;
            case "srs": format = EmbossAlignFormat.Srs; return true;
            case "srspair": format = EmbossAlignFormat.SrsPair; return true;
            case "tcoffee": format = EmbossAlignFormat.TCoffee; return true;
            default: return false;
        }
    }

    /// <summary>The EMBOSS name of a format (the <c>Align_format:</c> line).</summary>
    public static string FormatName(EmbossAlignFormat format) => format.ToString().ToLowerInvariant();

    /// <summary>All formats that need exactly two rows (the pair markup, SAM, score, match, T-Coffee ...).</summary>
    public static bool IsPairFormat(EmbossAlignFormat format) => format is >= EmbossAlignFormat.Markx0 and <= EmbossAlignFormat.Markx10
        or EmbossAlignFormat.Match or EmbossAlignFormat.Pair or EmbossAlignFormat.Sam or EmbossAlignFormat.Score or EmbossAlignFormat.SrsPair;

    /// <summary>Writes the alignments (one <c>ajAlignWrite</c>) and the tail, returning the file.</summary>
    public static string Format(EmbossAlignFormat format, IReadOnlyList<EmbossPairAlignment> alignments, EmbossAlignOptions? options = null)
    {
        var w = new EmbossAlignmentWriter(format, options);
        w.Write(alignments);
        return w.Close();
    }

    /// <summary>Writes the rows in a sequence-style format (<c>fasta</c>, <c>msf</c>, <c>clustal</c>, ...) and returns the file.</summary>
    public static string Format(EmbossAlignFormat format, IReadOnlyList<EmbossAlignedSequence> rows, EmbossAlignOptions? options = null)
    {
        var w = new EmbossAlignmentWriter(format, options);
        w.WriteRows(rows);
        return w.Close();
    }

    /// <summary>One <c>ajAlignWrite</c> of a single pairwise alignment.</summary>
    public EmbossAlignmentWriter Write(EmbossPairAlignment alignment) => Write([alignment]);

    /// <summary>One <c>ajAlignWrite</c> of the alignments (EMBOSS writes the alignments of one call in one pass; needle / water call it once per alignment).</summary>
    public EmbossAlignmentWriter Write(IReadOnlyList<EmbossPairAlignment> alignments)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        ThrowIfClosed();
        if (alignments.Count == 0)
            return this;
        var data = alignments.Select(a => new Pair(a, _o)).ToList();
        switch (_format)
        {
            case EmbossAlignFormat.Unknown:
            case EmbossAlignFormat.Multiple:
            case EmbossAlignFormat.Pair:
            case EmbossAlignFormat.Simple:
                WriteSimple(data);
                break;
            case EmbossAlignFormat.Srs:
                WriteSrs(data, false);
                break;
            case EmbossAlignFormat.SrsPair:
                WriteSrs(data, true);
                break;
            case EmbossAlignFormat.Markx0: WriteMarks(data, 0); break;
            case EmbossAlignFormat.Markx1: WriteMarks(data, 1); break;
            case EmbossAlignFormat.Markx2: WriteMarks(data, 2); break;
            case EmbossAlignFormat.Markx3: WriteMarks(data, 3); break;
            case EmbossAlignFormat.Markx10: WriteMarks(data, 10); break;
            case EmbossAlignFormat.Match: WriteMatch(data); break;
            case EmbossAlignFormat.Score: WriteScore(data); break;
            case EmbossAlignFormat.Sam: WriteSam(data); break;
            case EmbossAlignFormat.TCoffee: WriteTCoffee(data); break;
            default:
                // the sequence-style formats write the rows of the first alignment (alignWriteSeqformat(thys, 0, ...))
                _seqOnly = true;
                WriteSeqFormat(data[0].Rows);
                break;
        }

        return this;
    }

    /// <summary>
    /// One <c>ajAlignWrite</c> of an alignment of any number of rows (<c>aligncopy</c>): <c>unknown</c> / <c>multiple</c> / <c>simple</c> (the consensus markup under the rows for three or more), <c>srs</c> and the sequence-style formats.
    /// The pair-markup formats take at most two rows: the program refuses them with the message of the exception.
    /// </summary>
    /// <exception cref="ArgumentException">A pair-markup format of more than two rows.</exception>
    public EmbossAlignmentWriter Write(EmbossMultiAlignment alignment)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        ThrowIfClosed();
        int n = alignment.Rows.Count;
        if (n < 2)
            throw new ArgumentException("An alignment has at least two rows.", nameof(alignment));
        if (n > 2 && IsPairFormat(_format))
            throw new ArgumentException($"Alignment format {FormatName(_format)} specifies at most 2 sequences, alignment has {n}", nameof(alignment));
        var d = new Pair(alignment, _o);
        switch (_format)
        {
            case EmbossAlignFormat.Unknown:
            case EmbossAlignFormat.Multiple:
            case EmbossAlignFormat.Simple:
                WriteSimple([d]);
                break;
            case EmbossAlignFormat.Srs:
                WriteSrs([d], false);
                break;
            case EmbossAlignFormat.TCoffee:
                WriteTCoffee([d]);
                break;
            default:
                if (IsPairFormat(_format))
                    throw new NotSupportedException($"The format '{FormatName(_format)}' of {n} rows is not ported.");
                _seqOnly = true;
                WriteSeqFormat(d.Rows);
                break;
        }

        return this;
    }

    private void WriteMarks(List<Pair> data, int markx)
    {
        foreach (var d in data)
            WriteMark(d, markx);
    }

    /// <summary>Writes rows in a sequence-style format (an <c>ajSeqout</c> set: <c>fasta</c>, <c>a2m</c>, <c>msf</c>, <c>clustal</c>, <c>mega</c>, <c>meganon</c>, <c>nexus</c>, <c>nexusnon</c>, <c>phylip</c>, <c>phylipnon</c>, <c>selex</c>, <c>treecon</c>).</summary>
    public EmbossAlignmentWriter WriteRows(IReadOnlyList<EmbossAlignedSequence> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);
        ThrowIfClosed();
        if (IsPairFormat(_format) || _format is EmbossAlignFormat.Srs or EmbossAlignFormat.Simple or EmbossAlignFormat.Multiple or EmbossAlignFormat.Unknown or EmbossAlignFormat.TCoffee)
            throw new ArgumentException($"Format '{FormatName(_format)}' writes pairwise alignments (use Write).", nameof(rows));
        _seqOnly = true;
        WriteSeqFormat([.. rows.Select(r => Standardize(r))]);
        return this;
    }

    private int _taken, _unnamed;

    /// <summary>The text written since the last call (a program that writes the alignments and another output to one stream needs the order in which they were written).</summary>
    public string Take()
    {
        string text = _out.ToString(_taken, _out.Length - _taken);
        _taken = _out.Length;
        return text;
    }

    /// <summary>Writes the tail (<c>ajAlignClose</c>) and returns the whole file.</summary>
    public string Close()
    {
        if (!_closed)
        {
            _closed = true;
            if (!_seqOnly && HasHeader(_format))
                WriteTail();
        }

        return _out.ToString();
    }

    private void ThrowIfClosed()
    {
        if (_closed)
            throw new InvalidOperationException("The writer is closed.");
    }

    // ---------------------------------------------------------------- data ----

    private sealed class Pair
    {
        public readonly EmbossAlignedSequence[] Seq;
        public readonly string[] Row;
        public readonly int[] Start;
        public readonly int[] End;
        public readonly int[] Len;
        public readonly int[] Offset;
        public readonly int[] Offend;
        public readonly bool[] Rev;
        public readonly string DefaultMatrix = "";
        public readonly string ScoreText;
        public readonly string ExtraSubHeader;
        public readonly int LenAli;
        public readonly int NumId;
        public readonly int NumSim;
        public readonly int NumGap;
        public readonly string Marks;     // | : . and blank, the pair-markup of the simple format

        /// <summary>An alignment of three or more rows (<c>aligncopy</c>): the statistics and the markup come from <c>alignConsStats</c>.</summary>
        public Pair(EmbossMultiAlignment m, EmbossAlignOptions o)
        {
            Seq = [.. m.Rows.Select(Standardize)];
            int n = Seq.Length;
            Row = new string[n];
            Start = new int[n];
            Array.Fill(Start, 1);
            End = new int[n];
            Len = new int[n];
            Offset = new int[n];
            Offend = new int[n];
            Rev = new bool[n];
            for (int i = 0; i < n; i++)
            {
                Row[i] = Seq[i].Row;
                End[i] = Row[i].Length;
                Offset[i] = Seq[i].Offset;
                Offend[i] = Seq[i].Offend;
                Rev[i] = Seq[i].Reversed;
                Len[i] = Row[i].Length + Offset[i] + Offend[i] - Row[i].Count(c => GapChars.Contains(c));
            }

            if (Row.Any(r => r.Length != Row[0].Length))
                throw new ArgumentException("The aligned rows must have equal length.", nameof(m));
            ScoreText = "";
            ExtraSubHeader = "";
            bool nucleotide = IsNuc(Row[0]);
            var matrix = o.Matrix ?? SubstitutionMatrix.Get(nucleotide ? "EDNAFULL" : "EBLOSUM62");
            DefaultMatrix = o.Matrix is null ? MatrixFileName(matrix) : "";
            var st = EmbossAlignmentConsensus.Calculate(Row, [.. Seq.Select(x => (float)x.Weight)], matrix);
            LenAli = Row[0].Length;
            NumId = st.Identical;
            NumSim = st.Similar;
            NumGap = st.Gaps;
            Marks = EmbossAlignmentConsensus.Markup(st.Consensus, '|', ':', '.', ' ');
        }

        public Pair(EmbossPairAlignment a, EmbossAlignOptions o)
        {
            Seq = [Standardize(a.First), Standardize(a.Second)];
            Row = new string[2];
            Start = [1, 1];
            End = new int[2];
            Len = new int[2];
            Offset = new int[2];
            Offend = new int[2];
            Rev = new bool[2];
            for (int i = 0; i < 2; i++)
            {
                Row[i] = Seq[i].Row;
                End[i] = Row[i].Length;
                Offset[i] = Seq[i].Offset;
                Offend[i] = Seq[i].Offend;
                Rev[i] = Seq[i].Reversed;
                Len[i] = Row[i].Length + Offset[i] + Offend[i] - Row[i].Count(c => GapChars.Contains(c));
            }

            if (Row[0].Length != Row[1].Length)
                throw new ArgumentException("The two aligned rows must have equal length.", nameof(a));
            ScoreText = ScoreString(a.Score, o.IntegerScores);
            ExtraSubHeader = a.ExtraSubHeader;
            var matrix = o.Matrix ?? (IsNuc(Row[0]) && IsNuc(Row[1]) ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62);
            var al = new AlignmentResult(Row[0], Row[1], 0, AlignmentType.Global, 0, 0, 0, 0);
            var st = SequenceAligner.CalculateStatistics(al, matrix);
            LenAli = Row[0].Length;
            NumId = st.Matches;
            NumSim = st.Positives;
            NumGap = st.Gaps;
            Marks = SequenceAligner.EmbossColumnMarks(Row[0], Row[1], matrix);
        }

        public EmbossAlignedSequence[] Rows => Seq;
        public int SeqBegin(int i) => Rev[i] ? Len[i] - Offend[i] - Start[i] + 1 : Offset[i] + Start[i];
        public int SeqEnd(int i)
        {
            int iend = End[i] != 0 ? End[i] : Len[i];
            return Rev[i] ? Len[i] - Offend[i] - iend + 1 : Offset[i] + iend;
        }

        public int Increment(int i) => Rev[i] ? -1 : 1;

        public int GapBegin(int i)
        {
            int ret = Rev[i] ? Offend[i] : Offset[i];
            string cp = Row[i];
            int k = Start[i] - 1;
            int n = 0;
            while (k + n < cp.Length && GapChars.Contains(cp[k + n]))
                n++;
            return ret + n;
        }

        public int GapEnd(int i)
        {
            int ret = Rev[i] ? Offset[i] : Offend[i];
            string cp = Row[i];
            int pos = End[i];
            int n = 0;
            for (; pos > 0; pos--)
            {
                // the terminating NUL counts as a gap character (strchr finds it in the test string)
                bool gap = pos >= cp.Length || GapChars.Contains(cp[pos]);
                if (!gap)
                    return ret + n - 1;
                n++;
            }

            return ret + End[i];
        }
    }

    private static EmbossAlignedSequence Standardize(EmbossAlignedSequence s)
    {
        var sb = new StringBuilder(s.Row.Length);
        foreach (char c in s.Row)
            sb.Append(GapChars.Contains(c) ? '-' : c);
        return s with { Row = sb.ToString() };
    }

    private static bool InSet(string row, string set)
    {
        foreach (char c in row)
            if (!set.Contains(char.ToUpperInvariant(c)))
                return false;
        return true;
    }

    private static bool IsNuc(string row) => EmbossSequenceReader.IsNucleotide(row);

    private static bool IsProt(string row) => InSet(row, "ACDEFGHIKLMNPQRSTVWYBJOUXZ?*.~-");

    private static string Left(string s, int width) => (s.Length > width ? s[..width] : s).PadRight(width);

    private static string Right(object v, int width) => Convert.ToString(v, CultureInfo.InvariantCulture)!.PadLeft(width);

    private static string Float(double v, int decimals) => CFormat.Fixed(v, decimals);

    // ajAlignSetScoreR / ajAlignSetGapR: "%.3f" of a float, then up to two trailing zeros removed
    private static string ScoreString(double v, bool integer = false)
    {
        if (integer)
            return ((int)Math.Round(v, MidpointRounding.AwayFromZero)).ToString(CultureInfo.InvariantCulture);
        string t = CFormat.Fixed((float)v, 3);
        for (int i = 1; i < 3; i++)
        {
            if (t[^1] != '0')
                break;
            t = t[..^1];
        }

        return t;
    }

    private static bool HasHeader(EmbossAlignFormat f) => f is EmbossAlignFormat.Unknown or EmbossAlignFormat.Markx0 or EmbossAlignFormat.Markx1 or EmbossAlignFormat.Markx2 or EmbossAlignFormat.Markx3
        or EmbossAlignFormat.Markx10 or EmbossAlignFormat.Match or EmbossAlignFormat.Multiple or EmbossAlignFormat.Pair or EmbossAlignFormat.Simple or EmbossAlignFormat.Score
        or EmbossAlignFormat.Srs or EmbossAlignFormat.SrsPair;

    // ------------------------------------------------------------- dates ----

    // ajTimeRefTodayFmt("report"): "%a %e %b %Y %H:%M:%S"
    private string ReportDate() => EmbossDates.Report(_o.RunDate);

    // %D of ajTimeRefToday: "%d/%m/%y"
    private string ShortDate()
    {
        var d = _o.RunDate;
        return string.Create(CultureInfo.InvariantCulture, $"{d.Day:D2}/{d.Month:D2}/{d.Year % 100:D2}");
    }

    // ------------------------------------------------------------ header ----

    private string SeqName(EmbossAlignedSequence s) => _o.ShowUsa ? s.Usa : s.Name;

    private void WriteHeader(Pair d)
    {
        if (!HasHeader(_format))
            return;
        var o = _out;
        if (_count == 0)
        {
            o.Append("########################################\n");
            o.Append("# Program: ").Append(_o.Program).Append('\n');
            o.Append("# Rundate: ").Append(ReportDate()).Append('\n');
            o.Append("# Commandline: ").Append(_o.Program).Append('\n');
            if (_o.CommandLine.Count > 0)
                o.Append("#    ").Append(string.Join("\n#    ", _o.CommandLine)).Append('\n');
            o.Append("# Align_format: ").Append(FormatName(_format)).Append('\n');
            o.Append("# Report_file: ").Append(_o.ReportFile).Append('\n');
            o.Append("########################################\n");
        }

        o.Append("\n#=======================================\n#\n");
        o.Append("# Aligned_sequences: ").Append(d.Seq.Length).Append('\n');
        for (int i = 0; i < d.Seq.Length; i++)
        {
            var s = d.Seq[i];
            string extra = "";
            if (_o.ShowAccession)
                extra += " (" + s.Accession + ")";
            if (_o.ShowDescription)
                extra += " " + s.Description;
            o.Append("# ").Append(i + 1).Append(": ").Append(SeqName(s)).Append(extra).Append('\n');
        }

        string matrixName = PrintedMatrixName();
        if (matrixName.Length == 0)
            matrixName = d.DefaultMatrix;
        if (matrixName.Length > 0)
            o.Append("# Matrix: ").Append(matrixName).Append('\n');
        if (_o.GapPenalty is { } gp)
            o.Append("# Gap_penalty: ").Append(ScoreString(gp, _o.IntegerScores)).Append('\n');
        if (_o.ExtendPenalty is { } ep)
            o.Append("# Extend_penalty: ").Append(ScoreString(ep, _o.IntegerScores)).Append('\n');
        if (_format is not (EmbossAlignFormat.Match or EmbossAlignFormat.Score))
        {
            string sub = SubHeader(d);
            o.Append("#\n# ").Append(sub.Replace("\n", "\n# ")).Append("\n#\n");
        }
        else if (d.ExtraSubHeader.Length > 0)
        {
            // the plain-rows alignment of match / score has no standard sub-header, only the program's own (needleall -nobrief)
            o.Append("#\n# ").Append(d.ExtraSubHeader.Replace("\n", "\n# ")).Append("\n#\n");
        }

        o.Append("#=======================================\n\n");
        _count++;
    }

    // ajAlignSetSubStandard (+ the program's own sub-header appended before it by ajAlignSetSubHeaderApp: SubHeaderPre puts the standard text in front)
    private static string SubHeader(Pair d)
    {
        var sb = new StringBuilder();
        sb.Append("Length: ").Append(d.LenAli).Append('\n');
        if (d.LenAli > 0)
        {
            sb.Append("Identity:   ").Append(Right(d.NumId, 5)).Append('/').Append(d.LenAli).Append(" (").Append(Right(Pct(d.NumId, d.LenAli), 4)).Append("%)\n");
            sb.Append("Similarity: ").Append(Right(d.NumSim, 5)).Append('/').Append(d.LenAli).Append(" (").Append(Right(Pct(d.NumSim, d.LenAli), 4)).Append("%)\n");
            sb.Append("Gaps:       ").Append(Right(d.NumGap, 5)).Append('/').Append(d.LenAli).Append(" (").Append(Right(Pct(d.NumGap, d.LenAli), 4)).Append("%)\n");
        }

        if (d.ScoreText.Length > 0)
            sb.Append("Score: ").Append(d.ScoreText).Append('\n');
        if (d.ExtraSubHeader.Length > 0)
            sb.Append(d.ExtraSubHeader);
        return sb.ToString();
    }

    // (float) 100. * (float) n / (float) len, printed "%4.1f"
    private static string Pct(int n, int len)
    {
        float pct = 100f * n / len;
        return Float(pct, 1);
    }

    private string PrintedMatrixName()
    {
        if (_o.MatrixName is { } name)
            return name;
        return _o.Matrix is { } m ? MatrixFileName(m) : "";
    }

    private static string MatrixFileName(SubstitutionMatrix m)
    {
        if (SubstitutionMatrix.BuiltInNames.Contains(m.Name))
            return m.Name;
        string embl = "E" + m.Name;
        return SubstitutionMatrix.BuiltInNames.Contains(embl) ? embl : m.Name;
    }

    private void WriteTail()
    {
        _out.Append('\n').Append("#---------------------------------------\n").Append("#---------------------------------------\n");
    }

    // ------------------------------------------------------ simple / srs ----

    private void WriteSimple(List<Pair> data)
    {
        int iwidth = _o.Width != 0 ? _o.Width : 50;
        foreach (var d in data)
        {
            WriteHeader(d);
            string cons = d.Marks;
            int ilen = d.LenAli;
            int nseq = d.Seq.Length;
            int[] incs = [.. Enumerable.Range(0, nseq).Select(d.Increment)];
            int[] ipos = [.. Enumerable.Range(0, nseq).Select(k => d.SeqBegin(k) - incs[k])];
            for (int i = 0; i < ilen; i += iwidth)
            {
                string mrk = "";
                for (int iseq = 0; iseq < nseq; iseq++)
                {
                    int istart = i;
                    int iend = Math.Min(ilen - 1, istart + iwidth - 1);
                    string tmp = d.Row[iseq].Substring(istart, iend - istart + 1).Replace('.', '-');
                    mrk = cons.Substring(istart, iend - istart + 1);
                    int icnt = incs[iseq] * (tmp.Length - tmp.Count(c => c == '-') - tmp.Count(c => c == ' '));
                    if (nseq == 2 && iseq == 1)
                        _out.Append("                     ").Append(mrk).Append('\n');
                    if (tmp.Length > 0)
                        _out.Append(Left(SeqName(d.Seq[iseq]), 13)).Append(' ').Append(Right(ipos[iseq] + incs[iseq], 6)).Append(' ').Append(tmp).Append(' ').Append(Right(ipos[iseq] + icnt, 6)).Append('\n');
                    else
                        _out.Append(Left(SeqName(d.Seq[iseq]), 13)).Append('\n');
                    ipos[iseq] += icnt;
                }

                if (nseq > 2)
                    _out.Append("                     ").Append(mrk).Append('\n');
                _out.Append('\n');
            }
        }
    }

    private void WriteSrs(List<Pair> data, bool pairMode)
    {
        int iwidth = _o.Width != 0 ? _o.Width : 50;
        foreach (var d in data)
        {
            WriteHeader(d);
            string cons = d.Marks;
            int ilen = d.LenAli;
            int nseq = d.Seq.Length;
            int[] incs = [.. Enumerable.Range(0, nseq).Select(d.Increment)];
            int[] ipos = [.. Enumerable.Range(0, nseq).Select(k => d.SeqBegin(k) - incs[k])];
            for (int i = 0; i < ilen; i += iwidth)
            {
                for (int iseq = 0; iseq < nseq; iseq++)
                {
                    int istart = i;
                    int iend = Math.Min(ilen - 1, istart + iwidth - 1);
                    string tmp = d.Row[iseq].Substring(istart, iend - istart + 1).Replace('.', '-');
                    string mrk = cons.Substring(istart, iend - istart + 1);
                    int icnt = incs[iseq] * (tmp.Length - tmp.Count(c => c == '-') - tmp.Count(c => c == ' '));
                    if (pairMode && iseq == 1)
                        _out.Append("                     ").Append(mrk).Append('\n');
                    if (tmp.Length > 0)
                    {
                        if (icnt != 0)
                            _out.Append(Left(SeqName(d.Seq[iseq]), 13)).Append(' ').Append(Right(ipos[iseq] + incs[iseq], 6)).Append(' ').Append(tmp).Append(' ').Append(Right(ipos[iseq] + icnt, 6)).Append('\n');
                        else
                            _out.Append(Left(SeqName(d.Seq[iseq]), 13)).Append(' ').Append(Right(ipos[iseq], 6)).Append(' ').Append(tmp).Append(' ').Append(Right(ipos[iseq], 6)).Append('\n');
                    }
                    else
                    {
                        _out.Append(Left(SeqName(d.Seq[iseq]), 13)).Append('\n');
                    }

                    ipos[iseq] += icnt;
                }

                _out.Append('\n');
            }
        }
    }

    // --------------------------------------------------------------- markx ----

    private static string MarkLine(string cons, char id, char sim, char mis, char gap)
    {
        // alignSim over the | : . markup of the pair: identity -> id, similar -> sim, mismatch -> mis, anything else -> gap
        var sb = new StringBuilder(cons.Length);
        foreach (char c in cons)
            sb.Append(c switch { '|' => id, ':' => sim, '.' => mis, _ => gap });
        return sb.ToString();
    }

    private void WriteMark(Pair d, int markx)
    {
        int iwidth = _o.Width != 0 ? _o.Width : 50;
        _markCount++;
        int nseq = 2;
        int ilen = d.LenAli;
        WriteHeader(d);
        bool isNuc = IsNuc(d.Row[0]);
        float ident = (float)d.NumId / d.LenAli;
        string cons;
        if (markx == 1)
            cons = MarkLine(d.Marks, ' ', ' ', 'X', ' ');
        else if (markx == 2)
        {
            // alignConsSet(second row) then alignDiff against the first row: '.' where the residues agree
            var c = d.Row[1].ToCharArray();
            for (int i = 0; i < Math.Min(c.Length, d.Row[0].Length); i++)
            {
                char a = c[i];
                char b = d.Row[0][i];
                if (a is ' ' or '-' || b is ' ' or '-')
                    continue;
                if (char.ToUpperInvariant(a) == char.ToUpperInvariant(b))
                    c[i] = '.';
            }

            cons = new string(c);
        }
        else
        {
            cons = MarkLine(d.Marks, ':', '.', ' ', ' ');
        }

        int[] incs = [d.Increment(0), d.Increment(1)];
        int[] ipos = [d.SeqBegin(0) - incs[0], d.SeqBegin(1) - incs[1]];
        bool[] num = new bool[2];
        var o = _out;
        if (markx == 10)
        {
            o.Append(">>>").Append(d.Seq[0].Name).Append(", ").Append(d.Len[0]).Append(isNuc ? " nt" : " aa").Append(" vs ").Append(d.Seq[1].Name).Append(", ").Append(d.Len[1]).Append(isNuc ? " nt" : " aa").Append('\n');
            o.Append("; mp_name: EMBOSS\n");
            o.Append("; mp_ver: ").Append(_o.ProgramVersion).Append('\n');
            o.Append("; pg_name: ").Append(_o.Program).Append('\n');
            o.Append("; pg_ver: ").Append(_o.ProgramVersion).Append('\n');
            o.Append("; pg_matrix: ").Append(PrintedMatrixName()).Append('\n');
            o.Append("; pg_gap-pen: -").Append(ScoreString(_o.GapPenalty ?? 0, _o.IntegerScores)).Append(" -").Append(ScoreString(_o.ExtendPenalty ?? 0, _o.IntegerScores)).Append('\n');
            o.Append(">>#").Append(_markCount).Append('\n');
            o.Append("; sw_score: ").Append(d.ScoreText).Append('\n');
            o.Append("; sw_ident: ").Append(Float(ident, 3)).Append('\n');
            o.Append("; sw_overlap: ").Append(d.LenAli).Append('\n');
        }

        bool lastline = false;
        if (markx <= 2)
        {
            for (int i = 0; i < ilen; i += iwidth)
            {
                for (int iseq = 0; iseq < nseq; iseq++)
                {
                    int istart = i;
                    int iend = istart + iwidth - 1;
                    if (iend >= ilen - 1)
                    {
                        iend = ilen - 1;
                        lastline = true;
                    }

                    string tmp = d.Row[iseq].Substring(istart, iend - istart + 1).Replace('.', '-');
                    string mrkcons = cons.Substring(istart, iend - istart + 1);
                    int icnt = tmp.Length - tmp.Count(c => c == '-') - tmp.Count(c => c == ' ');
                    int increment = incs[iseq];
                    if (iseq == 0)
                    {
                        int jpos = ipos[iseq];
                        var tmpnum = new StringBuilder();
                        int jj = 7;
                        foreach (char c in tmp)
                        {
                            jj++;
                            if (c != '-' && c != ' ')
                            {
                                jpos += increment;
                                if (jpos % 10 == 0)
                                {
                                    tmpnum.Append(Right(jpos, jj));
                                    jj = 0;
                                    num[iseq] = true;
                                }
                            }
                        }

                        if (lastline && !num[iseq])
                        {
                            tmpnum.Append(Right(jpos, jj));
                            jj = 0;
                        }

                        if (jj != 0)
                            tmpnum.Append(' ', jj);
                        o.Append(tmpnum).Append('\n');
                    }

                    string name6 = Right(SeqName(d.Seq[iseq]) is { Length: > 6 } nm ? nm[..6] : SeqName(d.Seq[iseq]), 6);
                    if (markx == 2 && nseq == 2 && iseq == 1)
                    {
                        o.Append(name6).Append(' ').Append(mrkcons).Append('\n');
                    }
                    else
                    {
                        if (nseq == 2 && iseq == 1)
                            o.Append("       ").Append(mrkcons).Append('\n');
                        if (tmp.Length > 0)
                            o.Append(name6).Append(' ').Append(tmp).Append('\n');
                        else
                            o.Append(name6).Append('\n');
                        if (iseq + 1 == nseq)
                        {
                            int jpos = ipos[iseq];
                            var tmpnum = new StringBuilder();
                            int jj = 7;
                            int kk = 7;
                            foreach (char c in tmp)
                            {
                                jj++;
                                if (c != '-' && c != ' ')
                                {
                                    kk = jj;
                                    jpos += increment;
                                    if (jpos % 10 == 0)
                                    {
                                        tmpnum.Append(Right(jpos, jj));
                                        jj = 0;
                                        kk = 0;
                                        num[iseq] = true;
                                    }
                                }
                            }

                            if (lastline && !num[iseq])
                            {
                                tmpnum.Append(Right(jpos, kk));
                                jj -= kk;
                            }

                            if (jj != 0)
                                tmpnum.Append(' ', jj);
                            o.Append(tmpnum).Append('\n');
                        }

                        ipos[iseq] += increment > 0 ? icnt : -icnt;
                    }
                }

                o.Append('\n');
            }
        }
        else
        {
            for (int iseq = 0; iseq < nseq; iseq++)
            {
                string seq = d.Row[iseq];
                string tmp = seq.Replace('.', '-');
                int icnt = tmp.Length - tmp.Count(c => c == '-') - tmp.Count(c => c == ' ');
                if (markx == 3)
                {
                    o.Append('>').Append(SeqName(d.Seq[iseq])).Append(" ..\n");
                }
                else
                {
                    o.Append('>').Append(SeqName(d.Seq[iseq])).Append(" ..\n");
                    o.Append("; sq_len: ").Append(d.Len[iseq]).Append('\n');
                    o.Append(IsNuc(d.Row[iseq]) ? "; sq_type: D\n" : "; sq_type: p\n");
                    o.Append("; al_start: ").Append(d.Offset[iseq] + d.Start[iseq]).Append('\n');
                    o.Append("; al_stop: ").Append(d.Offset[iseq] + d.Start[iseq] + icnt - 1).Append('\n');
                    o.Append("; al_display_start: ").Append(d.Offset[iseq] + d.Start[iseq]).Append('\n');
                }

                for (int j = 0; j < ilen; j += iwidth)
                {
                    int jend = Math.Min(ilen - 1, j + iwidth - 1);
                    string part = seq.Substring(j, jend - j + 1);
                    if (part.Length > 0)
                        o.Append(part).Append('\n');
                }
            }
        }
    }

    // --------------------------------------------------- match / score / sam ----

    private void WriteMatch(List<Pair> data)
    {
        WriteHeader(data[0]);
        foreach (var d in data)
        {
            int len0 = d.End[0] - d.Start[0] + 1;
            char fwd0 = d.Rev[0] ? '-' : '+';
            char fwd1 = d.Rev[1] ? '-' : '+';
            int endgaps = d.GapEnd(0) + d.GapEnd(1);
            int begingaps = d.GapBegin(0) + d.GapBegin(1);
            string n1 = Left(SeqName(d.Seq[0]), 15);
            string n2 = Left(SeqName(d.Seq[1]), 15);
            if (_o.Global)
            {
                len0 = len0 - endgaps - begingaps;
                _out.Append(Right(len0, 6)).Append(' ').Append(n1).Append(' ').Append(fwd0).Append(' ')
                    .Append(Right(d.SeqBegin(0) + d.GapBegin(1), 8)).Append("..").Append(Left((d.SeqEnd(0) - endgaps - d.GapBegin(0)).ToString(CultureInfo.InvariantCulture), 8))
                    .Append(' ').Append(n2).Append(' ').Append(fwd1).Append(' ')
                    .Append(Right(d.SeqBegin(1) + d.GapBegin(0), 8)).Append("..").Append(d.SeqEnd(1) - endgaps - d.GapBegin(1)).Append('\n');
            }
            else
            {
                _out.Append(Right(len0, 6)).Append(' ').Append(n1).Append(' ').Append(fwd0).Append(' ')
                    .Append(Right(d.SeqBegin(0), 8)).Append("..").Append(Left(d.SeqEnd(0).ToString(CultureInfo.InvariantCulture), 8))
                    .Append(' ').Append(n2).Append(' ').Append(fwd1).Append(' ')
                    .Append(Right(d.SeqBegin(1), 8)).Append("..").Append(d.SeqEnd(1)).Append('\n');
            }
        }
    }

    private void WriteScore(List<Pair> data)
    {
        foreach (var d in data)
        {
            string n1 = SeqName(d.Seq[0]);
            string n2 = SeqName(d.Seq[1]);
            if (d.ScoreText.Length > 0)
                _out.Append(n1).Append(' ').Append(n2).Append(' ').Append(d.LenAli).Append(" (").Append(d.ScoreText).Append(")\n");
            else
                _out.Append(n1).Append(' ').Append(n2).Append(' ').Append(d.LenAli).Append('\n');
        }
    }

    private uint _samFlag;

    private void WriteSam(List<Pair> data)
    {
        if (_count++ == 0)
        {
            _out.Append("@HD\tVN:1.3\n");
            _out.Append("@PG\tID:").Append(_o.Program).Append("\tVN:").Append(_o.ProgramVersion).Append("\tCL:").Append(string.Join(' ', _o.CommandLine)).Append('\n');
        }

        foreach (var d in data)
        {
            int refpos = d.SeqBegin(0) + 1 - d.Increment(0);
            int iseq = 1;
            if (d.Rev[iseq])
                _samFlag |= 0x0010;
            int qrystart, qryend, refstart;
            if (_o.Global)
            {
                qrystart = d.SeqBegin(iseq) - d.Increment(iseq);
                qryend = d.SeqEnd(iseq);
                refstart = refpos - 1;
            }
            else
            {
                qrystart = 0;
                qryend = d.LenAli;
                refstart = 0;
            }

            var mseq = new StringBuilder();
            var cigar = new StringBuilder();
            int mism = Cigar(d, iseq, 0, qrystart, qryend, refstart, mseq, cigar);
            _out.Append(SeqName(d.Seq[iseq])).Append('\t').Append(_samFlag).Append('\t').Append(SeqName(d.Seq[0])).Append('\t').Append(refpos).Append('\t').Append(255).Append('\t')
                .Append(cigar).Append("\t*\t0\t0\t").Append(mseq).Append("\t*\tAS:i:").Append(d.ScoreText).Append("\tNM:i:").Append(mism).Append('\n');
        }
    }

    // alignCigar (ajalign.c), including its run-length quirks (the I/D counters of the "gap" branches)
    private static int Cigar(Pair d, int qi, int ri, int qstart, int qend, int rstart, StringBuilder qseq, StringBuilder cigar)
    {
        string q = d.Row[qi];
        string r = d.Row[ri];
        int m = 0, ins = 0, del = 0, mism = 0;
        int qseqlen = q.Length;
        int qoffset = d.Offset[qi];
        int qoffend = d.Offend[qi];
        if (d.Rev[qi] && qoffend > 0)
            cigar.Append(qoffend).Append('H');
        else if (qoffset > 0)
            cigar.Append(qoffset).Append('H');
        for (; qstart < qend; qstart++, rstart++)
        {
            char qc = qstart < q.Length ? q[qstart] : '\0';
            char rc = rstart < r.Length ? r[rstart] : '\0';
            if (GapChars.Contains(qc) || qc == '\0')
            {
                del++;
                if (m > 0)
                {
                    cigar.Append(m).Append('M');
                    m = 0;
                }

                if (ins > 0)
                {
                    cigar.Append(del).Append('I');
                    ins = 0;
                }
            }
            else if (GapChars.Contains(rc) || rc == '\0')
            {
                ins++;
                if (m > 0)
                {
                    cigar.Append(m).Append('M');
                    m = 0;
                }

                if (del > 0)
                {
                    cigar.Append(ins).Append('D');
                    del = 0;
                }

                qseq.Append(qc);
            }
            else
            {
                m++;
                qseq.Append(qc);
                if (char.ToUpperInvariant(qc) != char.ToUpperInvariant(rc))
                    mism++;
                if (del > 0)
                {
                    cigar.Append(del).Append('D');
                    del = 0;
                }

                if (ins > 0)
                {
                    cigar.Append(ins).Append('I');
                    ins = 0;
                }
            }
        }

        if (m > 0)
            cigar.Append(m).Append('M');
        else if (del > 0)
            cigar.Append(del).Append('D');
        else if (ins > 0)
            cigar.Append(ins).Append('I');
        if (d.Rev[qi] && qend > qseqlen)
            cigar.Append(qseqlen - qend).Append('H');
        else if (qend < qseqlen)
            cigar.Append(qseqlen - qend).Append('H');
        return mism;
    }

    private void WriteTCoffee(List<Pair> data)
    {
        _seqOnly = true;
        var first = data[0];
        int nseq = first.Row.Length;
        _out.Append(nseq).Append('\n');
        for (int i = 0; i < nseq; i++)
        {
            string ungapped = first.Row[i].Replace(" ", "").Replace("-", "");
            _out.Append(SeqName(first.Seq[i])).Append(' ').Append(ungapped.Length).Append(' ').Append(ungapped).Append('\n');
        }

        foreach (var d in data)
        {
            // alignConsStats is not called by this writer: NumId is what the program stored (stretcher) or 0 (needle, water, aligncopy)
            int pidentity = _o.StatisticsSetByProgram ? (int)(0.5 + (100.0 * ((float)d.NumId / d.LenAli))) : 0;
            // the program prints %S of a string it never set (aligncopy: no score, matrix or gap penalties) as <null>
            string matrix = PrintedMatrixName();
            _out.Append("! score=").Append(d.ScoreText.Length == 0 && nseq > 2 ? "<null>" : d.ScoreText).Append('\n');
            _out.Append("! matrix=").Append(matrix.Length == 0 && nseq > 2 ? "<null>" : matrix).Append('\n');
            _out.Append("! gapopen=").Append(_o.GapPenalty is null && nseq > 2 ? "<null>" : ScoreString(_o.GapPenalty ?? 0, _o.IntegerScores)).Append(" gapext=").Append(_o.ExtendPenalty is null && nseq > 2 ? "<null>" : ScoreString(_o.ExtendPenalty ?? 0, _o.IntegerScores)).Append('\n');
            for (int s1 = 0; s1 < nseq - 1; s1++)
            {
                for (int s2 = s1 + 1; s2 < nseq; s2++)
                {
                    _out.Append('#').Append(s1 + 1).Append(' ').Append(s2 + 1).Append('\n');
                    int n1 = 1, n2 = 1;
                    for (int i = 0; i < d.LenAli; i++)
                    {
                        if (first.Row[s1][i] >= 'A' && first.Row[s2][i] >= 'A')
                            _out.Append(n1).Append(' ').Append(n2).Append(' ').Append(pidentity).Append(" 1 0\n");
                        if (first.Row[s1][i] >= 'A')
                            n1++;
                        if (first.Row[s2][i] >= 'A')
                            n2++;
                    }
                }
            }
        }

        _out.Append("! SEQ_1_TO_N\n");
    }

    // ------------------------------------------------------ seqout formats ----

    private static string GapS(string s) => s.Replace(' ', '-').Replace('.', '-').Replace('~', '-');

    private static string MaskIdent(string s, string first, char mask)
    {
        var c = s.ToCharArray();
        for (int i = 0; i < Math.Min(c.Length, first.Length); i++)
            if (c[i] == first[i])
                c[i] = mask;
        return new string(c);
    }

    private static string FmtBlock(string s, int len)
    {
        var sb = new StringBuilder();
        for (int i = 0; i < s.Length; i++)
        {
            if (i > 0 && i % len == 0)
                sb.Append(' ');
            sb.Append(s[i]);
        }

        return sb.ToString();
    }

    private void WriteSeqFormat(IReadOnlyList<EmbossAlignedSequence> rows)
    {
        if (rows.Count == 0)
            return;
        // the sequence output the alignment is written through (ajSeqoutSetNameDefaultS): a row without a name is called EMBOSS_001, EMBOSS_002 ... in the order of the unnamed rows
        if (rows.Any(r => r.Name.Length == 0))
            rows = [.. rows.Select(r => r.Name.Length == 0 ? r with { Name = "EMBOSS_" + (++_unnamed).ToString("000", CultureInfo.InvariantCulture) } : r)];

        switch (_format)
        {
            case EmbossAlignFormat.Fasta:
            case EmbossAlignFormat.A2m:
                foreach (var r in rows)
                {
                    _out.Append('>').Append(r.Name);
                    if (r.Accession.Length > 0)
                        _out.Append(' ').Append(r.Accession);
                    if (r.Description.Length > 0)
                        _out.Append(' ').Append(r.Description);
                    _out.Append('\n');
                    for (int i = 0; i < r.Row.Length; i += 60)
                        _out.Append(r.Row, i, Math.Min(60, r.Row.Length - i)).Append('\n');
                }

                break;
            case EmbossAlignFormat.Msf: WriteMsf(rows); break;
            case EmbossAlignFormat.Clustal: WriteClustal(rows); break;
            case EmbossAlignFormat.Mega: WriteMega(rows, true); break;
            case EmbossAlignFormat.Meganon: WriteMega(rows, false); break;
            case EmbossAlignFormat.Nexus:
            case EmbossAlignFormat.Paup: WriteNexus(rows, true); break;
            case EmbossAlignFormat.Nexusnon:
            case EmbossAlignFormat.Paupnon: WriteNexus(rows, false); break;
            case EmbossAlignFormat.Phylip: WritePhylip(rows, true); break;
            case EmbossAlignFormat.Phylipnon: WritePhylip(rows, false); break;
            case EmbossAlignFormat.Selex: WriteSelex(rows); break;
            case EmbossAlignFormat.Treecon:
                _out.Append(rows.Max(r => r.Row.Length)).Append('\n');
                foreach (var r in rows)
                    _out.Append(r.Name).Append('\n').Append(GapS(r.Row)).Append('\n');
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(rows), _format, "Not a sequence-style format.");
        }
    }

    private static int GcgCheck(string seq)
    {
        long check = 0;
        long count = 0;
        foreach (char c in seq)
        {
            count++;
            check += count * char.ToUpperInvariant(c);
            if (count == 57)
                count = 0;
        }

        return (int)(check % 10000);
    }

    // seqGapSL(seq, '.', '~', ilen): every gap character -> '.', the leading and trailing run -> '~', padded to ilen with '~'
    private static string GapLenMsf(string s, int ilen)
    {
        var c = s.Replace(' ', '.').Replace('-', '.').ToCharArray();
        if (c.Length > 0 && GapTest.Contains(c[0]))
            for (int i = 0; i < c.Length && GapTest.Contains(c[i]); i++)
                c[i] = '~';
        if (c.Length > 0 && GapTest.Contains(c[^1]))
            for (int i = c.Length - 1; i > 0 && GapTest.Contains(c[i]); i--)
                c[i] = '~';
        string r = new(c);
        return r.Length < ilen ? r + new string('~', ilen - r.Length) : r;
    }

    private void WriteMsf(IReadOnlyList<EmbossAlignedSequence> rows)
    {
        int ilen = rows.Max(r => r.Row.Length);
        int maxnamelen = Math.Max(10, rows.Max(r => r.Name.Length));
        var seqs = rows.Select(r => GapLenMsf(r.Row, ilen)).ToArray();
        int checktot = 0;
        foreach (var s in seqs)
        {
            checktot += GcgCheck(s);
            checktot %= 10000;
        }

        bool prot = _o.SeqoutType is { } t ? t == 'P' : !IsNuc(seqs[0]) && IsProt(seqs[0]);
        var o = _out;
        if (prot)
        {
            o.Append("!!AA_MULTIPLE_ALIGNMENT 1.0\n\n");
            o.Append("  ").Append(_o.ReportFile).Append(" MSF:  ").Append(ilen).Append(" Type: P ").Append(ShortDate()).Append(" CompCheck: ").Append(Right(checktot, 4)).Append(" ..\n\n");
        }
        else
        {
            o.Append("!!NA_MULTIPLE_ALIGNMENT 1.0\n\n");
            o.Append("  ").Append(_o.ReportFile).Append(" MSF: ").Append(ilen).Append(" Type: N ").Append(ShortDate()).Append(" CompCheck: ").Append(Right(checktot, 4)).Append(" ..\n\n");
        }

        for (int i = 0; i < rows.Count; i++)
            o.Append("  Name: ").Append(rows[i].Name.PadRight(maxnamelen)).Append(" Len: ").Append(seqs[i].Length).Append("  Check: ").Append(Right(GcgCheck(seqs[i]), 4)).Append(" Weight: ").Append(Float(rows[i].Weight, 2)).Append('\n');
        o.Append("\n//\n\n");
        for (int ipos = 1; ipos <= ilen; ipos += 50)
        {
            int iend = Math.Min(ilen, ipos + 49);
            string sqbeg = ipos.ToString(CultureInfo.InvariantCulture);
            string sqend = iend.ToString(CultureInfo.InvariantCulture);
            if (iend == ilen)
            {
                uint igap = (uint)(iend - ipos - sqbeg.Length);
                if (igap >= (uint)sqend.Length)
                {
                    int w = (int)igap;
                    o.Append(new string(' ', maxnamelen - 1)).Append(' ').Append(' ').Append(sqbeg).Append(' ').Append(w >= 0 ? sqend.PadLeft(w) : sqend.PadRight(-w)).Append('\n');
                }
                else
                {
                    o.Append("           ").Append(sqbeg).Append('\n');
                }
            }
            else
            {
                o.Append("           ").Append(sqbeg.PadRight(25)).Append(sqend.PadLeft(25)).Append('\n');
            }

            for (int i = 0; i < rows.Count; i++)
                o.Append(rows[i].Name.PadRight(maxnamelen)).Append(' ').Append(seqs[i], ipos - 1, iend - ipos + 1).Append('\n');
            o.Append('\n');
        }
    }

    private void WriteClustal(IReadOnlyList<EmbossAlignedSequence> rows)
    {
        int ilen = rows.Max(r => r.Row.Length);
        var seqs = rows.Select(r => GapS(r.Row).PadRight(ilen, '-')).ToArray();
        _out.Append("CLUSTAL W (1.83) multiple sequence alignment\n");
        _out.Append("\n\n");
        int iwidth = 60;
        for (int ipos = 1; ipos <= ilen; ipos += 60)
        {
            int iend = ipos + 60 - 1;
            if (iend > ilen)
            {
                iend = ilen;
                iwidth = ilen - ipos + 1;
            }

            for (int i = 0; i < rows.Count; i++)
                _out.Append(Left(rows[i].Name, 15)).Append(' ').Append(seqs[i], ipos - 1, iend - ipos + 1).Append('\n');
            _out.Append(new string(' ', 15)).Append(' ').Append(new string(' ', iwidth)).Append('\n');
            if (iend < ilen)
                _out.Append('\n');
        }
    }

    private void WriteMega(IReadOnlyList<EmbossAlignedSequence> rows, bool interleaved)
    {
        int isize = rows.Count;
        int ilen = 0;
        bool onedesc = true;
        int namewidth = rows[0].Name.Length;
        if (isize == 1)
        {
            onedesc = false;
        }
        else
        {
            ilen = rows[0].Row.Length;
            if (rows[0].Description.Length == 0)
                onedesc = false;
            for (int i = 1; i < isize; i++)
            {
                namewidth = Math.Max(namewidth, rows[i].Name.Length);
                ilen = Math.Max(ilen, rows[i].Row.Length);
                if (onedesc && rows[i].Description != rows[0].Description)
                    onedesc = false;
            }
        }

        if (namewidth > 40)
            namewidth = 40;
        if (!interleaved)
        {
            for (int i = 0; i < isize; i++)
            {
                ilen = Math.Max(ilen, rows[i].Row.Length);
                if (rows[i].Description != rows[0].Description)
                    onedesc = false;
            }
        }

        // ajSeqIsProt(first sequence): a typed sequence (seqret types every sequence by its own residues) is protein unless it is a nucleotide sequence; an untyped one only needs protein-compatible symbols
        bool prot = _o.SeqoutType.HasValue ? !IsNuc(rows[0].Row) && IsProt(rows[0].Row) : IsProt(rows[0].Row);
        var o = _out;
        o.Append("#mega\n");
        o.Append(onedesc ? "!Title: " + rows[0].Description + ";\n" : "!Title: Written by EMBOSS " + ShortDate() + ";\n");
        o.Append("!Format\n");
        o.Append(prot ? "    DataType=Protein" : "    DataType=Nucleotide").Append(interleaved ? " DataFormat=Interleaved\n" : "\n");
        o.Append("    Identical=. Indel=- Missing=?\n");
        if (!prot)
            o.Append("    CodeTable=Standard\n");
        o.Append("    ;\n\n");
        string first = GapS(rows[0].Row);
        if (interleaved)
        {
            o.Append('\n').Append('\n');
            for (int ipos = 1; ipos <= ilen; ipos += 50)
            {
                int iend = Math.Min(ilen, ipos + 49);
                string firstPart = Slice(first, ipos - 1, iend - 1);
                o.Append('\n');
                for (int i = 0; i < isize; i++)
                {
                    string sseq = GapS(Slice(rows[i].Row, ipos - 1, iend - 1));
                    if (i > 0)
                        sseq = MaskIdent(sseq, firstPart, '.');
                    if (!onedesc && ipos == 1 && rows[i].Description.Length > 0)
                        o.Append('[').Append(rows[i].Description).Append("]\n");
                    o.Append('#').Append(Left(rows[i].Name, namewidth)).Append(' ').Append(sseq).Append('\n');
                }
            }

            o.Append('\n');
        }
        else
        {
            for (int i = 0; i < isize; i++)
            {
                string sseq = GapS(rows[i].Row);
                if (i > 0)
                    sseq = MaskIdent(sseq, first, '.');
                if (!onedesc && rows[i].Description.Length > 0)
                    o.Append('[').Append(rows[i].Description).Append("]\n");
                o.Append('#').Append(Left(rows[i].Name, namewidth)).Append('\n').Append(sseq).Append('\n');
            }
        }
    }

    // ajStrAssignSubS clamps the range into the string: a start behind the end yields the last character
    private static string Slice(string s, int from, int to)
    {
        if (from >= s.Length)
            return s.Length > 0 ? s[^1..] : "";
        return s.Substring(from, Math.Min(to, s.Length - 1) - from + 1);
    }

    private void WriteNexus(IReadOnlyList<EmbossAlignedSequence> rows, bool interleaved)
    {
        int isize = rows.Count;
        int ilen = rows.Max(r => r.Row.Length);
        var o = _out;
        o.Append("#NEXUS\n");
        o.Append("[TITLE: Written by EMBOSS ").Append(ShortDate()).Append("]\n\n");
        o.Append("begin data;\n");
        o.Append("dimensions ntax=").Append(isize).Append(" nchar=").Append(ilen).Append(";\n");
        string datatype = _o.SeqoutType == 'P' ? "datatype=protein missing=X" : "datatype=DNA missing=N";
        o.Append(interleaved ? "format interleave " : "format ").Append(datatype).Append(" gap=-;\n");
        o.Append('\n');
        o.Append("matrix\n");
        if (interleaved)
        {
            var seqs = rows.Select(r => r.Row.PadRight(ilen, '-')).ToArray();
            for (int ipos = 1; ipos <= ilen; ipos += 50)
            {
                int iend = Math.Min(ilen, ipos + 49);
                if (ipos > 1)
                    o.Append('\n');
                for (int i = 0; i < isize; i++)
                    o.Append(Left(rows[i].Name, 20)).Append(' ').Append(GapS(seqs[i].Substring(ipos - 1, iend - ipos + 1))).Append('\n');
            }
        }
        else
        {
            for (int i = 0; i < isize; i++)
                o.Append(rows[i].Name).Append('\n').Append(GapS(rows[i].Row)).Append('\n');
        }

        o.Append(";\n\n");
        o.Append("end;\n");
        o.Append("begin assumptions;\n");
        o.Append("options deftype=unord;\n");
        o.Append("end;\n");
    }

    private void WritePhylip(IReadOnlyList<EmbossAlignedSequence> rows, bool interleaved)
    {
        int isize = rows.Count;
        int ilen = rows.Max(r => r.Row.Length);
        var o = _out;
        o.Append(interleaved ? " " : "").Append(isize).Append(' ').Append(ilen).Append('\n');
        var seqs = rows.Select(r => r.Row.PadRight(ilen, '-')).ToArray();
        if (interleaved)
        {
            for (int ipos = 1; ipos <= ilen; ipos += 50)
            {
                int iend = Math.Min(ilen, ipos + 49);
                for (int i = 0; i < isize; i++)
                {
                    string sseq = FmtBlock(GapS(seqs[i].Substring(ipos - 1, iend - ipos + 1)), 10);
                    if (ipos == 1)
                        o.Append(Left(rows[i].Name, 10)).Append(sseq).Append('\n');
                    else
                        o.Append(new string(' ', 10)).Append(sseq).Append('\n');
                }

                if (iend < ilen)
                    o.Append('\n');
            }
        }
        else
        {
            int iend = 0;
            for (int n = 0; n < isize; n++)
            {
                for (int ipos = 1; ipos <= ilen; ipos += 50)
                {
                    iend = Math.Min(ilen, ipos + 49);
                    string sseq = FmtBlock(GapS(seqs[n].Substring(ipos - 1, iend - ipos + 1)), 10);
                    if (ipos == 1)
                        o.Append(Left(rows[n].Name, 10)).Append(sseq).Append('\n');
                    else
                        o.Append(new string(' ', 10)).Append(sseq).Append('\n');
                }

                if (iend < ilen)
                    o.Append('\n');
            }
        }
    }

    private void WriteSelex(IReadOnlyList<EmbossAlignedSequence> rows)
    {
        int n = rows.Count;
        int nlen = rows.Max(r => r.Name.Length);
        int slen = rows.Max(r => r.Row.Length);
        int len = rows[n - 1].Row.Length;     // EMBOSS reuses 'len' (the length of the LAST sequence) as the line-end test
        var names = rows.Select(r => r.Name.PadRight(nlen)).ToArray();
        var aseqs = rows.Select(r => r.Row.PadRight(slen)).ToArray();
        for (int i = 0; i < n; i++)
        {
            _out.Append("#=SQ ").Append(names[i]).Append(' ').Append(Float(rows[i].Weight, 2)).Append(" - - 0..0:0 ");
            _out.Append(rows[i].Description.Length > 0 ? rows[i].Description + "\n" : "-\n");
        }

        _out.Append('\n');
        for (int i = 0; i < slen; i += 50)
        {
            for (int j = 0; j < n; j++)
            {
                string p = aseqs[j][i..];
                if (i + 50 >= len)
                    _out.Append(names[j]).Append(' ').Append(p).Append('\n');
                else
                    _out.Append(names[j]).Append(' ').Append(Left(p, 50)).Append('\n');
            }

            if (i + 50 < len)
                _out.Append('\n');
        }
    }
}
