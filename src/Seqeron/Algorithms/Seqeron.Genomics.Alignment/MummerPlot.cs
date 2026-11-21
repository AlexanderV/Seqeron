using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;
using static Seqeron.Genomics.Alignment.MummerPerl;

namespace Seqeron.Genomics.Alignment;

/// <summary>The switches of <c>mummerplot</c> (MUMmer 3.23).</summary>
public sealed record MummerPlotOptions
{
    /// <summary><c>-b|breaklen</c>: highlight alignments with breakpoints further than this many nucleotides from the nearest sequence end.</summary>
    public int? BreakLength { get; init; }

    /// <summary><c>--[no]color</c>: colour plot lines by percent similarity (true) or turn all colour off (false).</summary>
    public bool? Color { get; init; }

    /// <summary><c>-c|--[no]coverage</c>: a reference coverage plot (default for tiling input).</summary>
    public bool? Coverage { get; init; }

    /// <summary><c>-f|--[no]filter</c>: only the "best" alignments of a delta file (<c>delta-filter -r -q</c>).</summary>
    public bool? Filter { get; init; }

    /// <summary><c>-l|--[no]layout</c>: lay a delta multiplot out in an intelligible fashion (needs <see cref="RFile"/> / <see cref="QFile"/>).</summary>
    public bool? Layout { get; init; }

    /// <summary><c>-p|prefix</c>: the prefix of the output files.</summary>
    public string Prefix { get; init; } = "out";

    /// <summary><c>-rv</c>: reverse video for x11 plots.</summary>
    public bool ReverseVideo { get; init; }

    /// <summary><c>-r|IdR</c>: plot one reference sequence id.</summary>
    public string? IdR { get; init; }

    /// <summary><c>-q|IdQ</c>: plot one query sequence id.</summary>
    public string? IdQ { get; init; }

    /// <summary><c>-R|Rfile</c>: an ordered set of reference sequences (a multi-FASTA file or a list of ids, lengths and directions); the name is the key of <c>files</c> of <see cref="MummerPlot.Run"/>.</summary>
    public string? RFile { get; init; }

    /// <summary><c>-Q|Qfile</c>: the same for the query sequences.</summary>
    public string? QFile { get; init; }

    /// <summary><c>-rport</c>: the port to send the reference id and position to on a mouse double click in the x11 plot window.</summary>
    public int? RPort { get; init; }

    /// <summary><c>-qport</c>: the same for the query.</summary>
    public int? QPort { get; init; }

    /// <summary><c>-s|size</c>: <c>small</c> (default), <c>medium</c> or <c>large</c>.</summary>
    public string Size { get; init; } = "small";

    /// <summary><c>--small</c> / <c>--medium</c> / <c>--large</c> (these win over <see cref="Size"/>; several of them warn).</summary>
    public bool Small { get; init; }

    /// <summary>See <see cref="Small"/>.</summary>
    public bool Medium { get; init; }

    /// <summary>See <see cref="Small"/>.</summary>
    public bool Large { get; init; }

    /// <summary><c>-S|--SNP</c>: highlight the SNP locations (<c>show-snps -H -T -l</c>, needs the sequence files named in the delta header).</summary>
    public bool? Snp { get; init; }

    /// <summary><c>-t|terminal</c>: <c>x11</c> (default), <c>postscript</c> or <c>png</c>.</summary>
    public string Terminal { get; init; } = "x11";

    /// <summary><c>--x11</c> / <c>--postscript</c> / <c>--png</c> (these win over <see cref="Terminal"/>; several of them warn).</summary>
    public bool X11 { get; init; }

    /// <summary>See <see cref="X11"/>.</summary>
    public bool Postscript { get; init; }

    /// <summary>See <see cref="X11"/>.</summary>
    public bool Png { get; init; }

    /// <summary><c>-title</c>: the gnuplot plot title.</summary>
    public string? Title { get; init; }

    /// <summary><c>-x|xrange</c>: the x range <c>[min:max]</c> (a comma is accepted for the colon).</summary>
    public string? XRange { get; init; }

    /// <summary><c>-y|yrange</c>: the y range.</summary>
    public string? YRange { get; init; }

    /// <summary><c>--fat</c>: lay the sequences out by their fattest alignment only (implies <see cref="Layout"/>).</summary>
    public bool Fat { get; init; }

    /// <summary>The file name of the match file (printed in the messages).</summary>
    public string MatchName { get; init; } = "out.delta";

    /// <summary>The status of <c>system ("gnuplot --version")</c>: 0 = gnuplot 4 or newer is installed (the script checks nothing else), -1 = gnuplot is not found, any other value = an outdated gnuplot.</summary>
    public int GnuplotStatus { get; init; }

    /// <summary>What <c>gnuplot --version</c> writes to standard output when it runs.</summary>
    public string GnuplotVersion { get; init; } = "gnuplot 6.0 patchlevel 0\n";
}

/// <summary>A <c>die</c> of a ported perl script: its message (already ending with a newline) goes to standard error and the script ends with a non-zero status.</summary>
/// <param name="message">The message.</param>
public sealed class MummerScriptDieException(string message) : Exception(message);

/// <summary>Output of <see cref="MummerPlot.Run"/>.</summary>
/// <param name="Files">The files the script writes by name: <c>.fplot</c>, <c>.rplot</c>, <c>.hplot</c>, <c>.gp</c> and the <c>.filter</c> delta file.</param>
/// <param name="Stdout">Standard output (the version line of gnuplot).</param>
/// <param name="Stderr">Standard error (the progress lines and warnings).</param>
/// <param name="ExitCode">The exit status (0, or the status of <c>die</c>: 2).</param>
/// <param name="GnuplotCommand">The command line the script runs to render the plot (null when it does not run gnuplot).</param>
public sealed record MummerPlotResult(IReadOnlyDictionary<string, string> Files, string Stdout, string Stderr, int ExitCode, string? GnuplotCommand);

/// <summary>
/// <c>mummerplot</c> of MUMmer 3.23 (perl script, Adam M Phillippy, 3.5): the gnuplot plot of the alignments of a <c>.delta</c> (nucmer, promer), <c>.cluster</c>, <c>.out</c> (mummer -c) or <c>.tiling</c> (show-tiling) file: the data files
/// <c>.fplot</c> / <c>.rplot</c> / <c>.hplot</c> (forward, reverse and highlighted hits sorted by similarity), the gnuplot script <c>.gp</c> and, with <c>-f</c> / <c>-l</c>, the filtered delta file; the script then runs gnuplot (the command line is
/// returned: rendering is not part of the script's text output).
/// </summary>
public static class MummerPlot
{
    /// <summary>Runs <c>mummerplot [options] match-file</c>.</summary>
    /// <param name="matchText">The text of the match file.</param>
    /// <param name="files">The files the script reads besides the match file, by name: the sequence files named in the first line of a delta file (they define the sequences, offsets and lengths of the plot), the <c>-R</c> / <c>-Q</c> files.</param>
    /// <param name="options">The switches (null = defaults).</param>
    public static MummerPlotResult Run(string matchText, IReadOnlyDictionary<string, string>? files = null, MummerPlotOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(matchText);
        return new Engine(options ?? new MummerPlotOptions(), matchText, files ?? new Dictionary<string, string>()).Execute();
    }

    private sealed class Align(string sR, string eR, string sQ, string eQ, double sim, string simText, string lenR, string lenQ, string idR, string idQ)
    {
        public string SR = sR, ER = eR, SQ = sQ, EQ = eQ, SimText = simText, LenR = lenR, LenQ = lenQ, IdR = idR, IdQ = idQ;
        public double Sim = sim;
    }

    /// <summary>An entry of the reference / query hash: offset, length, direction.</summary>
    private sealed class Seq(double off, double len, double dir)
    {
        public double? Off = off;
        public double Len = len, Dir = dir;
    }

    private sealed class Engine
    {
        private const string X11 = "x11", Ps = "postscript", Png = "png";
        private readonly MummerPlotOptions _o;
        private readonly string _match;
        private readonly IReadOnlyDictionary<string, string> _fs;
        private readonly StringBuilder _stdout = new(), _stderr = new();
        private readonly Dictionary<string, string> _files = [];

        // OPT_* variables
        private int? _breaklen;
        private bool? _color, _coverage, _filter, _layout, _snp;
        private string _terminal, _size;
        private string? _idRfile, _idQfile, _xrange, _yrange;
        private int _gpstatus;
        private string? _hfile, _dfile, _pfile;
        private string _ffile = string.Empty, _rfile = string.Empty, _gfile = string.Empty;
        private readonly bool _onlyFattest;
        private string _mfile, _current = string.Empty;
        private int? _rport, _qport;

        public Engine(MummerPlotOptions o, string match, IReadOnlyDictionary<string, string> files)
        {
            _o = o;
            _match = match;
            _fs = files;
            _terminal = o.Terminal;
            _size = o.Size;
            _mfile = o.MatchName;
            _breaklen = o.BreakLength;
            _color = o.Color;
            _coverage = o.Coverage;
            _filter = o.Filter;
            _layout = o.Layout;
            _snp = o.Snp;
            _idRfile = o.RFile;
            _idQfile = o.QFile;
            _xrange = o.XRange;
            _yrange = o.YRange;
            _rport = o.RPort;
            _qport = o.QPort;
            _onlyFattest = o.Fat;
        }

        public MummerPlotResult Execute()
        {
            int exit = 0;
            string? cmd = null;
            try
            {
                cmd = Main();
            }
            catch (MummerScriptDieException e)
            {
                _stderr.Append(e.Message);
                exit = 2;
            }

            return new MummerPlotResult(_files, _stdout.ToString(), _stderr.ToString(), exit, cmd);
        }

        private void Warn(string s) => _stderr.Append(s);

        private static MummerScriptDieException Die(string s) => new(s);

        private string? Main()
        {
            ParseOptions();
            string parseFunc = GetParseFunc();
            if (parseFunc != "delta" && ((_filter ?? false) || (_layout ?? false) || (_snp ?? false)))
            {
                Warn("WARNING: -f -l -S only work with delta input\n");
                _filter = null;
                _layout = null;
                _snp = null;
            }

            var refs = new Dictionary<string, Seq>();
            var qrys = new Dictionary<string, Seq>();
            var refOrder = new List<string>();
            var qryOrder = new List<string>();
            if (_o.IdR is not null)
            {
                refs[_o.IdR] = new Seq(0, 0, 1);
                refOrder.Add(_o.IdR);
            }
            else if (_idRfile is not null)
            {
                ParseIds(_idRfile, refs, refOrder);
            }

            if (_o.IdQ is not null)
            {
                qrys[_o.IdQ] = new Seq(0, 0, 1);
                qryOrder.Add(_o.IdQ);
            }
            else if (_idQfile is not null)
            {
                ParseIds(_idQfile, qrys, qryOrder);
            }

            string matchSource = _match;
            _current = _match;
            if ((_filter ?? false) || (_layout ?? false))
            {
                Warn("Writing filtered delta file " + _dfile + "\n");
                string filtered = MummerDeltaFilter.Run(_match, new MummerDeltaFilterOptions { ReferenceLis = true, QueryLis = true });
                _files[_dfile!] = filtered;
                if (_filter ?? false)
                {
                    matchSource = filtered;
                    _current = filtered;
                    _mfile = _dfile!;
                }
            }

            var aligns = new List<Align>();
            switch (parseFunc)
            {
                case "delta": ParseDelta(matchSource, aligns); break;
                case "cluster": ParseCluster(matchSource, aligns); break;
                case "mummer": ParseMummer(matchSource, aligns); break;
                default: ParseTiling(matchSource, aligns); break;
            }

            if (_layout ?? false)
            {
                if (refs.Count > 0 || qrys.Count > 0)
                {
                    LayoutIds(refs, refOrder, qrys, qryOrder);
                }
                else
                {
                    Warn("WARNING: --layout option only works with -R or -Q\n");
                    _layout = null;
                }
            }

            PlotData(aligns, refs, refOrder, qrys, qryOrder);
            WriteGp(refs, qrys);
            if (_gpstatus == -1)
                return null;
            bool listener = _gpstatus == 0 && _terminal == X11;
            if (listener)
                Warn("Forking mouse listener\n");
            string command = RunGp();
            if (listener)
            {
                // the forked mouse listener (ListenGP) connects its sockets first: nothing listens on the ports
                if (_rport is int rp && rp != 0)
                    Warn("WARNING: Could not connect to rport " + rp.ToString(CultureInfo.InvariantCulture) + "\n");
                if (_qport is int qp && qp != 0)
                    Warn("WARNING: Could not connect to qport " + qp.ToString(CultureInfo.InvariantCulture) + "\n");
            }

            return command;
        }

        // ------------------------------------------------------------------ options

        private void ParseOptions()
        {
            int cnt = 0;
            if (_o.Png)
            {
                _terminal = Png;
                cnt++;
            }

            if (_o.Postscript)
            {
                _terminal = Ps;
                cnt++;
            }

            if (_o.X11)
            {
                _terminal = X11;
                cnt++;
            }

            if (cnt > 1)
                Warn("WARNING: Multiple terminals not allowed, using '" + _terminal + "'\n");
            cnt = 0;
            if (_o.Large)
            {
                _size = "large";
                cnt++;
            }

            if (_o.Medium)
            {
                _size = "medium";
                cnt++;
            }

            if (_o.Small)
            {
                _size = "small";
                cnt++;
            }

            if (cnt > 1)
                Warn("WARNING: Multiple sizes now allowed, using '" + _size + "'\n");

            _gpstatus = _o.GnuplotStatus;
            if (_gpstatus != -1)
                _stdout.Append(_o.GnuplotVersion);
            if (_gpstatus == -1)
            {
                Warn("WARNING: Could not find gnuplot, plot will not be rendered\n");
            }
            else if (_gpstatus != 0)
            {
                Warn("WARNING: Using outdated gnuplot, use v4.0 for best results\n");
                if (_color ?? false)
                {
                    Warn("WARNING: Turning of --color option for compatibility\n");
                    _color = null;
                }

                if (_terminal == Png && _size != "small")
                {
                    Warn("WARNING: Turning of --size option for compatibility\n");
                    _size = "small";
                }
            }

            if (_terminal is not (X11 or Ps or Png))
                throw Die("ERROR: Invalid terminal type, " + _terminal + "\n");
            if (_size is not ("small" or "medium" or "large"))
                throw Die("ERROR: Invalid terminal size, " + _size + "\n");
            if (True(_xrange))
            {
                _xrange = _xrange!.Replace(',', ':');
                if (!Regex.IsMatch(_xrange, @"^\[[0-9]+:[0-9]+\]$"))
                    throw Die("ERROR: Invalid xrange format, " + _xrange + "\n");
            }

            if (True(_yrange))
            {
                _yrange = _yrange!.Replace(',', ':');
                if (!Regex.IsMatch(_yrange, @"^\[[0-9]+:[0-9]+\]$"))
                    throw Die("ERROR: Invalid yrange format, " + _yrange + "\n");
            }

            string prefix = _o.Prefix;
            _ffile = prefix + ".fplot";
            _rfile = prefix + ".rplot";
            if (_breaklen is not null || _snp is not null)
                _hfile = prefix + ".hplot";
            if (_onlyFattest)
                _layout = true;
            if ((_filter ?? false) || (_layout ?? false))
                _dfile = prefix + ".filter";
            _gfile = prefix + ".gp";
            if (_terminal == Ps)
                _pfile = prefix + ".ps";
            else if (_terminal == Png)
                _pfile = prefix + ".png";
            if (_idRfile is not null && !_fs.ContainsKey(_idRfile))
                throw Die("ERROR: Could not read " + _idRfile + ", No such file or directory\n");
            if (_idQfile is not null && !_fs.ContainsKey(_idQfile))
                throw Die("ERROR: Could not read " + _idQfile + ", No such file or directory\n");
            if ((_rport is not null || _qport is not null) && (_terminal != X11 || _gpstatus != 0))
            {
                Warn("WARNING: Port options available only for v4.0 X11 plots\n");
                _rport = null;
                _qport = null;
            }

            if (_color is not null && _hfile is not null)
            {
                Warn("WARNING: Turning off --color option so highlighting is visible\n");
                _color = null;
            }
        }

        private string GetParseFunc()
        {
            var f = new Lines(_match);
            string? line = f.Next();
            if (line is null)
                throw Die("ERROR: Could not read " + _mfile + ", File is empty\n");
            if (Regex.IsMatch(line, @"^>\S+ [0-9]+ bases"))
                return "tiling";
            if (Regex.IsMatch(line, @"^> \S+"))
                return "mummer";
            var m = Regex.Match(line, @"^(\S+) (\S+)");
            if (m.Success)
            {
                _idRfile ??= m.Groups[1].Value;
                _idQfile ??= m.Groups[2].Value;
                line = f.Next();
                if (line is not null && (Regex.IsMatch(line, @"^NUCMER$") || Regex.IsMatch(line, @"^PROMER$")))
                {
                    f.Next();
                    line = f.Next();
                    if (line is null)
                        return "delta";
                    if (Regex.IsMatch(line, @"^[0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+$"))
                        return "delta";
                    if (Regex.IsMatch(line, @"^[ \-][1-3] [ \-][1-3]$"))
                        return "cluster";
                }
            }

            throw Die("ERROR: Could not read " + _mfile + ", Unrecognized file type\n");
        }

        // ------------------------------------------------------------------ id files

        private void ParseIds(string file, Dictionary<string, Seq> href, List<string> order)
        {
            if (!_fs.TryGetValue(file, out string? text))
            {
                Warn("WARNING: Could not open " + file + ", No such file or directory\n");
                Warn("WARNING: Trouble closing " + file + ", Bad file descriptor\n");
                return;
            }

            Seq? aref = null;
            bool isfasta = false;
            double offset = 0;
            var f = new Lines(text);
            string? line;
            while ((line = f.Next()) is not null)
            {
                if (Regex.IsMatch(line, @"^\s*$"))
                    continue;
                var m = Regex.Match(line, @"^>(\S+)");
                if (m.Success)
                {
                    string id = m.Groups[1].Value;
                    if (href.ContainsKey(id))
                    {
                        Warn("WARNING: Duplicate sequence '" + id + "' ignored\n");
                        aref = null;
                        continue;
                    }

                    if (!isfasta)
                        isfasta = true;
                    if (aref is not null)
                        offset += aref.Len - 1;
                    aref = new Seq(offset, 0, 1);
                    href[id] = aref;
                    order.Add(id);
                    continue;
                }

                if (isfasta && Regex.IsMatch(line, @"^\S+$"))
                {
                    if (aref is not null)
                        aref.Len += line.Length - 1;
                    continue;
                }

                m = Regex.Match(line, @"^(\S+)\s+([0-9]+)\s+([+-]?)$");
                if (!isfasta && m.Success)
                {
                    string id = m.Groups[1].Value;
                    if (href.ContainsKey(id))
                    {
                        Warn("WARNING: Duplicate sequence '" + id + "' ignored\n");
                        aref = null;
                        continue;
                    }

                    double dir = m.Groups[3].Value == "-" ? -1 : 1;
                    aref = new Seq(offset, Num(m.Groups[2].Value), dir);
                    offset += Num(m.Groups[2].Value) - 1;
                    href[id] = aref;
                    order.Add(id);
                    continue;
                }

                Warn("WARNING: Could not parse " + file + "\n" + line);
                href.Clear();
                order.Clear();
                break;
            }
        }

        // ------------------------------------------------------------------ parsers

        private double _tot, _sim;

        private void ParseDelta(string text, List<Align> aligns)
        {
            Warn("Reading delta file " + _mfile + "\n");
            var f = new Lines(text);
            f.Next();
            string? second = f.Next();
            bool ispromer = second is not null && second.StartsWith("PROMER", StringComparison.Ordinal);
            Align? cur = null;
            string lenR = "", lenQ = "", idR = "", idQ = "";
            string? line;
            while ((line = f.Next()) is not null)
            {
                var m = Regex.Match(line, @"^(-?[0-9]+)$");
                if (m.Success)
                {
                    double v = Num(m.Groups[1].Value);
                    if (v < 0)
                    {
                        _tot++;
                    }
                    else if (v == 0)
                    {
                        if (_tot == 0)
                            throw Die("Illegal division by zero\n");
                        double sim = (_tot - _sim) / _tot * 100.0;
                        cur!.Sim = sim;
                        cur.SimText = Str(sim);
                        aligns.Add(cur);
                        _tot = 0;
                    }

                    continue;
                }

                m = Regex.Match(line, @"^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) [0-9]+ ([0-9]+) [0-9]+$");
                if (m.Success && _tot == 0)
                {
                    cur = new Align(m.Groups[1].Value, m.Groups[2].Value, m.Groups[3].Value, m.Groups[4].Value, 0, "0", lenR, lenQ, idR, idQ);
                    _tot = Math.Abs(Num(m.Groups[1].Value) - Num(m.Groups[2].Value)) + 1;
                    _sim = Num(m.Groups[5].Value);
                    if (ispromer)
                        _tot /= 3.0;
                    continue;
                }

                m = Regex.Match(line, @"^>(\S+) (\S+) ([0-9]+) ([0-9]+)$");
                if (m.Success)
                {
                    idR = m.Groups[1].Value;
                    idQ = m.Groups[2].Value;
                    lenR = m.Groups[3].Value;
                    lenQ = m.Groups[4].Value;
                    _tot = 0;
                    continue;
                }

                throw Die("ERROR: Could not parse " + _mfile + "\n" + line);
            }
        }

        private void ParseCluster(string text, List<Align> aligns)
        {
            Warn("Reading cluster file " + _mfile + "\n");
            var f = new Lines(text);
            f.Next();
            f.Next();
            double dR = 0, dQ = 0;
            string lenR = "", lenQ = "", idR = "", idQ = "";
            string? line;
            while ((line = f.Next()) is not null)
            {
                var m = Regex.Match(line, @"^\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+\S+\s+\S+$");
                if (m.Success)
                {
                    double len = Num(m.Groups[3].Value) - 1;
                    double s1 = Num(m.Groups[1].Value), s2 = Num(m.Groups[2].Value);
                    aligns.Add(new Align(m.Groups[1].Value, Str(s1 + (dR == 1 ? len : -len)), m.Groups[2].Value, Str(s2 + (dQ == 1 ? len : -len)), 100, "100", lenR, lenQ, idR, idQ));
                    continue;
                }

                if (Regex.IsMatch(line, @"^[ \-][1-3] [ \-][1-3]$"))
                {
                    dR = line.StartsWith('-') ? -1 : 1;
                    dQ = Regex.IsMatch(line, @"-[1-3]$") ? -1 : 1;
                    continue;
                }

                m = Regex.Match(line, @"^>(\S+) (\S+) ([0-9]+) ([0-9]+)$");
                if (m.Success)
                {
                    idR = m.Groups[1].Value;
                    idQ = m.Groups[2].Value;
                    lenR = m.Groups[3].Value;
                    lenQ = m.Groups[4].Value;
                    continue;
                }

                throw Die("ERROR: Could not parse " + _mfile + "\n" + line);
            }
        }

        private void ParseMummer(string text, List<Align> aligns)
        {
            Warn("Reading mummer file " + _mfile + " (use mummer -c)\n");
            var f = new Lines(text);
            double dQ = 0;
            string lenQ = "", idQ = "";
            string? line;
            while ((line = f.Next()) is not null)
            {
                var m = Regex.Match(line, @"^\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)$");
                if (m.Success)
                {
                    double len = Num(m.Groups[3].Value) - 1;
                    double s1 = Num(m.Groups[1].Value), s2 = Num(m.Groups[2].Value);
                    aligns.Add(new Align(m.Groups[1].Value, Str(s1 + len), m.Groups[2].Value, Str(s2 + (dQ == 1 ? len : -len)), 100, "100", "0", lenQ, "REF", idQ));
                    continue;
                }

                m = Regex.Match(line, @"^\s+(\S+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)$");
                if (m.Success)
                {
                    double len = Num(m.Groups[4].Value) - 1;
                    double s1 = Num(m.Groups[2].Value), s2 = Num(m.Groups[3].Value);
                    aligns.Add(new Align(m.Groups[2].Value, Str(s1 + len), m.Groups[3].Value, Str(s2 + (dQ == 1 ? len : -len)), 100, "100", "0", lenQ, m.Groups[1].Value, idQ));
                    continue;
                }

                m = Regex.Match(line, @"^> (\S+)");
                if (m.Success)
                {
                    idQ = m.Groups[1].Value;
                    dQ = Regex.IsMatch(line, @"^> \S+ Reverse") ? -1 : 1;
                    var l = Regex.Match(line, @"Len = ([0-9]+)");
                    lenQ = l.Success ? l.Groups[1].Value : "0";
                    continue;
                }

                throw Die("ERROR: Could not parse " + _mfile + "\n" + line);
            }
        }

        private void ParseTiling(string text, List<Align> aligns)
        {
            Warn("Reading tiling file " + _mfile + "\n");
            var f = new Lines(text);
            string lenR = "", idR = "";
            string? line;
            while ((line = f.Next()) is not null)
            {
                var m = Regex.Match(line, @"^(\S+)\s+\S+\s+\S+\s+([0-9]+)\s+\S+\s+(\S+)\s+([+-])\s+(\S+)$");
                if (m.Success)
                {
                    double len = Num(m.Groups[2].Value) - 1;
                    string start = m.Groups[1].Value;
                    string endText = Str(Num(start) + len);
                    string s2 = "1", e2 = "1";
                    if (m.Groups[4].Value == "-")
                        s2 = Str(1 + len);
                    else
                        e2 = Str(1 + len);
                    aligns.Add(new Align(start, endText, s2, e2, Num(m.Groups[3].Value), m.Groups[3].Value, lenR, m.Groups[2].Value, idR, m.Groups[5].Value));
                    continue;
                }

                m = Regex.Match(line, @"^>(\S+) ([0-9]+) bases$");
                if (m.Success)
                {
                    idR = m.Groups[1].Value;
                    lenR = m.Groups[2].Value;
                    continue;
                }

                throw Die("ERROR: Could not parse " + _mfile + "\n" + line);
            }

            _coverage ??= true;
        }

        // ------------------------------------------------------------------ layout

        private sealed class Chain(double len)
        {
            public bool Placed;
            public double Len = len;
            public readonly Dictionary<string, double[]> Links = [];
            public readonly List<string> Order = [];
        }

        private sealed class View(double[] v, bool query)
        {
            public double Slope { get => v[0]; set => v[0] = value; }

            public double Lo { get => query ? v[3] : v[1]; set { if (query) v[3] = value; else v[1] = value; } }

            public double Hi { get => query ? v[4] : v[2]; set { if (query) v[4] = value; else v[2] = value; } }

        }

        private void LayoutIds(Dictionary<string, Seq> refs, List<string> refOrder, Dictionary<string, Seq> qrys, List<string> qryOrder)
        {
            var rc = new Dictionary<string, Chain>();
            var qc = new Dictionary<string, Chain>();
            var rcOrder = new List<string>();
            var rl = new List<(string Id, double Slope)>();
            var ql = new List<(string Id, double Slope)>();
            string btab = MummerShowCoords.Run(_files[_dfile!], new MummerShowCoordsOptions { Btab = true, Now = new DateTime(2023, 11, 14, 22, 13, 20, DateTimeKind.Utc) }).Text;
            foreach (string raw in btab.Split('\n'))
            {
                if (raw.Length == 0)
                    continue;
                var al = raw.Split('\t');
                if (al.Length != 21)
                    throw Die("ERROR: Could not read show-coords pipe, invalid btab format\n");
                double sR = Num(al[8]), eR = Num(al[9]), sQ = Num(al[6]), eQ = Num(al[7]);
                double lenR = Num(al[18]), lenQ = Num(al[2]);
                string idR = al[5], idQ = al[0];
                if (!refs.ContainsKey(idR) || !qrys.ContainsKey(idQ))
                    continue;
                double dR = sR < eR ? 1 : -1, dQ = sQ < eQ ? 1 : -1;
                double slope = dR == dQ ? 1 : -1;
                double loR = dR == 1 ? sR : eR, hiR = dR == 1 ? eR : sR;
                double loQ = dQ == 1 ? sQ : eQ, hiQ = dQ == 1 ? eQ : sQ;
                if (_onlyFattest && qc.TryGetValue(idQ, out var qold) && qold.Links.Count > 0)
                {
                    string oldR = qold.Order[0];
                    var val = qold.Links[oldR];
                    if (val[2] - val[1] > hiR - loR)
                        continue;
                    rc[oldR].Links.Remove(idQ);
                    rc[oldR].Order.Remove(idQ);
                    qc.Remove(idQ);
                }

                if (!rc.TryGetValue(idR, out var rchain))
                {
                    rc[idR] = rchain = new Chain(lenR);
                    rcOrder.Add(idR);
                }

                if (!qc.TryGetValue(idQ, out var qchain))
                    qc[idQ] = qchain = new Chain(lenQ);
                if (!rchain.Links.TryGetValue(idQ, out var existing) || !qchain.Links.ContainsKey(idR) || hiR - loR > existing[2] - existing[1])
                {
                    var v = new[] { slope, loR, hiR, loQ, hiQ };
                    if (!rchain.Links.ContainsKey(idQ))
                        rchain.Order.Add(idQ);
                    if (!qchain.Links.ContainsKey(idR))
                        qchain.Order.Add(idR);
                    rchain.Links[idQ] = v;
                    qchain.Links[idR] = v;
                }
            }

            foreach (string idR in rcOrder.OrderByDescending(k => rc[k].Len))
                Span(idR, new Side(rc, rl), new Side(qc, ql), false);

            foreach (var kv in refs)
                kv.Value.Off = null;
            foreach (var kv in qrys)
                kv.Value.Off = null;
            double roff = 0;
            foreach (var r in rl)
            {
                refs[r.Id].Off = roff;
                refs[r.Id].Dir = r.Slope;
                roff += refs[r.Id].Len - 1;
            }

            foreach (string id in refOrder)
            {
                if (refs[id].Off is null)
                {
                    refs[id].Off = roff;
                    roff += refs[id].Len - 1;
                }
            }

            double qoff = 0;
            foreach (var q in ql)
            {
                qrys[q.Id].Off = qoff;
                qrys[q.Id].Dir = q.Slope;
                qoff += qrys[q.Id].Len - 1;
            }

            foreach (string id in qryOrder)
            {
                if (qrys[id].Off is null)
                {
                    qrys[id].Off = qoff;
                    qoff += qrys[id].Len - 1;
                }
            }
        }

        private sealed record Side(Dictionary<string, Chain> Chains, List<(string Id, double Slope)> Placed);

        private static void Span(string x, Side first, Side second, bool firstIsQuery)
        {
            var post = new List<string>();
            var xn = first.Chains[x];
            foreach (string y in xn.Order.OrderBy(k => new View(xn.Links[k], firstIsQuery).Lo))
            {
                var yn = second.Chains[y];
                if (yn.Placed)
                    continue;
                yn.Placed = true;
                double len = yn.Len;
                double slope = new View(xn.Links[y], firstIsQuery).Slope;
                if (slope == -1)
                {
                    foreach (string xx in yn.Order)
                    {
                        var v = new View(yn.Links[xx], !firstIsQuery);
                        v.Slope *= -1;
                        double loy = v.Lo, hiy = v.Hi;
                        v.Lo = len - hiy + 1;
                        v.Hi = len - loy + 1;
                    }
                }

                second.Placed.Add((y, slope));
                if (len > xn.Len)
                    Recurse(y, second, first, !firstIsQuery);
                else
                    post.Add(y);
            }

            foreach (string y in post)
                Recurse(y, second, first, !firstIsQuery);
        }

        private static void Recurse(string y, Side now, Side other, bool nowIsQuery) => Span(y, now, other, nowIsQuery);

        // ------------------------------------------------------------------ plot data

        private void PlotData(List<Align> aligns, Dictionary<string, Seq> refs, List<string> refOrder, Dictionary<string, Seq> qrys, List<string> qryOrder)
        {
            Warn("Writing plot files " + _ffile + ", " + _rfile + (_hfile is not null ? ", " + _hfile + "\n" : "\n"));
            var fwd = new StringBuilder("#-- forward hits sorted by %sim\n0 0 0\n0 0 0\n\n\n");
            var rev = new StringBuilder("#-- reverse hits sorted by %sim\n0 0 0\n0 0 0\n\n\n");
            var hlt = _hfile is not null ? new StringBuilder("#-- highlighted hits sorted by %sim\n0 0 0\n0 0 0\n\n\n") : null;
            bool isplotted = false, ismultiref = false, ismultiqry = false;
            string? pidR = null, pidQ = null;
            string plenR = string.Empty, plenQ = string.Empty;
            foreach (var al in aligns.OrderBy(x => x.Sim))
            {
                if (pidR is null)
                {
                    plenR = al.LenR;
                    plenQ = al.LenQ;
                    pidR = al.IdR;
                    pidQ = al.IdQ;
                }

                double refoff, reflen, refdir, qryoff, qrylen, qrydir;
                if (refs.Count > 0)
                {
                    if (!refs.TryGetValue(al.IdR, out var rs))
                        continue;
                    (refoff, reflen, refdir) = (rs.Off ?? 0, rs.Len, rs.Dir);
                }
                else
                {
                    (refoff, reflen, refdir) = (0, Num(al.LenR), 1);
                }

                if (qrys.Count > 0)
                {
                    if (!qrys.TryGetValue(al.IdQ, out var qs))
                        continue;
                    (qryoff, qrylen, qrydir) = (qs.Off ?? 0, qs.Len, qs.Dir);
                }
                else
                {
                    (qryoff, qrylen, qrydir) = (0, Num(al.LenQ), 1);
                }

                double sR = Num(al.SR), eR = Num(al.ER), sQ = Num(al.SQ), eQ = Num(al.EQ);
                if (refdir == -1)
                {
                    sR = reflen - sR + 1;
                    eR = reflen - eR + 1;
                }

                if (qrydir == -1)
                {
                    sQ = qrylen - sQ + 1;
                    eQ = qrylen - eQ + 1;
                }

                var targets = new List<StringBuilder>();
                if (_breaklen is int bl && ((sR - 1 > bl && sQ - 1 > bl && reflen - sR > bl && qrylen - sQ > bl) || (eR - 1 > bl && eQ - 1 > bl && reflen - eR > bl && qrylen - eQ > bl)))
                    targets.Add(hlt!);
                targets.Add((sR < eR) == (sQ < eQ) ? fwd : rev);
                sR += refoff;
                eR += refoff;
                sQ += qryoff;
                eQ += qryoff;
                foreach (var fh in targets)
                {
                    if (_coverage ?? false)
                        fh.Append(Str(sR)).Append(" 10 ").Append(al.SimText).Append('\n').Append(Str(eR)).Append(" 10 ").Append(al.SimText).Append("\n\n\n")
                            .Append(Str(sR)).Append(' ').Append(al.SimText).Append(" 0\n").Append(Str(eR)).Append(' ').Append(al.SimText).Append(" 0\n\n\n");
                    else
                        fh.Append(Str(sR)).Append(' ').Append(Str(sQ)).Append(' ').Append(al.SimText).Append('\n').Append(Str(eR)).Append(' ').Append(Str(eQ)).Append(' ').Append(al.SimText).Append("\n\n\n");
                }

                if (!ismultiref && al.IdR != pidR)
                    ismultiref = true;
                if (!ismultiqry && al.IdQ != pidQ)
                    ismultiqry = true;
                isplotted = true;
            }

            if (_snp is not null)
            {
                Warn("Determining SNPs from sequence and alignment data\n");
                SnpsInto(hlt!, refs, qrys);
            }

            _files[_ffile] = fwd.ToString();
            _files[_rfile] = rev.ToString();
            if (hlt is not null)
                _files[_hfile!] = hlt.ToString();

            if (refs.Count == 0)
            {
                if (ismultiref)
                    Warn("WARNING: Multiple ref sequences overlaid, try -R or -r\n");
                else if (pidR is not null)
                {
                    refs[pidR] = new Seq(0, Num(plenR), 1);
                    refOrder.Add(pidR);
                }
            }

            if (qrys.Count == 0)
            {
                if (ismultiqry && !(_coverage ?? false))
                    Warn("WARNING: Multiple qry sequences overlaid, try -Q, -q or -c\n");
                else if (pidQ is not null)
                {
                    qrys[pidQ] = new Seq(0, Num(plenQ), 1);
                    qryOrder.Add(pidQ);
                }
            }

            if (!isplotted)
                throw Die("ERROR: No alignment data to plot\n");
        }

        private void SnpsInto(StringBuilder hlt, Dictionary<string, Seq> refs, Dictionary<string, Seq> qrys)
        {
            var delta = new Lines(_current);
            string? names = delta.Next();
            var nm = names is null ? [] : Split(names);
            string refText = nm.Count > 0 && _fs.TryGetValue(nm[0], out var rt) ? rt : string.Empty;
            string qryText = nm.Count > 1 && _fs.TryGetValue(nm[1], out var qt) ? qt : string.Empty;
            string snps = MummerShowSnps.Run(_current, refText, qryText, new MummerShowSnpsOptions { PrintHeader = false, Tabular = true, ShowLength = true }).Text;
            foreach (string raw in snps.Split('\n'))
            {
                if (raw.Length == 0)
                    continue;
                var s = raw.Split('\t');
                if (s.Length != 14)
                    throw Die("ERROR: Could not read show-snps pipe, invalid format\n");
                double pR = Num(s[0]), pQ = Num(s[3]);
                string lenR = s[8], lenQ = s[9], idR = s[12], idQ = s[13];
                double refoff, reflen, refdir, qryoff, qrylen, qrydir;
                if (refs.Count > 0)
                {
                    if (!refs.TryGetValue(idR, out var rs))
                        continue;
                    (refoff, reflen, refdir) = (rs.Off ?? 0, rs.Len, rs.Dir);
                }
                else
                {
                    (refoff, reflen, refdir) = (0, Num(lenR), 1);
                }

                if (qrys.Count > 0)
                {
                    if (!qrys.TryGetValue(idQ, out var qs))
                        continue;
                    (qryoff, qrylen, qrydir) = (qs.Off ?? 0, qs.Len, qs.Dir);
                }
                else
                {
                    (qryoff, qrylen, qrydir) = (0, Num(lenQ), 1);
                }

                if (refdir == -1)
                    pR = reflen - pR + 1;
                if (qrydir == -1)
                    pQ = qrylen - pQ + 1;
                pR += refoff;
                pQ += qryoff;
                if (_coverage ?? false)
                    hlt.Append(Str(pR)).Append(" 10 0\n").Append(Str(pR)).Append(" 10 0\n\n\n");
                else
                    hlt.Append(Str(pR)).Append(' ').Append(Str(pQ)).Append(" 0\n").Append(Str(pR)).Append(' ').Append(Str(pQ)).Append(" 0\n\n\n");
            }
        }

        // ------------------------------------------------------------------ gnuplot script

        private void WriteGp(Dictionary<string, Seq> refs, Dictionary<string, Seq> qrys)
        {
            Warn("Writing gnuplot script " + _gfile + "\n");
            var g = new StringBuilder();
            double size = TerminalSize(_terminal, _size);
            bool covOrColor = (_coverage ?? false) || (_color ?? false);
            string pTerm, pSize;
            double ps, lw;
            switch (_terminal)
            {
                case X11:
                    pTerm = _gpstatus == 0 ? "x11 font \"Courier,8\"" : "x11";
                    ps = 1;
                    lw = covOrColor ? 3 : 2;
                    pSize = "set size 1,1";
                    break;
                case Ps:
                    pTerm = _color is false ? "postscript monochrome" : "postscript color";
                    pTerm += " solid \"Courier\" 8";
                    ps = 0.5;
                    lw = covOrColor ? 4 : 2;
                    pSize = (_coverage ?? false) ? "set size " + Str(1.0 * size) + "," + Str(0.5 * size) : "set size " + Str(1.0 * size) + "," + Str(1.0 * size);
                    break;
                default:
                    pTerm = _gpstatus == 0 ? "png tiny size " + Str(size) + "," + Str(size) : "png small";
                    if (_color is false)
                        pTerm += " xffffff x000000 x000000 x000000 x000000 x000000 x000000 x000000 x000000";
                    ps = 1;
                    lw = 3;
                    pSize = (_coverage ?? false) ? "set size 1,.375" : "set size 1,1";
                    break;
            }

            string pWith = covOrColor ? "w l" : "w lp";
            string pFormat = "set format \"%.0f\"";
            string pLs, pKey;
            if (_gpstatus == 0)
            {
                pLs = "set style line";
                pKey = "unset key";
                pFormat += "\nset mouse format \"%.0f\"\nset mouse mouseformat \"[%.0f, %.0f]\"";
            }
            else
            {
                pLs = "set linestyle";
                pKey = "set nokey";
            }

            var lt = _hfile is not null ? new[] { 0, 2, 2, 1 } : new[] { 0, 1, 3, 2 };
            var refk = refs.Keys.ToList();
            var qryk = qrys.Keys.ToList();
            string xlabel, ylabel;
            string xrange, yrange;
            int border = 0;
            g.Append("set terminal ").Append(pTerm).Append('\n');
            if (_pfile is not null)
                g.Append("set output \"").Append(_pfile).Append("\"\n");
            if (_o.Title is not null)
                g.Append("set title \"").Append(_o.Title).Append("\"\n");

            if (refk.Count == 1)
            {
                xlabel = refk[0];
                xrange = Str(refs[xlabel].Len);
            }
            else
            {
                double total = 0;
                g.Append("set xtics rotate ( \\\n");
                foreach (string k in refk.OrderBy(k => refs[k].Off ?? 0))
                {
                    total += refs[k].Len;
                    double tic = (refs[k].Off ?? 0) + 1;
                    g.Append(" \"").Append(refs[k].Dir == 1 ? "" : "*").Append(k).Append("\" ").Append(Str(tic)).Append(", \\\n");
                }

                g.Append(" \"\" ").Append(Str(total)).Append(" \\\n)\n");
                xlabel = "REF";
                xrange = Str(total);
            }

            if (xrange == "0")
                xrange = "*";

            if (_coverage ?? false)
            {
                ylabel = "%SIM";
                yrange = "110";
            }
            else if (qryk.Count == 1)
            {
                ylabel = qryk[0];
                yrange = Str(qrys[ylabel].Len);
            }
            else
            {
                double total = 0;
                g.Append("set ytics ( \\\n");
                foreach (string k in qryk.OrderBy(k => qrys[k].Off ?? 0))
                {
                    total += qrys[k].Len;
                    double tic = (qrys[k].Off ?? 0) + 1;
                    g.Append(" \"").Append(qrys[k].Dir == 1 ? "" : "*").Append(k).Append("\" ").Append(Str(tic)).Append(", \\\n");
                }

                g.Append(" \"\" ").Append(Str(total)).Append(" \\\n)\n");
                ylabel = "QRY";
                yrange = Str(total);
            }

            if (yrange == "0")
                yrange = "*";
            if (xrange != "*" && refk.Count == 1)
                border |= 10;
            if (yrange != "*" && qryk.Count == 1)
                border |= 5;
            if (_coverage ?? false)
                border |= 5;

            g.Append(pSize).Append('\n').Append("set grid\n").Append(pKey).Append('\n').Append("set border ").Append(border.ToString(CultureInfo.InvariantCulture)).Append('\n')
                .Append("set tics scale 0\n").Append("set xlabel \"").Append(xlabel).Append("\"\n").Append("set ylabel \"").Append(ylabel).Append("\"\n").Append(pFormat).Append('\n');
            g.Append(_xrange is not null ? "set xrange " + _xrange + "\n" : "set xrange [1:" + xrange + "]\n");
            g.Append(_yrange is not null ? "set yrange " + _yrange + "\n" : "set yrange [1:" + yrange + "]\n");
            if (_color ?? false)
            {
                g.Append("set zrange [0:100]\nset colorbox default\nset cblabel \"%similarity\"\nset cbrange [0:100]\nset cbtics 20\nset pm3d map\nset palette model RGB defined ( \\\n")
                    .Append("  0 \"#000000\", \\\n  4 \"#DD00DD\", \\\n  6 \"#0000DD\", \\\n  7 \"#00DDDD\", \\\n  8 \"#00DD00\", \\\n  9 \"#DDDD00\", \\\n 10 \"#DD0000\"  \\\n)\n");
            }

            for (int s = 1; s <= 3; s++)
            {
                string ss = pLs + " " + s.ToString(CultureInfo.InvariantCulture) + " ";
                ss += (_color ?? false) ? " palette" : " lt " + lt[s].ToString(CultureInfo.InvariantCulture);
                ss += " lw " + Str(lw);
                if (!(_coverage ?? false) || s == 3)
                    ss += " pt 6 ps " + Str(ps);
                g.Append(ss).Append('\n');
            }

            g.Append((_color ?? false) ? "splot \\\n" : "plot \\\n");
            g.Append(" \"").Append(_ffile).Append("\" title \"FWD\" ").Append(pWith).Append(" ls 1, \\\n")
                .Append(" \"").Append(_rfile).Append("\" title \"REV\" ").Append(pWith).Append(" ls 2")
                .Append(_hfile is null ? "\n" : ", \\\n \"" + _hfile + "\" title \"HLT\" w lp ls 3");
            if (_terminal == X11)
                g.Append("\nprint \"-- INTERACTIVE MODE --\"\nprint \"consult gnuplot docs for command list\"\nprint \"mouse 1: coords to clipboard\"\nprint \"mouse 2: mark on plot\"\nprint \"mouse 3: zoom box\"\nprint \"'h' for help in plot window\"\nprint \"enter to exit\"\npause -1\n");
            _files[_gfile] = g.ToString();
        }

        private static double TerminalSize(string terminal, string size)
        {
            int k = 2;
            if (size == "small")
                k = 0;
            else if (size == "medium")
                k = 1;
            double[] x11 = [500, 700, 900], ps = [1, 2, 3], png = [800, 1024, 1400];
            double[] table = png;
            if (terminal == X11)
                table = x11;
            else if (terminal == Ps)
                table = ps;
            return table[k];
        }

        private string RunGp()
        {
            Warn(_pfile is not null ? "Rendering plot " + _pfile + "\n" : "Rendering plot to screen\n");
            string cmd = "gnuplot";
            if (_terminal == X11)
            {
                double size = TerminalSize(X11, _size);
                cmd += " -geometry " + Str(size) + "x";
                if (_coverage ?? false)
                    size = Num(F0(size * .375));
                cmd += Str(size) + "+0+0 -title mummerplot";
                if (_color is false)
                    cmd += " -mono -xrm 'gnuplot*line1Dashes: 0' -xrm 'gnuplot*line2Dashes: 0' -xrm 'gnuplot*line3Dashes: 0'";
                if (_o.ReverseVideo)
                    cmd += " -rv -xrm 'gnuplot*background: black' -xrm 'gnuplot*textColor: white' -xrm 'gnuplot*borderColor: white' -xrm 'gnuplot*axisColor: white'";
            }

            return cmd + " " + _gfile;
        }
    }
}
