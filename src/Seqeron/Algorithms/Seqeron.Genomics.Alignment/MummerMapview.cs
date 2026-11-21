using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;
using static Seqeron.Genomics.Alignment.MummerPerl;

namespace Seqeron.Genomics.Alignment;

/// <summary>The switches of <c>mapview</c> (MUMmer 3.23).</summary>
public sealed record MummerMapviewOptions
{
    /// <summary><c>-d</c>: the largest distance in base pairs between linked matches (0 = the default 50 000).</summary>
    public int MaxDistance { get; init; }

    /// <summary><c>-f</c>: the output format <c>fig</c> (default), <c>pdf</c> or <c>ps</c> (the last two are produced by <c>fig2dev</c> from the FIG text, which the result returns together with the command line).</summary>
    public string Format { get; init; } = "fig";

    /// <summary><c>-m</c>: the magnification of <c>fig2dev</c> (0 = none).</summary>
    public double Magnification { get; init; }

    /// <summary><c>-n</c>: the number of output files that partition the reference (0 = the default 10).</summary>
    public int FileCount { get; init; }

    /// <summary><c>-p</c>: the prefix of the output files (empty = <c>PROMER_graph</c> / <c>NUCMER_graph</c>: the program name in the coords file).</summary>
    public string Prefix { get; init; } = string.Empty;

    /// <summary><c>-v</c>: list the output files and their ranges.</summary>
    public bool Verbose { get; init; }

    /// <summary><c>-x1</c> / <c>-x2</c>: the lower and upper coordinate bound of the display (0 = unset).</summary>
    public int X1 { get; init; }

    /// <summary>See <see cref="X1"/>.</summary>
    public int X2 { get; init; }

    /// <summary><c>-g</c>: the input is the output of <c>mgaps</c>; this is the reference sequence id (as in the first column of the GFF files).</summary>
    public string MgapsReference { get; init; } = string.Empty;

    /// <summary><c>-I</c>: display the names of the query sequences.</summary>
    public bool PrintQueryNames { get; init; }

    /// <summary><c>-Ir</c>: display the names of the reference genes.</summary>
    public bool PrintGeneNames { get; init; }

    /// <summary>The file name of the coords file as the script receives it (the <c>mgaps</c> input gets the derived file <c>&lt;name&gt;coords</c>).</summary>
    public string CoordsName { get; init; } = "coords";

    /// <summary>Name of the UTR file (used only to decide whether it equals the CDS file: the exon lists are then not split).</summary>
    public string UtrName { get; init; } = "utr.gff";

    /// <summary>Name of the CDS file.</summary>
    public string CdsName { get; init; } = "cds.gff";
}

/// <summary>Output of <see cref="MummerMapview.Run"/>.</summary>
/// <param name="Figs">The FIG files the script writes (by name; with <c>pdf</c> / <c>ps</c> the script writes them as <c>tmp&lt;prefix&gt;_n.fig</c>, converts and removes them).</param>
/// <param name="Commands">The <c>fig2dev</c> command lines of the script (empty for <c>fig</c>); their outputs are the files named in <see cref="Outputs"/>.</param>
/// <param name="Outputs">The names of the converted files (<c>pdf</c> / <c>ps</c>), in the order of the commands.</param>
/// <param name="Stdout">Standard output (the verbose lines, the reference lines of a multi-reference run, the messages of wrong options).</param>
/// <param name="Stderr">Standard error.</param>
/// <param name="ExitCode">The exit status (1 when the GFF files do not belong to the coords file).</param>
/// <param name="Side">The side file of a <c>-g</c> run (<c>&lt;coords&gt;coords</c>: the mgaps input in the format of show-coords).</param>
public sealed record MummerMapviewResult(IReadOnlyDictionary<string, string> Figs, IReadOnlyList<string> Commands, IReadOnlyList<string> Outputs, string Stdout, string Stderr, int ExitCode, IReadOnlyDictionary<string, string> Side);

/// <summary>
/// <c>mapview</c> of MUMmer 3.23 (perl script of TIGR, <c>scripts/mapview.pl</c>): turns the output of <c>show-coords -rl</c> (or of <c>mgaps</c>, <c>-g</c>) into FIG files that show the alignments along the reference
/// (percent identity on the vertical axis, matches of one query contig linked by a line of one of three colours), optionally with the genes of a UTR and a CDS file in GFF format (5' UTR, exons, 3' UTR, introns), partitioned into
/// <c>-n</c> files. The script is ported with the behaviour of Perl's loose typing it relies on (numeric strings, <c>split</c>, undefined values, <c>int</c>, <c>%.0f</c>).
/// </summary>
public static class MummerMapview
{
    /// <summary>Runs <c>mapview [options] coords [utr cds]</c>.</summary>
    /// <param name="coords">The text of the coords file.</param>
    /// <param name="utrGff">The UTR file in GFF format (null = none; with a UTR file and no CDS file the script stops: its genes cannot be placed).</param>
    /// <param name="cdsGff">The CDS file in GFF format.</param>
    /// <param name="options">The switches (null = defaults).</param>
    public static MummerMapviewResult Run(string coords, string? utrGff = null, string? cdsGff = null, MummerMapviewOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(coords);
        return new Engine(options ?? new MummerMapviewOptions(), coords, utrGff, cdsGff).Execute();
    }

    private sealed class Engine
    {
        private const int Yorig = 3000, YdistPid = 2000, Gap = 800;
        private const double Yscale = YdistPid / 50.0, Xscale = 14.5;
        private const string ErrGff = "\n  ERROR in the input files ! The reference seq ID can't be found in GFF files !\n      The first column in the GFF file should be the ID of the reference seq. \n      The alignments file should provide the same info in the column before the last one.\n      \n      Here are some example records for the GFF file:\n     \n      gnl|FlyBase|X\tDmel3\tinitial-exon\t2155\t2413\t.\t-\t.\tX_CG3038.1\n      gnl|FlyBase|X\tDmel3\tlast-exon\t1182\t2077\t.\t-\t.\tX_CG3038.1\n      ...\n      The fields are :\n      <seq_ID> <source> <exon type> <start> <end> <score> <strand> <frame> <gene_name>\n    ";

        private readonly MummerMapviewOptions _o;
        private string _alignText;
        private readonly string? _utr, _cds;
        private readonly bool _hasUtr, _hasCds;
        private readonly StringBuilder _stdout = new(), _stderr = new();
        private readonly Dictionary<string, StringBuilder> _figs = [];
        private readonly List<string> _figOrder = [];
        private readonly List<string> _commands = [], _outputs = [];
        private readonly Dictionary<string, string> _side = [];
        private StringBuilder _o_ = new();

        // globals of the script
        private double _matchDist, _noOutfiles, _x1win, _x2win;
        private bool _startfind, _endfind;
        private string _format = "fig", _prog = string.Empty, _procfile = string.Empty;
        private string? _mgaps;
        private int _indS1, _indE1, _indPidy, _indLenchr;
        private readonly Dictionary<string, string?> _hRefContigId = [];
        private int _mref;
        private double _lenrefseqs, _nobpinfile;
        private readonly List<string> _file = [];
        private int _fileno;
        private double _startcoord, _endcoord;
        private string _startText = string.Empty;   // $startcoord is undefined until a range sets it
        private List<string> _a = [];
        private string _tmpIdQrycontig = string.Empty;
        private int[] _linkcolors = [31, 14, 11];
        private int _linkcolor;
        private readonly Dictionary<string, (int Color, string Pos, string End)> _hQrycontig = [];
        private double _xlastPrint;
        private int _nrf;
        private string? _tmpname, _sign, _tmpgene;
        private string? _cds5, _cds3, _utr5, _utr3;
        private readonly Dictionary<string, string> _hcdsEnds = [], _utrEnds = [];
        private readonly Dictionary<string, string> _hContigGenes = [];
        private readonly Dictionary<string, List<string?>> _geneinfo = [];
        private readonly Dictionary<string, string> _utr5End = [];
        private int _testGffFormat;
        private List<string> _g1 = [];
        private string _tmpcontig = string.Empty;
        private double _startdrawX, _proclen, _len;
        private readonly Dictionary<string, string> _hLastOnLevel = [];
        private double _y;

        public Engine(MummerMapviewOptions o, string coords, string? utr, string? cds)
        {
            _o = o;
            _alignText = coords;
            _utr = utr;
            _cds = cds;
            _hasUtr = utr is not null;
            _hasCds = cds is not null;
        }

        public MummerMapviewResult Execute()
        {
            int exit = 0;
            try
            {
                Main();
            }
            catch (MummerExitException e)
            {
                exit = e.Code;
            }

            var figs = new Dictionary<string, string>();
            foreach (string n in _figOrder)
                figs[n] = _figs[n].ToString();
            return new MummerMapviewResult(figs, _commands, _outputs, _stdout.ToString(), _stderr.ToString(), exit, _side);
        }

        private void Out(string s) => _o_.Append(s);

        private void Open(string name)
        {
            _o_ = new StringBuilder();
            _figs[name] = _o_;
            if (!_figOrder.Contains(name))
                _figOrder.Add(name);
        }

        private static string Get(IReadOnlyList<string> a, int i) => At(a, i) ?? string.Empty;

        private void Main()
        {
            _matchDist = _o.MaxDistance;
            _noOutfiles = _o.FileCount;
            _x1win = _o.X1;
            _x2win = _o.X2;
            _mgaps = _o.MgapsReference.Length > 0 ? _o.MgapsReference : null;
            _format = _o.Format;
            // the checks of the coordinates only print (the script goes on)
            if (_x1win < 0 || _x2win < 0)
                _stdout.Append("ERROR2 : coords x1,x2 should be positive integers !!\n");
            if (_x1win > _x2win)
                _stdout.Append("ERROR3 : wrong range coords : x1 >= x2 !!!\n");
            if (_mgaps is not null)
                FormatMgaps();

            if (_format.Length == 0)
                _format = "fig";
            _startfind = !(_x1win != 0 && _x2win != 0);
            _endfind = false;
            if (_noOutfiles == 0)
                _noOutfiles = 10;
            if (_matchDist == 0)
                _matchDist = 50000;

            var f = new Lines(_alignText);
            f.Next();
            string? progLine = f.Next();
            _prog = progLine is null ? string.Empty : Chomp(progLine);
            f.Next();
            string? header = f.Next();
            var a = new List<string>();
            if (header is not null)
            {
                foreach (Match m in Regex.Matches(header, @"\s+(\||\[.+?\])"))
                    a.Add(m.Groups[1].Value);
            }

            for (int ind = 0; ind < a.Count; ind++)
            {
                if (a[ind] == "[S1]")
                    _indS1 = ind;
                else if (a[ind] == "[E1]")
                    _indE1 = ind;
                else if (a[ind] == "[% IDY]")
                    _indPidy = ind;
                else if (a[ind] == "[LEN R]")
                    _indLenchr = ind;
            }

            f.Next();
            _mref = -1;
            string? line;
            while ((line = f.Next()) is not null)
            {
                _a = Split(Chomp(line));
                string key = Get(_a, -2);
                if (!_hRefContigId.ContainsKey(key))
                {
                    _hRefContigId[key] = At(_a, _indLenchr);
                    _lenrefseqs += Num(At(_a, _indLenchr));
                    _mref++;
                }
            }

            _nobpinfile = Int(_lenrefseqs / _noOutfiles);

            int argc = 1 + (_hasUtr ? 1 : 0) + (_hasCds ? 1 : 0);
            if (argc > 1)
            {
                GetCdsEnds();
                GetUtrCdsInfo();
                TestOverlap();
            }
            else if (_mref == 0)
            {
                _fileno = (int)_noOutfiles;
                SetStart(0);
                _endcoord = 0;
                for (int i = 0; i < _fileno; i++)
                {
                    _endcoord = _startcoord + _nobpinfile - 1;
                    if (_endcoord > _lenrefseqs)
                        _endcoord = _lenrefseqs;
                    SetFile(_file, i, _startText + " " + Str(_endcoord));
                    SetStart(_endcoord + 1);
                }
            }

            if (_mref == 0)
            {
                for (int i = 0; i < _fileno; i++)
                {
                    _nrf = i;
                    SetOutputFname();
                    var se = SplitSpaces(_file[i]);
                    _startcoord = Num(At(se, 0));
                    _startText = Get(se, 0);
                    _endcoord = Num(At(se, 1));
                    Open(_procfile + Str(_nrf) + ".fig");
                    PrintHeader();
                    if (_o.Verbose && _format == "fig")
                        _stdout.Append(_procfile).Append(Str(_nrf)).Append(".fig\t range : ").Append(_startText).Append('\t').Append(Str(_endcoord)).Append(" \n");
                    double xs = 0;
                    double xe = Int((_endcoord - _startcoord + 1) / Xscale);
                    PrintGrid(xs, xe, _startText, _endcoord);
                    _tmpIdQrycontig = string.Empty;
                    _linkcolor = _linkcolors[0];

                    var g = new Lines(_alignText);
                    for (int k = 0; k < 5; k++)
                        g.Next();
                    while ((line = g.Next()) is not null)
                    {
                        _a = Split(Chomp(line));
                        if (Num(At(_a, _indS1)) > _endcoord)
                            break;
                        if (Num(At(_a, _indS1)) < _startcoord && Num(At(_a, _indE1)) > _startcoord)
                            Set(_a, _indS1, Str(_startcoord));
                        if (Num(At(_a, _indS1)) < _endcoord && _endcoord < Num(At(_a, _indE1)))
                            Set(_a, _indE1, Str(_endcoord));
                        if (Num(At(_a, _indS1)) >= _startcoord && Num(At(_a, _indE1)) <= _endcoord)
                        {
                            double x1 = Int((Num(At(_a, _indS1)) - _startcoord) / Xscale);
                            double x2 = Int((Num(At(_a, _indE1)) - _startcoord) / Xscale);
                            PrintAlign(x1, x2);
                        }
                    }

                    _hQrycontig.Clear();
                    if (_hasUtr)
                        PrintGenes();
                    PrintLegend();
                    ChangeFileFormat();
                }
            }
            else
            {
                SetOutputFname();
                _tmpIdQrycontig = string.Empty;
                _linkcolor = _linkcolors[0];
                _startdrawX = 0;
                _proclen = 0;
                bool first = true;
                _nrf = 0;
                var g = new Lines(_alignText);
                for (int k = 0; k < 5; k++)
                    g.Next();
                while ((line = g.Next()) is not null)
                {
                    _a = Split(Chomp(line));
                    if (Get(_a, -2) != _tmpcontig)
                    {
                        _hQrycontig.Clear();
                        _tmpcontig = Get(_a, -2);
                        if (first)
                        {
                            first = false;
                            _nrf++;
                            Open(_procfile + Str(_nrf) + ".fig");
                            PrintHeader();
                            if (_o.Verbose && _format == "fig")
                                _stdout.Append(_procfile).Append(Str(_nrf)).Append(".fig\n");
                            _len = Num(_hRefContigId.GetValueOrDefault(Get(_a, -2)));
                        }
                        else
                        {
                            _startdrawX += Int(_len / Xscale) + Gap;
                            _len = Num(_hRefContigId.GetValueOrDefault(Get(_a, -2)));
                            if (_proclen + _len > _nobpinfile && _proclen != 0)
                            {
                                PrintLegend();
                                ChangeFileFormat();
                                _nrf++;
                                Open(_procfile + Str(_nrf) + ".fig");
                                PrintHeader();
                                if (_o.Verbose && _format == "fig")
                                    _stdout.Append('\n').Append(_procfile).Append(Str(_nrf)).Append(".fig\n");
                                _proclen = 0;
                                _startdrawX = 0;
                            }
                        }

                        double xs = _startdrawX;
                        double xe = _startdrawX + Int(_len / Xscale);
                        PrintGrid(xs, xe, "0", _len);
                        _stdout.Append(Get(_a, -2)).Append('\t').Append(_hRefContigId.GetValueOrDefault(Get(_a, -2)) ?? string.Empty).Append('\n');
                        if (_hasUtr)
                            PrintGenesMr();
                        _proclen += _len;
                    }

                    double x1 = _startdrawX + Int(Num(At(_a, _indS1)) / Xscale);
                    double x2 = _startdrawX + Int(Num(At(_a, _indE1)) / Xscale);
                    PrintAlign(x1, x2);
                }

                PrintLegend();
                ChangeFileFormat();
            }
        }

        private void SetStart(double v)
        {
            _startcoord = v;
            _startText = Str(v);
        }

        private void SetOutputFname()
        {
            _procfile = _o.Prefix.Length == 0 ? _prog + "_graph_" : _o.Prefix + "_";
            if (_format != "fig")
                _procfile = "tmp" + _procfile;
        }

        private void PrintHeader() => Out("#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nMultiple\n-2\n1200 2\n");

        // ------------------------------------------------------------------ GFF files

        private static IEnumerable<string> GffLines(string? text)
        {
            if (text is null)
                yield break;
            var l = new Lines(text);
            string? line;
            while ((line = l.Next()) is not null)
                yield return Chomp(line);
        }

        private void GetCdsEnds()
        {
            _testGffFormat = 0;
            foreach (string line in GffLines(_cds))
            {
                if (!True(line))
                    continue;
                var a = Split(line);
                if (_hRefContigId.ContainsKey(Get(a, 0)))
                    _testGffFormat++;
                string genename = Get(a, 8);
                if (genename != (_tmpname ?? string.Empty))
                {
                    if (_sign == "+")
                        _hcdsEnds[_tmpname ?? string.Empty] = Joined(_cds5, _cds3);
                    else if (_sign == "-")
                        _hcdsEnds[_tmpname ?? string.Empty] = Joined(_cds3, _cds5);
                    _tmpname = genename;
                    _sign = Get(a, 6);
                }

                if (_sign == "-")
                    (a[3], a[4]) = (a[4], a[3]);
                string t = Get(a, 2);
                if (t == "single-exon")
                {
                    _cds5 = Get(a, 3);
                    _cds3 = Get(a, 4);
                }
                else if (t == "initial-exon")
                {
                    _cds5 = Get(a, 3);
                }
                else if (t == "last-exon")
                {
                    _cds3 = Get(a, 4);
                }
            }

            if (_sign == "+")
                _hcdsEnds[_tmpname ?? string.Empty] = Joined(_cds5, _cds3);
            else if (_sign == "-")
                _hcdsEnds[_tmpname ?? string.Empty] = Joined(_cds3, _cds5);
            TestFormatGff();
        }

        private static string Joined(string? x, string? y) => (x ?? string.Empty) + " " + (y ?? string.Empty);

        private void TestFormatGff()
        {
            if (_testGffFormat == 0)
            {
                _stderr.Append(ErrGff).Append(" \n");
                throw new MummerExitException(1);
            }
        }

        private List<string> lastUtrLine = [];

        private void GetUtrCdsInfo()
        {
            _testGffFormat = 0;
            foreach (string line in GffLines(_utr))
            {
                if (!True(line))
                    continue;
                var a = Split(line);
                lastUtrLine = a;
                string genename = Get(a, 8);
                if (_hRefContigId.ContainsKey(Get(a, 0)))
                    _testGffFormat++;
                if (genename != (_tmpname ?? string.Empty))
                {
                    if (_sign == "+")
                        _utrEnds[_tmpname ?? string.Empty] = Joined(_utr5, _utr3);
                    else if (_sign == "-")
                        _utrEnds[_tmpname ?? string.Empty] = Joined(_utr3, _utr5);
                    if (True(_tmpgene))
                        GetUtrCdsEnds();
                    _tmpgene = Get(a, 3) + " " + Get(a, 4);
                    _tmpname = genename;
                    _sign = Get(a, 6);
                    string c = Get(a, 0);
                    _hContigGenes[c] = _hContigGenes.GetValueOrDefault(c, string.Empty) + " " + Get(a, 8);
                }
                else
                {
                    string exon = Get(a, 3) + " " + Get(a, 4);
                    _tmpgene = _sign == "-" ? exon + ";" + _tmpgene : _tmpgene + ";" + exon;
                }

                if (_sign == "-")
                    (a[3], a[4]) = (a[4], a[3]);
                string t = Get(a, 2);
                if (t == "single-exon")
                {
                    _utr5 = Get(a, 3);
                    _utr3 = Get(a, 4);
                }
                else if (t == "initial-exon")
                {
                    _utr5 = Get(a, 3);
                }
                else if (t == "last-exon")
                {
                    _utr3 = Get(a, 4);
                }

                if (!_geneinfo.ContainsKey(Get(a, 8)))
                    _geneinfo[Get(a, 8)] = [ "1" ];
            }

            GetUtrCdsEnds();
            if (_sign == "+")
                _utrEnds[_tmpname ?? string.Empty] = Joined(_utr5, _utr3);
            else if (_sign == "-")
                _utrEnds[_tmpname ?? string.Empty] = Joined(_utr3, _utr5);
            string c0 = Get(lastUtrLine, 0);
            _hContigGenes[c0] = _hContigGenes.GetValueOrDefault(c0, string.Empty) + " " + Get(lastUtrLine, 8);
            TestFormatGff();
        }

        private void GetUtrCdsEnds()
        {
            var u5 = new StringBuilder();
            var ex = new StringBuilder();
            var u3 = new StringBuilder();
            string exText;
            if (_o.UtrName == _o.CdsName)
            {
                exText = _tmpgene ?? string.Empty;
            }
            else
            {
                var exons = SplitChar(_tmpgene, ';');
                var cds = Split(_hcdsEnds.GetValueOrDefault(_tmpname ?? string.Empty));
                for (int i = 0; i < exons.Count; i++)
                {
                    var coord = Split(exons[i]);
                    double c0 = Num(At(cds, 0)), c1 = Num(At(cds, 1)), k0 = Num(At(coord, 0)), k1 = Num(At(coord, 1));
                    if (c0 > k0)
                    {
                        if (c0 > k1)
                        {
                            u5.Append(Get(coord, 0)).Append(' ').Append(Get(coord, 1)).Append(';');
                        }
                        else
                        {
                            u5.Append(Get(coord, 0)).Append(' ').Append(Str(c0 - 1)).Append(';');
                            if (c1 < k1)
                            {
                                ex.Append(Get(cds, 0)).Append(' ').Append(Get(cds, 1)).Append(';');
                                u3.Append(Str(c1 + 1)).Append(' ').Append(Get(coord, 1)).Append(';');
                            }
                            else
                            {
                                ex.Append(Get(cds, 0)).Append(' ').Append(Get(coord, 1)).Append(';');
                            }
                        }
                    }
                    else if (c1 > k0)
                    {
                        if (c1 > k1)
                        {
                            ex.Append(Get(coord, 0)).Append(' ').Append(Get(coord, 1)).Append(';');
                        }
                        else
                        {
                            ex.Append(Get(coord, 0)).Append(' ').Append(Get(cds, 1)).Append(';');
                            u3.Append(Str(c1 + 1)).Append(' ').Append(Get(coord, 1)).Append(';');
                        }
                    }
                    else
                    {
                        u3.Append(Get(coord, 0)).Append(' ').Append(Get(coord, 1)).Append(';');
                    }
                }

                exText = Chop(ex.ToString());
            }

            string name = _tmpname ?? string.Empty;
            if (!_geneinfo.TryGetValue(name, out var info))
                _geneinfo[name] = info = [];
            while (info.Count < 5)
                info.Add(null);
            info[1] = _sign;
            info[2] = Chop(u5.ToString());
            info[3] = exText;
            info[4] = Chop(u3.ToString());
        }

        private static void SetFile(List<string> file, int index, string value)
        {
            while (file.Count <= index)
                file.Add(string.Empty);
            file[index] = value;
        }

        private static string Chop(string s) => s.Length > 0 ? s[..^1] : s;

        private void TestOverlap()
        {
            if (_mref == 0)
            {
                _fileno = 0;
                _endcoord = 0;
            }

            foreach (string kcontgid in _hContigGenes.Keys.OrderBy(k => k, StringComparer.Ordinal))
            {
                var all = SplitSpaces(_hContigGenes[kcontgid]);
                for (int i = 1; i < all.Count; i++)
                {
                    _g1 = Split(_utrEnds.GetValueOrDefault(all[i]));
                    _utr5End[all[i]] = Get(_g1, 0);
                    double a0 = Num(At(_g1, 0)), a1 = Num(At(_g1, 1));
                    for (int j = i + 1; j < all.Count; j++)
                    {
                        var g2 = Split(_utrEnds.GetValueOrDefault(all[j]));
                        double b0 = Num(At(g2, 0)), b1 = Num(At(g2, 1));
                        if ((b0 >= a0 && b0 <= a1) || (b1 >= a0 && b1 <= a1) || (a0 >= b0 && a0 <= b1))
                        {
                            double li = Num(Level(all[i])), lj = Num(Level(all[j]));
                            if (li == lj)
                                SetLevel(all[j], Str(li + 1));
                        }
                    }

                    if (!_endfind && _mref == 0)
                        SetTheRangeForEachFile();
                }
            }

            if (_mref == 0)
                SetFile(_file, _fileno++, _startText + " " + Str(_endcoord));
            _utrEnds.Clear();
        }

        private string? Level(string gene) => _geneinfo.TryGetValue(gene, out var l) && l.Count > 0 ? l[0] : null;

        private void SetLevel(string gene, string value)
        {
            if (!_geneinfo.TryGetValue(gene, out var l))
                _geneinfo[gene] = l = [];
            if (l.Count == 0)
                l.Add(value);
            else
                l[0] = value;
        }

        private void SetTheRangeForEachFile()
        {
            double currstart = Num(At(_g1, 0)), currend = Num(At(_g1, 1));
            if (!_startfind && (_x1win <= currstart || _x1win <= currend))
            {
                currstart = _x1win;
                _startfind = true;
            }

            if (_startfind && _x1win != 0 && _x2win != 0 && (_x2win <= currstart || _x2win <= currend))
            {
                currend = _x2win;
                _endfind = true;
            }

            if (_startfind)
            {
                if (_endcoord == 0)
                {
                    SetStart(_x1win != 0 ? _x1win : 0);
                    _endcoord = currend;
                }
                else if (currend > _endcoord)
                {
                    if (currend - _startcoord < _nobpinfile)
                    {
                        _endcoord = currend;
                    }
                    else
                    {
                        SetFile(_file, _fileno++, _startText + " " + Str(_endcoord));
                        SetStart(_endcoord + 1);
                        _endcoord = currend;
                    }
                }
            }
        }

        // ------------------------------------------------------------------ drawing

        private void PrintAlign(double x1, double x2)
        {
            if (Num(At(_a, _indPidy)) < 50)
                Set(_a, _indPidy, "50");
            Set(_a, _indPidy, Str(Int(Num(At(_a, _indPidy)))));
            string qid = Get(_a, _a.Count - 1);
            double y;
            if (_mgaps is not null)
            {
                y = Yorig + 250 + YdistPid - Yscale * 2;
                if (Regex.IsMatch(qid, "rev$"))
                    y -= 25 * Yscale;
            }
            else
            {
                y = Yorig + 250 + YdistPid - Yscale * (Num(At(_a, _indPidy)) - 50);
            }

            if (x1 == x2)
                x2++;
            if (qid == _tmpIdQrycontig)
            {
                PrintConnections(_hQrycontig[_tmpIdQrycontig].Pos, x1, y);
            }
            else
            {
                if (_o.PrintQueryNames && x1 - _xlastPrint > 400)
                {
                    Out("4 0 0 5 0 0 8 0.0000 4 90 270 ");
                    Out("\t" + F0(x1) + " " + F0(y) + " ");
                    Out(qid + "\\001\n");
                    _xlastPrint = x1;
                }

                if (_hQrycontig.TryGetValue(qid, out var seen) && Num(At(_a, _indS1)) - Num(seen.End) < _matchDist)
                {
                    _linkcolor = seen.Color;
                    PrintConnections(seen.Pos, x1, y);
                }
                else
                {
                    // unshift (@linkcolors, pop (@linkcolors))
                    _linkcolors = [_linkcolors[^1], .. _linkcolors[..^1]];
                    _linkcolor = _linkcolors[0];
                    _hQrycontig[qid] = (_linkcolor, seen.Pos, seen.End);
                }
            }

            _tmpIdQrycontig = qid;
            var cur = _hQrycontig[qid];
            _hQrycontig[qid] = (cur.Color, Str(x2) + " " + Str(y), Get(_a, _indE1));
            Out("2 1 0 2 4 0 40 0 -1 0.000 0 0 -1 0 0 2\n");
            Out("\t" + Str(x1) + " " + Str(y) + " " + Str(x2) + " " + Str(y) + "\n");
            Out("2 1 0 5 20 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
            Out("\t " + Str(x1) + " " + F0(Yorig + 150) + " " + Str(x2) + " " + F0(Yorig + 150) + "\n");
        }

        private void PrintConnections(string setc1, double setx2, double sety2)
        {
            var p = setc1.Split(' ');
            double setx1 = Num(p.Length > 0 ? p[0] : null), sety1 = Num(p.Length > 1 ? p[1] : null);
            if (_mgaps is not null)
            {
                if (setx1 > setx2)
                    (setx1, setx2) = (setx2, setx1);
                double dist = Int((setx2 - setx1) / 2);
                double xcenter = setx1 + dist;
                double heightArcUp, yoffcenter;
                if (setx2 - setx1 > 4000)
                {
                    heightArcUp = 20 * Yscale;
                    yoffcenter = Int((dist * dist + heightArcUp * heightArcUp) * (1 / (2 * heightArcUp))) - heightArcUp;
                }
                else
                {
                    heightArcUp = Int(0.447 * dist);
                    yoffcenter = 2 * heightArcUp;
                }

                Out("5 1 0 2 " + _linkcolor.ToString(CultureInfo.InvariantCulture) + " 0 50 0 -1 0.000 0 0 0 0 ");
                Out(Fixed3(xcenter) + " " + Fixed3(sety1 + yoffcenter) + " " + Str(setx1) + " " + Str(sety1) + " " + Str(xcenter) + " " + F0(sety1 - heightArcUp) + " " + Str(setx2) + " " + Str(sety1) + " \n");
            }
            else
            {
                Out("2 1 0 1 " + _linkcolor.ToString(CultureInfo.InvariantCulture) + " 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
                Out("\t" + setc1 + " " + Str(setx2) + " " + Str(sety2) + "\n");
            }
        }

        private static string Fixed3(double v) => Infrastructure.CFormat.Fixed(v, 3);

        private void DrawGene(string kname, double offset, bool multi)
        {
            double tmpx2 = 0;
            var info = _geneinfo[kname];
            for (int l = 2; l < 5; l++)
            {
                var c = SplitChar(l < info.Count ? info[l] : null, ';');
                int colorend = l == 2 ? 27 : 1;
                if (l == 3)
                    colorend = 2;
                if (info.Count > 1 && info[1] == "-")
                {
                    if (l == 2)
                        colorend = 1;
                    else if (l == 4)
                        colorend = 27;
                }

                for (int k = 0; k < c.Count; k++)
                {
                    var e = Split(c[k]);
                    double x1, x2;
                    if (multi)
                    {
                        x1 = _startdrawX + Int(Num(At(e, 0)) / Xscale);
                        x2 = _startdrawX + Int(Num(At(e, 1)) / Xscale);
                    }
                    else
                    {
                        x1 = Int((Num(At(e, 0)) - _startcoord) / Xscale);
                        x2 = Int((Num(At(e, 1)) - _startcoord) / Xscale);
                    }

                    if (x1 == x2)
                        x2++;
                    if (tmpx2 != 0 && x1 - tmpx2 > 1)
                    {
                        Out("2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
                        Out("\t " + Str(tmpx2) + " " + Str(_y) + " " + Str(x1) + " " + Str(_y) + "\n");
                    }

                    tmpx2 = x2;
                    Out("2 1 0 5 " + colorend.ToString(CultureInfo.InvariantCulture) + " 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
                    Out("\t " + Str(x1) + " " + Str(_y) + " " + Str(x2) + " " + Str(_y) + "\n");
                }
            }

            _ = offset;
        }

        private void PrintGenesMr()
        {
            _hLastOnLevel.Clear();
            var g = SplitSpaces(_hContigGenes.GetValueOrDefault(Get(_a, -2)));
            for (int i = 1; i < g.Count; i++)
            {
                string kname = g[i];
                _y = Yorig - 100 - 200 * Num(Level(kname));
                double xid = _startdrawX + Int(Num(_utr5End.GetValueOrDefault(kname)) / Xscale);
                if (_o.PrintGeneNames)
                    PrintIdGenes(kname, xid);
                DrawGene(kname, _startdrawX, true);
            }
        }

        private void PrintIdGenes(string kname, double xid)
        {
            double printidspace;
            string level = Level(kname) ?? string.Empty;
            if (_hLastOnLevel.TryGetValue(level, out var last))
                printidspace = Int((Num(_utr5End.GetValueOrDefault(kname)) - Num(_utr5End.GetValueOrDefault(last))) / Xscale);
            else
                printidspace = 601;
            if (printidspace > 600)
            {
                Out("4 0 0 5 0 0 6 0.0000 4 90 270 ");
                Out("\t" + F0(xid) + " " + F0(_y - 50) + " ");
                Out(kname + "\\001\n");
                _hLastOnLevel[level] = kname;
            }
        }

        private void PrintGenes()
        {
            _hLastOnLevel.Clear();
            // sort { $Utr5End{$a} <=> $Utr5End{$b} } keys %Utr5End: the genes in the order of their 5' ends (equal ends are in the hash order of perl: here the order of first appearance)
            foreach (string kname in _utr5End.Keys.OrderBy(k => Num(_utr5End[k])))
            {
                double u5 = Num(_utr5End[kname]);
                if (u5 > _startcoord && u5 < _endcoord)
                {
                    _y = Yorig - 100 - 200 * Num(Level(kname));
                    double xid = Int((u5 - _startcoord) / Xscale);
                    if (_o.PrintGeneNames)
                        PrintIdGenes(kname, xid);
                    DrawGene(kname, 0, false);
                }
            }
        }

        private void PrintGrid(double xs, double xe, string startcontg, double endcontg)
        {
            _xlastPrint = 0;
            Out("2 1 0 10 11 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
            Out("\t " + Str(xs) + " " + F0(Yorig + 50) + " " + Str(xe) + " " + F0(Yorig + 50) + "\n");
            for (int percent = 50; percent < 101; percent += 25)
            {
                Out("2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n");
                string yy = F0(Yorig + 250 + YdistPid - (percent - 50) * Yscale);
                Out("\t" + Str(xs) + " " + yy + " " + Str(xe) + " " + yy + "\n");
            }

            double increment = 10000 / Xscale;
            int noIncr = 0;
            double xmark = xs, xmarkFloat = xs;
            while (xmark < xe)
            {
                Out("2 1 0 1 0 7 60 0 -1 0.000 0 0 -1 0 0 2\n");
                Out("\t" + Str(xmark) + " " + F0(Yorig + YdistPid + 250) + " " + Str(xmark) + " " + F0(Yorig + YdistPid + 300) + "\n");
                Out("4 0 0 100 0 0 8 0.0000 4 135 405 ");
                Out("\t " + F0(xmark) + " " + F0(Yorig + YdistPid + 400));
                Out(" " + noIncr.ToString(CultureInfo.InvariantCulture) + "k\\001\n");
                noIncr += 10;
                xmarkFloat += increment;
                xmark = Int(xmarkFloat);
            }

            Out("4 0 0 50 0 0 14 0.0000 4 135 450 " + Str(xs) + " " + Str(Yorig) + " " + startcontg + "\\001\n");
            Out("4 0 0 50 0 0 14 0.0000 4 135 810 " + F0(xe - Str(xe).Length * 125) + " " + Str(Yorig) + " " + Str(endcontg) + "\\001\n");
            if (_mref != 0)
            {
                Out("4 0 0 5 0 0 8 0.0000 4 135 405 ");
                Out("\t" + F0(xs) + " " + F0(Yorig + 70) + " ");
                Out(Get(_a, -2) + "\\001\n");
            }

            if (_mgaps is null)
            {
                for (int percent = 50; percent < 101; percent += 25)
                {
                    Out("4 0 0 100 0 0 8 0.0000 4 135 405 ");
                    Out("\t" + F0(xs - 200) + " " + F0(Yorig + YdistPid + 250 - (percent - 50) * Yscale + 20));
                    Out(" " + percent.ToString(CultureInfo.InvariantCulture) + "%\\001\n");
                    Out("4 0 0 100 0 0 8 0.0000 4 135 405 ");
                    Out("\t" + F0(xe + 20) + " " + F0(Yorig + YdistPid + 250 - (percent - 50) * Yscale + 20));
                    Out(" " + percent.ToString(CultureInfo.InvariantCulture) + "%\\001\n");
                }
            }
            else
            {
                Out("4 0 0 100 0 0 7 1.5710 4 135 405 ");
                Out("\t" + F0(xs - 50) + " " + F0(Yorig + YdistPid + 250 - 5 * Yscale + 10));
                Out(" + qry strand\\001\n");
                Out("4 0 0 100 0 0 7 1.5710 4 135 405 ");
                Out("\t" + F0(xs - 50) + " " + F0(Yorig + YdistPid + 250 - 30 * Yscale + 10));
                Out(" - qry strand\\001\n");
            }
        }

        private void PrintLegend()
        {
            Out("4 0 0 100 0 0 8 0.0000 4 135 405 ");
            Out("\t" + F0(100) + " " + F0(Yorig + YdistPid + 1100) + " ");
            Out(" Legend \\001\n");
            double y = Yorig + YdistPid + 1300;
            void Block(int color, string text, double yy)
            {
                Out("2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
                Out("\t 70 " + Str(yy) + " 99 " + Str(yy) + "\n");
                Out("2 1 0 5 " + color.ToString(CultureInfo.InvariantCulture) + " 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
                Out("\t 100 " + Str(yy) + " 200 " + Str(yy) + "\n");
                Out("2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
                Out("\t 200 " + Str(yy) + " 230 " + Str(yy) + "\n");
                Out("4 0 0 100 0 0 8 0.0000 4 135 405 ");
                Out("\t" + F0(300) + " " + F0(yy + 30) + " ");
                Out(text + "\\001\n");
            }

            Block(27, " 5' utr ", y);
            y += 150;
            Block(2, " cds ", y);
            y += 150;
            Block(1, " 3' utr ", y);
            y += 150;
            Out("2 1 0 2 4 0 50 0 -1 0.000 0 0 -1 0 0 2\n");
            Out("\t100 " + Str(y) + " 200 " + Str(y) + "\n");
            Out("4 0 0 100 0 0 8 0.0000 4 135 405 ");
            Out("\t" + F0(300) + " " + F0(y + 30) + " ");
            Out(" match found by " + _prog + " \\001\n");
        }

        /// <summary>Closes the current FIG file; with <c>pdf</c> / <c>ps</c> the script converts it with <c>fig2dev</c> and removes it.</summary>
        private void ChangeFileFormat()
        {
            if (_format == "fig")
                return;
            string name = _procfile + Str(_nrf) + ".fig";
            var m = Regex.Match(_procfile, "^tmp(.+)");
            string outfile = m.Groups[1].Value + Str(_nrf) + "." + _format;
            string cmd = "fig2dev -L " + _format + " -x 100";
            if (_o.Magnification != 0)
                cmd += " -m " + Str(_o.Magnification);
            cmd += " -M " + name + " " + outfile;
            _commands.Add(cmd);
            _outputs.Add(outfile);
            if (_o.Verbose)
            {
                _stdout.Append(outfile);
                if (_mref != 0)
                    _stdout.Append('\n');
                else
                    _stdout.Append("\t range : ").Append(_startText).Append('\t').Append(Str(_endcoord)).Append(" \n");
            }
        }

        // ------------------------------------------------------------------ mgaps input

        private void FormatMgaps()
        {
            string tmpfile2 = _o.CoordsName + "coords";
            // get_ref_len
            bool firstrow = true;
            double maxlenref = 0;
            var l = new Lines(_alignText);
            string? line;
            while ((line = l.Next()) is not null)
            {
                var a = Split(Chomp(line));
                if (firstrow)
                {
                    firstrow = false;
                    if (!Get(a, 0).StartsWith('>'))
                    {
                        _stdout.Append("\nWrong file format for MGAPS file : ").Append(_o.CoordsName).Append(" ! \n");
                        throw new MummerExitException(0);
                    }
                }

                if (Get(a, 0).StartsWith('>'))
                    continue;
                if (Get(a, 0) != "#")
                {
                    double e1 = Num(At(a, 0)) + Num(At(a, 2));
                    maxlenref = maxlenref < e1 ? e1 : maxlenref;
                }
            }

            var rows = new List<(double Key, string Line)>();
            string idquery = string.Empty;
            l = new Lines(_alignText);
            while ((line = l.Next()) is not null)
            {
                var a = Split(Chomp(line));
                if (Get(a, 0).StartsWith('>'))
                {
                    idquery = Get(a, 1) + (Get(a, 2) == "Reverse" ? "_rev" : string.Empty);
                }
                else if (Get(a, 0) != "#")
                {
                    double e1 = Num(At(a, 0)) + Num(At(a, 2));
                    string text = Get(a, 0) + "\t" + Str(e1) + "\t|" + "\t-\t-\t|" + "\t-\t-\t|" + "\t-\t|" + "\t" + Str(maxlenref) + "\t-\t|" + " " + _mgaps + "\t" + idquery + "\n";
                    rows.Add((Num(Get(a, 0)), text));
                }
            }

            // sort -n -k 1: numeric on the first field, ties by the whole line
            rows.Sort((x, y) =>
            {
                int c = x.Key.CompareTo(y.Key);
                return c != 0 ? c : string.CompareOrdinal(x.Line, y.Line);
            });
            var sb = new StringBuilder();
            sb.Append(_o.CoordsName).Append('\n').Append("Mgaps\n\n")
                .Append("    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n")
                .Append("===============================================================================================================================\n");
            foreach (var r in rows)
                sb.Append(r.Line);
            _alignText = sb.ToString();
            _side[tmpfile2] = _alignText;
            _stderr.Append("\n**************************************** \n");
            _stderr.Append("New input file created : ").Append(tmpfile2).Append('\n');
            _stderr.Append("**************************************** \n\n");
        }
    }
}
