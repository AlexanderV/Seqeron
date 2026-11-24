#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: the script mafft as a program (arguments -> variables -> pipeline -> output)

    /// <summary>
    /// One run of <c>core/mafft.tmpl</c>. The variables have the names and the values of the script (strings where the script echoes them or passes them on as text); the control flow is the script's, in the script's order, so that the
    /// first error message and the exit status are the script's.
    /// </summary>
    internal sealed partial class MafftScriptRun
    {
        private const string Version = "v7.505 (2022/Apr/10)";

        private readonly List<string> _argv;
        private readonly IReadOnlyDictionary<string, string> _files;
        private readonly string? _stdin;
        private readonly StringBuilder _stderr = new();
        private readonly StringBuilder _stdout = new();
        private readonly Dictionary<string, string> _outFiles = new(StringComparer.Ordinal);
        private int _pos;

        /// <summary>Signals <c>exit n</c> of the script.</summary>
        private sealed class ExitException(int code) : Exception
        {
            public int Code { get; } = code;
        }

        internal MafftScriptRun(IReadOnlyList<string> argv, IReadOnlyDictionary<string, string> files, string? stdin)
        {
            _argv = argv.Select(MafftProgramIo.ToBytes).ToList();
            var f = new Dictionary<string, string>(StringComparer.Ordinal);
            foreach (var kv in files)
                f[kv.Key] = MafftProgramIo.ToBytes(kv.Value);
            _files = f;
            _stdin = stdin == null ? null : MafftProgramIo.ToBytes(stdin);
        }

        // ---------------------------------------------------------------- the shell

        private int Remaining => _argv.Count - _pos;
        private string Cur => _pos < _argv.Count ? _argv[_pos] : "";
        private void Shift()
        {
            if (_pos < _argv.Count)
                _pos++;
        }

        private void Err(string s) => _stderr.Append(s);

        private void Exit(int code) => throw new ExitException(code);

        /// <summary><c>expr STRING : "[0-9]"</c> (the pattern is anchored at the start of the string).</summary>
        private static bool StartsWithDigit(string s) => s.Length > 0 && s[0] >= '0' && s[0] <= '9';

        /// <summary><c>expr STRING : "[0-9\-]"</c>.</summary>
        private static bool StartsWithDigitOrDash(string s) => s.Length > 0 && (char.IsAsciiDigit(s[0]) || s[0] == '-' || s[0] == '\\');

        /// <summary>
        /// <c>expr "$1" - 0</c>: the integer, or null with the message <c>expr: non-integer argument</c> on standard error (the variable is then empty and the script fails later with messages that name its line numbers).
        /// </summary>
        private int? ExprInt(string s)
        {
            int i = 0;
            if (s.Length > 0 && s[0] == '-')
                i = 1;
            if (s.Length > i && s.Skip(i).All(char.IsAsciiDigit) && long.TryParse(s, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long v) && v >= int.MinValue && v <= int.MaxValue)
                return (int)v;
            Err("expr: non-integer argument\n");
            return null;
        }

        private int ExprIntOrFail(string s)
        {
            return ExprInt(s) ?? throw new NotSupportedException($"'{s}' is not an integer: the script's variable is then empty and bash goes on with error messages that carry the line numbers of the script and runs the programs with garbage (not ported).");
        }

        /// <summary><c>v=`expr "$1" - 0`</c> followed by <c>if ! expr "$1" : "[0-9]"</c>: the message and <c>exit</c> (status 0) when the argument does not start with a digit.</summary>
        private int ExprThenDigitCheck(string arg, string message)
        {
            int? v = ExprInt(arg);
            if (!StartsWithDigit(arg))
            {
                Err(message);
                Exit(0);
            }
            return v ?? ExprIntOrFail(arg);
        }

        /// <summary>awk <c>0.0 + "str"</c>: the numeric prefix of the string.</summary>
        private static double Num(string s) => MafftGlobals.CAtof(s);

        // ---------------------------------------------------------------- the variables of the script

        private string _outputfile = "";
        private int _namelength = -1;
        private bool _anysymbol;
        private string _parallelizationstrategy = "BAATARI2";
        private string _kappa = " ";
        private string _sbstmodel = " -b 62 ";
        private string _fmodel = " ";
        private string _nmodel = " ";
        private string _gexp = "0";
        private string _gop = "1.53";
        private string _gopdist = "1.53";
        private string _aof = "0.000";
        private int _cycle = 2;
        /// <summary>The number of tree builds as it stands when the script decides about <c>-T</c> (<c>--retree 0</c>), before the number of sequences resets it.</summary>
        private int _retreeUser = 2;
        private int _iterate;
        private bool _fft = true;
        private string _distance = "ktuples";
        private bool _forcefft;
        private string _memopt = " ";
        private string _weightopt = " ";
        private string _GGOP = "-6.00", _LGOP = "-6.00", _LEXP = "-0.000", _GEXP = "-0.000";
        private string _lgop = "-2.00", _lexp = "-0.100", _laof = "0.100";
        private string _pggop = "-2.00", _pgexp = "-0.10", _pgaof = "0.10";
        private string _rgop = "-1.530", _rgep = "-0.000";
        private string _seqtype = "  ";
        private string _weighti = "2.7", _weightr = "0.0", _weightm = "1.0";
        private bool _rnaalifold, _dafs, _mccaskill, _contrafold;
        private string _progressfile = "/dev/stderr";
        private string _anchorfile = "/dev/null";
        private string _anchoropt = "";
        private int _maxanchorseparation = 1000;
        private bool _sw;
        private string _algopt = "  ", _algoptit = "  ";
        private bool _pairspecified;
        private string _scorecalcopt = " ";
        private bool _coreout;
        private string _corethr = "0.5", _corewin = "100", _coreext = " ";
        private string _outputformat = "pir";
        private string _f2clext = "-N";
        private string _outorder = "input";
        private string _seed = "x", _seedtable = "x";
        private bool _auto;
        private int _groupsize = -1, _partsize = 50;
        private string _partdist = "ktuples";
        private string _partorderopt = " -x ";
        private bool _treeout, _nodeout, _distout;
        private string _treeinopt = "  ";
        private bool _treein;
        private string _seedfiles = "/dev/null";
        private string _seedtablefile = "/dev/null";
        private string _pdblist = "/dev/null", _ownlist = "/dev/null";
        private string _scorematrix = "/dev/null";
        private string _treeinfile = "/dev/null";
        private string _codonposfile = "/dev/null", _codonscorefile = "/dev/null";
        private string _rnascoremtx = " ";
        private string _laraparams = "/dev/null";
        private string _treealg = " -X 0.1 ";
        private string _sueff = "1.0";
        private string _maxambiguous = "1.0";
        private bool _dofilter;
        private string _scoreoutarg = " ";
        private int _numthreads, _numthreadsit = -1, _numthreadstb = -1;
        private int _randomseed;
        private string _addfile = "/dev/null";
        private string _addarg0 = " ", _addarg = " ", _addsinglearg = " ", _add2ndhalfarg = " ";
        private string _mapoutfile = "/dev/null";
        private int _fragment;
        private string _legacygapopt = " ";
        private string _mergetable = "/dev/null";
        private string _outnum = " ";
        private int _lastE = 5000, _lastM = 3;
        private string _lastSubopt = " ", _lastOnce = " ";
        private int _adjustdirection;
        private int _tuplesize = 6;
        private string _termgapopt = " -O ";
        private string _similarityoffset = "0.0";
        private string _unalignlevel = "0.0";
        private bool _unalignspecified;
        private string _spfactor = "100.0";
        private bool _shiftpenaltyspecified, _opdistspecified, _allowshift;
        private bool _newdash;
        private bool _exclude_ho;
        private string _fixthreshold = "0.0";
        private string _bunkatsuopt = " ";
        private int _npickup;
        private string _minimumweight = "0.00001";
        private string _usenaivepairscore = " ";
        private bool _oldgenafparam;
        private bool _sprigorous;
        private string _treeext = "none";
        private string _focusarg = " ";
        private string _lhlimit = " ";
        private string _oneiterationopt = "";
        private bool _mpi;
        private bool _debug;
        private int _er;
        private int _nadd;
        private string _infilename = "";

        /// <summary>The options whose binaries are not ported yet (the first one is named in the refusal).</summary>
        private readonly List<string> _portNeeded = new();

        private void NeedPort(string what)
        {
            if (!_portNeeded.Contains(what))
                _portNeeded.Add(what);
        }

        // ---------------------------------------------------------------- the run

        internal MafftProgramResult Execute()
        {
            int rc;
            try
            {
                rc = Script();
            }
            catch (ExitException e)
            {
                rc = e.Code;
            }
            var files = new Dictionary<string, string>(StringComparer.Ordinal);
            foreach (var kv in _outFiles)
                files[kv.Key] = MafftProgramIo.FromBytes(kv.Value);
            return new MafftProgramResult(rc, MafftProgramIo.FromBytes(_stdout.ToString()), MafftProgramIo.FromBytes(_stderr.ToString()), files);
        }

        /// <summary>Progress text (<c>1&gt;&gt;"$progressfile"</c>): standard error, nothing with <c>--quiet</c>, a file with <c>--progress</c>.</summary>
        private void Prog(string s)
        {
            if (_progressfile == "/dev/stderr")
                _stderr.Append(s);
            else if (_progressfile == "/dev/null")
                return;
            else
                _outFiles[_progressfile] = (_outFiles.TryGetValue(_progressfile, out var cur) ? cur : "") + s;
        }

        private int Script()
        {
            if (Remaining == 0)
                throw new NotSupportedException("mafft without arguments starts an interactive wizard that reads the terminal.");
            if (Cur == "--man")
                throw new NotSupportedException("--man shows the manual page of the installation (man mafft.1).");
            if (Cur == "--version")
            {
                Err(Version + "\n");
                return 0;
            }
            if (Cur is "--help" or "--info")
            {
                Shift();
                _er = 1;
            }
            ParseOptions();
            // the temporary directory (mktemp) is not modelled
            string infile = "";
            string addfileText = "";
            if (Remaining == 1)
            {
                string input = Cur;
                if (Readable(input) || input == "-")
                {
                    if (!Readable(_addfile))
                    {
                        Err(MafftProgram.ProgramPath + ": Cannot open " + _addfile + ".\n\n");
                        return 1;
                    }
                    string src = input == "-" ? (_stdin ?? "") : ReadFile(input);
                    infile = MafftProgramIo.CarriageReturnsToNewlines(src) + "\n";
                    addfileText = MafftProgramIo.DropEmptyLines(MafftProgramIo.CarriageReturnsToNewlines(ReadFile(_addfile)));
                    if (_dofilter)
                        addfileText = FilterAmbiguous(addfileText);
                    infile += addfileText;
                    _infilename = input;
                    foreach (var f in new[] { _scorematrix, _mergetable, _treeinfile, _codonposfile, _codonscorefile, _seedtablefile, _laraparams, _pdblist, _ownlist, _anchorfile })
                        _ = f;
                    var seedList = _seedfiles.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    for (int s = 1; s < seedList.Length; s++)
                    {
                        if (!Readable(seedList[s]))
                        {
                            Err(MafftProgram.ProgramPath + ": Cannot open " + seedList[s] + ".\n\n");
                            return 1;
                        }
                    }
                }
                else
                {
                    Err(MafftProgram.ProgramPath + ": Cannot open " + input + ".\n\n");
                    _er = 1;
                }
            }
            else
            {
                _er = 1;
            }
            return Rest(infile, addfileText);
        }

        /// <summary>
        /// <c>filter -m maxambiguous</c> (<c>filter.c</c>) on the <c>--add</c> file: the sequences whose share of letters outside the usual alphabet is above the limit are removed, runs of the unknown letter are shortened to one, the rest is written
        /// on one line (gaps removed). The program writes "maxunusual = " even for an empty file.
        /// </summary>
        private string FilterAmbiguous(string text)
        {
            double maxunusual = MafftGlobals.CAtof(_maxambiguous);
            Prog("maxunusual = " + MafftGlobals.F(maxunusual, 0, 6) + "\n");
            var records = MafftProgramIo.Records(text);
            var seqs = new List<string>();
            foreach (var r in records)
            {
                string raw = r.Raw;
                if (raw.IndexOfAny(new[] { '=', '<', '>' }) >= 0)
                {
                    Prog("\nCharacters '= < >' can be used only in the title lines in the --anysymbol or --text mode.\n\n");
                    Exit(1);
                }
                seqs.Add(string.Concat(raw.Where(ch => ch is not ('\n' or ' ' or '\r'))));
            }
            bool protein = _seqtype == "-P" || (_seqtype != "-D" && !MafftEngine.GuessNucleotide(seqs, new MafftOptions()));
            string usual = protein ? "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-" : "ATGCUatgcu-";
            char unknown = protein ? 'X' : 'n';
            var sb = new StringBuilder();
            int kept = 0;
            for (int i = 0; i < records.Count; i++)
            {
                string t = seqs[i].Replace("-", "", StringComparison.Ordinal);
                int count = 0;
                foreach (char ch in t)
                    if (usual.IndexOf(ch) < 0)
                        count++;
                double frac = t.Length == 0 ? double.NaN : (double)count / t.Length;
                if (!(frac <= maxunusual))
                    continue;
                var o = new StringBuilder();
                bool status = false;
                foreach (char ch in t)
                {
                    if (char.ToUpperInvariant(ch) != char.ToUpperInvariant(unknown))
                    {
                        o.Append(ch);
                        status = false;
                    }
                    else if (!status)
                    {
                        o.Append(unknown);
                        status = true;
                    }
                }
                sb.Append('>').Append(records[i].Name).Append('\n').Append(o).Append('\n');
                kept++;
            }
            if (kept < records.Count)
                Prog("\n\nRemoved " + (records.Count - kept).ToString(CultureInfo.InvariantCulture) + " sequence(s) where the frequency of ambiguous " + (protein ? "amino acids" : "bases") + " > " + MafftGlobals.F(maxunusual, 5, 3) + "\n\n\n");
            return sb.ToString();
        }

        private bool Readable(string name) => name == "/dev/null" || _files.ContainsKey(name);

        private string ReadFile(string name) => name == "/dev/null" ? "" : _files[name];

        // ---------------------------------------------------------------- the option loop (mafft.tmpl 238 - 1036)

        private void ParseOptions()
        {
            while (Remaining > 1)
            {
                string o = Cur;
                switch (o)
                {
                    case "--auto": _auto = true; break;
                    case "--anysymbol":
                    case "--preservecase": _anysymbol = true; break;
                    case "--clustalout": _outputformat = "clustal"; break;
                    case "--phylipout": _outputformat = "phylip"; break;
                    case "--reorder": _outorder = "aligned"; _partorderopt = " "; break;
                    case "--inputorder": _outorder = "input"; _partorderopt = " -x "; break;
                    case "--unweight": _weightopt = " -u "; break;
                    case "--termgappenalty": _termgapopt = " "; break;
                    case "--alga": _algopt = " "; _algoptit = " "; break;
                    case "--algq":
                        Err("\n--algq is no longer supported!\n\n");
                        Exit(1);
                        break;
                    case "--namelength":
                        Shift();
                        _namelength = ExprThenDigitCheck(Cur, "Specify the length of name in clustal format output!\n");
                        break;
                    case "--groupsize":
                        Shift();
                        _groupsize = ExprThenDigitCheck(Cur, "Specify groupsize!\n");
                        break;
                    case "--partsize":
                        Shift();
                        _partsize = ExprThenDigitCheck(Cur, "Specify partsize!\n");
                        break;
                    case "--parttree": _distance = "parttree"; _partdist = "ktuples"; break;
                    case "--dpparttree": _distance = "parttree"; _partdist = "localalign"; break;
                    case "--fastaparttree": _distance = "parttree"; _partdist = "fasta"; break;
                    case "--treeout": _treeout = true; break;
                    case "--nodeout": _nodeout = true; _treeout = true; break;
                    case "--distout": _distout = true; break;
                    case "--fastswpair": _distance = "fasta"; _pairspecified = true; _sw = true; break;
                    case "--fastapair": _distance = "fasta"; _pairspecified = true; _sw = false; break;
                    case "--averagelinkage": _treealg = " -X 1.0 "; _sueff = "1.0"; break;
                    case "--minimumlinkage": _treealg = " -X 0.0 "; _sueff = "0.0"; break;
                    case "--mixedlinkage": Shift(); _sueff = Cur; _treealg = " -X " + Cur; break;
                    case "--maxambiguous": Shift(); _maxambiguous = Cur; _dofilter = true; break;
                    case "--codonpos":
                        Shift();
                        _codonposfile = Cur;
                        if (!_files.ContainsKey(_codonposfile))
                        {
                            Err("Cannot open " + _codonposfile + "\n\n");
                            Exit(0);
                        }
                        NeedPort("--codonpos");
                        break;
                    case "--codonscore":
                        Shift();
                        _codonscorefile = Cur;
                        if (!_files.ContainsKey(_codonscorefile))
                        {
                            Err("Cannot open " + _codonscorefile + "\n\n");
                            Exit(0);
                        }
                        NeedPort("--codonscore");
                        break;
                    case "--noscore": _scorecalcopt = " -Z "; break;
                    case "--6mermultipair": _distance = "ktuplesmulti"; _tuplesize = 6; _pairspecified = true; break;
                    case "--10mermultipair": _distance = "ktuplesmulti"; _tuplesize = 10; _pairspecified = true; break;
                    case "--6merpair": _distance = "ktuples"; _tuplesize = 6; _pairspecified = true; break;
                    case "--10merpair": _distance = "ktuples"; _tuplesize = 10; _pairspecified = true; break;
                    case "--blastpair": _distance = "blast"; _pairspecified = true; break;
                    case "--lastmultipair": _distance = "lastmulti"; _pairspecified = true; break;
                    case "--globalpair": _distance = "global"; _pairspecified = true; break;
                    case "--shortlongpair":
                    case "--longshortpair":
                        _distance = "local"; _usenaivepairscore = "-Z"; _laof = "0.0"; _lexp = "0.0"; _pgaof = "0.0"; _pgexp = "0.0"; _pairspecified = true;
                        break;
                    case "--localpair": _distance = "local"; _pairspecified = true; break;
                    case "--lastpair": _distance = "last"; _pairspecified = true; break;
                    case "--multipair": _distance = "multi"; _pairspecified = true; break;
                    case "--hybridpair": _distance = "hybrid"; _pairspecified = true; break;
                    case "--scarnapair": _distance = "scarna"; _pairspecified = true; break;
                    case "--dafspair": _distance = "dafs"; _pairspecified = true; break;
                    case "--larapair": _distance = "lara"; _pairspecified = true; break;
                    case "--slarapair": _distance = "slara"; _pairspecified = true; break;
                    case "--foldalignpair":
                    case "--foldalignlocalpair": _distance = "foldalignlocal"; _pairspecified = true; break;
                    case "--foldalignglobalpair": _distance = "foldalignglobal"; _pairspecified = true; break;
                    case "--globalgenafpair":
                        _distance = "globalgenaf"; _pairspecified = true;
                        Err("\n--globalgenaf is no longer supported!\n\n");
                        Exit(1);
                        break;
                    case "--localgenafpair":
                    case "--genafpair": _distance = "localgenaf"; _pairspecified = true; break;
                    case "--oldgenafpair": _distance = "localgenaf"; _pairspecified = true; _oldgenafparam = true; break;
                    case "--memsave": _memopt = " -M -B "; break;
                    case "--nomemsave": _memopt = " -N "; break;
                    case "--nuc": _seqtype = "-D"; break;
                    case "--amino": _seqtype = "-P"; break;
                    case "--fft": _fft = true; _forcefft = true; break;
                    case "--nofft": _fft = false; break;
                    case "--quiet": _progressfile = "/dev/null"; break;
                    case "--debug": _debug = true; break;
                    case "--coreext": _coreext = " -c "; break;
                    case "--core": _coreout = true; break;
                    case "--adjustdirection": _adjustdirection = 1; break;
                    case "--adjustdirectionaccurately": _adjustdirection = 2; break;
                    case "--oneiteration": _oneiterationopt = " -r "; break;
                    case "--progress":
                        Shift();
                        _progressfile = Cur;
                        if (!(Cur.StartsWith('/') || (Cur.Length >= 2 && char.IsAsciiLetter(Cur[0]) && Cur[1] == ':')))
                        {
                            Err("Specify a progress file name with the absolute path!\n");
                            Exit(0);
                        }
                        break;
                    case "--out": Shift(); _outputfile = Cur; break;
                    case "--skipanchorsremoterthan":
                        Shift();
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify maximum gap length between anchors.\n");
                            Exit(0);
                        }
                        _maxanchorseparation = ExprIntOrFail(Cur);
                        break;
                    case "--anchors":
                        Shift();
                        _anchorfile = Cur;
                        _anchoropt = " -l ";
                        if (!_files.ContainsKey(_anchorfile))
                        {
                            Err("Cannot open " + _anchorfile + "\n\n");
                            Exit(0);
                        }
                        NeedPort("--anchors");
                        break;
                    case "--thread":
                        Shift();
                        if (!StartsWithDigitOrDash(Cur))
                        {
                            Err("Specify the number of threads.  Or, use --thread -1\n");
                            Exit(0);
                        }
                        _numthreads = ExprIntOrFail(Cur);
                        break;
                    case "--threadtb":
                        Shift();
                        if (!StartsWithDigitOrDash(Cur))
                        {
                            Err("Check the argument after --threadtb, the number of threads for the progressive step.\n");
                            Exit(0);
                        }
                        _numthreadstb = ExprIntOrFail(Cur);
                        break;
                    case "--threadit":
                        Shift();
                        if (!StartsWithDigitOrDash(Cur))
                        {
                            Err("Check the argument after --threadit, the number of threads for the iterative step.\n");
                            Exit(0);
                        }
                        _numthreadsit = ExprIntOrFail(Cur);
                        break;
                    case "--last_subopt": _lastSubopt = "-S"; break;
                    case "--last_once": _lastOnce = "-U"; break;
                    case "--last_m": Shift(); _lastM = ExprIntOrFail(Cur); break;
                    case "--last_e": Shift(); _lastE = ExprIntOrFail(Cur); break;
                    case "--randomseed": Shift(); _randomseed = ExprIntOrFail(Cur); break;
                    case "--bestfirst": _parallelizationstrategy = "BESTFIRST"; break;
                    case "--adhoc0": _parallelizationstrategy = "BAATARI0"; break;
                    case "--adhoc1": _parallelizationstrategy = "BAATARI1"; break;
                    case "--adhoc2":
                    case "--simplehillclimbing": _parallelizationstrategy = "BAATARI2"; break;
                    case "--scoreout": _scoreoutarg = "-S -B"; break;
                    case "--outnum": _outnum = "-n"; break;
                    case "--leavegappyregion":
                    case "--legacygappenalty": _legacygapopt = "-L"; break;
                    case "--merge":
                        Shift();
                        _mergetable = Cur;
                        if (!_files.ContainsKey(_mergetable))
                        {
                            Err("Cannot open " + _mergetable + "\n\n");
                            Exit(0);
                        }
                        break;
                    case "--addprofile": Shift(); _addarg0 = "-I"; _addfile = Cur; break;
                    case "--add": Shift(); _addarg0 = "-K -I"; _addfile = Cur; break;
                    case "--addfragments": Shift(); _addarg0 = "-K -I"; _addfile = Cur; _fragment = 1; break;
                    case "--addfull": Shift(); _addarg0 = "-K -I"; _addfile = Cur; _fragment = -1; break;
                    case "--addlong": Shift(); _addarg0 = "-K -I"; _addfile = Cur; _fragment = -2; break;
                    case "--addtotop": Shift(); _addarg0 = "-K -I"; _addfile = Cur; _fragment = -3; break;
                    case "--smoothing": _add2ndhalfarg += " -p "; break;
                    case "--keeplength": _add2ndhalfarg += " -Y "; break;
                    case "--compactmapout": _add2ndhalfarg += " -z -Y "; break;
                    case "--compactmapoutfile": Shift(); _add2ndhalfarg += " -z -Y "; _mapoutfile = Cur; break;
                    case "--mapout": _add2ndhalfarg += " -Z -Y "; break;
                    case "--mapoutfile": Shift(); _add2ndhalfarg += " -Z -Y "; _mapoutfile = Cur; break;
                    case "--maxiterate":
                        Shift();
                        _iterate = ExprThenDigitCheck(Cur, "Specify the number of iterations!\n");
                        break;
                    case "--retree":
                        Shift();
                        _cycle = ExprThenDigitCheck(Cur, "Specify the number of tree rebuilding!\n");
                        break;
                    case "--aamatrix":
                        Shift();
                        _f2clext = "-N";
                        _sbstmodel = " -b -1 ";
                        _scorematrix = Cur;
                        if (!_files.ContainsKey(_scorematrix))
                        {
                            Err("Cannot open " + _scorematrix + "\n\n");
                            Exit(0);
                        }
                        break;
                    case "--textmatrix":
                        Shift();
                        _f2clext = "-E"; _seqtype = "-P"; _fft = false; _sbstmodel = " -b -2 -a ";
                        _scorematrix = Cur;
                        if (!_files.ContainsKey(_scorematrix))
                        {
                            Err("Cannot open " + _scorematrix + "\n\n");
                            Exit(0);
                        }
                        NeedPort("--textmatrix");
                        break;
                    case "--text":
                        _f2clext = "-E"; _seqtype = "-P"; _fft = false; _sbstmodel = " -b -2 -a ";
                        NeedPort("--text");
                        break;
                    case "--treein":
                        Shift();
                        _treeinopt = " -U "; _treein = true; _treeinfile = Cur;
                        if (!_files.ContainsKey(_treeinfile))
                        {
                            Err("Cannot open " + _treeinfile + "\n\n");
                            Exit(0);
                        }
                        break;
                    case "--pileup": _treeext = "pileup"; break;
                    case "--randomchain": _treeext = "randomchain"; break;
                    case "--topin":
                        Shift();
                        Err("The --topin option has been disabled.\nThere was a bug in version < 6.530.\nThis bug has not yet been fixed.\n");
                        Exit(1);
                        break;
                    case "--mpi": _mpi = true; break;
                    case "--large":
                    case "--memsavetree": _treeext = "memsavetree"; break;
                    case "--memsavetreex": _treeext = "memsavetreex"; break;
                    case "--stepadd": _treeext = "stepadd"; break;
                    case "--youngestlinkage": _treeext = "youngestlinkage"; break;
                    case "--kappa":
                        Shift();
                        _kappa = " -k " + Cur + " ";
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify kappa value!\n");
                            Exit(0);
                        }
                        break;
                    case "--fmodel": _fmodel = " -a "; break;
                    case "--nwildcard": _nmodel = " -: "; break;
                    case "--nzero": _nmodel = "  "; break;
                    case "--jtt":
                    case "--kimura": Shift(); _f2clext = "-N"; _sbstmodel = " -j " + Cur; break;
                    case "--tm": Shift(); _f2clext = "-N"; _sbstmodel = " -m " + Cur; break;
                    case "--bl":
                        Shift();
                        _f2clext = "-N";
                        _sbstmodel = " -b " + Cur;
                        if (!StartsWithDigit(Cur))
                        {
                            Err("blosum " + Cur + "?\n");
                            Exit(0);
                        }
                        break;
                    case "--weighti":
                        Shift();
                        _weighti = Cur;
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify weighti value!\n");
                            Exit(0);
                        }
                        break;
                    case "--weightr":
                        Shift();
                        _weightr = Cur;
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify weightr value!\n");
                            Exit(0);
                        }
                        break;
                    case "--weightm":
                        Shift();
                        _weightm = Cur;
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify weightm value!\n");
                            Exit(0);
                        }
                        break;
                    case "--rnaalifold": _rnaalifold = true; break;
                    case "--mccaskill": _mccaskill = true; _contrafold = false; _dafs = false; break;
                    case "--contrafold": _mccaskill = false; _contrafold = true; _dafs = false; break;
                    case "--dafs": _mccaskill = false; _contrafold = false; _dafs = true; break;
                    case "--ribosum": _rnascoremtx = " -s "; break;
                    case "--op":
                        Shift();
                        _gop = Cur;
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify op!\n");
                            Exit(0);
                        }
                        break;
                    case "--opdist":
                        Shift();
                        _gopdist = Cur;
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify opdist!\n");
                            Exit(0);
                        }
                        _opdistspecified = true;
                        break;
                    case "--allowshift": _allowshift = true; break;
                    case "--shiftpenalty":
                        Shift();
                        _spfactor = Cur;
                        if (!StartsWithDigit(Cur))
                        {
                            Err("Specify sf!\n");
                            Exit(0);
                        }
                        _shiftpenaltyspecified = true;
                        break;
                    case "--exp":
                        Shift();
                        _gexp = MafftScript.AwkPrint(-1.0 * Num(Cur));
                        if (!StartsWithDigitOrDash(_gexp))
                        {
                            Err("\nSpecify a number for exp, like --exp 0.1\n'" + Cur + "' cannot be interpreted as a number..\n\n");
                            Exit(0);
                        }
                        break;
                    case "--ep":
                        Shift();
                        _aof = MafftScript.AwkPrint(-1.0 * Num(Cur));
                        if (!StartsWithDigitOrDash(_aof))
                        {
                            Err("\nSpecify a number for ep, like --ep 0.1\n'" + Cur + "' cannot be interpreted as a number..\n\n");
                            Exit(0);
                        }
                        break;
                    case "--rop": Shift(); _rgop = Cur; break;
                    case "--rep": Shift(); _rgep = Cur; break;
                    case "--lop": Shift(); _lgop = Cur; break;
                    case "--LOP": Shift(); _LGOP = Cur; break;
                    case "--lep": Shift(); _laof = Cur; break;
                    case "--lexp": Shift(); _lexp = Cur; break;
                    case "--LEXP": Shift(); _LEXP = Cur; break;
                    case "--GEXP": Shift(); _GEXP = Cur; break;
                    case "--GOP": Shift(); _GGOP = Cur; break;
                    case "--gop": Shift(); _pggop = Cur; break;
                    case "--gep": Shift(); _pgaof = Cur; break;
                    case "--gexp": Shift(); _pgexp = Cur; break;
                    case "--laraparams": Shift(); _laraparams = Cur; break;
                    case "--corethr": Shift(); _corethr = Cur; break;
                    case "--corewin": Shift(); _corewin = Cur; break;
                    case "--strdir": Shift(); break;
                    case "--pdbidlist":
                        Err("--pdbidlist is temporarily unavailable, 2018/Dec.\n\n");
                        Exit(0);
                        break;
                    case "--pdbfilelist":
                        Err("--pdbfilelist is temporarily unavailable, 2018/Dec.\n\n");
                        Exit(0);
                        break;
                    case "--dash":
                        _newdash = true;
                        if (_distance != "local" && _distance != "localgenaf")
                            _distance = "global";
                        if (_iterate == 0)
                            _iterate = 3;
                        break;
                    case "--dashserver": Shift(); break;
                    case "--originalseqonly": break;
                    case "--excludehomologs": _exclude_ho = true; break;
                    case "--seedtable": Shift(); _seedtable = "y"; _seedtablefile = Cur; NeedPort("--seedtable"); break;
                    case "--seed": Shift(); _seed = "m"; _seedfiles += " " + Cur; NeedPort("--seed"); break;
                    case "--minimumweight": Shift(); _minimumweight = Cur; break;
                    case "--similaritylevel": Shift(); _similarityoffset = Cur; break;
                    case "--unalignlevel": Shift(); _unalignlevel = Cur; _unalignspecified = true; break;
                    case "--skipiterate": Shift(); _fixthreshold = Cur; break;
                    case "--bunkatsunashi": _bunkatsuopt = " -B "; break;
                    case "--sp": _sprigorous = true; break;
                    case "--focus": _focusarg = " -= "; break;
                    case "--lhlimit": Shift(); _lhlimit = " -q " + Cur + " "; break;
                    case "--sparsepickup": Shift(); _npickup = ExprIntOrFail(Cur); break;
                    default:
                        Err("Unknown option:  " + o + "\n");
                        _er = 1;
                        break;
                }
                Shift();
            }
        }
    }

    #endregion
}
