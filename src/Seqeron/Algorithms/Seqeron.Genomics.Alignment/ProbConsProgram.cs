using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The meaning of a ProbCons 1.12 command line (<see cref="ProbConsProgram.Parse"/>; <c>Main.cc</c> <c>ParseParams</c>).
/// </summary>
/// <param name="Options">Alignment parameters (<c>-c -ir -pre -a -go -ge -co -v</c>).</param>
/// <param name="InputFiles">The non-option arguments.</param>
/// <param name="Train">The <c>-t FILE</c> parameter output file (null = no training).</param>
/// <param name="TrainEmissions"><c>-e</c>.</param>
/// <param name="ParameterFile"><c>-p FILE</c> (null = built-in parameters).</param>
/// <param name="AnnotationFile"><c>-annot FILE</c>.</param>
/// <param name="AllPairs"><c>-pairs</c> (also set by <c>-viterbi</c>).</param>
/// <param name="Viterbi"><c>-viterbi</c>.</param>
/// <param name="ClustalW"><c>-clustalw</c>.</param>
public sealed record ProbConsCommand(ProbConsOptions Options, IReadOnlyList<string> InputFiles, string? Train, bool TrainEmissions, string? ParameterFile, string? AnnotationFile,
    bool AllPairs, bool Viterbi, bool ClustalW);

/// <summary>What a run of <see cref="ProbConsProgram.Run"/> produced.</summary>
/// <param name="ExitCode">0 = success, 1 = the program stopped with an error.</param>
/// <param name="Stdout">The alignment (MFA or Clustal format) or nothing.</param>
/// <param name="Stderr">Everything the program writes to the error stream: heading, parameter set, "Loading sequence file", progress, tree, errors.</param>
/// <param name="Files">Files the program wrote: pairwise alignments (<c>name1-name2.fasta</c> / <c>.aln</c>), the annotation file, the trained parameter file.</param>
public sealed record ProbConsProgramResult(int ExitCode, string Stdout, string Stderr, IReadOnlyDictionary<string, string> Files);

/// <summary>
/// ProbCons 1.12 as a command line (<c>probcons [-clustalw] [-c N] [-ir N] [-pre N] [-pairs] [-viterbi] [-v] [-annot F] [-t F] [-e] [-p F] [-a] [-go X] [-ge Y] [-co C] file...</c>): the argument
/// parser of <c>Main.cc</c> (exact option spellings, ranges and error texts, <c>strtol</c> / <c>strtod</c> number reading), the MFA / MSF / Clustal readers of <c>MultiSequence.h</c>
/// (<c>LoadMFA</c>, <c>ParseMSF</c>), the MFA and Clustal writers (<c>WriteMFA</c>, <c>WriteALN</c> with the conservation line) and the program's standard-error text. Input and output files are text
/// given as a dictionary (a missing input name is the program's "Could not open file" error).
/// </summary>
public static class ProbConsProgram
{
    private const string Version = "1.12";

    // ---------- number reading ----------

    /// <summary><c>strtol (data, &amp;end, 0)</c>: optional blanks and sign, <c>0x</c> hexadecimal, leading <c>0</c> octal, else decimal; trailing characters are ignored (<c>GetInteger</c>).</summary>
    internal static bool GetInteger(string data, out int value)
    {
        value = 0;
        int i = 0;
        while (i < data.Length && (data[i] == ' ' || (data[i] >= '\t' && data[i] <= '\r'))) i++;
        bool neg = false;
        if (i < data.Length && (data[i] == '+' || data[i] == '-'))
        {
            neg = data[i] == '-';
            i++;
        }

        int bas = 10;
        if (i + 2 < data.Length && data[i] == '0' && (data[i + 1] == 'x' || data[i + 1] == 'X') && Uri.IsHexDigit(data[i + 2]))
        {
            bas = 16;
            i += 2;
        }
        else if (i < data.Length && data[i] == '0')
        {
            bas = 8;
        }

        System.Numerics.BigInteger v = 0;
        int digits = 0;
        while (i < data.Length)
        {
            int d = DigitValue(data[i]);
            if (d >= bas) break;
            v = v * bas + d;
            digits++;
            i++;
        }

        if (digits == 0)
            return false; // no conversion performed: data == endPtr and the value is 0
        if (neg) v = -v;
        if (v < int.MinValue || v > int.MaxValue) return false;
        value = (int)v;
        return true;
    }

    private static int DigitValue(char c)
    {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        return 99;
    }

    /// <summary><c>strtod</c> of the longest decimal prefix (<c>GetFloat</c>); false when nothing was converted.</summary>
    internal static bool GetFloat(string data, out float value)
    {
        value = 0;
        var m = System.Text.RegularExpressions.Regex.Match(data, @"^[ \t\n\v\f\r]*[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?");
        if (!m.Success) return false;
        double d = double.Parse(m.Value.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture);
        if (double.IsInfinity(d)) return false; // errno = ERANGE and |value| >= 1e6: rejected; a large finite value is accepted

        value = (float)d;
        return true;
    }

    // ---------- parser ----------

    private const string Heading = "\nPROBCONS version " + Version + " - align multiple protein sequences and print to standard output\nWritten by Chuong Do\n\n";

    private const string Usage =
        "PROBCONS comes with ABSOLUTELY NO WARRANTY.  This is free software, and\nyou are welcome to redistribute it under certain conditions.  See the\nfile COPYING.txt for details.\n\n" +
        "Usage:\n       probcons [OPTION]... [MFAFILE]...\n\nDescription:\n       Align sequences in MFAFILE(s) and print result to standard output\n\n" +
        "       -clustalw\n              use CLUSTALW output format instead of MFA\n\n" +
        "       -c, --consistency REPS\n              use 0 <= REPS <= 5 (default: 2) passes of consistency transformation\n\n" +
        "       -ir, --iterative-refinement REPS\n              use 0 <= REPS <= 1000 (default: 100) passes of iterative-refinement\n\n" +
        "       -pre, --pre-training REPS\n              use 0 <= REPS <= 20 (default: 0) rounds of pretraining\n\n" +
        "       -pairs\n              generate all-pairs pairwise alignments\n\n" +
        "       -viterbi\n              use Viterbi algorithm to generate all pairs (automatically enables -pairs)\n\n" +
        "       -v, --verbose\n              report progress while aligning (default: off)\n\n" +
        "       -annot FILENAME\n              write annotation for multiple alignment to FILENAME\n\n" +
        "       -t, --train FILENAME\n              compute EM transition probabilities, store in FILENAME (default: no training)\n\n" +
        "       -e, --emissions\n              also reestimate emission probabilities (default: off)\n\n" +
        "       -p, --paramfile FILENAME\n              read parameters from FILENAME (default: )\n\n" +
        "       -a, --alignment-order\n              print sequences in alignment order rather than input order (default: off)\n\n";

    /// <summary>The text of the "Usage" screen (what the program prints to the error stream when called without arguments, after the heading).</summary>
    public static string UsageText => Usage;

    /// <summary>Parses the arguments like <c>ParseParams</c>.</summary>
    /// <exception cref="ArgumentException">An unrecognised option, a missing / invalid / out-of-range option value, <c>-e</c> without <c>-t</c>; the message is the program's error text.</exception>
    public static ProbConsCommand Parse(IReadOnlyList<string> arguments)
    {
        ArgumentNullException.ThrowIfNull(arguments);
        return ParseCore(arguments);
    }

    private static ProbConsCommand ParseCore(IReadOnlyList<string> argv)
    {
        var o = new ProbConsOptions();
        var inputs = new List<string>();
        string? train = null, paramFile = null, annot = null;
        bool trainEmissions = false, allPairs = false, viterbi = false, clustal = false;
        int i = -1;
        while (++i < argv.Count)
        {
            string a = argv[i];
            if (a.Length > 0 && a[0] == '-')
            {
                switch (a)
                {
                    case "-t":
                    case "--train":
                        if (i < argv.Count - 1) train = argv[++i];
                        else throw new ArgumentException("ERROR: Filename expected for option " + a);
                        break;
                    case "-e":
                    case "--emissions":
                        trainEmissions = true;
                        break;
                    case "-p":
                    case "--paramfile":
                        if (i < argv.Count - 1) paramFile = argv[++i];
                        else throw new ArgumentException("ERROR: Filename expected for option " + a);
                        break;
                    case "-c":
                    case "--consistency":
                        o = o with { ConsistencyReps = ReadInt(argv, ref i, 0, 5) };
                        break;
                    case "-ir":
                    case "--iterative-refinement":
                        o = o with { IterativeRefinementReps = ReadInt(argv, ref i, 0, 1000) };
                        break;
                    case "-pre":
                    case "--pre-training":
                        o = o with { PreTrainingReps = ReadInt(argv, ref i, 0, 20) };
                        break;
                    case "-go":
                    case "--gap-open":
                        o = o with { GapOpenPenalty = ReadFloat(argv, ref i, float.NegativeInfinity, 0, "must not be positive") };
                        break;
                    case "-ge":
                    case "--gap-extension":
                        o = o with { GapContinuePenalty = ReadFloat(argv, ref i, float.NegativeInfinity, 0, "must not be positive") };
                        break;
                    case "-pairs":
                        allPairs = true;
                        break;
                    case "-viterbi":
                        allPairs = true;
                        viterbi = true;
                        break;
                    case "-annot":
                        if (i < argv.Count - 1) annot = argv[++i];
                        else throw new ArgumentException("ERROR: FILENAME expected for option " + a);
                        break;
                    case "-clustalw":
                        clustal = true;
                        break;
                    case "-co":
                    case "--cutoff":
                        o = o with { Cutoff = ReadFloat(argv, ref i, 0, 1, "must be between 0 and 1") };
                        break;
                    case "-v":
                    case "--verbose":
                        o = o with { Verbose = true };
                        break;
                    case "-a":
                    case "--alignment-order":
                        o = o with { AlignmentOrder = true };
                        break;
                    default:
                        throw new ArgumentException("ERROR: Unrecognized option: " + a);
                }
            }
            else
            {
                inputs.Add(a);
            }
        }

        if (trainEmissions && train == null)
            throw new ArgumentException("ERROR: Training emissions (-e) requires training (-t)");
        o = o with { Annotate = annot != null };
        return new ProbConsCommand(o, inputs, train, trainEmissions, paramFile, annot, allPairs, viterbi, clustal);
    }

    private static int ReadInt(IReadOnlyList<string> argv, ref int i, int min, int max)
    {
        string opt = argv[i];
        if (i >= argv.Count - 1) throw new ArgumentException("ERROR: Integer expected for option " + opt);
        string text = argv[++i];
        if (!GetInteger(text, out int v)) throw new ArgumentException("ERROR: Invalid integer following option " + opt + ": " + text);
        if (v < min || v > max) throw new ArgumentException("ERROR: For option " + opt + ", integer must be between " + min + " and " + max + ".");
        return v;
    }

    private static float ReadFloat(IReadOnlyList<string> argv, ref int i, float min, float max, string rangeText)
    {
        string opt = argv[i];
        if (i >= argv.Count - 1) throw new ArgumentException("ERROR: Floating-point value expected for option " + opt);
        string text = argv[++i];
        if (!GetFloat(text, out float v)) throw new ArgumentException("ERROR: Invalid floating-point value following option " + opt + ": " + text);
        bool bad = max == 0 && float.IsNegativeInfinity(min) ? v > 0 : v < min || v > max;
        if (bad) throw new ArgumentException("ERROR: For option " + opt + ", floating-point value " + rangeText + ".");
        return v;
    }

    // ---------- sequence input ----------

    /// <summary>The <c>FileBuffer</c> of the source over a string: <c>eof</c> only turns true after a read past the end.</summary>
    private sealed class Buffer
    {
        private readonly string _s;
        private int _pos;
        private bool _isEof;

        public Buffer(string s) => _s = s;

        public bool Eof => _isEof;

        public bool Get(out char ch)
        {
            ch = '\0';
            if (_isEof) return false;
            if (_pos == _s.Length)
            {
                _isEof = true;
                return false;
            }

            ch = _s[_pos++];
            return true;
        }

        public void UnGet()
        {
            _pos--;
            _isEof = false;
        }

        public string GetLine()
        {
            var sb = new StringBuilder();
            while (Get(out char ch) && ch != '\n') sb.Append(ch);
            return sb.ToString();
        }
    }

    private sealed class Seq
    {
        public string Header = string.Empty;
        public string Residues = string.Empty;
        public bool Valid;
    }

    private static bool IsSpace(char c) => c == ' ' || (c >= '\t' && c <= '\r');

    /// <summary><c>Sequence (FileBuffer&amp;, stripGaps)</c>.</summary>
    private static Seq ReadSequence(Buffer infile, bool stripGaps)
    {
        var seq = new Seq { Header = "~" };
        string header = "~";
        while (!infile.Eof)
        {
            header = infile.GetLine();
            if (header.Length != 0) break;
        }

        seq.Header = header;
        if (header.Length > 0 && header[0] == '>')
        {
            header = header[1..];
            int s = 0;
            while (s < header.Length && IsSpace(header[s])) s++;
            header = header[s..];
            int e = header.Length;
            while (e > 0 && IsSpace(header[e - 1])) e--;
            header = header[..e];
            seq.Header = header;
            var data = new StringBuilder();
            while (infile.Get(out char ch))
            {
                if (ch == '>')
                {
                    infile.UnGet();
                    break;
                }

                if (IsSpace(ch)) continue;
                if (ch == '.') ch = '-';
                if (stripGaps && ch == '-') continue;
                if (!((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') || ch == '*' || ch == '-'))
                    throw new ArgumentException("ERROR: Unknown character encountered: " + ch);
                data.Append(ch);
            }

            seq.Residues = data.ToString();
            seq.Valid = seq.Residues.Length > 0;
        }

        return seq;
    }

    /// <summary><c>MultiSequence::LoadMFA</c>: appends the sequences of one file to <paramref name="names"/> / <paramref name="residues"/>.</summary>
    private static void LoadMfa(string text, List<string> names, List<string> residues)
    {
        var infile = new Buffer(text);
        while (true)
        {
            int index = names.Count;
            var seq = ReadSequence(infile, true);
            if (!seq.Valid)
            {
                if (index == 0)
                {
                    string header = seq.Header;
                    if (header.Length > 0 && header[0] != '>')
                    {
                        ParseMsf(infile, header, names, residues);
                        break;
                    }
                }

                break;
            }

            names.Add(seq.Header);
            residues.Add(seq.Residues);
        }

        if (names.Count == 0) throw new ArgumentException("ERROR: No sequences read.");
    }

    private static List<string> Tokens(string line)
    {
        var l = new List<string>();
        int i = 0;
        while (i < line.Length)
        {
            while (i < line.Length && IsSpace(line[i])) i++;
            int s = i;
            while (i < line.Length && !IsSpace(line[i])) i++;
            if (i > s) l.Add(line[s..i]);
        }

        return l;
    }

    /// <summary><c>MultiSequence::ParseMSF</c> (MSF and Clustal input, with the program's quirks: gaps are kept in the record (the program then reads past its sequence arrays: Segmentation fault in 2 of 5 probes, mis-aligned rows otherwise; this port aligns the gap-free sequences, see <c>Load</c> of the engine), letters are upper-cased, a last line without a newline is lost).</summary>
    private static void ParseMsf(Buffer infile, string header, List<string> names, List<string> residues)
    {
        var seqNames = new List<string>();
        var seqData = new List<StringBuilder>();
        bool missingHeader = false, clustalW = false;
        while (!infile.Eof && !header.Contains("..", StringComparison.Ordinal))
        {
            if (header.StartsWith("CLUSTAL", StringComparison.Ordinal) || header.StartsWith("PROBCONS", StringComparison.Ordinal))
            {
                clustalW = true;
                break;
            }

            header = infile.GetLine();
            if (header.Contains("//", StringComparison.Ordinal))
            {
                missingHeader = true;
                break;
            }
        }

        bool valid = true;
        while (valid)
        {
            header = infile.GetLine();
            if (infile.Eof) break;
            var toks = Tokens(header);
            if (toks.Count == 0) continue;
            string word = toks[0];
            if (clustalW && !IsSpace(header[0]) && !seqNames.Contains(word))
            {
                seqNames.Add(word);
                seqData.Add(new StringBuilder());
            }

            if (word == "Name:")
            {
                if (toks.Count > 1)
                {
                    seqNames.Add(toks[1]);
                    seqData.Add(new StringBuilder());
                }
                else
                {
                    valid = false;
                }
            }
            else if (seqNames.Contains(word))
            {
                AppendMsfChars(toks, seqData[seqNames.IndexOf(word)]);
            }
            else if (missingHeader)
            {
                seqNames.Add(word);
                seqData.Add(new StringBuilder());
                AppendMsfChars(toks, seqData[^1]);
            }
        }

        if (seqNames.Count == 0) throw new ArgumentException("ERROR: No sequences read!");
        for (int i = 0; i < seqNames.Count; i++)
        {
            if (seqData[i].Length == 0) throw new ArgumentException("ERROR: Sequence of zero length!");
            names.Add(seqNames[i]);
            residues.Add(seqData[i].ToString());
        }
    }

    // `while (in >> ch)` after the name: the characters of the remaining tokens (the name token was consumed by `in >> word`)
    private static void AppendMsfChars(List<string> toks, StringBuilder target)
    {
        for (int t = 1; t < toks.Count; t++)
        {
            foreach (char c0 in toks[t])
            {
                char ch = c0;
                if (ch >= 'a' && ch <= 'z') ch = (char)(ch - 'a' + 'A');
                if (ch == '.') ch = '-';
                if (!((ch >= 'A' && ch <= 'Z') || ch == '*' || ch == '-'))
                    throw new ArgumentException("ERROR: Unknown character encountered: " + ch);
                target.Append(ch);
            }
        }
    }

    // ---------- writers ----------

    /// <summary><c>MultiSequence::WriteMFA</c>: <c>&gt;header</c> and 60 residues per line.</summary>
    internal static string WriteMfa(IReadOnlyList<string> headers, IReadOnlyList<string> rows)
    {
        var sb = new StringBuilder();
        for (int i = 0; i < rows.Count; i++)
        {
            sb.Append('>').Append(headers[i]).Append('\n');
            string r = rows[i];
            int ct = 1;
            for (; ct <= r.Length; ct++)
            {
                sb.Append(r[ct - 1]);
                if (ct % 60 == 0) sb.Append('\n');
            }

            if ((ct - 1) % 60 != 0) sb.Append('\n');
        }

        return sb.ToString();
    }

    private static string FirstWord(string header)
    {
        int s = 0;
        while (s < header.Length && IsSpace(header[s])) s++;
        int e = s;
        while (e < header.Length && !IsSpace(header[e])) e++;
        return header[s..e];
    }

    /// <summary><c>MultiSequence::WriteALN</c>: the Clustal-like format with the program's conservation line (<c>GetAnnotationChar</c>).</summary>
    internal static string WriteAln(IReadOnlyList<string> headers, IReadOnlyList<string> rows)
    {
        var sb = new StringBuilder();
        sb.Append("PROBCONS version ").Append(Version).Append(" multiple sequence alignment\n");
        int n = rows.Count;
        var names = headers.Select(FirstWord).ToArray();
        int longest = 0;
        foreach (var nm in names) longest = Math.Max(longest, nm.Length);
        longest += 4;
        int written = 0;
        bool allDone = false;
        const int cols = 60;
        while (!allDone)
        {
            sb.Append('\n');
            allDone = true;
            for (int i = 0; i < n; i++)
            {
                if (written < rows[i].Length)
                {
                    sb.Append(names[i]).Append(' ', Math.Max(0, longest - names[i].Length));
                    for (int j = 0; j < cols; j++)
                    {
                        if (written + j < rows[i].Length) sb.Append(rows[i][written + j]);
                        else break;
                    }

                    sb.Append('\n');
                    if (written + cols < rows[i].Length) allDone = false;
                }
            }

            sb.Append(' ', longest);
            for (int j = 0; j < cols; j++)
            {
                var column = new List<char>();
                for (int i = 0; i < n; i++)
                    if (written + j < rows[i].Length) column.Add(rows[i][written + j]);
                if (column.Count > 0) sb.Append(AnnotationChar(column));
            }

            sb.Append('\n');
            written += cols;
        }

        return sb.ToString();
    }

    private static readonly (string Letters, char Mark)[] AnnotationGroups =
    {
        ("STA", ':'), ("NEQK", ':'), ("NHQK", ':'), ("NDEQ", ':'), ("QHRK", ':'), ("MILV", ':'), ("MILF", ':'), ("HY", ':'), ("FYW", ':'),
        ("CSA", '.'), ("ATV", '.'), ("SAG", '.'), ("STNK", '.'), ("STPA", '.'), ("SGND", '.'), ("SNDEQK", '.'), ("NDEQHK", '.'), ("NEHQRK", '.'), ("FVLIM", '.'), ("HFY", '.'),
    };

    private static char AnnotationChar(List<char> column)
    {
        var counts = new int[256];
        int all = column.Count;
        foreach (char c in column) counts[char.ToUpperInvariant(c) & 0xFF]++;
        all -= counts['-'];
        if (all == 1) return ' ';
        for (int i = 0; i < 256; i++)
            if (i != '-' && counts[i] == all) return '*';
        foreach (var (letters, mark) in AnnotationGroups)
        {
            int sum = 0;
            foreach (char l in letters) sum += counts[l];
            if (sum == all) return mark;
        }

        return ' ';
    }

    // ---------- run ----------

    /// <summary>Runs the program on text "files". See <see cref="ProbConsProgramResult"/>.</summary>
    /// <param name="arguments">The command-line arguments (without the program name).</param>
    /// <param name="files">Input files (name -> text): sequence files and the <c>-p</c> parameter file.</param>
    public static ProbConsProgramResult Run(IReadOnlyList<string> arguments, IReadOnlyDictionary<string, string>? files = null)
    {
        ArgumentNullException.ThrowIfNull(arguments);
        files ??= new Dictionary<string, string>();
        var err = new StringBuilder(Heading);
        var outFiles = new Dictionary<string, string>();
        try
        {
            if (arguments.Count < 1)
            {
                err.Append(Usage);
                return new ProbConsProgramResult(1, string.Empty, err.ToString(), outFiles);
            }

            var cmd = ParseCore(arguments);
            var p = ProbConsParameters.Default;
            if (cmd.ParameterFile != null)
            {
                if (!files.TryGetValue(cmd.ParameterFile, out var ptext))
                    throw new ArgumentException("ERROR: Unable to read parameter file: " + cmd.ParameterFile);
                try
                {
                    p = ProbConsParameters.Parse(ptext);
                }
                catch (ArgumentException e)
                {
                    throw new ArgumentException("ERROR: " + e.Message + ": " + cmd.ParameterFile);
                }
            }

            err.Append(p.Describe("Using parameter set:"));
            var o = cmd.Options with { Parameters = p };

            if (cmd.Train != null)
            {
                // the program loads and aligns the files one at a time: "Loading ..." and the progress lines of each pass alternate (a failing file stops the run after the earlier passes)
                var loaded = new List<ProbConsSet>();
                string? failedFile = null;
                ArgumentException? failure = null;
                foreach (var f in cmd.InputFiles)
                {
                    try
                    {
                        var names = new List<string>();
                        var residues = new List<string>();
                        LoadInto(f, files, names, residues);
                        loaded.Add(new ProbConsSet(residues, names));
                    }
                    catch (ArgumentException e)
                    {
                        failedFile = f;
                        failure = e;
                        break;
                    }
                }

                var result = SequenceAligner.ProbConsTrain(loaded, cmd.TrainEmissions, o);
                for (int i = 0; i < loaded.Count; i++)
                    err.Append("Loading sequence file: ").Append(cmd.InputFiles[i]).Append('\n').Append(result.SetLogs[i]);
                if (failure != null)
                {
                    err.Append("Loading sequence file: ").Append(failedFile).Append('\n');
                    throw failure;
                }

                var trained = result.Parameters;
                err.Append(trained.Describe("Trained parameter set:"));
                outFiles[cmd.Train] = trained.ToFileText();
                return new ProbConsProgramResult(0, string.Empty, err.ToString(), outFiles);
            }

            var allNames = new List<string>();
            var allResidues = new List<string>();
            foreach (var f in cmd.InputFiles)
            {
                err.Append("Loading sequence file: ").Append(f).Append('\n');
                LoadInto(f, files, allNames, allResidues);
            }

            if (allResidues.Count == 0)
            {
                // the Debian patch of Main.cc: DoAlign (): "Zero sequences delivered to function DoAlign()" and no output, exit status 0
                const string zero = "Zero sequences delivered to function DoAlign()\n";
                for (int ct = 0; ct < o.PreTrainingReps; ct++)
                    err.Append(zero).Append(p.Describe("Recomputed parameter set:"));
                err.Append(zero);
                return new ProbConsProgramResult(0, string.Empty, err.ToString(), outFiles);
            }

            o = o with { Names = allNames };
            if (cmd.AllPairs)
            {
                var res = SequenceAligner.ProbConsPairs(allResidues, cmd.Viterbi, o);
                err.Append(res.Log);
                foreach (var pa in res.Pairs)
                {
                    string h1 = allNames[pa.First], h2 = allNames[pa.Second];
                    string name = h1 + "-" + h2 + (cmd.ClustalW ? ".aln" : ".fasta");
                    var hs = new[] { h1, h2 };
                    var rows = new[] { pa.FirstRow, pa.SecondRow };
                    outFiles[name] = cmd.ClustalW ? WriteAln(hs, rows) : WriteMfa(hs, rows);
                }

                return new ProbConsProgramResult(0, string.Empty, err.ToString(), outFiles);
            }

            var run = SequenceAligner.ProbConsRun(allResidues, o);
            err.Append(run.Log);
            var outHeaders = run.Order.Select(i => allNames[i]).ToArray();
            var outRows = run.Alignment.AlignedSequences.ToArray();
            string stdout = cmd.ClustalW ? WriteAln(outHeaders, outRows) : WriteMfa(outHeaders, outRows);
            if (cmd.AnnotationFile != null && run.Annotation != null)
            {
                var sb = new StringBuilder();
                foreach (int v in run.Annotation) sb.Append(v.ToString(CultureInfo.InvariantCulture).PadLeft(4)).Append('\n');
                outFiles[cmd.AnnotationFile] = sb.ToString();
            }

            return new ProbConsProgramResult(0, stdout, err.ToString(), outFiles);
        }
        catch (ArgumentException e)
        {
            err.Append(e.Message.StartsWith("ERROR: ", StringComparison.Ordinal) ? e.Message : "ERROR: " + e.Message).Append('\n');
            return new ProbConsProgramResult(1, string.Empty, err.ToString(), outFiles);
        }
    }

    private static void LoadInto(string file, IReadOnlyDictionary<string, string> files, List<string> names, List<string> residues)
    {
        if (!files.TryGetValue(file, out var text))
            throw new ArgumentException("ERROR: Could not open file '" + file + "' for reading.");
        LoadMfa(text, names, residues);
    }
}
