using System.Globalization;
using System.Text;

#pragma warning disable S108, S3358

namespace Seqeron.Genomics.Alignment;

/// <summary>The outcome of <see cref="Minimap2Program.Run"/>: what the program would write to the standard output and to the standard error, and its exit status.</summary>
/// <param name="StdOut">The standard output (PAF or SAM).</param>
/// <param name="StdErr">The error and warning lines the program writes (the progress log <c>[M::...]</c> is not reproduced).</param>
/// <param name="ExitCode">0, or 1 when the program stops with an error.</param>
public sealed record Minimap2RunResult(string StdOut, string StdErr, int ExitCode)
{
    /// <summary>The files the program writes (<c>-o FILE</c>: the output goes to the file, the standard output stays empty).</summary>
    public IReadOnlyDictionary<string, string> Files { get; init; } = new Dictionary<string, string>();
}

/// <summary>
/// <c>minimap2 2.26</c> as a program (<c>main.c</c>, <c>options.c</c>, the readers of <c>bseq.c</c> / <c>kseq.h</c> and the pipeline of <c>map.c</c>): the command-line options with their getopt (<c>ketopt</c>) semantics (permutation of the arguments, abbreviated long options, attached
/// values), <c>-x</c> applied first and cumulatively, <c>mm_check_opt</c> with the program's error texts and exit status, the reference and the query files (FASTA / FASTQ text as <c>kseq</c> reads it, with its quirks; one or two query files = paired reads, an interleaved
/// file in fragment mode), <c>--junc-bed</c> / <c>--alt</c> files, the index in parts (<c>-I</c>: the reader takes sequences in mini-batches of <c>min(50 Mbases, -I)</c> and a part ends when the bases read exceed <c>-I</c>; every part is mapped on its own with the
/// <c>mid_occ</c> of the first part, the SAM header has no <c>@SQ</c> lines), the <c>--split-prefix</c> merge of the hits of all parts (<c>merge_hits</c>), the SAM header (with the permuted command line in <c>@PG</c>), the verbosity of warnings and errors (<c>-v</c>) and the PAF / SAM
/// output of <see cref="Minimap2Chaining"/>. The options that only select threads (<c>-t</c>), the mini-batch of the query reader (<c>-K</c>: the reads of one fragment are never split below 1 Mbase), the index dump (<c>-d</c>), the debug log on the standard error
/// (<c>--print-qname --print-seeds --print-chains --print-aln-seq --no-kalloc --cap-kalloc</c>), <c>--bucket-bits</c> and <c>-2</c> have no effect on the standard output; <c>-o FILE</c> moves the output to <see cref="Minimap2RunResult.Files"/>. The progress log (<c>[M::...]</c>)
/// is not reproduced. The program has no <c>--pe-ori</c> option in 2.26: the orientation of a pair is the preset's (<c>-x sr</c>: forward-reverse) or forward-forward.
/// </summary>
public static class Minimap2Program
{
    private const string Version = "2.26-r1175";

    private static readonly (string Name, int HasArg, int Val)[] LongOptions =
    [
        ("bucket-bits", 1, 300), ("mb-size", 1, 'K'), ("seed", 1, 302), ("no-kalloc", 0, 303), ("print-qname", 0, 304), ("no-self", 0, 'D'), ("print-seeds", 0, 306), ("max-chain-skip", 1, 307),
        ("min-dp-len", 1, 308), ("print-aln-seq", 0, 309), ("splice", 0, 310), ("cost-non-gt-ag", 1, 'C'), ("no-long-join", 0, 312), ("sr", 0, 313), ("frag", 1, 314), ("secondary", 1, 315),
        ("cs", 2, 316), ("end-bonus", 1, 317), ("no-pairing", 0, 318), ("splice-flank", 1, 319), ("idx-no-seq", 0, 320), ("end-seed-pen", 1, 321), ("for-only", 0, 322), ("rev-only", 0, 323),
        ("heap-sort", 1, 324), ("all-chain", 0, 'P'), ("dual", 1, 326), ("max-clip-ratio", 1, 327), ("min-occ-floor", 1, 328), ("MD", 0, 329), ("lj-min-ratio", 1, 330), ("score-N", 1, 331),
        ("eqx", 0, 332), ("paf-no-hit", 0, 333), ("split-prefix", 1, 334), ("no-end-flt", 0, 335), ("hard-mask-level", 0, 336), ("cap-sw-mem", 1, 337), ("max-qlen", 1, 338), ("max-chain-iter", 1, 339),
        ("junc-bed", 1, 340), ("junc-bonus", 1, 341), ("sam-hit-only", 0, 342), ("chain-gap-scale", 1, 343), ("alt", 1, 344), ("alt-drop", 1, 345), ("mask-len", 1, 346), ("rmq", 2, 347),
        ("qstrand", 0, 348), ("cap-kalloc", 1, 349), ("q-occ-frac", 1, 350), ("chain-skip-scale", 1, 351), ("print-chains", 0, 352), ("no-hash-name", 0, 353), ("secondary-seq", 0, 354),
        ("help", 0, 'h'), ("max-intron-len", 1, 'G'), ("version", 0, 'V'), ("min-count", 1, 'n'), ("min-chain-score", 1, 'm'), ("mask-level", 1, 'M'), ("min-dp-score", 1, 's'), ("sam", 0, 'a'),
    ];

    private const string ShortOptions = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:yYPo:e:U:J:";

    // ---------------------------------------------------------------------------------------------------- C library helpers

    private static (long Value, int Used) Strtol(string s, int start)
    {
        int i = start;
        while (i < s.Length && char.IsWhiteSpace(s[i]))
            i++;
        int numStart = i;
        bool neg = false;
        if (i < s.Length && (s[i] == '+' || s[i] == '-'))
            neg = s[i++] == '-';
        int digits = i;
        long v = 0;
        while (i < s.Length && char.IsAsciiDigit(s[i]))
        {
            v = unchecked((v * 10) + (s[i] - '0'));
            i++;
        }

        if (i == digits)
            return (0, start);
        _ = numStart;
        return (neg ? -v : v, i);
    }

    private static int Atoi(string s) => (int)Strtol(s, 0).Value;

    // strtod: the longest prefix that is a decimal floating-point number (no hex / inf / nan needed by the options)
    private static (double Value, int Used) Strtod(string s, int start)
    {
        int i = start;
        while (i < s.Length && char.IsWhiteSpace(s[i]))
            i++;
        int b = i;
        if (i < s.Length && (s[i] == '+' || s[i] == '-'))
            i++;
        int d0 = i;
        while (i < s.Length && char.IsAsciiDigit(s[i]))
            i++;
        int intDigits = i - d0;
        int fracDigits = 0;
        if (i < s.Length && s[i] == '.')
        {
            int j = i + 1;
            while (j < s.Length && char.IsAsciiDigit(s[j]))
                j++;
            fracDigits = j - i - 1;
            if (intDigits > 0 || fracDigits > 0)
                i = j;
        }

        if (intDigits == 0 && fracDigits == 0)
            return (0, start);
        if (i < s.Length && (s[i] == 'e' || s[i] == 'E'))
        {
            int j = i + 1;
            if (j < s.Length && (s[j] == '+' || s[j] == '-'))
                j++;
            int e0 = j;
            while (j < s.Length && char.IsAsciiDigit(s[j]))
                j++;
            if (j > e0)
                i = j;
        }

        return (double.Parse(s.AsSpan(b, i - b), NumberStyles.Float, CultureInfo.InvariantCulture), i);
    }

    private static double Atof(string s) => Strtod(s, 0).Value;

    // mm_parse_num2: strtod with the suffixes G / M / K (either case)
    private static (long Value, int Used) ParseNum2(string s, int start)
    {
        var (x, p) = Strtod(s, start);
        if (p < s.Length && (s[p] == 'G' || s[p] == 'g'))
        {
            x *= 1e9;
            p++;
        }
        else if (p < s.Length && (s[p] == 'M' || s[p] == 'm'))
        {
            x *= 1e6;
            p++;
        }
        else if (p < s.Length && (s[p] == 'K' || s[p] == 'k'))
        {
            x *= 1e3;
            p++;
        }

        return ((long)(x + .499), p);
    }

    private static long ParseNum(string s) => ParseNum2(s, 0).Value;

    // ---------------------------------------------------------------------------------------------------- ketopt

    private sealed class KetOpt
    {
        internal int Ind = 1, Opt;
        internal string? Arg;
        private int i = 1, pos, nArgs;

        private static void Permute(string[] argv, int j, int n)
        {
            string p = argv[j];
            int k;
            for (k = 0; k < n; ++k)
                argv[j - k] = argv[j - k - 1];
            argv[j - k] = p;
        }

        internal int Next(string[] argv)
        {
            int argc = argv.Length, opt = -1, i0, j;
            while (i < argc && (argv[i][0] != '-' || argv[i].Length == 1))
            {
                ++i;
                ++nArgs;
            }

            Arg = null;
            i0 = i;
            if (i >= argc || argv[i][0] != '-' || argv[i].Length == 1)
            {
                Ind = i - nArgs;
                return -1;
            }

            if (argv[i].Length > 1 && argv[i][1] == '-')
            {
                if (argv[i].Length == 2)
                {
                    Permute(argv, i, nArgs);
                    ++i;
                    Ind = i - nArgs;
                    return -1;
                }

                Opt = 0;
                opt = '?';
                pos = -1;
                string a = argv[i];
                for (j = 2; j < a.Length && a[j] != '='; ++j)
                {
                }

                int nExact = 0, nPartial = 0, oExact = -1, oPartial = -1;
                for (int k = 0; k < LongOptions.Length; ++k)
                {
                    string nm = LongOptions[k].Name;
                    if (string.CompareOrdinal(a, 2, nm, 0, j - 2) == 0 && nm.Length >= j - 2)
                    {
                        if (nm.Length == j - 2)
                        {
                            ++nExact;
                            oExact = k;
                        }
                        else
                        {
                            ++nPartial;
                            oPartial = k;
                        }
                    }
                }

                if (nExact > 1 || (nExact == 0 && nPartial > 1))
                    return '?';
                int o = nExact == 1 ? oExact : nPartial == 1 ? oPartial : -1;
                if (o >= 0)
                {
                    Opt = opt = LongOptions[o].Val;
                    if (j < a.Length && a[j] == '=')
                        Arg = a[(j + 1)..];
                    if (LongOptions[o].HasArg == 1 && j >= a.Length)
                    {
                        if (i < argc - 1)
                            Arg = argv[++i];
                        else
                            opt = ':';
                    }
                }
            }
            else
            {
                if (pos == 0)
                    pos = 1;
                opt = Opt = argv[i][pos++];
                int p = ShortOptions.IndexOf((char)opt);
                if (opt == ':' || p < 0)
                {
                    opt = '?';
                }
                else if (p + 1 < ShortOptions.Length && ShortOptions[p + 1] == ':')
                {
                    if (pos >= argv[i].Length)
                    {
                        if (i < argc - 1)
                            Arg = argv[++i];
                        else
                            opt = ':';
                    }
                    else
                    {
                        Arg = argv[i][pos..];
                    }

                    pos = -1;
                }
            }

            if (pos < 0 || pos >= argv[i].Length)
            {
                ++i;
                pos = 0;
                if (nArgs > 0)
                    for (j = i0; j < i; ++j)
                        Permute(argv, j, nArgs);
            }

            Ind = i - nArgs;
            return opt;
        }
    }

    // ---------------------------------------------------------------------------------------------------- sequence files (bseq.c / kseq.h)

    // kseq.h (shipped with minimap2 v2.26) on an in-memory text: ks_getc, ks_getuntil / ks_getuntil2 and kseq_read
    private sealed class KSeq(string text)
    {
        private int pos;
        private int lastChar;
        internal string Name = "", Seq = "", Qual = "", Comment = "";

        private int GetC() => pos < text.Length ? text[pos++] : -1;

        // ks_getuntil2: the characters up to a delimiter (a line feed, or any white space), appended to str; returns the new length of str, or -1 at the end of the text with nothing read; dret = the delimiter (0 when none was found)
        private int GetUntil(bool line, StringBuilder str, out int dret, bool append)
        {
            dret = 0;
            if (!append)
                str.Clear();
            if (pos >= text.Length)
                return -1;
            bool gotany = false;
            while (pos < text.Length)
            {
                char ch = text[pos];
                bool isDelim = line ? ch == '\n' : ch is ' ' or '\t' or '\n' or '\v' or '\f' or '\r';
                if (isDelim)
                {
                    gotany = true;
                    dret = ch;
                    pos++;
                    break;
                }

                str.Append(ch);
                pos++;
            }

            if (!gotany && str.Length == 0)
                return -1;
            if (line && str.Length > 1 && str[^1] == '\r')
                str.Length--;
            return str.Length;
        }

        // kseq_read: >= 0 the sequence length, -1 at the end, -2 for a malformed FASTQ record
        internal int Read()
        {
            int c;
            if (lastChar == 0)
            {
                while ((c = GetC()) >= 0 && c != '>' && c != '@')
                {
                }

                if (c < 0)
                    return c;
                lastChar = c;
            }

            var name = new StringBuilder();
            var comment = new StringBuilder();
            var seq = new StringBuilder();
            var qual = new StringBuilder();
            if (GetUntil(false, name, out c, false) < 0)
                return -1;
            if (c != '\n')
                GetUntil(true, comment, out _, false);
            while ((c = GetC()) >= 0 && c != '>' && c != '+' && c != '@')
            {
                if (c == '\n')
                    continue;
                seq.Append((char)c);
                GetUntil(true, seq, out _, true);
            }

            if (c is '>' or '@')
                lastChar = c;
            Name = name.ToString();
            Comment = comment.ToString();
            Seq = seq.ToString();
            Qual = "";
            if (c != '+')
                return seq.Length;
            while ((c = GetC()) >= 0 && c != '\n')
            {
            }

            if (c == -1)
                return -2;
            while (GetUntil(true, qual, out _, true) >= 0 && qual.Length < seq.Length)
            {
            }

            lastChar = 0;                               // we have not come to the next header line
            Qual = qual.ToString();
            if (seq.Length != qual.Length)
                return -2;
            return seq.Length;
        }
    }

    /// <summary>The records of a FASTA / FASTQ text as <c>kseq</c> reads them (the name is the first word of the header, the comment the rest of the line; multi-line sequences; a <c>+</c> line starts the quality of a FASTQ record; a malformed FASTQ record ends the
    /// reading as the program's reader does) with <c>U</c> converted to <c>T</c> as <c>mm_bseq</c> does.</summary>
    public static IReadOnlyList<MinimapQuery> ReadSequences(string text) => ReadSequences(text, out _);

    internal static IReadOnlyList<MinimapQuery> ReadSequences(string text, out string? failedAfter)
    {
        ArgumentNullException.ThrowIfNull(text);
        failedAfter = null;
        var records = new List<MinimapQuery>();
        var ks = new KSeq(text);
        int ret;
        while ((ret = ks.Read()) >= 0)
            records.Add(new MinimapQuery(ks.Name, ConvertU(ks.Seq), ks.Qual.Length > 0 ? ks.Qual : null, ks.Comment.Length > 0 ? ks.Comment : null));
        if (ret < -1)
            failedAfter = records.Count > 0 ? records[^1].Name : "";
        return records;
    }

    private static string ConvertU(string s)
    {
        if (s.IndexOfAny(['u', 'U']) < 0)
            return s;
        var c = s.ToCharArray();
        for (int i = 0; i < c.Length; i++)
            if (c[i] is 'u' or 'U')
                c[i]--;
        return new string(c);
    }

    // ---------------------------------------------------------------------------------------------------- the program

    private sealed class Parsed
    {
        internal MinimapMapOptions Map = new();
        internal MinimapOutputOptions Out = new();
        internal int OldBestN = -1;
        internal string? Rg, JuncBed, AltList;
        internal bool Help;
        internal StringBuilder WarningText = new();
        internal string Warnings => WarningText.ToString();
        internal bool SplitPrefix;
        internal string? OutFile;
        internal bool IdxNoSeq, DumpIndex;
        internal int Threads = 3;
        internal long BatchSize = 8_000_000_000L;
        internal int Verbose = 3;
    }

    private static void YesNo(Parsed p, string arg, string option, Func<bool, MinimapMapOptions> apply)
    {
        // yes_or_no(): yes / y and no / n; anything else is a warning and leaves the flag alone
        if (arg is "yes" or "y")
            p.Map = apply(true);
        else if (arg is "no" or "n")
            p.Map = apply(false);
        else
            p.WarningText.Append("[WARNING]\u001b[1;31m option '--" + option + "' only accepts 'yes' or 'no'.\u001b[0m\n");
    }

    /// <summary>
    /// Run <c>minimap2</c> with the given command line.
    /// </summary>
    /// <param name="args">The arguments after the program name: options and file names (any order: the options are permuted before the files as getopt does).</param>
    /// <param name="files">The text of the files named on the command line (the reference, the query files, <c>--junc-bed</c> and <c>--alt</c> files).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>The standard output, the error text and the exit status.</returns>
    public static Minimap2RunResult Run(IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files, CancellationToken cancellationToken = default)
    {
        var p = new Parsed();
        var r = RunCore(args, files, p, cancellationToken);
        if (p.OutFile is { } f && f != "-")
            return r with { StdOut = "", Files = new Dictionary<string, string> { [f] = r.StdOut } };
        return r;
    }

    private static Minimap2RunResult RunCore(IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files, Parsed p, CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(args);
        ArgumentNullException.ThrowIfNull(files);
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        var argv = new string[args.Count + 1];
        argv[0] = "minimap2";
        for (int i = 0; i < args.Count; i++)
            argv[i + 1] = args[i] ?? throw new ArgumentException("A command-line argument is null.", nameof(args));
        string err = "";
        // the first pass applies -x (and finds errors)
        var ko = new KetOpt();
        int c;
        string? preset = null;
        while ((c = ko.Next(argv)) >= 0)
        {
            if (c == 'x')
            {
                try
                {
                    p.Map = MinimapPresets.ApplyOn(p.Map, ko.Arg);
                    preset = ko.Arg;
                }
                catch (ArgumentException)
                {
                    return new Minimap2RunResult("", "[ERROR] unknown preset '" + ko.Arg + "'\n", 1);
                }
            }
            else if (c == ':')
            {
                return new Minimap2RunResult("", "[ERROR] missing option argument\n", 1);
            }
            else if (c == '?')
            {
                return new Minimap2RunResult("", "[ERROR] unknown option in \"" + argv[Math.Max(0, ko.Ind - 1)] + "\"\n", 1);
            }
        }

        _ = preset;
        ko = new KetOpt();
        int nPos0 = -1;
        while ((c = ko.Next(argv)) >= 0)
        {
            string a = ko.Arg ?? "";
            var m = p.Map;
            switch (c)
            {
                case 'w': p.Map = m with { Seed = m.Seed with { WindowSize = Atoi(a) } }; break;
                case 'k': p.Map = m with { Seed = m.Seed with { KmerSize = Atoi(a) } }; break;
                case 'H': p.Map = m with { Seed = m.Seed with { HomopolymerCompressed = true } }; break;
                case 'v': p.Verbose = Atoi(a); break;
                case 'o':
                    if (a != "-")
                        p.OutFile = a;
                    break;
                case 'd': p.DumpIndex = true; break;
                case 't': p.Threads = Atoi(a); break;
                case 'K': case '2': case 'x': break;
                case 320: p.IdxNoSeq = true; break;
                case 'I': p.BatchSize = ParseNum(a); break;
                case 'g': p.Map = m with { MaxGap = (int)ParseNum(a) }; break;
                case 'G': p.Map = MinimapPresets.WithMaxIntronLength(m, (int)ParseNum(a)); break;
                case 'F': p.Map = m with { MaxFragmentLength = (int)ParseNum(a) }; break;
                case 'N': p.OldBestN = m.BestN; p.Map = m with { BestN = Atoi(a) }; break;
                case 'p': p.Map = m with { PrimaryRatio = (float)Atof(a) }; break;
                case 'M': p.Map = m with { MaskLevel = (float)Atof(a) }; break;
                case 'c': p.Map = m with { Cigar = true }; p.Out = p.Out with { CigarTag = true }; break;
                case 'D': p.Map = m with { Seed = m.Seed with { SkipSelf = true } }; break;
                case 'P': p.Map = m with { AllChains = true }; break;
                case 'X': p.Map = m with { AllChains = true, NoLongJoin = true, Seed = m.Seed with { SkipSelf = true, SkipDual = true } }; break;
                case 'a': p.Map = m with { Cigar = true }; p.Out = p.Out with { Sam = true }; break;
                case 'Q': p.Out = p.Out with { NoQuality = true }; break;
                case 'Y': p.Out = p.Out with { SoftClip = true }; break;
                case 'L': p.Out = p.Out with { LongCigar = true }; break;
                case 'y': p.Out = p.Out with { CopyComment = true }; break;
                case 'T': p.Map = m with { Seed = m.Seed with { SdustThreshold = Atoi(a) } }; break;
                case 'n': p.Map = m with { MinCount = Atoi(a) }; break;
                case 'm': p.Map = m with { MinChainScore = Atoi(a) }; break;
                case 'A': p.Map = m with { Scoring = m.Scoring with { Match = Atoi(a) } }; break;
                case 'B': p.Map = m with { Scoring = m.Scoring with { Mismatch = Atoi(a) } }; break;
                case 's': p.Map = m with { Scoring = m.Scoring with { MinDpMax = Atoi(a) } }; break;
                case 'C': p.Map = m with { Scoring = m.Scoring with { NonCanonical = Atoi(a) } }; break;
                case 'e': p.Map = m with { Seed = m.Seed with { OccDist = (int)ParseNum(a) } }; break;
                case 'R': p.Rg = a; break;
                case 'h': p.Help = true; break;
                case 'J':
                    {
                        int t = Atoi(a);
                        if (t == 0)
                            p.Map = m with { SpliceOld = true };
                        else if (t == 1)
                            p.Map = m with { SpliceOld = false };
                        break;
                    }

                case 302: p.Map = m with { HashSeed = Atoi(a) }; break;
                case 306: case 309: p.Threads = 1; break;
                case 303: case 304: case 300: case 349: case 352: break;
                case 307: p.Map = m with { MaxChainSkip = Atoi(a) }; break;
                case 339: p.Map = m with { MaxChainIterations = Atoi(a) }; break;
                case 308: p.Map = m with { MinKswLength = Atoi(a) }; break;
                case 310: p.Map = m with { Splice = true }; break;
                case 312: p.Map = m with { NoLongJoin = true }; break;
                case 313: p.Map = m with { ShortRead = true }; break;
                case 317: p.Map = m with { Scoring = m.Scoring with { EndBonus = Atoi(a) } }; break;
                case 318: p.Map = m with { NoPairing = true }; break;
                case 321: p.Map = m with { AnchorExtensionShift = Atoi(a) }; break;
                case 322: p.Map = m with { Seed = m.Seed with { ForwardOnly = true } }; break;
                case 323: p.Map = m with { Seed = m.Seed with { ReverseOnly = true } }; break;
                case 327: p.Map = m with { MaxClipRatio = (float)Atof(a) }; break;
                case 328: p.Map = m with { Seed = m.Seed with { MinMidOcc = Atoi(a) } }; break;
                case 329: p.Out = p.Out with { Md = true }; break;
                case 331: p.Map = m with { Scoring = m.Scoring with { Ambiguous = Atoi(a) } }; break;
                case 332: p.Map = m with { Eqx = true }; break;
                case 333: p.Out = p.Out with { PafNoHit = true }; break;
                case 334: p.SplitPrefix = true; break;
                case 335: p.Map = m with { NoEndFilter = true }; break;
                case 336: p.Map = m with { HardMaskLevel = true }; break;
                case 337: p.Map = m with { Scoring = m.Scoring with { MaxSmithWatermanMatrix = ParseNum(a) > int.MaxValue ? int.MaxValue : (int)ParseNum(a) } }; break;
                case 338: p.Map = m with { MaxQueryLength = (int)ParseNum(a) }; break;
                case 340: p.JuncBed = a; break;
                case 341: p.Map = m with { Scoring = m.Scoring with { JunctionBonus = Atoi(a) } }; break;
                case 342: p.Out = p.Out with { SamHitOnly = true }; break;
                case 343: p.Map = m with { ChainGapScale = (float)Atof(a) }; break;
                case 351: p.Map = m with { ChainSkipScale = (float)Atof(a) }; break;
                case 344: p.AltList = a; break;
                case 345: p.Map = m with { AltDrop = (float)Atof(a) }; break;
                case 346: p.Map = m with { MaskLength = (int)ParseNum(a) }; break;
                case 348: p.Map = m with { Seed = m.Seed with { QueryStrand = true }, NoInversion = true }; break;
                case 350: p.Map = m with { Seed = m.Seed with { QueryOccurrenceFraction = (float)Atof(a) } }; break;
                case 353: p.Map = m with { NoHashName = true }; break;
                case 354: p.Out = p.Out with { SecondarySeq = true }; break;
                case 330: p.WarningText.Append("[WARNING] \u001b[1;31m --lj-min-ratio has been deprecated.\u001b[0m\n"); break;
                case 314: YesNo(p, a, "frag", v => m with { FragmentMode = v }); break;
                case 315: YesNo(p, a, "secondary", v => m with { OutputSecondary = v }); break;      // --secondary=yes clears MM_F_NO_PRINT_2ND
                case 316:
                    p.Map = m with { Cigar = true };
                    if (ko.Arg is null || a == "short")
                        p.Out = p.Out with { Cs = MinimapCsMode.Short };
                    else if (a == "long")
                        p.Out = p.Out with { Cs = MinimapCsMode.Long };
                    else if (a == "none")
                        p.Out = p.Out with { Cs = MinimapCsMode.None };
                    else
                        p.Out = p.Out with { Cs = MinimapCsMode.Short };
                    if (a is not ("short" or "long" or "none") && ko.Arg is not null && p.Verbose >= 2)
                        p.WarningText.Append("[WARNING]\u001b[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\u001b[0m\n");
                    break;
                case 319: YesNo(p, a, "splice-flank", v => m with { SpliceFlank = v }); break;
                case 324: YesNo(p, a, "heap-sort", v => m with { Seed = m.Seed with { HeapSort = v } }); break;
                case 326: YesNo(p, a, "dual", v => m with { Seed = m.Seed with { SkipDual = !v } }); break;                                       // --dual=yes clears MM_F_NO_DUAL
                case 347:
                    if (ko.Arg is not null)
                        YesNo(p, a, "rmq", v => m with { UseRmq = v });
                    else
                        p.Map = m with { UseRmq = true };
                    break;
                case 'S':
                    p.Map = m with { Cigar = true };
                    p.Out = p.Out with { Cs = MinimapCsMode.Long };
                    if (p.Verbose >= 2)
                        p.WarningText.Append("[WARNING]\u001b[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\u001b[0m\n");
                    break;
                case 'V': return new Minimap2RunResult(Version + "\n", p.Warnings, 0);              // puts(MM_VERSION); return 0 at once: the options after it are not read
                case 'r':
                    {
                        var (bw, used) = ParseNum2(a, 0);
                        long bwl = m.BandwidthLong;
                        if (used < a.Length && a[used] == ',')
                            bwl = ParseNum2(a, used + 1).Value;
                        p.Map = m with { Bandwidth = (int)bw, BandwidthLong = (int)bwl };
                        break;
                    }

                case 'U':
                    {
                        var (lo, used) = Strtol(a, 0);
                        long hi = m.Seed.MaxMidOcc;
                        if (used < a.Length && a[used] == ',')
                            hi = Strtol(a, used + 1).Value;
                        p.Map = m with { Seed = m.Seed with { MinMidOcc = (int)lo, MaxMidOcc = (int)hi } };
                        break;
                    }

                case 'f':
                    {
                        var (x, used) = Strtod(a, 0);
                        var seed = m.Seed;
                        if (x < 1.0)
                            seed = seed with { MidOccFraction = (float)x, MidOcc = 0 };
                        else
                            seed = seed with { MidOcc = (int)(x + .499) };
                        if (used < a.Length && a[used] == ',')
                            seed = seed with { MaxOcc = (int)(Strtod(a, used + 1).Value + .499) };
                        p.Map = m with { Seed = seed };
                        break;
                    }

                case 'u':
                    {
                        if (a.Length > 0 && a[0] == 'b')
                            p.Map = m with { SpliceForward = true, SpliceReverse = true };
                        else if (a.Length > 0 && a[0] == 'f')
                            p.Map = m with { SpliceForward = true, SpliceReverse = false };
                        else if (a.Length > 0 && a[0] == 'r')
                            p.Map = m with { SpliceReverse = true, SpliceForward = false };
                        else if (a.Length > 0 && a[0] == 'n')
                            p.Map = m with { SpliceForward = false, SpliceReverse = false };
                        else
                            return new Minimap2RunResult("", "[ERROR]\u001b[1;31m unrecognized cDNA direction\u001b[0m\n", 1);
                        break;
                    }

                case 'z':
                    {
                        var (z, used) = Strtol(a, 0);
                        long zi = z;
                        if (used < a.Length && a[used] == ',')
                            zi = Strtol(a, used + 1).Value;
                        p.Map = m with { Scoring = m.Scoring with { ZDrop = (int)z, ZDropInverted = (int)zi } };
                        break;
                    }

                case 'O':
                    {
                        var (q, used) = Strtol(a, 0);
                        long q2 = q;
                        if (used < a.Length && a[used] == ',')
                            q2 = Strtol(a, used + 1).Value;
                        p.Map = m with { Scoring = m.Scoring with { GapOpen = (int)q, GapOpen2 = (int)q2 } };
                        break;
                    }

                case 'E':
                    {
                        var (e, used) = Strtol(a, 0);
                        long e2 = e;
                        if (used < a.Length && a[used] == ',')
                            e2 = Strtol(a, used + 1).Value;
                        p.Map = m with { Scoring = m.Scoring with { GapExtend = (int)e, GapExtend2 = (int)e2 } };
                        break;
                    }
            }

            nPos0 = ko.Ind;
        }

        if (nPos0 < 0)
            nPos0 = ko.Ind;
        int ind = ko.Ind;
        var o = p.Map;
        if (o.Splice && o.FragmentMode)
            return new Minimap2RunResult("", "[ERROR]\u001b[1;31m --splice and --frag should not be specified at the same time.\u001b[0m\n", 1);
        // mm_check_opt
        if (o.BestN == 0 && p.Verbose >= 2 && CheckOptions(p) is null)
            p.WarningText.Append("[WARNING]\u001b[1;31m '-N 0' reduces mapping accuracy. Please use '--secondary=no' instead.\u001b[0m\n");
        string? checkError = CheckOptions(p);
        if (checkError is not null)
            return new Minimap2RunResult("", p.Warnings + checkError, 1);
        if (o.BestN == 0)
        {
            p.WarningText.Append("[WARNING]\u001b[1;31m changed '-N 0' to '-N " + p.OldBestN.ToString(CultureInfo.InvariantCulture) + " --secondary=no'.\u001b[0m\n");
            o = o with { BestN = p.OldBestN, OutputSecondary = false };
            p.Map = o;
        }

        var positional = new List<string>();
        for (int i = ind; i < argv.Length; i++)
            positional.Add(argv[i]);
        if (positional.Count == 0 || p.Help)
            return new Minimap2RunResult(p.Help ? Usage(o, p.Threads) : "", p.Help ? p.Warnings : p.Warnings + Usage(o, p.Threads), p.Help ? 0 : 1);
        if (o.ShortRead && positional.Count > 3)
            return new Minimap2RunResult("", p.Warnings + "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n", 1);
        if (!files.TryGetValue(positional[0], out string? refText))
            return new Minimap2RunResult("", p.Warnings + "[ERROR] failed to open file '" + positional[0] + "': No such file or directory\n", 1);
        if (positional.Count < 2 && !p.DumpIndex)
            return new Minimap2RunResult("", p.Warnings + "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n", 1);
        _ = err;
        _ = nPos0;
        return Execute(p, argv, positional, refText, files);
    }

    private static string? CheckOptions(Parsed p)
    {
        var o = p.Map;
        var sc = o.Scoring;
        string Red(string m) => p.Verbose >= 1 ? "[ERROR]\u001b[1;31m " + m + "\u001b[0m\n" : "";
        if (o.Bandwidth > o.BandwidthLong)
            return Red("with '-rNUM1,NUM2', NUM1 (" + o.Bandwidth + ") can't be larger than NUM2 (" + o.BandwidthLong + ")");
        if (o.UseRmq && (o.ShortRead || o.Splice))
            return Red("--rmq doesn't work with --sr or --splice");
        if (p.SplitPrefix && (p.Out.Cs != MinimapCsMode.None || p.Out.Md))
            return Red("--cs or --MD doesn't work with --split-prefix");
        if (o.Seed.KmerSize <= 0 || o.Seed.WindowSize <= 0)
            return Red("-k and -w must be positive");
        if (o.BestN < 0)
            return Red("-N must be no less than 0");
        if (o.PrimaryRatio < 0.0f || o.PrimaryRatio > 1.0f)
            return Red("-p must be within 0 and 1 (including 0 and 1)");
        if (o.Seed.ForwardOnly && o.Seed.ReverseOnly)
            return Red("--for-only and --rev-only can't be applied at the same time");
        if (sc.GapExtend <= 0 || sc.GapOpen <= 0)
            return Red("-O and -E must be positive");
        if ((sc.GapOpen != sc.GapOpen2 || sc.GapExtend != sc.GapExtend2) && !(sc.GapExtend > sc.GapExtend2 && sc.GapOpen + sc.GapExtend < sc.GapOpen2 + sc.GapExtend2))
            return Red("dual gap penalties violating E1>E2 and O1+E1<O2+E2");
        if (sc.GapOpen + sc.GapExtend + sc.GapOpen2 + sc.GapExtend2 > 127)
            return Red("scoring system violating ({-O}+{-E})+({-O2}+{-E2}) <= 127");
        if (sc.ZDrop < sc.ZDropInverted)
            return Red("Z-drop should not be less than inversion-Z-drop");
        if (!o.OutputSecondary && o.AllChains)
            return Red("-X/-P and --secondary=no can't be applied at the same time");
        if (o.Seed.QueryStrand && (p.Out.Sam || o.Splice || o.FragmentMode || o.Seed.HomopolymerCompressed))
            return Red("--qstrand doesn't work with -a, -H, --frag or --splice");
        return null;
    }

    private static string Usage(MinimapMapOptions o, int threads)
    {
        var sb = new StringBuilder();
        string G(double v) => CFormat.General(v, 6);
        sb.Append("Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
        sb.Append("Options:\n  Indexing:\n");
        sb.Append("    -H           use homopolymer-compressed k-mer (preferrable for PacBio)\n");
        sb.Append("    -k INT       k-mer size (no larger than 28) [").Append(o.Seed.KmerSize).Append("]\n");
        sb.Append("    -w INT       minimizer window size [").Append(o.Seed.WindowSize).Append("]\n");
        sb.Append("    -I NUM       split index for every ~NUM input bases [8G]\n");
        sb.Append("    -d FILE      dump index to FILE []\n  Mapping:\n");
        sb.Append("    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [").Append(G(o.Seed.MidOccFraction)).Append("]\n");
        sb.Append("    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [").Append(o.MaxGap).Append("]\n");
        sb.Append("    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
        sb.Append("    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
        sb.Append("    -r NUM[,NUM] chaining/alignment bandwidth and long-join bandwidth [").Append(o.Bandwidth).Append(',').Append(o.BandwidthLong).Append("]\n");
        sb.Append("    -n INT       minimal number of minimizers on a chain [").Append(o.MinCount).Append("]\n");
        sb.Append("    -m INT       minimal chaining score (matching bases minus log gap penalty) [").Append(o.MinChainScore).Append("]\n");
        sb.Append("    -X           skip self and dual mappings (for the all-vs-all mode)\n");
        sb.Append("    -p FLOAT     min secondary-to-primary score ratio [").Append(G(o.PrimaryRatio)).Append("]\n");
        sb.Append("    -N INT       retain at most INT secondary alignments [").Append(o.BestN).Append("]\n  Alignment:\n");
        var sc = o.Scoring;
        sb.Append("    -A INT       matching score [").Append(sc.Match).Append("]\n");
        sb.Append("    -B INT       mismatch penalty (larger value for lower divergence) [").Append(sc.Mismatch).Append("]\n");
        sb.Append("    -O INT[,INT] gap open penalty [").Append(sc.GapOpen).Append(',').Append(sc.GapOpen2).Append("]\n");
        sb.Append("    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [").Append(sc.GapExtend).Append(',').Append(sc.GapExtend2).Append("]\n");
        sb.Append("    -z INT[,INT] Z-drop score and inversion Z-drop score [").Append(sc.ZDrop).Append(',').Append(sc.ZDropInverted).Append("]\n");
        sb.Append("    -s INT       minimal peak DP alignment score [").Append(sc.MinDpMax ?? 80).Append("]\n");
        sb.Append("    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
        sb.Append("    -J INT       splice mode. 0: original minimap2 model; 1: miniprot model [1]\n  Input/Output:\n");
        sb.Append("    -a           output in the SAM format (PAF by default)\n");
        sb.Append("    -o FILE      output alignments to FILE [stdout]\n");
        sb.Append("    -L           write CIGAR with >65535 ops at the CG tag\n");
        sb.Append("    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
        sb.Append("    -c           output CIGAR in PAF\n");
        sb.Append("    --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]\n");
        sb.Append("    --MD         output the MD tag\n");
        sb.Append("    --eqx        write =/X CIGAR operators\n");
        sb.Append("    -Y           use soft clipping for supplementary alignments\n");
        sb.Append("    -t INT       number of threads [").Append(threads).Append("]\n");
        sb.Append("    -K NUM       minibatch size for mapping [500M]\n");
        sb.Append("    --version    show version number\n  Preset:\n");
        sb.Append("    -x STR       preset (always applied before other options; see minimap2.1 for details) []\n");
        sb.Append("                 - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping\n");
        sb.Append("                 - map-hifi - PacBio HiFi reads vs reference mapping\n");
        sb.Append("                 - ava-pb/ava-ont - PacBio/Nanopore read overlap\n");
        sb.Append("                 - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence\n");
        sb.Append("                 - splice/splice:hq - long-read/Pacbio-CCS spliced alignment\n");
        sb.Append("                 - sr - genomic short-read mapping\n");
        sb.Append("\nSee `man ./minimap2.1' for detailed description of these and other advanced command-line options.\n");
        return sb.ToString();
    }

    // mm_idx_gen: the index is built in parts; a part ends when the bases read so far exceed -I (default 8G), and the reader takes sequences in mini-batches of min(50 Mbases, -I) (the last sequence of a mini-batch is the one that
    // reaches the size)
    private static List<List<MinimapQuery>> SplitReference(IReadOnlyList<MinimapQuery> records, long batchSize)
    {
        var parts = new List<List<MinimapQuery>>();
        long mini = Math.Min(50_000_000L, batchSize);
        int pos = 0;
        while (pos < records.Count)
        {
            var part = new List<MinimapQuery>();
            long sumLen = 0;
            while (sumLen <= batchSize)
            {
                long size = 0;
                int before = part.Count;
                while (pos < records.Count)
                {
                    var r = records[pos++];
                    part.Add(r);
                    size += r.Sequence.Length;
                    if (size >= mini)
                        break;
                }

                if (part.Count == before)
                    break;
                for (int i = before; i < part.Count; i++)
                    sumLen += part[i].Sequence.Length;
            }

            parts.Add(part);
        }

        return parts;
    }

    // one entry of the temporary file of --split-prefix: the hits of one read (mm_err_fwrite of n_reg, rep_len, frag_gap, the hits and their extra data)
    private sealed record SplitEntry(List<Minimap2Chaining.Reg> Regs, int RepLen, int FragGap, int MidOcc);

    private static Minimap2RunResult Execute(Parsed p, string[] argv, List<string> positional, string refText, IReadOnlyDictionary<string, string> files)
    {
        var o = p.Map;
        var oo = p.Out;
        var sb = new StringBuilder();
        var reference = ReadSequences(refText, out string? refFailed);
        string RefWarning(string? after) => after is null ? "" : after.Length > 0 ? "[WARNING]\u001b[1;31m failed to parse the FASTA/FASTQ record next to '" + after + "'. Continue anyway.\u001b[0m\n" : "[WARNING]\u001b[1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.\u001b[0m\n";
        var parts = SplitReference(reference, p.BatchSize);
        if (parts.Count == 0)
            return new Minimap2RunResult("", p.Warnings, 0);
        if (p.IdxNoSeq && o.Cigar && positional.Count > 1)
            return new Minimap2RunResult("", p.Warnings + "[ERROR] the prebuilt index doesn't contain sequences.\n", 1);
        // the program computes the output with every secondary hit present and filters them when printing (--secondary=no)
        // mm_mapopt_update: the transcript-strand options switch the splice mode on
        var mapOptions = o with { OutputSecondary = true, SplitPrefix = p.SplitPrefix, Splice = o.Splice || o.SpliceForward || o.SpliceReverse };
        oo = oo with { NoSecondary = !o.OutputSecondary, ReadGroup = p.Rg };
        string[] cl = argv.Skip(1).ToArray();
        string?[] allNames = [.. reference.Select(r => (string?)r.Name)];
        int[] allLengths = [.. reference.Select(r => r.Sequence.Length)];
        if (oo.Sam)
        {
            var rgCheck = Minimap2Chaining.ParseReadGroupChecked(p.Rg, out string? rgError);
            if (rgCheck is null)
            {
                // the header is printed without the read group line, then the program stops
                sb.Append(Minimap2Chaining.FormatSamHeader(allNames, allLengths, p.Rg, cl, parts.Count == 1 && !p.SplitPrefix));
                return new Minimap2RunResult(sb.ToString(), p.Warnings + (p.Verbose >= 1 ? rgError + "\n" : ""), 1);
            }

            if (parts.Count == 1 && !p.SplitPrefix)
                sb.Append(Minimap2Chaining.FormatSamHeader(allNames, allLengths, p.Rg, cl));
            else
                sb.Append(Minimap2Chaining.FormatSamHeader(allNames, allLengths, p.Rg, cl, false));
            if (parts.Count > 1 && !p.SplitPrefix && p.Verbose >= 2 && positional.Count > 1)
                p.WarningText.Append("[WARNING]\u001b[1;31m For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.\u001b[0m\n");
        }

        if (positional.Count == 1)
            return new Minimap2RunResult(sb.ToString(), p.Warnings, 0);
        string? junc = null, alt = null;
        if (p.JuncBed is not null)
            files.TryGetValue(p.JuncBed, out junc);
        if (p.AltList is not null)
            files.TryGetValue(p.AltList, out alt);
        var queryFiles = positional.Skip(1).ToList();
        // a query file that cannot be opened stops the program after the files before it (fragment mode: before anything is mapped)
        string? missing = queryFiles.FirstOrDefault(f => !files.ContainsKey(f));
        if (missing is not null)
        {
            int firstBad = queryFiles.IndexOf(missing);
            queryFiles = o.FragmentMode ? [] : queryFiles.Take(firstBad).ToList();
        }
        bool withQual = oo.Sam && !oo.NoQuality;
        // the fragments of the mapping pass: every query file on its own, or in fragment mode all files together
        List<List<MinimapQuery>> Group(List<List<MinimapQuery>> perFile, bool fragMode, bool pairFiles)
        {
            var result = new List<List<MinimapQuery>>();
            if (!fragMode)
            {
                foreach (var f in perFile)
                    foreach (var r in f)
                        result.Add([r]);
                return result;
            }

            var flat = new List<MinimapQuery>();
            if (pairFiles)
            {
                int n = perFile.Min(x => x.Count);
                for (int k = 0; k < n; k++)
                    foreach (var x in perFile)
                        flat.Add(x[k]);
            }
            else
            {
                flat.AddRange(perFile[0]);
            }

            for (int i = 1, j = 0; i <= flat.Count; ++i)
            {
                if (i == flat.Count || !Minimap2Chaining.QnameSame(flat[i - 1].Name, flat[i].Name))
                {
                    result.Add(flat.GetRange(j, i - j));
                    j = i;
                }
            }

            return result;
        }

        List<MinimapQuery> Prepare(IReadOnlyList<MinimapQuery> reads) => [.. reads.Select(q => withQual || q.Quality is null ? q : q with { Quality = null })];
        var queryFailed = new List<string?>();
        var perFileReads = queryFiles.Select(f =>
        {
            var reads = ReadSequences(files[f], out string? failed);
            queryFailed.Add(failed);
            return Prepare(reads);
        }).ToList();
        if (perFileReads.Count == 0)
            perFileReads.Add([]);
        // map_file_frag: one file at a time unless the fragment mode is on; n files = n segments in the fragment mode
        var mapFragments = !o.FragmentMode ? Group(perFileReads, false, false) : Group(perFileReads, true, perFileReads.Count > 1);
        var perPartEntries = new List<List<SplitEntry>>();
        int midOcc = o.Seed.MidOcc;
        var ridShift = new List<int>();
        int shift = 0;
        foreach (var part in parts)
        {
            AlignmentCancellation.ThrowIfRequested();
            string[] targets = [.. part.Select(r => r.Sequence)];
            string?[] names = [.. part.Select(r => (string?)r.Name)];
            if (part == parts[0])
                p.WarningText.Append(RefWarning(refFailed));
            if (p.Verbose >= 2)
                foreach (var r in part)
                    if (r.Sequence.Length == 0)
                        p.WarningText.Append("[WARNING] the length database sequence '").Append(r.Name).Append("' is 0\n");
            foreach (string? failed in queryFailed)
                p.WarningText.Append(RefWarning(failed));
            int[] lengths = [.. targets.Select(t => t.Length)];
            // mm_mapopt_update sets mid_occ from the first part; the later parts use that value
            var partOptions = mapOptions with { Seed = mapOptions.Seed with { MidOcc = midOcc } };
            var mi = Minimap2Chaining.CreateIndex(targets, names, partOptions, junc, alt);
            var seeds = Minimap2Chaining.CreateSeedIndex(targets, partOptions);
            if (midOcc <= 0)
                midOcc = Minimap2Chaining.SeedIndexMidOcc(seeds);
            ridShift.Add(shift);
            shift += targets.Length;
            var entries = new List<SplitEntry>();
            foreach (var frag in mapFragments)
            {
                AlignmentCancellation.ThrowIfRequested();
                var fr = Minimap2Chaining.MapFragmentRegsWithIndex(mi, seeds, [.. frag.Select(q => q.Sequence)], partOptions, names, frag[0].Name, [.. frag.Select(q => q.Name)]);
                if (p.SplitPrefix)
                {
                    for (int i = 0; i < frag.Count; i++)
                        entries.Add(new SplitEntry(fr.Regs[i], fr.RepLens[i], fr.FragGap, fr.MidOcc));
                }
                else
                {
                    var result = Minimap2Chaining.ToResultWithIndex(fr, partOptions);
                    sb.Append(Minimap2Chaining.FormatFragment(frag, result, names, lengths, targets, oo, o.Seed.QueryStrand));
                }
            }

            perPartEntries.Add(entries);
            if (missing is not null)
                break;
        }

        if (p.SplitPrefix)
        {
            // mm_split_merge: the reads are read again (n_fp = the number of query files: more than one file means fragments of one read of every file) and merged with the stored hits in the order of storage
            if (oo.Sam)
                for (int i = 0; i < allNames.Length; i++)
                    sb.Append("@SQ\tSN:").Append(allNames[i]).Append("\tLN:").Append(allLengths[i].ToString(CultureInfo.InvariantCulture)).Append('\n');
            var mergeFragments = Group(perFileReads, o.FragmentMode || perFileReads.Count > 1, perFileReads.Count > 1);
            int k = 0;
            string[] noTargets = [];
            foreach (var frag in mergeFragments)
            {
                int[] qlens = [.. frag.Select(q => q.Sequence.Length)];
                var partFr = new List<Minimap2Chaining.FragRegs>();
                foreach (var entries in perPartEntries)
                {
                    var regs = new List<Minimap2Chaining.Reg>[frag.Count];
                    var reps = new int[frag.Count];
                    int gap = 0, mid = 0;
                    for (int i = 0; i < frag.Count; i++)
                    {
                        var e = entries[k + i];
                        regs[i] = e.Regs;
                        reps[i] = e.RepLen;
                        gap = e.FragGap;
                        mid = e.MidOcc;
                    }

                    partFr.Add(new Minimap2Chaining.FragRegs { Regs = regs, RepLens = reps, FragGap = gap, MidOcc = mid });
                }

                k += frag.Count;
                var merged = Minimap2Chaining.MergeParts(partFr, ridShift, qlens, mapOptions);
                sb.Append(Minimap2Chaining.FormatFragment(frag, Minimap2Chaining.ToResultWithIndex(merged, mapOptions), allNames, allLengths, noTargets, oo, o.Seed.QueryStrand));
            }
        }

        if (missing is not null)
            return new Minimap2RunResult(sb.ToString(), p.Warnings + (p.Verbose >= 1 ? "ERROR: failed to open file '" + missing + "': No such file or directory\n" : "") + "ERROR: failed to map the query file\n", 1);
        return new Minimap2RunResult(sb.ToString(), p.Warnings, 0);
    }
}
