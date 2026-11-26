// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: LastalArguments.cc, LastdbArguments.cc, getoptUtil.hh), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Last;

/// <summary>The exit of a program with a status and no further message (<c>--help</c>, <c>--version</c>).</summary>
public sealed class LastProgramExitException(string output) : Exception("exit")
{
    public string Output { get; } = output;
}

/// <summary>The scan of an argument vector as GNU <c>getopt</c> does it (options may follow the other arguments; a cluster of flags; an option argument attached or separate).</summary>
internal static class LastGetopt
{
    /// <summary>Returns (option, argument) pairs in order and the remaining arguments; errors are reported as the C library does and surface as an "option error".</summary>
    public static (List<(char Option, string? Argument)> Options, List<string> Remaining, string? Error) Scan(IReadOnlyList<string> args, string optstring, string programName)
    {
        var options = new List<(char, string?)>();
        var rest = new List<string>();
        var err = new StringBuilder();
        bool error = false;
        int i = 0;
        while (i < args.Count)
        {
            string a = args[i++];
            if (a == "--")
            {
                for (int j = i; j < args.Count; j++)
                    rest.Add(args[j]);
                break;
            }

            if (a.Length < 2 || a[0] != '-')
            {
                rest.Add(a);
                continue;
            }

            for (int k = 1; k < a.Length; k++)
            {
                char c = a[k];
                int idx = c == ':' ? -1 : optstring.IndexOf(c);
                if (idx < 0)
                {
                    err.Append(programName).Append(": invalid option -- '").Append(c).Append("'\n");
                    options.Add(('?', null));
                    error = true;
                    continue;
                }

                bool takesArgument = idx + 1 < optstring.Length && optstring[idx + 1] == ':';
                if (!takesArgument)
                {
                    options.Add((c, null));
                    continue;
                }

                if (k + 1 < a.Length)
                {
                    options.Add((c, a[(k + 1)..]));
                }
                else if (i < args.Count)
                {
                    options.Add((c, args[i++]));
                }
                else
                {
                    err.Append(programName).Append(": option requires an argument -- '").Append(c).Append("'\n");
                    options.Add(('?', null));
                    error = true;
                }

                break;
            }
        }

        return (options, rest, error ? err.ToString() : null);
    }
}

/// <summary>The command line of <c>lastal</c> (the options that the port implements; the other options of the program stop with a message that says so).</summary>
internal sealed class LastalArguments
{
    public const string OptionString = "hVvf:r:q:p:X:a:b:A:B:c:F:x:y:z:d:e:D:E:s:S:MT:m:l:L:n:N:C:K:k:W:i:P:R:u:w:t:g:G:j:J:Q:";

    public string ProgramName = "lastal";
    public char OutputFormat = 'm';
    public int OutputType = 3;
    public int ScoreType = -1;
    public int Strand = -1;
    public bool IsQueryStrandMatrix;
    public int Globality;
    public bool IsKeepLowercase = true;
    public int TantanSetting = -1;
    public int MaskLowercase = -1;
    public double MaxEvalue = -1;
    public double QueryLettersPerRandomAlignment = 1e6;
    public double MinScoreGapped = -1;
    public int MinScoreGapless = -1;
    public int MatchScore = -1;
    public int MismatchCost = -1;
    public int GapPairCost = -1;
    public string MatrixFile = "";
    public int AmbiguousLetterOpt;
    public int MaxDropGapped = 100;
    public char MaxDropGappedSuffix = '%';
    public int MaxDropGapless = -1;
    public int MaxDropFinal = 100;
    public char MaxDropFinalSuffix = '%';
    public long MinHitDepth = 1;
    public long MaxHitDepth = -1;
    public long OneHitMultiplicity = 10;
    public long MaxGaplessAlignmentsPerQueryPosition;
    public long MaxAlignmentsPerQueryStrand = -1;
    public long CullingLimitForGaplessAlignments;
    public long CullingLimitForFinalAlignments = -1;
    public long QueryStep = 1;
    public long MinimizerWindow;
    public long BatchSize;
    public int NumOfThreads = 1;
    public long MaxRepeatDistance = 1000;
    public double Temperature = -1;
    public double Gamma = 1;
    public int Verbosity;
    public List<int> DelOpenCosts = [], DelGrowCosts = [], InsOpenCosts = [], InsGrowCosts = [];
    public string LastdbName = "";
    public List<string> Inputs = [];
    public string? Message;

    private static void BadOption(char opt, string arg) => throw new LastProgramException("bad option value: -" + opt + " " + arg);

    private static long ParseLong(string s)
    {
        string t = s.Trim();
        if (t.Length == 0 || !long.TryParse(t, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long v))
            throw new LastProgramException("can't interpret: " + s);
        return v;
    }

    private static double ParseDouble(string s)
    {
        string t = s.Trim();
        if (t.Length == 0 || !double.TryParse(t, NumberStyles.Float, CultureInfo.InvariantCulture, out double v))
            throw new LastProgramException("can't interpret: " + s);
        return v;
    }

    private static List<int> ParseIntList(string s)
    {
        var list = new List<int>();
        foreach (string w in s.Split(','))
        {
            if (!int.TryParse(w.Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int v))
                throw new LastProgramException("bad option value: " + s);
            list.Add(v);
        }

        return list;
    }

    private static int ParseScoreAndSuffix(string s, out char suffix)
    {
        int i = 0;
        while (i < s.Length && (char.IsDigit(s[i]) || (i == 0 && (s[i] == '-' || s[i] == '+'))))
            i++;
        if (i == 0 || !int.TryParse(s[..i], NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int x) || x < 0 || s.Length - i > 1)
            throw new LastProgramException("bad value: " + s);
        suffix = i < s.Length ? s[i] : '\0';
        return x;
    }

    private static char ParseOutputFormat(string text)
    {
        string s = text.ToLowerInvariant();
        if (s is "tab" or "0")
            return 't';
        if (s is "maf" or "1")
            return 'm';
        if (s == "blasttab")
            return 'b';
        if (s == "blasttab+")
            return 'B';
        return '\0';
    }

    public void FromArgs(IReadOnlyList<string> argv)
    {
        ProgramName = argv[0];
        var scanArgs = argv.Skip(1).ToList();
        if (scanArgs.Count > 0 && scanArgs[0] == "--help")
            scanArgs[0] = "-h";
        else if (scanArgs.Count > 0 && scanArgs[0] == "--version")
            scanArgs[0] = "-V";
        var (options, rest, error) = LastGetopt.Scan(scanArgs, OptionString, ProgramName);
        if (error is not null)
            Message = error;
        foreach (var (c, optarg) in options)
        {
            string arg = optarg ?? "";
            switch (c)
            {
                case 'h':
                    throw new LastProgramExitException(Help(ProgramName));
                case 'V':
                    throw new LastProgramExitException("lastal 1256\n");
                case 'v':
                    Verbosity++;
                    break;
                case 'f':
                    OutputFormat = ParseOutputFormat(arg);
                    if (OutputFormat == '\0')
                        BadOption(c, arg);
                    break;
                case 'r':
                    MatchScore = (int)ParseLong(arg);
                    if (MatchScore <= 0)
                        BadOption(c, arg);
                    break;
                case 'q':
                    MismatchCost = (int)ParseLong(arg);
                    if (MismatchCost < 0)
                        BadOption(c, arg);
                    break;
                case 'p':
                    MatrixFile = arg;
                    break;
                case 'X':
                    AmbiguousLetterOpt = (int)ParseLong(arg);
                    if (AmbiguousLetterOpt < 0 || AmbiguousLetterOpt > 3)
                        BadOption(c, arg);
                    break;
                case 'a':
                    DelOpenCosts = ParseIntList(arg);
                    break;
                case 'b':
                    DelGrowCosts = ParseIntList(arg);
                    break;
                case 'A':
                    InsOpenCosts = ParseIntList(arg);
                    break;
                case 'B':
                    InsGrowCosts = ParseIntList(arg);
                    break;
                case 'c':
                    GapPairCost = (int)ParseLong(arg);
                    if (GapPairCost <= 0)
                        BadOption(c, arg);
                    break;
                case 'x':
                    MaxDropGapped = ParseScoreAndSuffix(arg, out char sx);
                    MaxDropGappedSuffix = sx == '\0' ? '\0' : sx;
                    break;
                case 'y':
                    MaxDropGapless = (int)ParseLong(arg);
                    if (MaxDropGapless < 0)
                        BadOption(c, arg);
                    break;
                case 'z':
                    MaxDropFinal = ParseScoreAndSuffix(arg, out char sz);
                    MaxDropFinalSuffix = sz;
                    break;
                case 'd':
                    MinScoreGapless = (int)ParseLong(arg);
                    if (MinScoreGapless < 0)
                        BadOption(c, arg);
                    break;
                case 'e':
                    MinScoreGapped = ParseDouble(arg);
                    if (MinScoreGapped < 0)
                        BadOption(c, arg);
                    break;
                case 'D':
                    QueryLettersPerRandomAlignment = ParseDouble(arg);
                    if (QueryLettersPerRandomAlignment <= 0)
                        BadOption(c, arg);
                    break;
                case 'E':
                    MaxEvalue = ParseDouble(arg);
                    if (MaxEvalue <= 0)
                        BadOption(c, arg);
                    break;
                case 'm':
                    OneHitMultiplicity = ParseLong(arg);
                    break;
                case 'l':
                    MinHitDepth = ParseLong(arg);
                    break;
                case 'L':
                    MaxHitDepth = ParseLong(arg);
                    break;
                case 'k':
                    QueryStep = ParseLong(arg);
                    if (QueryStep <= 0)
                        BadOption(c, arg);
                    break;
                case 'W':
                    MinimizerWindow = ParseLong(arg);
                    break;
                case 's':
                    Strand = (int)ParseLong(arg);
                    if (Strand < 0 || Strand > 2)
                        BadOption(c, arg);
                    break;
                case 'S':
                    IsQueryStrandMatrix = ParseLong(arg) != 0;
                    break;
                case 'T':
                    Globality = (int)ParseLong(arg);
                    if (Globality < 0 || Globality > 1)
                        BadOption(c, arg);
                    break;
                case 'n':
                    MaxGaplessAlignmentsPerQueryPosition = ParseLong(arg);
                    if (MaxGaplessAlignmentsPerQueryPosition <= 0)
                        BadOption(c, arg);
                    break;
                case 'N':
                    MaxAlignmentsPerQueryStrand = ParseLong(arg);
                    break;
                case 'C':
                    CullingLimitForGaplessAlignments = ParseLong(arg);
                    break;
                case 'K':
                    CullingLimitForFinalAlignments = ParseLong(arg);
                    break;
                case 'i':
                    BatchSize = ParseSize(arg);
                    if (BatchSize <= 0)
                        BadOption(c, arg);
                    break;
                case 'P':
                    NumOfThreads = (int)ParseLong(arg);
                    break;
                case 'R':
                    if (arg.Length < 2 || arg[0] < '0' || arg[0] > '1' || arg[1] < '0' || arg[1] > '3' || arg.Length > 2)
                        BadOption(c, arg);
                    IsKeepLowercase = arg[0] == '1';
                    TantanSetting = arg[1] - '0';
                    break;
                case 'u':
                    MaskLowercase = (int)ParseLong(arg);
                    if (MaskLowercase < 0 || MaskLowercase > 3)
                        BadOption(c, arg);
                    break;
                case 'w':
                    MaxRepeatDistance = ParseLong(arg);
                    break;
                case 't':
                    Temperature = ParseDouble(arg);
                    if (Temperature <= 0)
                        BadOption(c, arg);
                    break;
                case 'g':
                    Gamma = ParseDouble(arg);
                    if (Gamma <= 0)
                        BadOption(c, arg);
                    break;
                case 'j':
                    OutputType = (int)ParseLong(arg);
                    if (OutputType < 0 || OutputType > 7)
                        BadOption(c, arg);
                    break;
                case 'J':
                    ScoreType = (int)ParseLong(arg);
                    if (ScoreType < 0 || ScoreType > 1)
                        BadOption(c, arg);
                    break;
                case 'F':
                case 'G':
                case 'M':
                case 'Q':
                    throw new LastProgramException("option -" + c + " is not ported (DNA-versus-protein alignment, greedy extension and other input formats)");
                case '?':
                    throw new LastProgramException("bad option");
            }
        }

        if (rest.Count == 0)
            throw new LastProgramException("please give me a database name and sequence file(s)\n\n" + Usage(ProgramName));
        LastdbName = rest[0];
        Inputs = rest.Skip(1).ToList();
    }

    private static long ParseSize(string s)
    {
        var m = System.Text.RegularExpressions.Regex.Match(s.Trim(), @"^(-?\d+)\s*([KMGTP]?)$");
        if (!m.Success)
            throw new LastProgramException("can't interpret: " + s);
        long x = long.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture);
        int i = m.Groups[2].Value.Length == 0 ? 0 : "KMGTP".IndexOf(m.Groups[2].Value[0]) + 1;
        while (i-- > 0)
            x *= 1024;
        return x;
    }

    public static string Usage(string program) =>
        "Usage: " + program + " [options] lastdb-name fasta-sequence-file(s)\nFind and align similar sequences.\n\nCosmetic options:\n-h, --help: show all options and their default settings, and exit\n-V, --version: show version information, and exit\n-v: be verbose: write messages about what lastal is doing\n-f: output format: TAB, MAF, BlastTab, BlastTab+ (default=MAF)";

    public static string Help(string program) =>
        Usage(program) + "\n\nE-value options (default settings):\n-D: query letters per random alignment (1e+06)\n-E: maximum expected alignments per square giga (1e+18/D/refSize/numOfStrands)\n\nScore options (default settings):\n-r: match score   (2 if -M, else  6 if 1<=Q<=4, else 1 if DNA)\n-q: mismatch cost (3 if -M, else 18 if 1<=Q<=4, else 1 if DNA)\n-p: match/mismatch score matrix (protein-protein: BL62, DNA-protein: BL80)\n-X: N/X is ambiguous in: 0=neither sequence, 1=reference, 2=query, 3=both (0)\n-a: gap existence cost (DNA: 7, protein: 11, 1<=Q<=4: 21)\n-b: gap extension cost (DNA: 1, protein:  2, 1<=Q<=4:  9)\n-A: insertion existence cost (a)\n-B: insertion extension cost (b)\n-c: unaligned residue pair cost (off)\n-F: frameshift cost(s) (off)\n-x: maximum score drop for preliminary gapped alignments (z)\n-y: maximum score drop for gapless alignments (min[t*10, x])\n-z: maximum score drop for final gapped alignments (e-1)\n-d: minimum score for gapless alignments (min[e, 2500/n query letters per hit])\n-e: minimum score for gapped alignments\n\nInitial-match options (default settings):\n-m: maximum initial matches per query position (10)\n-l: minimum length for initial matches (1)\n-L: maximum length for initial matches (infinity)\n-k: use initial matches starting at every k-th position in each query (1)\n-W: use \"minimum\" positions in sliding windows of W consecutive positions\n\nMiscellaneous options (default settings):\n-s: strand: 0=reverse, 1=forward, 2=both (2 for DNA, 1 for protein)\n-S: score matrix applies to forward strand of: 0=reference, 1=query (0)\n-K: omit alignments whose query range lies in >= K others with > score (off)\n-C: omit gapless alignments in >= C others with > score-per-length (off)\n-P: number of parallel threads (1)\n-i: query batch size (8 KiB, unless there is > 1 thread or lastdb volume)\n-M: find minimum-difference alignments (faster but cruder)\n-T: type of alignment: 0=local, 1=overlap (0)\n-n: maximum gapless alignments per query position (infinity if m=0, else m)\n-N: stop after the first N alignments per query strand\n-R: lowercase & simple-sequence options (the same as was used by lastdb)\n-u: mask lowercase during extensions: 0=never, 1=gapless,\n    2=gapless+postmask, 3=always (2 if lastdb -c and Q!=pssm, else 0)\n-w: suppress repeats inside exact matches, offset by <= this distance (1000)\n-G: genetic code (1)\n-t: 'temperature' for calculating probabilities (1/lambda)\n-g: 'gamma' parameter for gamma-centroid and LAMA (1)\n-j: output type: 0=match counts, 1=gapless, 2=redundant gapped, 3=gapped,\n                 4=column ambiguity estimates, 5=gamma-centroid, 6=LAMA,\n                 7=expected counts (3)\n-J: score type: 0=ordinary, 1=full (1 for new-style frameshifts, else 0)\n-Q: input format: fastx, keep, sanger, solexa, illumina, prb, pssm\n                  (default=fasta)\n";

    /// <summary><c>setDefaultsFromAlphabet</c> for a DNA database.</summary>
    public void SetDefaultsFromAlphabet(bool isKeepRefLowercase, int refTantanSetting, bool isCaseSensitiveSeeds, int realNumOfThreads)
    {
        if (Strand < 0)
            Strand = 2;
        if (MatchScore < 0)
            MatchScore = 1;
        if (MismatchCost < 0)
            MismatchCost = 1;
        if (DelOpenCosts.Count == 0)
            DelOpenCosts = [7];
        if (DelGrowCosts.Count == 0)
            DelGrowCosts = [1];
        if (InsOpenCosts.Count == 0)
            InsOpenCosts = [.. DelOpenCosts];
        if (InsGrowCosts.Count == 0)
            InsGrowCosts = [.. DelGrowCosts];
        if (TantanSetting < 0)
        {
            IsKeepLowercase = isKeepRefLowercase;
            TantanSetting = refTantanSetting < 3 ? refTantanSetting : 1;
        }

        if (MaskLowercase < 0)
            MaskLowercase = isCaseSensitiveSeeds ? 2 : 0;
        if (BatchSize == 0)
        {
            if (realNumOfThreads == 1)
                BatchSize = 0x2000;
            else if (OutputType == 0)
                BatchSize = 0x1000000;
            else
                BatchSize = 0x800000;
        }

        if (MaxGaplessAlignmentsPerQueryPosition == 0)
            MaxGaplessAlignmentsPerQueryPosition = OneHitMultiplicity > 0 ? OneHitMultiplicity : -1;
        if (DelOpenCosts.Count != DelGrowCosts.Count || InsOpenCosts.Count != InsGrowCosts.Count)
            throw new LastProgramException("bad gap costs");
        if (DelOpenCosts.Count > 1 || InsOpenCosts.Count > 1)
            throw new LastProgramException("piecewise linear gap costs not implemented");
        if (ScoreType < 0)
            ScoreType = 0;
    }

    /// <summary><c>setDefaultsFromMatrix</c>.</summary>
    public void SetDefaultsFromMatrix(double lambda, double minScore, double maxEvalueDefault)
    {
        if (MaxEvalue < 0)
            MaxEvalue = maxEvalueDefault;
        if (OutputType < 2 && MinScoreGapped < 0)
            MinScoreGapped = MinScoreGapless;
        if (MinScoreGapped < 0)
        {
            if (OutputType > 0 && minScore < 0)
                throw new LastProgramException("can't calculate E-values.\nTo proceed without E-values, set a score threshold with option -e.");
            MinScoreGapped = minScore;
        }

        if (MinScoreGapped > int.MaxValue)
            throw new LastProgramException("the alignment score threshold is too big");
        if (OutputType < 2)
            MinScoreGapless = (int)Math.Ceiling(MinScoreGapped);
        if (Temperature < 0)
            Temperature = 1 / lambda;
        int defaultMaxScoreDrop = (int)Math.Max(Math.Ceiling(MinScoreGapped) - 1, 0.0);
        defaultMaxScoreDrop = Math.Max(defaultMaxScoreDrop, MaxDropGapless);
        if (MaxDropFinalSuffix is '%' or '\0')
            MaxDropFinal = Percent(defaultMaxScoreDrop, MaxDropFinal, MaxDropFinalSuffix);
        else
            throw new LastProgramException("option -z with a suffix other than % is not ported");
        if (MaxDropGappedSuffix is '%' or '\0')
            MaxDropGapped = Percent(MaxDropFinal, MaxDropGapped, MaxDropGappedSuffix);
        else
            throw new LastProgramException("option -x with a suffix other than % is not ported");
        if (MaxDropGapless < 0)
        {
            MaxDropGapless = Temperature < 0 ? 0 : (int)(10.0 * Temperature + 0.5);
            MaxDropGapless = Math.Min(MaxDropGapless, MaxDropGapped);
        }
    }

    private static int Percent(int val, int percentage, char suffix)
    {
        // the default suffix is '%'; a value given without a suffix is an absolute score (the constructor sets suffix '%' only for the defaults)
        if (suffix == '\0')
            return percentage;
        return percentage == 100 ? val : val * percentage / 100;
    }

    private static string FormatG(double x)
    {
        string s = x.ToString("G6", CultureInfo.InvariantCulture);
        int e = s.IndexOf('E');
        if (e >= 0)
        {
            string mant = s[..e];
            string exp = s[(e + 1)..];
            char sign = exp[0] == '-' ? '-' : '+';
            exp = exp.TrimStart('+', '-').PadLeft(2, '0');
            s = mant + "e" + sign + exp;
        }

        return s;
    }

    /// <summary><c>LastalArguments::writeCommented</c>.</summary>
    public string WriteCommented()
    {
        string List(List<int> v) => string.Join(",", v.Select(x => x.ToString(CultureInfo.InvariantCulture)));
        var sb = new StringBuilder();
        sb.Append("# a=").Append(List(DelOpenCosts)).Append(" b=").Append(List(DelGrowCosts)).Append(" A=").Append(List(InsOpenCosts)).Append(" B=").Append(List(InsGrowCosts));
        if (GapPairCost > 0)
            sb.Append(" c=").Append(GapPairCost.ToString(CultureInfo.InvariantCulture));
        sb.Append(" e=").Append(FormatG(MinScoreGapped)).Append(" d=").Append(MinScoreGapless.ToString(CultureInfo.InvariantCulture));
        sb.Append(" x=").Append(MaxDropGapped.ToString(CultureInfo.InvariantCulture)).Append(" y=").Append(MaxDropGapless.ToString(CultureInfo.InvariantCulture)).Append(" z=").Append(MaxDropFinal.ToString(CultureInfo.InvariantCulture));
        sb.Append(" D=").Append(FormatG(QueryLettersPerRandomAlignment));
        if (MaxEvalue >= 0)
            sb.Append(" E=").Append(FormatG(MaxEvalue));
        sb.Append('\n');
        sb.Append("# R=").Append(IsKeepLowercase ? 1 : 0).Append(TantanSetting).Append(" u=").Append(MaskLowercase).Append(" s=").Append(Strand).Append(" S=").Append(IsQueryStrandMatrix ? 1 : 0).Append(" M=0 T=").Append(Globality);
        sb.Append(" m=").Append(OneHitMultiplicity.ToString(CultureInfo.InvariantCulture)).Append(" l=").Append(MinHitDepth.ToString(CultureInfo.InvariantCulture));
        if (MaxHitDepth + 1 > 0)
            sb.Append(" L=").Append(MaxHitDepth.ToString(CultureInfo.InvariantCulture));
        sb.Append(" n=").Append(MaxGaplessAlignmentsPerQueryPosition.ToString(CultureInfo.InvariantCulture));
        if (MaxAlignmentsPerQueryStrand + 1 > 0)
            sb.Append(" N=").Append(MaxAlignmentsPerQueryStrand.ToString(CultureInfo.InvariantCulture));
        if (CullingLimitForGaplessAlignments != 0)
            sb.Append(" C=").Append(CullingLimitForGaplessAlignments.ToString(CultureInfo.InvariantCulture));
        if (CullingLimitForFinalAlignments + 1 > 0)
            sb.Append(" K=").Append(CullingLimitForFinalAlignments.ToString(CultureInfo.InvariantCulture));
        sb.Append(" k=").Append(QueryStep.ToString(CultureInfo.InvariantCulture));
        if (MinimizerWindow > 1)
            sb.Append(" W=").Append(MinimizerWindow.ToString(CultureInfo.InvariantCulture));
        sb.Append(" w=").Append(MaxRepeatDistance.ToString(CultureInfo.InvariantCulture)).Append(" t=").Append(FormatG(Temperature));
        if (OutputType > 4 && OutputType < 7)
            sb.Append(" g=").Append(FormatG(Gamma));
        sb.Append(" j=").Append(OutputType).Append(" Q=0\n");
        sb.Append("# ").Append(LastdbName).Append('\n');
        return sb.ToString();
    }
}
