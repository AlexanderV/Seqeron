// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

/// <summary>What a ClustalW command line asks the program to do (the else-if chain at the end of <c>CommandLineParser::parseParams</c>).</summary>
public enum ClustalWVerb
{
    /// <summary><c>-ALIGN</c> (also the default when nothing else is asked): full multiple alignment (<c>Clustal::align</c>).</summary>
    Align = 0,

    /// <summary><c>-PROFILE1 -PROFILE2 -PROFILE</c> (<c>Clustal::profileAlign</c>).</summary>
    Profile = 1,

    /// <summary><c>-PROFILE1 -PROFILE2 -SEQUENCES</c> (<c>Clustal::sequencesAlignToProfile</c>).</summary>
    Sequences = 2,

    /// <summary><c>-CONVERT</c>: the input sequences in another file format.</summary>
    Convert = 3,

    /// <summary><c>-USETREE=file</c> without <c>-ALIGN</c>: alignment along a user guide tree (<c>Clustal::doAlignUseOldTree</c>).</summary>
    AlignUsingTree = 4,

    /// <summary><c>-NEWTREE=file</c> without <c>-ALIGN</c>: only the guide tree (<c>Clustal::doGuideTreeOnly</c>).</summary>
    GuideTree = 5,

    /// <summary><c>-TREE</c>: neighbour-joining tree of the alignment.</summary>
    Tree = 6,

    /// <summary><c>-BOOTSTRAP</c>: bootstrapped neighbour-joining tree.</summary>
    Bootstrap = 7,
}

/// <summary>The meaning of a ClustalW command line (<see cref="ClustalWCommandLine.Parse"/>): verb, file names, and the options for DNA and for protein input.</summary>
/// <param name="Verb">What to do.</param>
/// <param name="Type"><c>-TYPE</c> when given, null = decided from the input (<see cref="ClustalWAligner.GuessClustalWSequenceType"/>).</param>
/// <param name="Dna">The options when the sequences are DNA.</param>
/// <param name="Protein">The options when the sequences are protein.</param>
/// <param name="Run">Guide-tree, row-order, output and structure options (the guide trees' file contents are filled from the files given to the parser).</param>
/// <param name="InFile"><c>-INFILE</c> (or the first argument).</param>
/// <param name="Profile1"><c>-PROFILE1</c>.</param>
/// <param name="Profile2"><c>-PROFILE2</c>.</param>
/// <param name="OutFile"><c>-OUTFILE</c>.</param>
/// <param name="NewTreeFile">The name given to <c>-NEWTREE</c> / <c>-NEWTREE1</c> / <c>-NEWTREE2</c> (<c>-NEWTREE1/2</c>: the program writes <c>&lt;profile&gt;.dnd</c> regardless).</param>
/// <param name="Messages">What the program prints to the error stream while reading the line (ignored values, warnings).</param>
public sealed record ClustalWCommand(
    ClustalWVerb Verb,
    ClustalWSequenceType? Type,
    ClustalWOptions Dna,
    ClustalWProteinOptions Protein,
    ClustalWRunOptions Run,
    string? InFile,
    string? Profile1,
    string? Profile2,
    string? OutFile,
    string? NewTreeFile,
    IReadOnlyList<string> Messages);

/// <summary>
/// The command line of ClustalW 2.1 (<c>CommandLineParser::checkParam</c> / <c>parseParams</c> / <c>setOptionalParam</c>): option names are case-insensitive, <c>-</c> or <c>/</c> starts an option, values follow <c>=</c>;
/// a value of an enumerated option is matched by its first distinguishing letters (<c>CommandLineParser::findMatch</c>); numbers are read like <c>sscanf</c>. Every option the program has is accepted;
/// the effect of each is stored in the options of <see cref="ClustalWCommand"/> (alignment parameters for DNA and for protein, tree / output / structure parameters in <see cref="ClustalWRunOptions"/>).
/// </summary>
public static class ClustalWCommandLine
{
    private static readonly string[] BoolFlags =
    [
        "help", "fullhelp", "quiet", "check", "options", "align", "bootstrap", "tree", "quicktree", "convert", "interactive", "batch", "profile", "sequences", "negative", "noweights", "endgaps", "nopgap",
        "nohgap", "novgap", "kimura", "tossgaps", "nosecstr1", "nosecstr2", "pim",
    ];

    private static readonly string[] ValueFlags =
    [
        "infile", "profile1", "profile2", "newtree", "usetree", "newtree1", "usetree1", "newtree2", "usetree2", "iteration", "type", "matrix", "dnamatrix", "gapopen", "gapext", "hgapresidues", "maxdiv", "gapdist",
        "pwmatrix", "pwdnamatrix", "pwgapopen", "pwgapext", "ktuple", "window", "pairgap", "topdiags", "score", "transweight", "seed", "bootlabels", "debug", "output", "outputtree", "outfile", "outorder", "case",
        "seqnos", "seqno_range", "range", "secstrout", "helixgap", "strandgap", "loopgap", "terminalgap", "helixendin", "helixendout", "strandendin", "strandendout", "numiter", "clustering", "maxseqlen", "stats",
    ];

    /// <summary>
    /// Reads a ClustalW command line. <paramref name="files"/> gives the content of the files the line names (<c>-USETREE=x</c>, <c>-MATRIX=x</c>, <c>-DNAMATRIX=x</c>, <c>-PWMATRIX=x</c>, <c>-PWDNAMATRIX=x</c>): file name to text.
    /// </summary>
    /// <param name="arguments">The arguments (<c>-ALIGN</c>, <c>-GAPOPEN=10</c>, ...); a first argument without <c>-</c> or <c>/</c> is the input file.</param>
    /// <param name="files">Contents of named files.</param>
    /// <exception cref="ArgumentException">The program's refusal: unknown option, missing or malformed value.</exception>
    public static ClustalWCommand Parse(IReadOnlyList<string> arguments, IReadOnlyDictionary<string, string>? files = null)
    {
        ArgumentNullException.ThrowIfNull(arguments);
        var messages = new List<string>();
        var flags = new Dictionary<string, string>(StringComparer.Ordinal);
        string? inFile = null;
        for (int i = 0; i < arguments.Count; i++)
        {
            string a = arguments[i] ?? "";
            if (i == 0 && a.Length > 0 && a[0] is not ('-' or '/'))
            {
                inFile = a;
                continue;
            }

            string body = new string(a.Skip(1).Where(c => c >= ' ' && c < 127).ToArray());
            string name = body, value = "";
            int eq = body.IndexOf('=', StringComparison.Ordinal);
            if (eq >= 0)
            {
                name = body[..eq];
                value = body[(eq + 1)..];
            }

            name = name.ToLowerInvariant();
            bool known = BoolFlags.Contains(name) || ValueFlags.Contains(name);
            if (!known)
                throw new ArgumentException($"Error: unknown option -{name}", nameof(arguments));
            if (ValueFlags.Contains(name) && value.Length == 0)
                throw new ArgumentException($"Error: parameter required for -{name}", nameof(arguments));
            // a value that is not a file name is lower-cased
            if (value.Length > 0 && !(name is "infile" or "profile1" or "profile2" or "newtree" or "usetree" or "newtree1" or "usetree1" or "newtree2" or "usetree2" or "matrix" or "dnamatrix" or "pwmatrix" or "pwdnamatrix" or "outfile" or "stats"))
                value = value.ToLowerInvariant();
            flags[name] = value;
        }

        if (flags.TryGetValue("infile", out var inf))
            inFile = inf.Length > 0 ? inf : throw new ArgumentException("Bad sequence file name", nameof(arguments));

        var dna = new ClustalWOptions();
        var protein = new ClustalWProteinOptions();
        var run = new ClustalWRunOptions();
        var output = new ClustalWOutputOptions();
        var structure = new ClustalWStructureOptions();
        ClustalWSequenceType? type = null;

        if (flags.TryGetValue("type", out var typeValue))
        {
            int t = FindMatch(typeValue, ["protein", "dna"]);
            if (t == 0)
                type = ClustalWSequenceType.Protein;
            else if (t == 1)
                type = ClustalWSequenceType.Dna;
            else
                messages.Add("Unknown sequence type " + typeValue);
        }

        // ---- numbers ----
        if (flags.TryGetValue("numiter", out var numiter))
        {
            int n = ScanInt(numiter, "numiter", "int", 0);
            if (n <= 0)
                throw new ArgumentException("Cannot use a negative value for number of iterations.", nameof(arguments));
            dna = dna with { NumIterations = n };
            protein = protein with { NumIterations = n };
        }

        if (flags.TryGetValue("iteration", out var iteration) && iteration.Length > 0)
        {
            int t = FindMatch(iteration, ["tree", "alignment", "none"]);
            var mode = t switch { 0 => ClustalWIteration.Tree, 1 => ClustalWIteration.Alignment, _ => ClustalWIteration.None };
            if (t < 0 || t > 2)
                messages.Add("Unknown option for iteration. Setting to NONE");
            dna = dna with { Iteration = mode };
            protein = protein with { Iteration = mode };
        }

        if (flags.TryGetValue("clustering", out var clustering) && clustering.Length > 0)
        {
            int t = FindMatch(clustering, ["nj", "upgma"]);
            var mode = t == 1 ? ClustalWClustering.Upgma : ClustalWClustering.NeighborJoining;
            if (t is < 0 or > 1)
                messages.Add("Unknown option for clustering algorithm. Using default");
            dna = dna with { Clustering = mode };
            protein = protein with { Clustering = mode };
        }

        if (flags.TryGetValue("score", out var score) && score.Length > 0)
        {
            int t = FindMatch(score, ["percent", "absolute"]);
            if (t is 0 or 1)
            {
                dna = dna with { Quick = (dna.Quick ?? new ClustalWQuickOptions()) with { Percent = t == 0 } };
                protein = protein with { Pairwise = protein.Pairwise with { Quick = protein.Pairwise.Quick with { Percent = t == 0 } } };
            }
            else
            {
                messages.Add("Unknown SCORE type: " + score);
            }
        }

        if (flags.TryGetValue("output", out var outputValue) && outputValue.Length > 0)
        {
            int t = FindMatch(outputValue, ["gcg", "gde", "pir", "phylip", "nexus", "fasta", "clustal"]);
            if (t < 0 || t > 6)
                throw new ArgumentException("Unknown OUTPUT type: " + outputValue, nameof(arguments));
            output = output with
            {
                Format = t switch
                {
                    0 => ClustalWOutputFormat.Gcg,
                    1 => ClustalWOutputFormat.Gde,
                    2 => ClustalWOutputFormat.Pir,
                    3 => ClustalWOutputFormat.Phylip,
                    4 => ClustalWOutputFormat.Nexus,
                    5 => ClustalWOutputFormat.Fasta,
                    _ => ClustalWOutputFormat.Clustal,
                },
            };
        }

        // ---- the fast pairwise stage: ktuple, pairgap, topdiags, window (DNA: ktuple <= 4, protein: <= 2; the others only above the ktuple in force) ----
        int ktupDna = dna.Quick?.KTuple ?? 2, ktupProt = protein.Pairwise.Quick.KTuple;
        var quickDna = dna.Quick ?? new ClustalWQuickOptions();
        var quickProt = protein.Pairwise.Quick;
        if (flags.TryGetValue("ktuple", out var ktuple))
        {
            int k = ScanInt(ktuple, "ktuple", "integer", 0);
            if (k > 0)
            {
                if (k <= 4)
                {
                    quickDna = quickDna with { KTuple = k, PairGap = k + 4 };
                    ktupDna = k;
                }
                else
                {
                    messages.Add($"WARNING: Ignoring invalid ktuple of {k} (must be <=4)");
                }

                if (k <= 2)
                {
                    quickProt = quickProt with { KTuple = k, PairGap = k + 3 };
                    ktupProt = k;
                }
                else
                {
                    messages.Add($"WARNING: Ignoring invalid ktuple of {k} (must be <=2)");
                }
            }
        }

        if (flags.TryGetValue("pairgap", out var pairgap))
        {
            int g = ScanInt(pairgap, "pairgap", "integer", 0);
            if (g > 0 && g > ktupDna)
                quickDna = quickDna with { PairGap = g };
            if (g > 0 && g > ktupProt)
                quickProt = quickProt with { PairGap = g };
        }

        if (flags.TryGetValue("topdiags", out var topdiags))
        {
            int g = ScanInt(topdiags, "topdiags", "integer", 0);
            if (g > 0 && g > ktupDna)
                quickDna = quickDna with { TopDiagonals = g };
            if (g > 0 && g > ktupProt)
                quickProt = quickProt with { TopDiagonals = g };
        }

        if (flags.TryGetValue("window", out var window))
        {
            int g = ScanInt(window, "window", "integer", 0);
            if (g > 0 && g > ktupDna)
                quickDna = quickDna with { Window = g };
            if (g > 0 && g > ktupProt)
                quickProt = quickProt with { Window = g };
        }

        dna = dna with { Quick = quickDna };
        protein = protein with { Pairwise = protein.Pairwise with { Quick = quickProt } };

        if (flags.ContainsKey("negative"))
        {
            dna = dna with { NegativeMatrix = true };
            protein = protein with { NegativeMatrix = true };
        }

        if (flags.ContainsKey("noweights"))
        {
            dna = dna with { NoWeights = true };
            protein = protein with { NoWeights = true };
        }

        if (flags.ContainsKey("quicktree"))
        {
            dna = dna with { PairwiseMode = ClustalWPairwiseMode.Quick };
            protein = protein with { PairwiseMode = ClustalWPairwiseMode.Quick };
        }

        // ---- matrices ----
        if (flags.TryGetValue("pwmatrix", out var pwm) && pwm.Length > 0)
        {
            string lower = pwm.ToLowerInvariant();
            if (TryProteinPairwiseMatrix(lower, out var pm))
                protein = protein with { Pairwise = protein.Pairwise with { Matrix = pm } };
            else
                protein = protein with { Pairwise = protein.Pairwise with { UserMatrix = ParseMatrixFile(files, pwm) } };
        }

        if (flags.TryGetValue("matrix", out var mat) && mat.Length > 0)
        {
            string lower = mat.ToLowerInvariant();
            if (TryProfileMatrix(lower, out var m))
            {
                protein = protein with { Matrix = m };
            }
            else
            {
                var parsed = ParseMatrixSeries(files, mat);
                protein = protein with { UserMatrix = parsed.Matrix, UserMatrixSeries = parsed.Series };
            }
        }

        if (flags.TryGetValue("pwdnamatrix", out var pwd) && pwd.Length > 0)
        {
            string lower = pwd.ToLowerInvariant();
            if (lower is "iub" or "clustalw")
                dna = dna with { PairwiseDnaMatrix = lower == "iub" ? ClustalWDnaMatrix.Iub : ClustalWDnaMatrix.Clustalw };
            else
                dna = dna with { PairwiseUserDnaMatrix = ParseMatrixFile(files, pwd), UserDnaMatrixUse = ClustalWUserDnaMatrixUse.CommandLine };
        }

        if (flags.TryGetValue("dnamatrix", out var dm) && dm.Length > 0)
        {
            string lower = dm.ToLowerInvariant();
            if (lower is "iub" or "clustalw")
                dna = dna with { DnaMatrix = lower == "iub" ? ClustalWDnaMatrix.Iub : ClustalWDnaMatrix.Clustalw };
            else
                dna = dna with { UserDnaMatrix = ParseMatrixFile(files, dm), UserDnaMatrixUse = ClustalWUserDnaMatrixUse.CommandLine };
        }

        if (flags.TryGetValue("maxdiv", out var maxdiv))
        {
            int v = ScanInt(maxdiv, "maxdiv", "integer", 0);
            if (v >= 0)
            {
                dna = dna with { DivergenceCutoff = v };
                protein = protein with { DivergenceCutoff = v };
            }
        }

        if (flags.TryGetValue("gapdist", out var gapdist))
        {
            int v = ScanInt(gapdist, "gapdist", "integer", 0);
            if (v >= 0)
                protein = protein with { GapSeparationDistance = v };
        }

        if (flags.TryGetValue("case", out var caseValue) && caseValue.Length > 0)
        {
            int t = FindMatch(caseValue, ["lower", "upper"]);
            if (t is 0 or 1)
                output = output with { Lowercase = t == 0 };
            else
                messages.Add("Unknown case " + caseValue);
        }

        if (flags.TryGetValue("seqnos", out var seqnos) && seqnos.Length > 0)
        {
            int t = FindMatch(seqnos, ["off", "on"]);
            if (t is 0 or 1)
                output = output with { SequenceNumbers = t == 1 };
            else
                messages.Add("Unknown SEQNO option " + seqnos);
        }

        if (flags.TryGetValue("seqno_range", out var seqnoRange) && seqnoRange.Length > 0)
        {
            int t = FindMatch(seqnoRange, ["off", "on"]);
            if (t is 0 or 1)
                output = output with { SequenceNumberRange = t == 1 };
            else
                messages.Add("Unknown Sequence range  option " + seqnoRange);
        }

        if (flags.TryGetValue("range", out var range) && range.Length > 0)
        {
            var m = Regex.Match(range, @"^\s*([+-]?\d+)[ :,\-]\s*([+-]?\d+)");
            if (m.Success)
                output = output with { Range = (int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture), int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture)) };
            else
                messages.Add("setRange:  Syntax Error: Cannot set range, should be from:to ");
        }

        if (flags.TryGetValue("gapopen", out var gapopen))
        {
            double v = ScanFloat(gapopen, "gapopen", "real number");
            if (v >= 0.0)
            {
                dna = dna with { GapOpen = v };
                protein = protein with { GapOpen = v };
            }
        }

        if (flags.TryGetValue("gapext", out var gapext))
        {
            double v = ScanFloat(gapext, "gapext", "real number");
            if (v >= 0.0)
            {
                dna = dna with { GapExtend = v };
                protein = protein with { GapExtend = v };
            }
        }

        if (flags.TryGetValue("transweight", out var tw))
            dna = dna with { TransitionWeight = (float)ScanFloat(tw, "transweight", "real number") };

        if (flags.TryGetValue("pwgapopen", out var pgo))
        {
            double v = ScanFloat(pgo, "pwgapopen", "real number");
            if (v >= 0.0)
            {
                dna = dna with { PairwiseGapOpen = v };
                protein = protein with { Pairwise = protein.Pairwise with { GapOpen = v } };
            }
        }

        if (flags.TryGetValue("pwgapext", out var pge))
        {
            double v = ScanFloat(pge, "pwgapext", "real number");
            if (v >= 0.0)
            {
                dna = dna with { PairwiseGapExtend = v };
                protein = protein with { Pairwise = protein.Pairwise with { GapExtend = v } };
            }
        }

        if (flags.TryGetValue("outorder", out var outorder) && outorder.Length > 0)
        {
            int t = FindMatch(outorder, ["input", "aligned"]);
            if (t is 0 or 1)
                run = run with { OutputOrder = t == 0 ? ClustalWOutputOrder.Input : ClustalWOutputOrder.Aligned };
            else
                messages.Add("Unknown OUTPUT ORDER type " + outorder);
        }

        // -ENDGAPS only ever writes false (the program's default is false as well)
        if (flags.ContainsKey("nopgap"))
            protein = protein with { NoPreferredGapPenalties = true };
        if (flags.ContainsKey("nohgap"))
            protein = protein with { NoHydrophilicGapPenalties = true };
        if (flags.ContainsKey("novgap"))
            protein = protein with { VariableRegionGapPenalties = true };
        if (flags.TryGetValue("hgapresidues", out var hgap))
            protein = protein with { HydrophilicResidues = hgap.ToUpperInvariant() };   // UserParameters::setHydResidues upper-cases

        // ---- secondary-structure parameters ----
        if (flags.ContainsKey("nosecstr1"))
            structure = structure with { UseStructure1 = false };
        if (flags.ContainsKey("nosecstr2"))
            structure = structure with { UseStructure2 = false };
        if (flags.TryGetValue("secstrout", out var secstrout) && secstrout.Length > 0)
        {
            int t = FindMatch(secstrout, ["structure", "mask", "both", "none"]);
            if (t is >= 0 and <= 3)
                output = output with { Structure = (ClustalWStructureOutput)t };
            else
                messages.Add("Unknown case " + secstrout);
        }

        structure = ApplyStructureNumber(flags, "helixgap", 1, 9, structure, (s, v) => s with { HelixPenalty = v });
        structure = ApplyStructureNumber(flags, "strandgap", 1, 9, structure, (s, v) => s with { StrandPenalty = v });
        structure = ApplyStructureNumber(flags, "loopgap", 1, 9, structure, (s, v) => s with { LoopPenalty = v });
        structure = ApplyStructureNumber(flags, "terminalgap", 1, 9, structure, (s, v) => s with { TerminalPenalty = v });
        structure = ApplyStructureNumber(flags, "helixendin", 0, 3, structure, (s, v) => s with { HelixEndInside = v });
        structure = ApplyStructureNumber(flags, "helixendout", 0, 3, structure, (s, v) => s with { HelixEndOutside = v });
        structure = ApplyStructureNumber(flags, "strandendin", 0, 3, structure, (s, v) => s with { StrandEndInside = v });
        structure = ApplyStructureNumber(flags, "strandendout", 0, 3, structure, (s, v) => s with { StrandEndOutside = v });

        // ---- verbs, guide-tree files ----
        string? profile1 = flags.GetValueOrDefault("profile1"), profile2 = flags.GetValueOrDefault("profile2");
        if (profile2 != null && profile1 == null)
            throw new ArgumentException("Only 1 profile file (profile 2) specified.", nameof(arguments));
        bool sequencesMode = flags.ContainsKey("sequences");   // -PROFILE and -SEQUENCES both set profileType in this order: -SEQUENCES wins
        bool profileCommand = profile2 != null;

        string? newTree = null;
        bool doUseTree = false, doGuideTreeOnly = false;
        if (flags.TryGetValue("usetree", out var useTree))
        {
            RequireInput(inFile, "Cannot align sequences.  No input file");
            if (useTree.Length == 0)
                throw new ArgumentException("Cannot align sequences.  No tree file specified", nameof(arguments));
            run = run with { UseTree = ReadFile(files, useTree) };
            doUseTree = true;
        }

        if (flags.TryGetValue("newtree", out var nt))
        {
            RequireInput(inFile, "Cannot align sequences.  No input file");
            if (nt.Length == 0)
                throw new ArgumentException("Cannot align sequences.  No tree file specified", nameof(arguments));
            newTree = nt;
            doGuideTreeOnly = true;
        }

        foreach (var (flag, set) in new (string, Action<string>)[]
        {
            ("usetree1", v => { run = run with { UseTree1 = ReadFile(files, v) }; doUseTree = true; }),
            ("usetree2", v => { run = run with { UseTree2 = ReadFile(files, v) }; doUseTree = true; }),
            ("newtree1", v => { newTree = v; run = run with { GuideTreeOnly = true }; }),
            ("newtree2", v => { newTree = v; run = run with { GuideTreeOnly = true }; }),
        })
        {
            if (!flags.TryGetValue(flag, out var v))
                continue;
            if (profile1 == null)
                throw new ArgumentException("Cannot align profiles.  No input file", nameof(arguments));
            if (profileCommand && sequencesMode)
                throw new ArgumentException($"Error: invalid option -{flag}", nameof(arguments));
            if (v.Length == 0 && flag.StartsWith("use", StringComparison.Ordinal))
                throw new ArgumentException("Cannot align profiles.  No tree file specified", nameof(arguments));
            set(v);
        }

        bool doTree = flags.ContainsKey("tree"), doBootstrap = flags.ContainsKey("bootstrap"), doConvert = flags.ContainsKey("convert"), doAlign = flags.ContainsKey("align");
        if (!doTree && !doBootstrap && inFile != null && !profileCommand && !doUseTree && !doGuideTreeOnly && !doConvert)
            doAlign = true;
        if (doTree && inFile == null)
            throw new ArgumentException("Cannot draw tree.  No input alignment file", nameof(arguments));
        if (doBootstrap && inFile == null)
            throw new ArgumentException("Cannot bootstrap tree. No input alignment file", nameof(arguments));
        if (doAlign && inFile == null)
            throw new ArgumentException("Cannot align sequences.  No input file", nameof(arguments));

        // the else-if chain of CommandLineParser::parseParams
        ClustalWVerb verb;
        if (profileCommand)
            verb = sequencesMode ? ClustalWVerb.Sequences : ClustalWVerb.Profile;
        else if (doAlign)
            verb = ClustalWVerb.Align;
        else if (doConvert)
            verb = ClustalWVerb.Convert;
        else if (doUseTree)
            verb = ClustalWVerb.AlignUsingTree;
        else if (doGuideTreeOnly)
            verb = ClustalWVerb.GuideTree;
        else if (doTree)
            verb = ClustalWVerb.Tree;
        else
            verb = doBootstrap ? ClustalWVerb.Bootstrap : ClustalWVerb.Align;
        if (verb == ClustalWVerb.GuideTree)
            run = run with { GuideTreeOnly = true };
        run = run with { Output = output, Structure = structure };
        return new ClustalWCommand(verb, type, dna, protein, run, inFile, profile1, profile2, flags.GetValueOrDefault("outfile"), newTree, messages);
    }

    private static void RequireInput(string? inFile, string message)
    {
        if (inFile == null)
            throw new ArgumentException(message, nameof(inFile));
    }

    private static string ReadFile(IReadOnlyDictionary<string, string>? files, string name)
    {
        if (files != null && files.TryGetValue(name, out var text))
            return text;
        throw new ArgumentException($"Cannot open file [{name}]: its text is not among the files given.", nameof(files));
    }

    private static ClustalWUserMatrix ParseMatrixFile(IReadOnlyDictionary<string, string>? files, string name)
    {
        try
        {
            return ClustalWUserMatrix.Parse(ReadFile(files, name), name);
        }
        catch (FormatException e)
        {
            throw new ArgumentException(e.Message, nameof(files), e);
        }
    }

    private static ClustalWMatrixFile ParseMatrixSeries(IReadOnlyDictionary<string, string>? files, string name)
    {
        try
        {
            return ClustalWMatrixFile.Parse(ReadFile(files, name), n => ReadFile(files, n), name);
        }
        catch (FormatException e)
        {
            throw new ArgumentException(e.Message, nameof(files), e);
        }
    }

    private static bool TryProteinPairwiseMatrix(string name, out ClustalWPairwiseMatrix matrix)
    {
        (bool ok, matrix) = name switch
        {
            "blosum" => (true, ClustalWPairwiseMatrix.Blosum),
            "pam" => (true, ClustalWPairwiseMatrix.Pam),
            "gonnet" => (true, ClustalWPairwiseMatrix.Gonnet),
            "id" => (true, ClustalWPairwiseMatrix.Identity),
            _ => (false, ClustalWPairwiseMatrix.Gonnet),
        };
        return ok;
    }

    private static bool TryProfileMatrix(string name, out ClustalWProfileMatrix matrix)
    {
        (bool ok, matrix) = name switch
        {
            "blosum" => (true, ClustalWProfileMatrix.Blosum),
            "pam" => (true, ClustalWProfileMatrix.Pam),
            "gonnet" => (true, ClustalWProfileMatrix.Gonnet),
            "id" => (true, ClustalWProfileMatrix.Identity),
            _ => (false, ClustalWProfileMatrix.Gonnet),
        };
        return ok;
    }

    private static ClustalWStructureOptions ApplyStructureNumber(IReadOnlyDictionary<string, string> flags, string flag, int min, int max, ClustalWStructureOptions s, Func<ClustalWStructureOptions, int, ClustalWStructureOptions> set)
    {
        if (!flags.TryGetValue(flag, out var text))
            return s;
        int v = ScanInt(text, flag, "integer", 0);
        return v >= min && v <= max ? set(s, v) : s;
    }

    /// <summary><c>sscanf(text, "%d", &amp;v)</c>: leading blanks, a sign and digits; the rest is ignored; no digits is the program's "Bad option" refusal.</summary>
    private static int ScanInt(string text, string option, string expected, int fallback)
    {
        if (text.Length == 0)
            return fallback;
        var m = Regex.Match(text, @"^\s*([+-]?\d+)");
        if (!m.Success)
            throw new ArgumentException($"Bad option: -{option} (expected {expected})", nameof(text));
        return long.TryParse(m.Groups[1].Value, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long v) ? (int)Math.Clamp(v, int.MinValue, int.MaxValue) : fallback;
    }

    /// <summary><c>sscanf(text, "%f", &amp;v)</c> (a single-precision value, widened).</summary>
    private static double ScanFloat(string text, string option, string expected)
    {
        if (text.Length == 0)
            return 0.0;
        var m = Regex.Match(text, @"^\s*([+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?)");
        if (!m.Success)
            throw new ArgumentException($"Bad option: -{option} (expected {expected})", nameof(text));
        return (double)float.Parse(m.Groups[1].Value, NumberStyles.Float, CultureInfo.InvariantCulture);
    }

    /// <summary><c>CommandLineParser::findMatch</c>: the index of the list entry that the characters of <paramref name="probe"/> single out position by position (-1 when none).</summary>
    internal static int FindMatch(string probe, IReadOnlyList<string> list)
    {
        int match = 0;
        for (int i = 0; i < probe.Length; i++)
        {
            int count = 0;
            for (int j = 0; j < list.Count; j++)
            {
                char c = i < list[j].Length ? list[j][i] : '\0';
                if (probe[i] == c)
                {
                    match = j;
                    count++;
                }
            }

            if (count == 0)
                return -1;
            if (count == 1)
                return match;
        }

        return -1;
    }
}
