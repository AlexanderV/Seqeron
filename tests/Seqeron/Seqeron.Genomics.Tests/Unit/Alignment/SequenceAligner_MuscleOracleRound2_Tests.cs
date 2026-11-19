using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10c part a2, B09 F70 / F76): round-2 oracle of the MUSCLE 3.8.1551 command lines and options of
/// <see cref="SequenceAligner.MuscleRun"/>, <see cref="SequenceAligner.MuscleProfileAlign"/>, <see cref="SequenceAligner.MuscleRefine"/>, <see cref="SequenceAligner.MuscleScoreAlignment"/>,
/// <see cref="SequenceAligner.MuscleProfileScore"/>, <see cref="SequenceAligner.MuscleGuideTree"/>, <see cref="SequenceAligner.MuscleRefineWindows"/> and <see cref="SequenceAligner.MuscleLocalAlign"/>
/// against MUSCLE v3.8.1551 (the unmodified Debian muscle3 3.8.1551-2 source, x86_64, Debian 12; the Docker package docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2 run by the user,
/// archive returned/A1-10c-muscle-r2-out.zip, converted by scripts/alignment/generate_muscle3_r2_oracle.py). 630 of the 638 sets (the 8 sets of -matrix blosum62nostar, whose console output the
/// harness could not decode, are re-run in round 3): every option of the algorithm (iterations 3-16 and the default 8, anchors,
/// objective scores, the six weightings, distances, clusterings, roots, -brenner, -matrix, -usetree, the tree / distance / weight files, -clusteronly, the compiled-out options) and the
/// command lines -profile, -refine, -spscore, -maketree, -ppscore, -refinew, -sw, -profdb. Every set must give MUSCLE's rows in MUSCLE's written order, MUSCLE's files, MUSCLE's printed score
/// - or, where the program refuses (31 runs: error exit, crash), the matching refusal of the port.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MuscleOracleRound2_Tests
{
    internal sealed record Row(string Id, string Kind, string Type, string Options, string Data, string Rc, string[] Names, string Rows, string Tree1, string Tree2,
        string Dist1, string Dist2, string Weights, string UsedTree, string Stdout, string MatrixText);

    private static string Unesc(string s) => s.Contains('\\')
        ? System.Text.RegularExpressions.Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", var c => c })
        : s;

    internal static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.muscle3_r2_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], f[3], f[4], f[5], f[6].Length == 0 ? [] : f[6].Split(','), Unesc(f[7]), Unesc(f[8]), Unesc(f[9]), Unesc(f[10]), Unesc(f[11]), Unesc(f[12]),
                Unesc(f[13]), Unesc(f[14]), Unesc(f[15])));
        }

        return list;
    });

    public static IEnumerable<TestCaseData> Sets(params string[] kinds) =>
        Oracle.Value.Where(r => kinds.Contains(r.Kind)).Select(r => new TestCaseData(r.Id).SetArgDisplayNames(r.Id));

    public static IEnumerable<TestCaseData> AlignSets() => Sets("align");
    public static IEnumerable<TestCaseData> ModeSets() => Sets("profile", "refine", "spscore", "maketree", "ppscore", "refinew", "sw", "profdb");

    private static double D(string v) => double.Parse(v, CultureInfo.InvariantCulture);

    private static MuscleClusterMethod Cluster(string v) => v switch
    {
        "upgmb" => MuscleClusterMethod.Upgmb,
        "upgma" => MuscleClusterMethod.Upgma,
        "upgmamin" => MuscleClusterMethod.UpgmaMin,
        "upgmamax" => MuscleClusterMethod.UpgmaMax,
        "neighborjoining" => MuscleClusterMethod.NeighborJoining,
        _ => throw new ArgumentException(v),
    };

    private static MuscleRoot Root(string v) => v switch
    {
        "pseudo" => MuscleRoot.Pseudo,
        "midlongestspan" => MuscleRoot.MidLongestSpan,
        "minavgleafdist" => MuscleRoot.MinAvgLeafDist,
        _ => throw new ArgumentException(v),
    };

    private static MuscleDistance Dist(string v) => v switch
    {
        "kmer6_6" => MuscleDistance.Kmer6_6,
        "kmer4_6" => MuscleDistance.Kmer4_6,
        "kmer20_3" => MuscleDistance.Kmer20_3,
        "kmer20_4" => MuscleDistance.Kmer20_4,
        "kbit20_3" => MuscleDistance.Kbit20_3,
        "pctidkimura" => MuscleDistance.PctIdKimura,
        "pctidlog" => MuscleDistance.PctIdLog,
        "pwkimura" => MuscleDistance.PwKimura,
        "pwscoredist" => MuscleDistance.PwScoreDist,
        "scoredist" => MuscleDistance.ScoreDist,
        "edit" => MuscleDistance.Edit,
        _ => throw new ArgumentException(v),
    };

    private static MuscleSeqWeight Weight(string v) => v switch
    {
        "clustalw" => MuscleSeqWeight.ClustalW,
        "none" => MuscleSeqWeight.None,
        "henikoff" => MuscleSeqWeight.Henikoff,
        "henikoffpb" => MuscleSeqWeight.HenikoffPb,
        "gsc" => MuscleSeqWeight.Gsc,
        "threeway" => MuscleSeqWeight.ThreeWay,
        _ => throw new ArgumentException(v),
    };

    /// <summary>The option string of the Docker package (<c>key=value</c>, comma separated; the key is the program's switch) as <see cref="MuscleOptions"/>; <c>-maxiters</c> absent = the program's default 8.</summary>
    internal static MuscleOptions Parse(Row r, IReadOnlyList<string>? names)
    {
        var o = new MuscleOptions
        {
            SequenceType = r.Type switch { "protein" => MuscleSequenceType.Protein, "dna" => MuscleSequenceType.Dna, "rna" => MuscleSequenceType.Rna, _ => MuscleSequenceType.Auto },
            MaxIterations = 8,
            TreeOrder = true,
            SequenceNames = names,
        };
        if (r.Options.Length == 0)
            return o;
        foreach (var kv in r.Options.Split(','))
        {
            var p = kv.Split('=', 2);
            var v = p.Length > 1 ? p[1] : "";
            o = p[0] switch
            {
                "maxiters" => o with { MaxIterations = int.Parse(v) },
                "maxtrees" => o with { MaxTrees = int.Parse(v) },
                "gapopen" => o with { GapOpen = D(v) },
                "gapextend" => o with { GapExtend = D(v) },
                "center" => o with { Center = D(v) },
                "sueff" => o with { Sueff = float.Parse(v, CultureInfo.InvariantCulture) },
                "distance1" => o with { Distance1 = Dist(v) },
                "distance2" => o with { Distance2 = Dist(v) },
                "cluster1" => o with { Cluster1 = Cluster(v) },
                "cluster2" => o with { Cluster2 = Cluster(v) },
                "root1" => o with { Root1 = Root(v) },
                "root2" => o with { Root2 = Root(v) },
                "weight" => o with { Weight1 = Weight(v), Weight2 = Weight(v) },
                "weight1" => o with { Weight1 = Weight(v) },
                "weight2" => o with { Weight2 = Weight(v) },
                "objscore" => o with { ObjectiveScore = v switch { "spm" => MuscleObjectiveScore.Spm, "sp" => MuscleObjectiveScore.Sp, "dp" => MuscleObjectiveScore.Dp, "xp" => MuscleObjectiveScore.Xp, "ps" => MuscleObjectiveScore.Ps, "spf" => MuscleObjectiveScore.Spf, _ => throw new ArgumentException(kv) } },
                "termgaps" => o with { TermGaps = v switch { "half" => MuscleTermGaps.Half, "full" => MuscleTermGaps.Full, "ext" => MuscleTermGaps.Ext, _ => throw new ArgumentException(kv) } },
                "noanchors" => o with { Anchors = false },
                "anchorspacing" => o with { AnchorSpacing = int.Parse(v) },
                "smoothwindow" => o with { SmoothWindow = int.Parse(v) },
                "smoothscoreceil" => o with { SmoothScoreCeil = D(v) },
                "minbestcolscore" => o with { MinBestColScore = D(v) },
                "minsmoothscore" => o with { MinSmoothScore = D(v) },
                "brenner" => o with { Brenner = true },
                "matrix" => o with { UserMatrix = r.MatrixText },
                "usetree" => o with { UserTree = r.UsedTree },
                "clusteronly" => o with { ClusterOnly = true },
                "computeweights" => o with { ComputeWeights = true },
                "le" => o with { ScoreModel = MuscleScoreModel.LogExpectation },
                "sp" => o with { ScoreModel = MuscleScoreModel.SumOfPairs },
                "sv" => o with { ScoreModel = MuscleScoreModel.SumOfPairsVtml },
                "spn" => o with { ScoreModel = MuscleScoreModel.SumOfPairsNucleotide },
                "diags" => o with { Diags1 = true, Diags2 = true },
                "diags1" => o with { Diags1 = true },
                "diags2" => o with { Diags2 = true },
                "diaglength" => o with { DiagLength = int.Parse(v) },
                "diagbreak" => o with { DiagBreak = int.Parse(v) },
                "diagmargin" => o with { DiagMargin = int.Parse(v) },
                "hydro" => o with { Hydro = int.Parse(v) },
                "hydrofactor" => o with { HydroFactor = D(v) },
                "tomhydro" => o with { TomHydro = true },
                "gapopen2" => o with { GapOpen2 = D(v) },
                "gapextend2" => o with { GapExtend2 = D(v) },
                "gapambig" => o with { GapAmbig = D(v) },
                "maxmb" => o with { MaxMb = int.Parse(v) },
                "maxhours" => o with { MaxHours = D(v) },
                "refinewindow" => o with { RefineWindow = int.Parse(v) },
                "stable" => o with { Stable = true },
                "tree1" or "tree2" or "distmx1" or "distmx2" or "core" or "nocore" or "group" or "dimer" or "verbose" => o,   // files are always produced by the port; the others have no effect
                _ => throw new ArgumentException(kv),
            };
        }

        return o;
    }

    private static string[] Split(string s) => s.Length == 0 ? [] : s.Split(';');

    /// <summary>
    /// The 8 sets of the package that are NOT in the oracle: <c>-matrix blosum62nostar</c> (R330, R335, R336, R342, R344, R347, R348, R350), whose console output the harness could not decode
    /// (RUNNER_ERROR); they are re-run in round 3 (A2-1, handoff/A1-10c-muscle-r3). Naming them makes a row dropped silently anywhere else (or a regenerated TSV with another gap) fail this test;
    /// when round 3 returns, this list becomes empty and the count 638.
    /// </summary>
    private static readonly string[] MissingUntilRound3 = ["R330", "R335", "R336", "R342", "R344", "R347", "R348", "R350"];

    [Test]
    public void OracleHas630Sets_TheEightMissingIdsAreExactlyTheBlosum62NoStarOnes()
    {
        // the package has the 638 sets R0 ... R637; the oracle has every one of them except the eight round-3 sets
        var expected = Enumerable.Range(0, 638).Select(i => "R" + i).Except(MissingUntilRound3).ToHashSet();
        var oracleIds = Oracle.Value.Select(r => r.Id).ToHashSet();
        Assert.That(expected.Except(oracleIds), Is.Empty, "a set of the package is missing from the oracle although it is not one of the 8 round-3 sets");
        Assert.That(oracleIds.Except(expected), Is.Empty, "the oracle has a row that is not a set of the package (or one of the 8 round-3 sets came back: empty MissingUntilRound3)");
    }

    [Test]
    public void OracleHas630Sets_AllKindsPresent()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Oracle.Value, Has.Count.EqualTo(630));
            Assert.That(Oracle.Value.Count(r => r.Kind == "align"), Is.EqualTo(441));
            Assert.That(Oracle.Value.Count(r => r.Kind == "profile"), Is.EqualTo(36));
            Assert.That(Oracle.Value.Count(r => r.Kind == "refine"), Is.EqualTo(35));
            Assert.That(Oracle.Value.Count(r => r.Kind == "spscore"), Is.EqualTo(22));
            Assert.That(Oracle.Value.Count(r => r.Kind == "maketree"), Is.EqualTo(23));
            Assert.That(Oracle.Value.Count(r => r.Kind == "ppscore"), Is.EqualTo(27));
            Assert.That(Oracle.Value.Count(r => r.Kind == "refinew"), Is.EqualTo(18));
            Assert.That(Oracle.Value.Count(r => r.Kind == "sw"), Is.EqualTo(23));
            Assert.That(Oracle.Value.Count(r => r.Kind == "profdb"), Is.EqualTo(5));
        });
    }

    /// <summary>
    /// The refusal of the port that models a failed run of the program: a crash (segmentation fault, exit status 139 = -11) and "-stable not supported" -> <see cref="NotSupportedException"/>;
    /// an invalid distance for the stage ("Invalid DISTANCE_n") -> <see cref="ArgumentException"/>; the other messages of the program's own Quit() ("Clust::SetLeafCount(1)", "g_MuscleWeights = 0")
    /// -> <see cref="InvalidOperationException"/>.
    /// </summary>
    private static void AssertRefused(Row r, TestDelegate run)
    {
        string why = $"{r.Id} {r.Kind} {r.Type} {r.Options}: oracle rc={r.Rc} {r.Stdout}";
        if (r.Rc == "FAILED-11" || r.Stdout.Contains("not supported in this version"))
            Assert.Throws<NotSupportedException>(run, why);
        else if (r.Stdout.Contains("nvalid DISTANCE"))
            Assert.Throws<ArgumentException>(run, why);
        else
            Assert.Throws<InvalidOperationException>(run, why);
    }

    [TestCaseSource(nameof(AlignSets))]
    public void Align_RowsWrittenOrderAndFiles_EqualMuscle3(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var seqs = Split(r.Data);
        var names = Enumerable.Range(0, seqs.Length).Select(i => "s" + i).ToArray();
        var o = Parse(r, names);
        if (r.Rc != "0")
        {
            AssertRefused(r, () => SequenceAligner.MuscleRun(seqs, o));
            return;
        }

        var res = SequenceAligner.MuscleRun(seqs, o);
        Assert.Multiple(() =>
        {
            Assert.That(res.Alignment.AlignedSequences.ToArray(), Is.EqualTo(Split(r.Rows)), $"{id} {r.Type} {r.Options}: rows");
            Assert.That(res.OutputOrder.Select(i => names[i]).ToArray(), Is.EqualTo(r.Names), $"{id}: written order");
            if (r.Tree1.Length > 0)
                Assert.That(res.Tree1, Is.EqualTo(r.Tree1), $"{id}: -tree1");
            if (r.Tree2.Length > 0)
                Assert.That(res.Tree2, Is.EqualTo(r.Tree2), $"{id}: -tree2");
            if (r.Dist1.Length > 0)
                Assert.That(res.DistanceMatrix1, Is.EqualTo(r.Dist1), $"{id}: -distmx1");
            if (r.Dist2.Length > 0)
                Assert.That(res.DistanceMatrix2, Is.EqualTo(r.Dist2), $"{id}: -distmx2");
            if (r.Weights.Length > 0)
                Assert.That(res.WeightsFile, Is.EqualTo(r.Weights), $"{id}: -computeweights");
        });
    }

    /// <summary>The names the program uses in the profile modes: 'a0..' / 'b0..' (FASTA labels of the two input files).</summary>
    private static string[] Labelled(string prefix, int n) => Enumerable.Range(0, n).Select(i => prefix + i).ToArray();

    [TestCaseSource(nameof(ModeSets))]
    public void CommandLine_EqualsMuscle3(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var parts = r.Data.Split('|');
        var a = Split(parts[0]);
        var b = parts.Length > 1 ? Split(parts[1]) : [];
        string why = $"{id} {r.Kind} {r.Type} {r.Options}";
        var o = Parse(r, Labelled("s", a.Length));
        TestDelegate run = r.Kind switch
        {
            "profile" => () => SequenceAligner.MuscleProfileAlign(a, b, o),
            "refine" => () => SequenceAligner.MuscleRefine(a, o),
            "spscore" => () => SequenceAligner.MuscleScoreAlignment(a, o),
            "maketree" => () => SequenceAligner.MuscleGuideTree(a, o),
            "ppscore" => () => SequenceAligner.MuscleProfileScore(a, b, o),
            "refinew" => () => SequenceAligner.MuscleRefineWindows(a, o),
            "sw" => () => SequenceAligner.MuscleLocalAlign(a, b, o),
            "profdb" => () => SequenceAligner.MuscleProfileDatabase(a, b, o),
            _ => throw new ArgumentException(r.Kind),
        };
        if (r.Rc != "0")
        {
            AssertRefused(r, run);
            return;
        }

        switch (r.Kind)
        {
            case "profile":
                Assert.That(SequenceAligner.MuscleProfileAlign(a, b, o), Is.EqualTo(Split(r.Rows)), why);
                Assert.That(r.Names, Is.EqualTo(Labelled("a", a.Length).Concat(Labelled("b", b.Length))), why);
                break;
            case "refine":
                Assert.That(SequenceAligner.MuscleRefine(a, o), Is.EqualTo(Split(r.Rows)), why);
                break;
            case "spscore":
                Assert.That(Printed(SequenceAligner.MuscleScoreAlignment(a, o)), Is.EqualTo(PrintedOf(r.Stdout, "SP=")), why);
                break;
            case "maketree":
                Assert.That(SequenceAligner.MuscleGuideTree(a, o), Is.EqualTo(r.Rows), why);
                break;
            case "ppscore":
                Assert.That(Printed(SequenceAligner.MuscleProfileScore(a, b, o)), Is.EqualTo(PrintedOf(r.Stdout, "Score=")), why);
                break;
            case "refinew":
                Assert.That(SequenceAligner.MuscleRefineWindows(a, o), Is.EqualTo(Split(r.Rows)), why);
                break;
            case "sw":
                Assert.That(SequenceAligner.MuscleLocalAlign(a, b, o), Is.EqualTo(Split(r.Rows)), why);
                break;
            default:
                Assert.Fail(why);
                break;
        }
    }

    /// <summary>"%.4g" of the program.</summary>
    private static string Printed(double v) => SequenceAligner.MuscleEngine.FormatG(v, 4);

    private static string PrintedOf(string stdout, string key)
    {
        int i = stdout.IndexOf(key, StringComparison.Ordinal);
        Assert.That(i, Is.GreaterThanOrEqualTo(0), "no " + key + " in the program output: " + stdout);
        string rest = stdout[(i + key.Length)..];
        int end = rest.IndexOfAny(['\n', '\r', ' ', ';']);
        return end < 0 ? rest : rest[..end];
    }
}
