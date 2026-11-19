using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14, B09 F102): ClustalW 2.1 <c>-ITERATION=TREE | ALIGNMENT</c> / <c>-NUMITER</c> (<c>Iteration::removeFirstIterate</c> / <c>iterationOnTreeNode</c>,
/// <c>ObjectiveScore</c>, <c>MSA::calcPairwiseForProfileAlign</c> / <c>doProfileAlign</c>, <c>TreeInterface::getWeightsForProfileAlignNJ</c>) against the installed <c>clustalw</c> 2.1.
/// Oracle (TestData/Alignment/clustalw_iteration_oracle.tsv, scripts/alignment/generate_clustalw_iteration_oracle.py): 700 sets (350 DNA, 350 protein; 3-14 sequences up to ~140 residues)
/// with <c>-ITERATION=tree|alignment</c>, <c>-NUMITER</c> and the usual option variants; the FASTA rows AND the "FINAL score: n" lines (the ObjectiveScore of the best re-alignment of every
/// removeFirstIterate call) must be identical.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWIteration_Tests
{
    private sealed record Row(string Id, string Type, string[] Seqs, string[] Rows, string Options, long[] Finals);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_iteration_oracle.tsv";
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
            list.Add(new Row(f[0], f[1], f[2].Split(';'), f[3].Split(';'), f[4], f[5].Length == 0 ? [] : f[5].Split(',').Select(long.Parse).ToArray()));
        }

        return list;
    });

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static ClustalWProfileMatrix ProfileMatrix(string name) => name switch
    {
        "GONNET" => ClustalWProfileMatrix.Gonnet,
        "BLOSUM" => ClustalWProfileMatrix.Blosum,
        "PAM" => ClustalWProfileMatrix.Pam,
        "ID" => ClustalWProfileMatrix.Identity,
        _ => throw new ArgumentException(name),
    };

    private static ClustalWIteration Mode(string word) => word switch
    {
        "tree" => ClustalWIteration.Tree,
        "alignment" => ClustalWIteration.Alignment,
        _ => throw new ArgumentException(word),
    };

    private static ClustalWOptions DnaOptions(Row r, bool withIteration = true)
    {
        var o = new ClustalWOptions();
        foreach (var kv in r.Options.Length == 0 ? [] : r.Options.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "iteration" when withIteration: o = o with { Iteration = Mode(p[1]) }; break;
                case "numiter" when withIteration: o = o with { NumIterations = int.Parse(p[1]) }; break;
                case "iteration" or "numiter": break;
                case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "noweights": o = o with { NoWeights = true }; break;
                case "negative": o = o with { NegativeMatrix = true }; break;
                case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                case "transweight": o = o with { TransitionWeight = float.Parse(p[1], CultureInfo.InvariantCulture) }; break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static ClustalWProteinOptions ProteinOptions(Row r, bool withIteration = true)
    {
        var o = new ClustalWProteinOptions();
        foreach (var kv in r.Options.Length == 0 ? [] : r.Options.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "iteration" when withIteration: o = o with { Iteration = Mode(p[1]) }; break;
                case "numiter" when withIteration: o = o with { NumIterations = int.Parse(p[1]) }; break;
                case "iteration" or "numiter": break;
                case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "noweights": o = o with { NoWeights = true }; break;
                case "negative": o = o with { NegativeMatrix = true }; break;
                case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                case "matrix": o = o with { Matrix = ProfileMatrix(p[1]) }; break;
                case "nopgap": o = o with { NoPreferredGapPenalties = true }; break;
                case "nohgap": o = o with { NoHydrophilicGapPenalties = true }; break;
                case "gapdist": o = o with { GapSeparationDistance = int.Parse(p[1]) }; break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static (string[] Rows, IReadOnlyList<long> Finals) Run(Row r) => r.Type == "dna"
        ? SequenceAligner.MultipleAlignClustalWWithIterationScores(r.Seqs, DnaOptions(r))
        : SequenceAligner.MultipleAlignClustalWWithIterationScores(r.Seqs, ProteinOptions(r));

    private static string[] RunWithoutIteration(Row r) => r.Type == "dna"
        ? SequenceAligner.MultipleAlignClustalW(r.Seqs, DnaOptions(r, withIteration: false)).AlignedSequences.ToArray()
        : SequenceAligner.MultipleAlignClustalW(r.Seqs, ProteinOptions(r, withIteration: false)).AlignedSequences.ToArray();

    [Test]
    public void Oracle_Has700Sets_OfBothTypesAndBothModes()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Count(r => r.Type == "dna"), Is.EqualTo(350));
            Assert.That(rows.Count(r => r.Type == "protein"), Is.EqualTo(350));
            Assert.That(rows.Count(r => r.Options.Contains("iteration=tree")), Is.GreaterThan(250));
            Assert.That(rows.Count(r => r.Options.Contains("iteration=alignment")), Is.GreaterThan(250));
        });
    }

    [Test]
    public void MultipleAlignClustalW_Iteration_Rows_And_FinalScores_EqualClustalW21_OnAllOracleSets()
    {
        int badRows = 0, badScores = 0, changed = 0;
        var first = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var (rows, finals) = Run(r);
            if (!rows.SequenceEqual(r.Rows))
            {
                badRows++;
                if (first.Count < 3)
                    first.Add($"set {r.Id} ({r.Type}, {r.Options}): rows differ");
            }

            if (!finals.SequenceEqual(r.Finals))
            {
                badScores++;
                if (first.Count < 3)
                    first.Add($"set {r.Id} ({r.Type}, {r.Options}): FINAL scores expected [{string.Join(',', r.Finals)}], got [{string.Join(',', finals)}]");
            }

            if (!r.Rows.SequenceEqual(RunWithoutIteration(r)))
                changed++;
        }

        Assert.Multiple(() =>
        {
            Assert.That(badRows, Is.Zero, string.Join("\n", first));
            Assert.That(badScores, Is.Zero, string.Join("\n", first));
            Assert.That(changed, Is.GreaterThan(300), "the iteration must change the rows of a large share of the oracle sets (else the test would not exercise it)");
        });
    }

    [Test]
    public void Default_IsNoIteration_AndThreeIterations()
    {
        Assert.Multiple(() =>
        {
            Assert.That(new ClustalWOptions().Iteration, Is.EqualTo(ClustalWIteration.None));
            Assert.That(new ClustalWOptions().NumIterations, Is.EqualTo(3));
            Assert.That(new ClustalWProteinOptions().Iteration, Is.EqualTo(ClustalWIteration.None));
            Assert.That(new ClustalWProteinOptions().NumIterations, Is.EqualTo(3));
        });
    }

    [Test]
    public void Options_AreValidated()
    {
        string[] seqs = ["ACGTACGTAA", "ACGTCCGTAA", "ACGAACGTAA"];
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions { Iteration = ClustalWIteration.Alignment, NumIterations = 0 }));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWProteinOptions { Iteration = ClustalWIteration.Tree, NumIterations = -1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions { Iteration = (ClustalWIteration)9 }));
    }

    [Test]
    public void TwoSequences_AreNeverIterated()
    {
        string[] seqs = ["ACGTACGTAAGG", "ACGTCCGTAAGG"];
        var plain = SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions()).AlignedSequences;
        foreach (var mode in new[] { ClustalWIteration.Tree, ClustalWIteration.Alignment })
            Assert.That(SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions { Iteration = mode }).AlignedSequences, Is.EqualTo(plain));
    }

    [Test]
    public void UpgmaClustering_WithIteration_KeepsEverySequenceIntact()
    {
        // -CLUSTERING=UPGMA -ITERATION (getWeightsForProfileAlignUPGMA) is ported from the source; the installed program crashes on every UPGMA run, so no row oracle exists yet
        // (docs/Validation/review-2026-09/handoff/A2-14-upgma-iteration): this only locks the integrity of the result.
        var rnd = new Random(5);
        for (int t = 0; t < 60; t++)
        {
            string template = new(Enumerable.Range(0, rnd.Next(20, 60)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var seqs = Enumerable.Range(0, rnd.Next(3, 8))
                .Select(_ => new string(template.Where(_ => rnd.NextDouble() > 0.1).Select(c => rnd.NextDouble() < 0.1 ? "ACGT"[rnd.Next(4)] : c).ToArray()))
                .Where(s => s.Length > 0).ToArray();
            if (seqs.Length < 3)
                continue;
            foreach (var mode in new[] { ClustalWIteration.Alignment, ClustalWIteration.Tree })
            {
                var rows = SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions { Clustering = ClustalWClustering.Upgma, Iteration = mode }).AlignedSequences;
                Assert.That(rows.Select(r => r.Replace("-", "")), Is.EqualTo(seqs));
                Assert.That(rows.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
            }
        }
    }
}
