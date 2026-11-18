using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10c, B09 F65 / F69): row-for-row oracle of <see cref="SequenceAligner.MultipleAlignMuscle(IEnumerable{string}, MuscleOptions?, MatrixScoring?)"/>
/// against MUSCLE v3.8.1551 (built unmodified from the Debian muscle3 3.8.1551-2 source, x86_64, Debian 12; the Docker package docs/Validation/review-2026-09/handoff/A1-10c-muscle,
/// run by the user, archive returned/A1-10c-muscle-out.zip, scripts/alignment/generate_muscle3_oracle.py). 420 sets (270 protein, 110 DNA, 40 RNA) with the option variants
/// -maxiters 1|2, -le/-sp/-sv, -gapopen/-gapextend/-center, -cluster1/-cluster2, -sueff, -distance1, -maxtrees 0, leading M, wildcards, lower case, duplicates, short sequences.
/// Every set must give MUSCLE's rows (input order) and MUSCLE's output order.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MuscleOracle_Tests
{
    internal sealed record Row(string Id, string Type, string Options, string[] Seqs, int[] Order, string[]? Rows);

    internal static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.muscle3_oracle.tsv";
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
            bool failed = f[5] == "FAILED";
            list.Add(new Row(f[0], f[1], f[2], f[3].Split(';'), failed ? [] : f[4].Split(',').Select(int.Parse).ToArray(), failed ? null : f[5].Split(';')));
        }

        return list;
    });

    public static IEnumerable<TestCaseData> Sets() => Oracle.Value.Select(r => new TestCaseData(r.Id).SetArgDisplayNames(r.Id));

    /// <summary>The oracle option string of the Docker package (<c>key=value</c>, comma separated) as <see cref="MuscleOptions"/>.</summary>
    internal static MuscleOptions Parse(string type, string options)
    {
        var o = new MuscleOptions
        {
            SequenceType = type switch { "protein" => MuscleSequenceType.Protein, "dna" => MuscleSequenceType.Dna, "rna" => MuscleSequenceType.Rna, _ => MuscleSequenceType.Auto },
        };
        if (options.Length == 0)
            return o;
        foreach (var kv in options.Split(','))
        {
            var p = kv.Split('=', 2);
            var v = p[1];
            o = p[0] switch
            {
                "maxiters" => o with { MaxIterations = int.Parse(v) },
                "maxtrees" => o with { MaxTrees = int.Parse(v) },
                "gapopen" => o with { GapOpen = double.Parse(v, System.Globalization.CultureInfo.InvariantCulture) },
                "gapextend" => o with { GapExtend = double.Parse(v, System.Globalization.CultureInfo.InvariantCulture) },
                "center" => o with { Center = double.Parse(v, System.Globalization.CultureInfo.InvariantCulture) },
                "sueff" => o with { Sueff = float.Parse(v, System.Globalization.CultureInfo.InvariantCulture) },
                "distance1" => o with { Distance1 = v switch { "kmer6_6" => MuscleDistance.Kmer6_6, "kmer4_6" => MuscleDistance.Kmer4_6, _ => throw new ArgumentException(kv) } },
                "cluster1" => o with { Cluster1 = Cluster(v) },
                "cluster2" => o with { Cluster2 = Cluster(v) },
                "model" => o with { ScoreModel = v switch { "le" => MuscleScoreModel.LogExpectation, "sp" => MuscleScoreModel.SumOfPairs, "sv" => MuscleScoreModel.SumOfPairsVtml, _ => throw new ArgumentException(kv) } },
                _ => throw new ArgumentException(kv),
            };
        }

        return o;
    }

    private static MuscleClusterMethod Cluster(string v) => v switch
    {
        "upgmb" => MuscleClusterMethod.Upgmb,
        "upgma" => MuscleClusterMethod.Upgma,
        "upgmamin" => MuscleClusterMethod.UpgmaMin,
        "upgmamax" => MuscleClusterMethod.UpgmaMax,
        _ => throw new ArgumentException(v),
    };

    [Test]
    public void OracleHas420Sets_AllRun()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Oracle.Value, Has.Count.EqualTo(420));
            Assert.That(Oracle.Value.Count(r => r.Type == "protein"), Is.EqualTo(270));
            Assert.That(Oracle.Value.Count(r => r.Type == "dna"), Is.EqualTo(110));
            Assert.That(Oracle.Value.Count(r => r.Type == "rna"), Is.EqualTo(40));
            Assert.That(Oracle.Value.All(r => r.Rows is not null), Is.True, "MUSCLE failed on a set");
        });
    }

    [TestCaseSource(nameof(Sets))]
    public void Rows_EqualMuscle3(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var o = Parse(r.Type, r.Options);
        var rows = SequenceAligner.MultipleAlignMuscle(r.Seqs, o).AlignedSequences.ToArray();
        Assert.That(rows, Is.EqualTo(r.Rows), $"{id} {r.Type} {r.Options}");
    }

    [TestCaseSource(nameof(Sets))]
    public void OutputOrder_EqualsMuscle3(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var o = Parse(r.Type, r.Options);
        var d = SequenceAligner.MultipleAlignMuscleDetailed(r.Seqs, o);
        Assert.That(d.TreeOrder, Is.EqualTo(r.Order), $"{id} {r.Type} {r.Options}");
    }
}
