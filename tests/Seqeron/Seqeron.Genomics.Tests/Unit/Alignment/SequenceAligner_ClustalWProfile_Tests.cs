using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14, B09 F103): ClustalW 2.1 profile alignment commands, <c>-PROFILE1=a -PROFILE2=b -PROFILE</c> (<see cref="SequenceAligner.AlignClustalWProfiles(IEnumerable{string}, IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/>,
/// <c>Clustal::profileAlign</c>) and <c>-SEQUENCES</c> (<see cref="SequenceAligner.AlignClustalWSequencesToProfile(IEnumerable{string}, IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/>,
/// <c>Clustal::sequencesAlignToProfile</c> / <c>MSA::seqsAlignToProfile</c>) against the installed <c>clustalw</c> 2.1.
/// Oracle (TestData/Alignment/clustalw_profile_oracle.tsv, scripts/alignment/generate_clustalw_profile_oracle.py): 1 200 sets (300 each of profile / sequences x DNA / protein; the profiles are
/// <c>clustalw -ALIGN</c> alignments with gaps, 1-13 sequences) with -GAPOPEN -GAPEXT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE -PWGAPOPEN -PWGAPEXT -PWMATRIX -TRANSWEIGHT -MATRIX -NOPGAP -NOHGAP -GAPDIST;
/// every output row must be identical.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWProfile_Tests
{
    private sealed record Row(string Id, string Mode, string Type, string[] Profile1, string[] Profile2, string[] Rows, string Options);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_profile_oracle.tsv";
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
            list.Add(new Row(f[0], f[1], f[2], f[3].Split(';'), f[4].Split(';'), f[5].Split(';'), f[6]));
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

    private static ClustalWPairwiseMatrix PairwiseMatrix(string name) => name switch
    {
        "GONNET" => ClustalWPairwiseMatrix.Gonnet,
        "BLOSUM" => ClustalWPairwiseMatrix.Blosum,
        "PAM" => ClustalWPairwiseMatrix.Pam,
        "ID" => ClustalWPairwiseMatrix.Identity,
        _ => throw new ArgumentException(name),
    };

    private static string[] Align(Row r)
    {
        var opts = r.Options.Length == 0 ? [] : r.Options.Split(',');
        if (r.Type == "dna")
        {
            var o = new ClustalWOptions();
            foreach (var kv in opts)
            {
                var p = kv.Split('=');
                switch (p[0])
                {
                    case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                    case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                    case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                    case "noweights": o = o with { NoWeights = true }; break;
                    case "negative": o = o with { NegativeMatrix = true }; break;
                    case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                    case "transweight": o = o with { TransitionWeight = float.Parse(p[1], CultureInfo.InvariantCulture) }; break;
                    case "pwgapopen": o = o with { PairwiseGapOpen = D(p[1]) }; break;
                    case "pwgapext": o = o with { PairwiseGapExtend = D(p[1]) }; break;
                    default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
                }
            }

            return (r.Mode == "profile"
                ? SequenceAligner.AlignClustalWProfiles(r.Profile1, r.Profile2, o)
                : SequenceAligner.AlignClustalWSequencesToProfile(r.Profile1, r.Profile2, o)).AlignedSequences.ToArray();
        }
        else
        {
            var o = new ClustalWProteinOptions();
            foreach (var kv in opts)
            {
                var p = kv.Split('=');
                switch (p[0])
                {
                    case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                    case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                    case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                    case "noweights": o = o with { NoWeights = true }; break;
                    case "negative": o = o with { NegativeMatrix = true }; break;
                    case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                    case "matrix": o = o with { Matrix = ProfileMatrix(p[1]) }; break;
                    case "pwmatrix": o = o with { Pairwise = o.Pairwise with { Matrix = PairwiseMatrix(p[1]) } }; break;
                    case "pwgapopen": o = o with { Pairwise = o.Pairwise with { GapOpen = D(p[1]) } }; break;
                    case "pwgapext": o = o with { Pairwise = o.Pairwise with { GapExtend = D(p[1]) } }; break;
                    case "nopgap": o = o with { NoPreferredGapPenalties = true }; break;
                    case "nohgap": o = o with { NoHydrophilicGapPenalties = true }; break;
                    case "gapdist": o = o with { GapSeparationDistance = int.Parse(p[1]) }; break;
                    default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
                }
            }

            return (r.Mode == "profile"
                ? SequenceAligner.AlignClustalWProfiles(r.Profile1, r.Profile2, o)
                : SequenceAligner.AlignClustalWSequencesToProfile(r.Profile1, r.Profile2, o)).AlignedSequences.ToArray();
        }
    }

    [Test]
    public void Oracle_Has300SetsOfEachModeAndType()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Count(r => r.Mode == "profile" && r.Type == "dna"), Is.EqualTo(300));
            Assert.That(rows.Count(r => r.Mode == "profile" && r.Type == "protein"), Is.EqualTo(300));
            Assert.That(rows.Count(r => r.Mode == "sequences" && r.Type == "dna"), Is.EqualTo(300));
            Assert.That(rows.Count(r => r.Mode == "sequences" && r.Type == "protein"), Is.EqualTo(300));
        });
    }

    [Test]
    public void ProfileAndSequencesAlignment_Rows_EqualClustalW21_OnAllOracleSets()
    {
        int bad = 0;
        var first = new List<string>();
        foreach (var r in Oracle.Value)
        {
            string[] got = Align(r);
            if (!got.SequenceEqual(r.Rows))
            {
                bad++;
                if (first.Count < 3)
                    first.Add($"set {r.Id} ({r.Mode}, {r.Type}, {r.Options}): expected {string.Join('/', r.Rows)}, got {string.Join('/', got)}");
            }
        }

        Assert.That(bad, Is.Zero, string.Join("\n", first));
    }

    [Test]
    public void Result_KeepsBothInputsInOrder_AndAddsOnlyGaps()
    {
        string[] p1 = ["ACGTACGT-ACGT", "ACGTACGTAACGT"];
        string[] p2 = ["ACGTTCGTACGT", "ACGTTCGTACCT"];
        var r = SequenceAligner.AlignClustalWProfiles(p1, p2);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Has.Length.EqualTo(4));
            Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(p1.Concat(p2).Select(x => x.Replace("-", ""))));
            Assert.That(r.AlignedSequences.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
        });
        var s = SequenceAligner.AlignClustalWSequencesToProfile(p1, p2);
        Assert.That(s.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(p1.Concat(p2).Select(x => x.Replace("-", ""))));
    }

    [Test]
    public void InvalidInput_IsRefused()
    {
        Assert.Throws<ArgumentException>(() => SequenceAligner.AlignClustalWProfiles(Array.Empty<string>(), ["ACGT"]));
        Assert.Throws<ArgumentException>(() => SequenceAligner.AlignClustalWProfiles(["ACGT"], Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => SequenceAligner.AlignClustalWProfiles(["ACGT"], ["AC*T"]));
        Assert.Throws<ArgumentException>(() => SequenceAligner.AlignClustalWSequencesToProfile(["ACGT"], Array.Empty<string>()));
        Assert.Throws<NotSupportedException>(() => SequenceAligner.AlignClustalWProfiles(["ACGT"], ["ACGA"], new ClustalWOptions { Iteration = ClustalWIteration.Tree }));
        Assert.Throws<NotSupportedException>(() => SequenceAligner.AlignClustalWSequencesToProfile(["ACGT"], ["ACGA"], new ClustalWProteinOptions { Iteration = ClustalWIteration.Tree }));
    }
}
