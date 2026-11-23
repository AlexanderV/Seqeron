using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit R14-10, B09 F303): <see cref="PairwiseAligner"/> against Biopython 1.88 <c>PairwiseAligner</c> on 698 recorded cases: <c>mode = "fogsaa"</c> (the FOGSAA port gives Biopython's alignment itself, coordinates and
/// score, and its warnings), the <c>wildcard</c> letter, scalar and matrix scoring, linear / affine / end-specific / directional gap scores, both strands, the algorithm name. Global and local alignments are the same SET of optimal alignments
/// (the library enumerates in its own order); a configuration with end gap scores that are neither the interior score nor zero gives one optimal alignment of the set. <c>scripts/alignment/generate_pairwise_aligner_oracle.py</c>.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class PairwiseAligner_Oracle_Tests
{
    private static readonly Lazy<List<JsonObject>> Records = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.pairwise_aligner_oracle.jsonl";
        using var stream = TestDataResource.Open(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        var list = new List<JsonObject>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length > 0)
                list.Add((JsonObject)JsonNode.Parse(line)!);
        }

        return list;
    });

    private static IEnumerable<TestCaseData> Cases(string mode) =>
        Records.Value.Where(r => r["config"]!["mode"]!.GetValue<string>() == mode).Select(r => new TestCaseData(r.ToJsonString()).SetName($"{mode}_{r["id"]}_{r["config"]!["kind"]}"));

    private static IEnumerable<TestCaseData> GlobalCases() => Cases("global");

    private static IEnumerable<TestCaseData> LocalCases() => Cases("local");

    private static IEnumerable<TestCaseData> FogsaaCases() => Cases("fogsaa");

    private static PairwiseAligner Make(JsonNode config)
    {
        var aligner = new PairwiseAligner
        {
            Mode = config["mode"]!.GetValue<string>() switch
            {
                "local" => PairwiseAlignerMode.Local,
                "fogsaa" => PairwiseAlignerMode.Fogsaa,
                _ => PairwiseAlignerMode.Global,
            },
        };
        foreach (var (key, value) in config["props"]!.AsObject())
        {
            double v = value!.GetValue<double>();
            if (key == "match_score")
                aligner.MatchScore = v;
            else if (key == "mismatch_score")
                aligner.MismatchScore = v;
            else
                aligner.SetGapScore(key, v);
        }

        if (config["matrix"] is JsonNode m)
            aligner.SubstitutionMatrix = SubstitutionMatrix.Get(m.GetValue<string>());
        if (config["wildcard"] is JsonNode w)
            aligner.Wildcard = w.GetValue<string>()[0];
        return aligner;
    }

    private static string Key(long[][] coordinates) => string.Join("|", coordinates.Select(r => string.Join(",", r)));

    private static string Key(JsonNode n) => string.Join("|", n["coordinates"]!.AsArray().Select(r => string.Join(",", r!.AsArray().Select(x => x!.GetValue<long>()))));

    [TestCaseSource(nameof(GlobalCases))]
    public void Global_ScoreAlgorithmAndAlignments_EqualBiopython(string json) => Check(json, exactSet: true);

    [TestCaseSource(nameof(LocalCases))]
    public void Local_ScoreAlgorithmAndAlignments_EqualBiopython(string json) => Check(json, exactSet: true);

    [TestCaseSource(nameof(FogsaaCases))]
    public void Fogsaa_AlignmentScoreAndWarnings_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var aligner = Make(rec["config"]!);
        string a = rec["a"]!.GetValue<string>(), b = rec["b"]!.GetValue<string>();
        char strand = rec["strand"]!.GetValue<string>()[0];
        Assert.That(aligner.Algorithm, Is.EqualTo(rec["algorithm"]!.GetValue<string>()));
        double score = aligner.Score(a, b, strand);
        Assert.That(score, Is.EqualTo(rec["score"]!.GetValue<double>()).Within(1e-9));
        var alignments = aligner.Align(a, b, strand);
        Assert.That(alignments.Score, Is.EqualTo(rec["score"]!.GetValue<double>()).Within(1e-9));
        Assert.That(alignments, Has.Count.EqualTo(1));
        Assert.That(Key(alignments[0].Coordinates), Is.EqualTo(Key(rec["alignments"]![0]!)));
        Assert.That(aligner.Warnings.OrderBy(x => x, StringComparer.Ordinal), Is.EqualTo(rec["warnings"]!.AsArray().Select(x => x!.GetValue<string>()).OrderBy(x => x, StringComparer.Ordinal)));
    }

    private static void Check(string json, bool exactSet)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var aligner = Make(rec["config"]!);
        string a = rec["a"]!.GetValue<string>(), b = rec["b"]!.GetValue<string>();
        char strand = rec["strand"]!.GetValue<string>()[0];
        Assert.That(aligner.Algorithm, Is.EqualTo(rec["algorithm"]!.GetValue<string>()));
        double score = aligner.Score(a, b, strand);
        // aligner.score() and the score of the alignments agree, except for strand '-' with different left and right end gap scores (a defect of Biopython 1.88: its score() disagrees with the alignments it returns); the alignments decide
        double reported = rec["score"]!.GetValue<double>();
        double alignmentScore = rec["alignments"]!.AsArray().Count > 0 ? rec["alignments"]![0]!["score"]!.GetValue<double>() : reported;
        if (Math.Abs(alignmentScore - reported) > 1e-9)
        {
            Assert.That(strand, Is.EqualTo('-'), "only strand '-' shows the discrepancy");
            reported = alignmentScore;
        }

        Assert.That(score, Is.EqualTo(reported).Within(1e-9), "score");
        var alignments = aligner.Align(a, b, strand);
        var expected = rec["alignments"]!.AsArray().Select(x => Key(x!)).ToHashSet();
        var actual = alignments.Select(x => Key(x.Coordinates)).ToList();
        string kind = rec["config"]!["kind"]!.GetValue<string>();
        bool local = rec["config"]!["mode"]!.GetValue<string>() == "local";
        if (exactSet && (local || kind != "end_scores"))
        {
            Assert.That(actual.ToHashSet(), Is.EquivalentTo(expected), "the set of optimal alignments");
            Assert.That(actual, Has.Count.EqualTo(rec["count"]!.GetValue<int>()), "number of alignments");
        }
        else
        {
            Assert.That(actual, Is.Not.Empty);
            Assert.That(expected, Does.Contain(actual[0]), "one optimal alignment of the set");
        }
    }

    [Test]
    public void Oracle_HasTheAnnouncedInventory()
    {
        Assert.That(Records.Value, Has.Count.EqualTo(698));
        Assert.That(Records.Value.Count(r => r["config"]!["mode"]!.GetValue<string>() == "fogsaa"), Is.GreaterThan(200));
        Assert.That(Records.Value.Count(r => r["config"]!["wildcard"] is not null), Is.GreaterThan(20));
        Assert.That(Records.Value.Count(r => r["strand"]!.GetValue<string>() == "-"), Is.GreaterThan(100));
        // Biopython 1.88: score() differs from the score of its own alignments for strand '-' with asymmetric end gap scores
        Assert.That(Records.Value.Count(r => r["alignments"]!.AsArray().Count > 0 && Math.Abs(r["alignments"]![0]!["score"]!.GetValue<double>() - r["score"]!.GetValue<double>()) > 1e-9), Is.EqualTo(3));
    }

    [Test]
    public void Wildcard_ScoresZeroAgainstEverything_AndPresetsAreBiopythons()
    {
        var aligner = new PairwiseAligner { MatchScore = 2, MismatchScore = -1, Wildcard = 'N' };
        aligner.SetGapScore("gap_score", -2);
        Assert.That(aligner.Score("ACNGT", "ACAGT"), Is.EqualTo(8));
        Assert.That(new PairwiseAligner("blastn").Algorithm, Is.EqualTo("Gotoh global alignment algorithm"));
        Assert.That(new PairwiseAligner("megablast").Algorithm, Is.EqualTo("Needleman-Wunsch"));
        Assert.That(new PairwiseAligner("blastp").OpenInternalInsertionScore, Is.EqualTo(-12.0));
        Assert.Throws<ArgumentException>(() => new PairwiseAligner("nosuch"));
    }
}
