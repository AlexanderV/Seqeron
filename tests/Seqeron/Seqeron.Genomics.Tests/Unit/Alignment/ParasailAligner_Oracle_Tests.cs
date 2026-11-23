using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit R14-11, B09 F304): <see cref="ParasailAligner"/> against parasail 1.3.4 on 520 recorded pairs x 6 result kinds (plain, _stats, _table, _rowcol, _stats_table, _stats_rowcol) over
/// <c>nw</c>, <c>sw</c> and the eleven semi-global variants (<c>sg</c>, <c>sg_qb</c>, <c>sg_qe</c>, <c>sg_db</c>, <c>sg_de</c>, <c>sg_qx</c>, <c>sg_dx</c>, <c>sg_qb_db</c>, <c>sg_qb_de</c>, <c>sg_qe_db</c>, <c>sg_qe_de</c>): score, <c>end_query</c>, <c>end_ref</c>,
/// matches / similar / length, the tables, the last row and column — equal on every value. The program runs with O(min(m, n)) cells unless a table is asked for; <c>scripts/alignment/generate_parasail_oracle.py</c>.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class ParasailAligner_Oracle_Tests
{
    private static readonly Lazy<List<JsonObject>> Records = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.parasail_oracle.jsonl";
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

    private static IEnumerable<TestCaseData> Cases() => Records.Value.Select(r => new TestCaseData(r.ToJsonString()).SetName($"{r["id"]}_{r["fn"]}"));

    private static ParasailMatrix MatrixOf(string name)
    {
        if (name == "blosum62")
            return ParasailMatrix.From(SubstitutionMatrix.Get("BLOSUM62"));
        var p = name.Split(':');
        return ParasailMatrix.Create("ACGT", int.Parse(p[1]), int.Parse(p[2]));
    }

    private static (ParasailMode Mode, ParasailEndGaps Ends) ModeOf(string fn)
    {
        if (fn == "nw")
            return (ParasailMode.Global, default);
        if (fn == "sw")
            return (ParasailMode.Local, default);
        Assert.That(ParasailEndGaps.TryParse(fn, out var flags), Is.True, fn);
        return (ParasailMode.SemiGlobal, flags);
    }

    private static int[][] Rows2(JsonNode n) => n.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<int>()).ToArray()).ToArray();

    private static int[] Row1(JsonNode n) => n.AsArray().Select(x => x!.GetValue<int>()).ToArray();

    private static void AssertTable(int[,]? actual, JsonNode expected, string label)
    {
        var e = Rows2(expected);
        Assert.That(actual, Is.Not.Null, label);
        Assert.That(actual!.GetLength(0), Is.EqualTo(e.Length), label);
        for (int i = 0; i < e.Length; i++)
            for (int j = 0; j < e[i].Length; j++)
                Assert.That(actual[i, j], Is.EqualTo(e[i][j]), $"{label}[{i},{j}]");
    }

    [TestCaseSource(nameof(Cases))]
    public void AllResultKinds_EqualParasail(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        string s1 = rec["s1"]!.GetValue<string>(), s2 = rec["s2"]!.GetValue<string>();
        int open = rec["open"]!.GetValue<int>(), gap = rec["gap"]!.GetValue<int>();
        var matrix = MatrixOf(rec["matrix"]!.GetValue<string>());
        var (mode, ends) = ModeOf(rec["fn"]!.GetValue<string>());
        foreach (var (variant, expected) in rec["results"]!.AsObject())
        {
            var output = new ParasailOutput(variant.Contains("stats"), variant.Contains("table"), variant.Contains("rowcol"));
            var r = ParasailAligner.Run(s1, s2, open, gap, matrix, mode, ends, output);
            string label = $"{rec["id"]} {rec["fn"]} {variant}";
            Assert.That(r.Score, Is.EqualTo(expected!["score"]!.GetValue<int>()), label + " score");
            Assert.That(r.EndQuery, Is.EqualTo(expected["end_query"]!.GetValue<int>()), label + " end_query");
            Assert.That(r.EndReference, Is.EqualTo(expected["end_ref"]!.GetValue<int>()), label + " end_ref");
            if (output.Statistics)
            {
                Assert.That(r.Matches, Is.EqualTo(expected["matches"]!.GetValue<int>()), label + " matches");
                Assert.That(r.Similar, Is.EqualTo(expected["similar"]!.GetValue<int>()), label + " similar");
                Assert.That(r.Length, Is.EqualTo(expected["length"]!.GetValue<int>()), label + " length");
            }
            else
            {
                Assert.That(r.Matches, Is.Null);
            }

            if (output.Table)
            {
                AssertTable(r.Tables!.Score, expected["score_table"]!, label + " score_table");
                if (output.Statistics)
                {
                    AssertTable(r.Tables.Matches, expected["matches_table"]!, label + " matches_table");
                    AssertTable(r.Tables.Similar, expected["similar_table"]!, label + " similar_table");
                    AssertTable(r.Tables.Length, expected["length_table"]!, label + " length_table");
                }
            }

            if (output.RowCol)
            {
                Assert.That(r.RowCol!.ScoreRow, Is.EqualTo(Row1(expected["score_row"]!)), label + " score_row");
                Assert.That(r.RowCol.ScoreColumn, Is.EqualTo(Row1(expected["score_col"]!)), label + " score_col");
                if (output.Statistics)
                {
                    Assert.That(r.RowCol.MatchesRow, Is.EqualTo(Row1(expected["matches_row"]!)), label + " matches_row");
                    Assert.That(r.RowCol.MatchesColumn, Is.EqualTo(Row1(expected["matches_col"]!)), label + " matches_col");
                    Assert.That(r.RowCol.SimilarRow, Is.EqualTo(Row1(expected["similar_row"]!)), label + " similar_row");
                    Assert.That(r.RowCol.SimilarColumn, Is.EqualTo(Row1(expected["similar_col"]!)), label + " similar_col");
                    Assert.That(r.RowCol.LengthRow, Is.EqualTo(Row1(expected["length_row"]!)), label + " length_row");
                    Assert.That(r.RowCol.LengthColumn, Is.EqualTo(Row1(expected["length_col"]!)), label + " length_col");
                }
            }
        }
    }

    [Test]
    public void Oracle_HasTheAnnouncedInventory()
    {
        Assert.That(Records.Value, Has.Count.EqualTo(520));
        Assert.That(Records.Value.Select(r => r["fn"]!.GetValue<string>()).Distinct().Count(), Is.EqualTo(13));
        Assert.That(Records.Value.All(r => r["results"]!.AsObject().Count == 6), Is.True);
        Assert.That(Records.Value.Count(r => r["matrix"]!.GetValue<string>() == "blosum62"), Is.GreaterThan(150));
    }

    [Test]
    public void Profile_ScoresManyReferences_LikeTheDirectCalls_AndMemoryIsInTheShorterSequence()
    {
        var matrix = ParasailMatrix.Create("ACGT", 2, -3);
        var profile = new ParasailProfile("ACGTACGTTAGC", matrix);
        string[] references = ["ACGTACGTTAGC", "TTTTACGTACGTAAAA", "GGGG", "ACGTAC"];
        var batch = profile.ScoreBatch(references, 5, 2, ParasailMode.SemiGlobal, ParasailEndGaps.All, new ParasailOutput(Statistics: true));
        var parallel = profile.ScoreBatch(references, 5, 2, ParasailMode.SemiGlobal, ParasailEndGaps.All, new ParasailOutput(Statistics: true), parallel: true);
        for (int k = 0; k < references.Length; k++)
        {
            var direct = ParasailAligner.Sg("sg", "ACGTACGTTAGC", references[k], 5, 2, matrix, new ParasailOutput(Statistics: true));
            Assert.That(batch[k], Is.EqualTo(direct));
            Assert.That(parallel[k], Is.EqualTo(batch[k]));
        }

        // a long reference against a short query: the outer loop runs over the long one, the arrays hold the short one
        string longReference = string.Concat(Enumerable.Repeat("ACGTTGCA", 20000));
        var big = ParasailAligner.Nw("ACGTTGCAAC", longReference, 5, 2, matrix);
        Assert.That(big.EndReference, Is.EqualTo(longReference.Length - 1));
        Assert.Throws<ArgumentException>(() => ParasailAligner.Nw("", "ACGT", 5, 2, matrix));
        Assert.Throws<ArgumentException>(() => ParasailAligner.Nw("ACGT", "ACGT", -1, 2, matrix));
        Assert.Throws<ArgumentException>(() => ParasailAligner.Nw("ACGU", "ACGT", 5, 2, matrix));
        Assert.Throws<ArgumentException>(() => ParasailAligner.Sg("sg_zz", "A", "A", 5, 2, matrix));
    }
}
