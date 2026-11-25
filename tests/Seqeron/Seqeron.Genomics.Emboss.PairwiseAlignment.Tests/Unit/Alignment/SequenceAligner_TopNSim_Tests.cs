// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.PairwiseAlignment;
using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A2-28d, B09 F145): <see cref="TopNTieOrder.Sim"/> = EMBOSS 6.6.0 <c>matcher</c> (Huang-Miller SIM, <c>matcher.c</c> ported line by line):
/// the complete alignment lists (score, coordinates, both rows) of <c>matcher -alternatives K</c> in the program's own tie order, including the 7 rows of
/// <c>topn_oracle.tsv</c> whose score list differs from the Waterman-Eggert order. Oracle: <c>matcher_sim_oracle.tsv</c> (scripts/alignment/generate_matcher_sim_oracle.py).
/// </summary>
[TestFixture]
public class SequenceAligner_TopNSim_Tests
{
    private sealed record Row(int Line, string Kind, int GapOpen, int GapExtend, string Seq1, string Seq2, int K, string Alignments);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.matcher_sim_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine();
        int line = 1;
        string? text;
        while ((text = reader.ReadLine()) is not null)
        {
            line++;
            if (text.Length == 0 || text[0] == '#') continue;
            var f = text.Split('\t');
            rows.Add(new Row(line, f[0], int.Parse(f[2], CultureInfo.InvariantCulture), int.Parse(f[3], CultureInfo.InvariantCulture), f[4], f[5], int.Parse(f[6], CultureInfo.InvariantCulture), f[7]));
        }
        return rows;
    }

    private static string Key(AlignmentResult r) =>
        string.Create(CultureInfo.InvariantCulture, $"{r.Score}|{r.StartPosition1}|{r.StartPosition2}|{r.EndPosition1}|{r.EndPosition2}|{r.AlignedSequence1}|{r.AlignedSequence2}");

    /// <summary>matcher -gapopen G -gapextend E = GapOpen -(G - E), GapExtend -E.</summary>
    private static IReadOnlyList<AlignmentResult> Sim(Row r, TopNTieOrder order = TopNTieOrder.Sim)
    {
        if (r.Kind == "dna")
            return EmbossPairwisePrograms.LocalAlignAffineTopN(r.Seq1, r.Seq2, new ScoringMatrix(5, -4, -(r.GapOpen - r.GapExtend), -r.GapExtend), r.K, 0, order);
        var ms = new MatrixScoring(SubstitutionMatrix.Blosum62, -(r.GapOpen - r.GapExtend), -r.GapExtend);
        return EmbossPairwisePrograms.LocalAlignAffineTopN(r.Seq1, r.Seq2, ms, r.K, 0, order);
    }

    [Test]
    public void Sim_EqualsEmbossMatcher_OnEveryRow_AlignmentByAlignment()
    {
        int equal = 0, abort = 0;
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var actual = Sim(r);
            if (r.Alignments == "ABORT")
            {
                // the program dies ("The number of alignments computed is too large") exactly when fewer than K local alignments exist
                abort++;
                if (actual.Count >= r.K) failures.Add($"line {r.Line}: program aborts but the port returned {actual.Count} >= K={r.K} alignments");
                continue;
            }
            var expected = r.Alignments.Split(';');
            var got = actual.Select(Key).ToArray();
            if (expected.SequenceEqual(got)) equal++;
            else failures.Add($"line {r.Line} {r.Kind} G={r.GapOpen} E={r.GapExtend} {r.Seq1} {r.Seq2} K={r.K}\n   expected {string.Join(';', expected)}\n   actual   {string.Join(';', got)}");
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
        Assert.That(equal, Is.EqualTo(Oracle.Value.Count - abort));
        Assert.That(Oracle.Value.Count, Is.GreaterThanOrEqualTo(1000));
        Assert.That(abort, Is.GreaterThan(0));
    }

    /// <summary>
    /// The 7 rows of topn_oracle.tsv (lines 46, 487, 502, 787, 793, 865, 1348) whose matcher score list differs from Waterman-Eggert order are exact in SIM order,
    /// and Waterman-Eggert order still differs on them (the two orders are really different programs).
    /// </summary>
    private static readonly int[] MatcherTieRows = [46, 487, 502, 787, 793, 865, 1348];

    [Test]
    public void Sim_ReproducesTheSevenTieRows_WhereWatermanEggertOrderDiffers()
    {
        var lines = MatcherTieRows;
        var topn = ReadTopN();
        foreach (int line in lines)
        {
            var f = topn[line - 2];
            var row = Oracle.Value.Single(x => x.Seq1 == f[7] && x.Seq2 == f[8] && x.K == f[12].Split(',').Length && x.GapExtend == (int)-double.Parse(f[6], CultureInfo.InvariantCulture)
                                              && x.GapOpen == (int)-(double.Parse(f[5], CultureInfo.InvariantCulture) + double.Parse(f[6], CultureInfo.InvariantCulture)));
            var sim = Sim(row).Select(a => a.Score.ToString(CultureInfo.InvariantCulture)).ToArray();
            var emboss = f[12].Split(',');
            Assert.That(sim, Is.EqualTo(emboss), $"topn_oracle line {line}");
            var we = Sim(row, TopNTieOrder.WatermanEggert).Select(a => a.Score.ToString(CultureInfo.InvariantCulture)).ToArray();
            Assert.That(we, Is.Not.EqualTo(emboss), $"topn_oracle line {line} should differ in Waterman-Eggert order");
        }
    }

    private static List<string[]> ReadTopN()
    {
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream("Seqeron.Genomics.Tests.TestData.Alignment.topn_oracle.tsv")!;
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<string[]>();
        string? s;
        while ((s = reader.ReadLine()) is not null)
            list.Add(s.Split('\t'));
        return list;
    }

    [Test]
    public void Sim_EveryReportedList_IsALegalWatermanEggertSequence()
    {
        // each alignment's score is the optimum of the problem with the earlier aligned pairs removed: scores never increase and no pair is reused
        foreach (var r in Oracle.Value.Where(x => x.Alignments != "ABORT").Take(300))
        {
            var list = Sim(r);
            var used = new HashSet<(int, int)>();
            int previous = int.MaxValue;
            foreach (var a in list)
            {
                Assert.That(a.Score, Is.LessThanOrEqualTo(previous), $"line {r.Line}");
                previous = a.Score;
                foreach (var p in SequenceAligner.AlignedPairs(a))
                    Assert.That(used.Add(p), Is.True, $"line {r.Line}: pair {p} reused");
            }
        }
    }

    [Test]
    public void Sim_FirstAlignment_HasTheOptimalScore_AndEmptyInputsGiveNothing()
    {
        foreach (var r in Oracle.Value.Where(x => x.Alignments != "ABORT").Take(200))
        {
            var sim = Sim(r with { K = 1 });
            var we = Sim(r with { K = 1 }, TopNTieOrder.WatermanEggert);
            Assert.That(sim, Has.Count.EqualTo(1));
            Assert.That(sim[0].Score, Is.EqualTo(we[0].Score), $"line {r.Line}");
        }
        Assert.That(EmbossPairwisePrograms.LocalAlignAffineTopN("", "ACGT", null, 3, 0, TopNTieOrder.Sim), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.LocalAlignAffineTopN("ACGT", "ACGT", null, 0, 0, TopNTieOrder.Sim));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.LocalAlignAffineTopN("ACGT", "ACGT", null, 1, 0, (TopNTieOrder)9));
    }

    [Test]
    public void Sim_MinScore_StopsTheList()
    {
        var r = Oracle.Value.First(x => x.Alignments != "ABORT" && x.K >= 3 && x.Alignments.Split(';').Length >= 3);
        var all = Sim(r);
        int cut = all[1].Score;
        var some = r.Kind == "dna"
            ? EmbossPairwisePrograms.LocalAlignAffineTopN(r.Seq1, r.Seq2, new ScoringMatrix(5, -4, -(r.GapOpen - r.GapExtend), -r.GapExtend), r.K, cut, TopNTieOrder.Sim)
            : EmbossPairwisePrograms.LocalAlignAffineTopN(r.Seq1, r.Seq2, new MatrixScoring(SubstitutionMatrix.Blosum62, -(r.GapOpen - r.GapExtend), -r.GapExtend), r.K, cut, TopNTieOrder.Sim);
        Assert.That(some.Select(Key), Is.EqualTo(all.TakeWhile(a => a.Score >= cut).Select(Key)));
    }
}
