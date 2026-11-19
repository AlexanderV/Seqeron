using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// multiple_alignment_statistics / alignment_distance_matrix — delegation to <c>MultipleAlignmentStatistics</c> (B09 F43).
/// Expected values: EMBOSS 6.6.0 infoalign (-refseq 1) / distmat (uncorrected, x100 two decimals) and the ClustalW 2.1 conservation line.
/// </summary>
[TestFixture]
public class MultipleAlignmentStatisticsToolTests
{
    private static readonly string[] Dna = { "ACGT-ACGTAC", "ACGTTACG-AC", "AGGT-ACGTAC", "A-GTTTCGTAA" };
    private static readonly string[] Protein = { "MKVSLAKAG-HST", "MKVALAWAGQHSA", "MRVGLADAGNHNV", "MKVCLAEAG-HTA" };

    [Test]
    public void MultipleAlignmentStatistics_Dna_InfoalignAndDistmatValues()
    {
        var r = AnalysisTools.MultipleAlignmentStatistics(Dna, "nucleotide", referenceSequence: 1);
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceCount, Is.EqualTo(4));
            Assert.That(r.AlignmentLength, Is.EqualTo(11));
            // infoalign -refseq 1: s2  10 11 1 1 9 0 1 18.181818
            var s2 = r.Sequences[1];
            Assert.That((s2.SequenceLength, s2.AlignmentLength, s2.GapRuns, s2.GapLength, s2.Identities, s2.Similarities, s2.Differences), Is.EqualTo((10, 11, 1, 1, 9, 0, 1)));
            Assert.That(s2.PercentChange, Is.EqualTo(18.181818).Within(1e-5));
            Assert.That(r.PairwiseIdentity[0][2], Is.EqualTo(90.0).Within(1e-9)); // distmat 10.00
            Assert.That(r.Columns[0].Conservation, Is.EqualTo("*"));
            Assert.That(r.Columns[1].Counts["C"], Is.EqualTo(2));
            Assert.That(r.GapFreeColumns, Is.EqualTo(8));
        });
    }

    [Test]
    public void MultipleAlignmentStatistics_Protein_ClustalwConservationLine_AndConsensusReference()
    {
        var r = AnalysisTools.MultipleAlignmentStatistics(Protein, "protein");
        Assert.Multiple(() =>
        {
            Assert.That(r.ConservationLine.TrimEnd(), Is.EqualTo("*:* ** ** *.."));
            Assert.That(r.FullyConservedColumns, Is.EqualTo(7));
            Assert.That(r.Sequences, Has.Length.EqualTo(4));
            Assert.That(r.Sequences.All(s => s.Identities + s.Similarities + s.Differences == s.SequenceLength), Is.True);
        });
    }

    [Test]
    public void AlignmentDistanceMatrix_Uncorrected_And_JukesCantor()
    {
        var u = AnalysisTools.AlignmentDistanceMatrix(Dna);
        var jc = AnalysisTools.AlignmentDistanceMatrix(Dna, "jukes-cantor", "nucleotide");
        Assert.Multiple(() =>
        {
            Assert.That(u.Method, Is.EqualTo("uncorrected"));
            Assert.That(u.Matrix[0][2]!.Value * 100, Is.EqualTo(10.0).Within(0.005), "distmat prints 10.00 (single precision)");
            Assert.That(u.Matrix[1][2]!.Value * 100, Is.EqualTo(11.11).Within(0.005));
            Assert.That(jc.Matrix[0][2]!.Value, Is.EqualTo(-0.75 * Math.Log(1 - 0.1 / 0.75)).Within(1e-5), "distmat computes in single precision");
            Assert.That(u.Matrix[3][3], Is.EqualTo(0.0));
        });
    }

    [Test]
    public void AlignmentDistanceMatrix_AllDistmatCorrections_PrintedValuesOfEmboss6_6_0()
    {
        // distmat 6.6.0 -nucmethod 1..5 on the same alignment, matrix[0][2] and matrix[0][3] x100 (two decimals)
        var expected = new (string Method, double D02, double D03)[]
        {
            ("jukes-cantor", 10.73, 26.35), ("kimura", 10.85, 27.26), ("tamura", 10.85, 27.34), ("tajima-nei", 11.36, 29.26), ("jin-nei", 11.81, 34.29),
        };
        foreach (var (method, d02, d03) in expected)
        {
            var r = AnalysisTools.AlignmentDistanceMatrix(Dna, method);
            Assert.That(r.Method, Is.EqualTo(method));
            Assert.That(r.Matrix[0][2]!.Value * 100, Is.EqualTo(d02).Within(0.0051), method);
            Assert.That(r.Matrix[0][3]!.Value * 100, Is.EqualTo(d03).Within(0.0051), method);
        }

        Assert.That(AnalysisTools.AlignmentDistanceMatrix(Dna, "tajima-nei").Matrix[0][1], Is.Null, "distmat prints nan");
        Assert.Throws<NotSupportedException>(() => AnalysisTools.AlignmentDistanceMatrix(Dna, "tajima-nei", position: 12));
    }

    [Test]
    public void Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(Dna, "rna"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MultipleAlignmentStatistics(Dna, referenceSequence: 5));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(new[] { "ACG", "AC" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlignmentDistanceMatrix(Dna, "nope"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlignmentDistanceMatrix(Dna, "tamura", "protein"), "nucleotide-only method on protein");
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlignmentDistanceMatrix(new[] { "ACGT" }));
    }
}
