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
            Assert.That(u.Matrix[0][2]!.Value * 100, Is.EqualTo(10.0).Within(1e-9));
            Assert.That(u.Matrix[1][2]!.Value * 100, Is.EqualTo(11.11).Within(0.005));
            Assert.That(jc.Matrix[0][2]!.Value, Is.EqualTo(-0.75 * Math.Log(1 - 0.1 / 0.75)).Within(1e-9));
            Assert.That(u.Matrix[3][3], Is.EqualTo(0.0));
        });
    }

    [Test]
    public void Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(Dna, "rna"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MultipleAlignmentStatistics(Dna, referenceSequence: 5));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultipleAlignmentStatistics(new[] { "ACG", "AC" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlignmentDistanceMatrix(Dna, "kimura"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlignmentDistanceMatrix(new[] { "ACGT" }));
    }
}
