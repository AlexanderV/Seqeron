using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// oligo_analysis <c>sequenceType</c> (RSAT -seqtype prot|other) → <c>MotifFinder.AnalyzeOligoStrings</c>;
/// shared_motifs_significance <c>degenerate</c> (RSAT -oneN / -onedeg with -return mseq,proba) →
/// <c>MotifFinder.FindSharedMotifs(…, pseudoFrequency, OligoDegeneracy)</c>.
/// Expected values: RSAT oligo-analysis oracle runs locked in MotifFinder_OligoSequenceType_Tests and
/// MotifFinder_SharedMotifsDegenerate_Tests.
/// </summary>
[TestFixture]
public class OligoSequenceTypeAndDegenerateSharedTests
{
    private const double Rel = 1e-10;

    [Test]
    public void OligoAnalysis_Protein_EqualsRsat()
    {
        // -seqtype prot -l 2 on MKLLVAAGLLKLMKXLLA* / mkllvqqKLLAA: n = 26, NPO 400, 14 tested.
        var r = AnalysisTools.OligoAnalysis("MKLLVAAGLLKLMKXLLA*", 2, 1, extraSequences: new[] { "mkllvqqKLLAA" },
            zscore: true, sequenceType: "protein");
        var ll = r.Motifs.Single(m => m.Sequence == "LL");
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceType, Is.EqualTo("protein"));
            Assert.That(r.AlphabetSize, Is.EqualTo(20));
            Assert.That(r.PossibleOligos, Is.EqualTo(400));
            Assert.That(r.TotalOccurrences, Is.EqualTo(26));
            Assert.That(r.TestedPatterns, Is.EqualTo(14));
            Assert.That(r.Strands, Is.EqualTo("single"));
            Assert.That(ll.OccurrenceProbability, Is.EqualTo(0.31619335865366638).Within(Rel).Percent);
            Assert.That(ll.ZScore, Is.EqualTo(0.63864701205606378).Within(Rel).Percent);
            Assert.That(ll.ReverseComplement, Is.Null);
        });
    }

    [Test]
    public void OligoAnalysis_Other_EqualsRsat()
    {
        // -seqtype other -bg equi -l 2: alphabet 13, NPO 169.
        var r = AnalysisTools.OligoAnalysis("Hello, World! hello world. AbC abc", 2, 1, "equiprobable", sequenceType: "other");
        Assert.Multiple(() =>
        {
            Assert.That(r.AlphabetSize, Is.EqualTo(13));
            Assert.That(r.PossibleOligos, Is.EqualTo(169));
            Assert.That(r.Motifs.Single(m => m.Sequence == "ll").OccurrenceProbability,
                Is.EqualTo(0.01194994022319412).Within(Rel).Percent);
        });
    }

    [Test]
    public void OligoAnalysis_DnaDefault_Unchanged()
    {
        var a = AnalysisTools.OligoAnalysis("ATGCATGCATGCAAATTTGGG", 3, 2);
        var b = AnalysisTools.OligoAnalysis("ATGCATGCATGCAAATTTGGG", 3, 2, sequenceType: "dna");
        Assert.Multiple(() =>
        {
            Assert.That(b.Motifs.Select(m => (m.Sequence, m.OccurrenceProbability)), Is.EqualTo(a.Motifs.Select(m => (m.Sequence, m.OccurrenceProbability))));
            Assert.That(b.SequenceType, Is.Null);
        });
    }

    [Test]
    public void OligoAnalysis_ResidueAlphabet_InvalidOptions_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("MKLLV", 2, 1, sequenceType: "rna"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("MKLLV", 2, 1, strands: "both", sequenceType: "protein"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("MKLLV", 2, 1, degenerate: "oneN", sequenceType: "protein"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("MKLLV", 2, 1, "bernoulli", residueFrequencies: new[] { 0.25, 0.25, 0.25, 0.25 }, sequenceType: "protein"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("MKLLV", 2, 1, calibrationTable: "ac\t1\t1\t1\n", sequenceType: "other"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("  ", 2, 1, sequenceType: "other"));
        });
    }

    [Test]
    public void SharedMotifsSignificance_Degenerate_EqualsRsatOracle()
    {
        string[] d3 = { "ACGTACGGATCC", "ATGCATGAAC", "ACGATGTT" };
        // oracle -l 3 -1str -onedeg -lth mseq 1: AYG mseq 3, ms_P 0.01270294085234148, NPO 528.
        var one = AnalysisTools.SharedMotifsSignificance(d3, 3, 1, degenerate: "onedeg");
        // oracle -l 3 -2str -oneN -lth mseq 2: ANG|CNT ms_P 0.28055070694211726, NPO 24.
        var both = AnalysisTools.SharedMotifsSignificance(d3, 3, 2, strands: "both", degenerate: "oneN");
        Assert.Multiple(() =>
        {
            Assert.That(one.Degenerate, Is.EqualTo("onedeg"));
            Assert.That(one.PossibleOligos, Is.EqualTo(528));
            var ayg = one.Motifs.Single(m => m.Sequence == "AYG");
            Assert.That(ayg.SequenceIndices, Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(ayg.MatchingSequenceProbability, Is.EqualTo(0.01270294085234148).Within(Rel).Percent);
            Assert.That(both.PossibleOligos, Is.EqualTo(24));
            var ang = both.Motifs.Single(m => m.Sequence == "ANG");
            Assert.That(ang.ReverseComplement, Is.EqualTo("CNT"));
            Assert.That(ang.MatchingSequenceProbability, Is.EqualTo(0.28055070694211726).Within(Rel).Percent);
            Assert.That(AnalysisTools.SharedMotifsSignificance(d3, 3, 2).Degenerate, Is.Null);
            Assert.Throws<ArgumentException>(() => AnalysisTools.SharedMotifsSignificance(d3, 3, 2, degenerate: "two"));
        });
    }
}
