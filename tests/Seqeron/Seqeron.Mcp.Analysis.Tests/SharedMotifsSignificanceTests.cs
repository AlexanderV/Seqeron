using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// shared_motifs_significance — delegation to the RSAT overload
/// <c>MotifFinder.FindSharedMotifs(seqs, k, minSequences, OligoBackgroundModel, OligoStrandMode)</c>.
/// Expected values: RSAT oligo-analysis mseq outputs locked in MotifFinder_OligoAnalysis_Tests (ms.fa).
/// </summary>
[TestFixture]
public class SharedMotifsSignificanceTests
{
    private static readonly string[] Ms = { "ACGTACGTTAGC", "TTACGTAGCAAC", "GGTAGCACGTTT", "CATTTTACG" };
    private const double Rel = 1e-10;

    [Test]
    public void SharedMotifsSignificance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.SharedMotifsSignificance(Ms, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SharedMotifsSignificance(null!, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SharedMotifsSignificance(Array.Empty<string>(), 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SharedMotifsSignificance(new[] { "ACGT", "" }, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.SharedMotifsSignificance(Ms, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.SharedMotifsSignificance(Ms, 4, 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SharedMotifsSignificance(Ms, 4, background: "uniform"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SharedMotifsSignificance(Ms, 4, strands: "minus"));
    }

    [Test]
    public void SharedMotifsSignificance_InputBernoulli_SingleStrand_EqualsRsat()
    {
        // oligo-analysis -l 4 -1str -return occ,mseq,proba -lth mseq 2
        var r = AnalysisTools.SharedMotifsSignificance(Ms, 4, 2);
        var acgt = r.Motifs.Single(m => m.Sequence == "ACGT");
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceCount, Is.EqualTo(4));
            Assert.That(r.PossiblePositions, Is.EqualTo(33));
            Assert.That(r.PossibleOligos, Is.EqualTo(256));
            Assert.That(acgt.SequenceIndices, Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(acgt.ExpectedFrequency, Is.EqualTo(0.0037555250723974999).Within(Rel).Percent);
            Assert.That(acgt.ExpectedMatchingSequences, Is.EqualTo(0.12225827585979898).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceProbability, Is.EqualTo(0.00011159466135309564).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceEValue!.Value, Is.EqualTo(0.028568233306392483).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceSignificance!.Value, Is.EqualTo(1.5441166160753867).Within(1e-11));
        });
    }

    [Test]
    public void SharedMotifsSignificance_BothStrands_EqualsRsat_AndLibrary()
    {
        var r = AnalysisTools.SharedMotifsSignificance(Ms, 4, 2, strands: "both");
        var cgta = r.Motifs.Single(m => m.Sequence == "CGTA");
        Assert.Multiple(() =>
        {
            Assert.That(r.Motifs, Has.Length.EqualTo(7));
            Assert.That(r.PossibleOligos, Is.EqualTo(136));
            Assert.That(r.Strands, Is.EqualTo("both"));
            Assert.That(cgta.ReverseComplement, Is.EqualTo("TACG"));
            Assert.That(cgta.SequenceIndices, Is.EqualTo(new[] { 0, 1, 3 }));
            Assert.That(cgta.MatchingSequenceProbability, Is.EqualTo(0.00087319490302470882).Within(Rel).Percent);
            Assert.That(cgta.MatchingSequenceSignificance!.Value, Is.EqualTo(0.9253498996301337).Within(1e-11));
        });

        var lib = global::Seqeron.Genomics.Analysis.MotifFinder.FindSharedMotifs(
            Ms.Select(s => new global::Seqeron.Genomics.Core.DnaSequence(s)), 4, 2,
            global::Seqeron.Genomics.Analysis.OligoBackgroundModel.BernoulliFromInput,
            global::Seqeron.Genomics.Analysis.OligoStrandMode.Both);
        Assert.That(r.Motifs.Select(m => (m.Sequence, m.ReverseComplement, m.MatchingSequenceProbability)),
            Is.EqualTo(lib.Motifs.Select(m => (m.Sequence, m.ReverseComplement, m.MatchingSequenceProbability))));
    }
}
