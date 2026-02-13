using VerifyNUnit;

namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot/approval tests for complex algorithm outputs.
/// Instead of verifying individual fields, the entire result is serialized and compared
/// against a committed .verified.txt golden master file.
///
/// Any change to the output requires explicit developer approval (updating the snapshot).
/// This catches unintentional behavioral regressions in complex outputs.
///
/// First run: creates .verified.txt files that must be committed to git.
/// Subsequent runs: compares output against .verified.txt and fails on diff.
/// </summary>
[TestFixture]
[Category("Snapshot")]
public class AlignmentSnapshotTests
{
    [Test]
    public Task GlobalAlign_KnownPair_MatchesSnapshot()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGT");

        return Verifier.Verify(new
        {
            result.AlignedSequence1,
            result.AlignedSequence2,
            result.Score,
            result.AlignmentType
        });
    }

    [Test]
    public Task LocalAlign_OverlappingMotif_MatchesSnapshot()
    {
        var result = SequenceAligner.LocalAlign("AAAAACGTACGTAAAAA", "CGTACGT");

        return Verifier.Verify(new
        {
            result.AlignedSequence1,
            result.AlignedSequence2,
            result.Score,
            result.StartPosition1,
            result.StartPosition2
        });
    }
}

/// <summary>
/// Snapshot tests for disorder prediction outputs.
/// </summary>
[TestFixture]
[Category("Snapshot")]
public class DisorderSnapshotTests
{
    [Test]
    public Task PredictDisorder_KnownSequence_MatchesSnapshot()
    {
        // W×10 (ordered) + P×20 (disordered) + W×10 (ordered)
        string sequence = new string('W', 10) + new string('P', 20) + new string('W', 10);
        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        return Verifier.Verify(new
        {
            result.OverallDisorderContent,
            result.MeanDisorderScore,
            Regions = result.DisorderedRegions.Select(r => new
            {
                r.Start,
                r.End,
                r.RegionType,
                r.Confidence,
                r.MeanScore
            }).ToList()
        });
    }

    [Test]
    public Task PredictDisorder_AllOrdered_MatchesSnapshot()
    {
        string sequence = new string('W', 30);
        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        return Verifier.Verify(new
        {
            result.OverallDisorderContent,
            result.MeanDisorderScore,
            RegionCount = result.DisorderedRegions.Count
        });
    }
}
