namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for splice site prediction.
///
/// Test Units: SPLICE-DONOR-001, SPLICE-ACCEPTOR-001, SPLICE-PREDICT-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Annotation")]
public class SplicingSnapshotTests
{
    // Exon–intron–exon with canonical GT...AG
    private const string TestSequence =
        "ATGATGAAAGCCGCCATGGCG" +
        "GTAAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCAG" +
        "GCGATGAAAGCCGCCATGGCG";

    [Test]
    public Task FindDonorSites_Snapshot()
    {
        var donors = SpliceSitePredictor.FindDonorSites(TestSequence, minScore: 0.1)
            .Select(d => new { d.Position, d.Motif, d.Score, d.Confidence })
            .ToList();

        return Verify(new { Donors = donors });
    }

    [Test]
    public Task FindAcceptorSites_Snapshot()
    {
        var acceptors = SpliceSitePredictor.FindAcceptorSites(TestSequence, minScore: 0.1)
            .Select(a => new { a.Position, a.Motif, a.Score, a.Confidence })
            .ToList();

        return Verify(new { Acceptors = acceptors });
    }

    [Test]
    public Task PredictGeneStructure_Snapshot()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(TestSequence,
            minExonLength: 10, minIntronLength: 30);

        return Verify(new
        {
            ExonCount = structure.Exons.Count,
            IntronCount = structure.Introns.Count,
            SplicedSequenceLength = structure.SplicedSequence.Length,
            structure.OverallScore,
            Exons = structure.Exons.Select(e => new { e.Start, e.End, e.Length, e.Type }).ToList()
        });
    }
}
