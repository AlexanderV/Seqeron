namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for CRISPR design algorithms.
///
/// Test Units: CRISPR-PAM-001, CRISPR-GUIDE-001, CRISPR-OFF-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("MolTools")]
public class CrisprSnapshotTests
{
    private const string TargetSequence =
        "ACGTACGTACGTACGTACGTGGGTTTTTTTTTTTTTTTTTTTTGGATCGATCGATCGATCG";

    [Test]
    public Task FindPamSites_SpCas9_Snapshot()
    {
        var dna = new DnaSequence(TargetSequence);
        var sites = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9)
            .Select(s => new { s.Position, s.PamSequence, s.TargetSequence, s.IsForwardStrand })
            .ToList();

        return Verify(new { Sites = sites });
    }

    [Test]
    public Task DesignGuideRnas_Snapshot()
    {
        var longSeq = new DnaSequence(string.Concat(Enumerable.Repeat("ACGTACGTACGTACGTACGT", 5)) + "GGG");
        var guides = CrisprDesigner.DesignGuideRnas(longSeq, 0, longSeq.Length - 1)
            .Take(3)
            .Select(g => new { g.Sequence, g.Position, g.GcContent, g.Score, g.HasPolyT })
            .ToList();

        return Verify(new { Guides = guides });
    }

    [Test]
    public Task EvaluateGuideRna_Snapshot()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var evaluation = CrisprDesigner.EvaluateGuideRna(guide);

        return Verify(new
        {
            evaluation.Sequence,
            evaluation.GcContent,
            evaluation.SeedGcContent,
            evaluation.HasPolyT,
            evaluation.SelfComplementarityScore,
            evaluation.Score
        });
    }
}
