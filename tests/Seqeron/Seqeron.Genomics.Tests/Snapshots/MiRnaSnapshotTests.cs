namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for miRNA analysis.
///
/// Test Units: MIRNA-SEED-001, MIRNA-TARGET-001, MIRNA-PRECURSOR-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Annotation")]
public class MiRnaSnapshotTests
{
    [Test]
    public Task CreateMiRna_Snapshot()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("hsa-miR-21-5p", "UAGCUUAUCAGACUGAUGUUGA");

        return Verify(new
        {
            mirna.Name,
            mirna.Sequence,
            mirna.SeedSequence,
            mirna.SeedStart,
            mirna.SeedEnd
        });
    }

    [Test]
    public Task FindTargetSites_Snapshot()
    {
        string mRna = "AUGCCAUUUUAGCUUAUCAGACAACUAUGAAUCCAAUUAGCUUAUCAGACAACUAUUU";
        var mirna = MiRnaAnalyzer.CreateMiRna("miR-test", "UAGCUUAUCAGACUGAUGUUGA");
        var sites = MiRnaAnalyzer.FindTargetSites(mRna, mirna, minScore: 0.1)
            .Select(s => new { s.Start, s.End, s.TargetSequence, s.Type, s.Score })
            .ToList();

        return Verify(new { Sites = sites });
    }

    [Test]
    public Task AlignMiRnaToTarget_Snapshot()
    {
        var duplex = MiRnaAnalyzer.AlignMiRnaToTarget(
            "UAGCUUAUCAGACUGAUGUUGA",
            "UCAACAUCAGUCUGAUAAGCUA");

        return Verify(new
        {
            duplex.MiRnaSequence,
            duplex.TargetSequence,
            duplex.AlignmentString,
            duplex.Matches,
            duplex.Mismatches,
            duplex.GUWobbles
        });
    }
}
