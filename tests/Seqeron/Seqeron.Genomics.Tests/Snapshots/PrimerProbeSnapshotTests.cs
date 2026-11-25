namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for primer and probe design.
///
/// Test Units: PRIMER-TM-001, PRIMER-DESIGN-001, PROBE-DESIGN-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("MolTools")]
public class PrimerProbeSnapshotTests
{
    [Test]
    public Task MeltingTemperature_Snapshot()
    {
        var primers = new[] { "ATGCATGCATGC", "GCGCGCGCGCGC", "ATATATATATATAT", "ACGTACGTACGTACGT" };
        var results = primers.Select(p => new
        {
            Primer = p,
            Tm = PrimerDesigner.CalculateMeltingTemperature(p),
            TmWithSalt = PrimerDesigner.CalculateMeltingTemperatureWithSalt(p),
            GcContent = PrimerDesigner.CalculateGcContent(p)
        }).ToList();

        return Verify(new { Results = results });
    }
}
