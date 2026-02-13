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

    [Test]
    public Task DesignPrimers_Snapshot()
    {
        string template = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 20));
        var dna = new DnaSequence(template);
        var result = PrimerDesigner.DesignPrimers(dna, targetStart: 80, targetEnd: 160);

        return Verify(new
        {
            result.IsValid,
            result.ProductSize,
            result.Message,
            Forward = result.Forward != null ? new { result.Forward.Sequence, result.Forward.GcContent, result.Forward.MeltingTemperature, result.Forward.Score } : null,
            Reverse = result.Reverse != null ? new { result.Reverse.Sequence, result.Reverse.GcContent, result.Reverse.MeltingTemperature, result.Reverse.Score } : null
        });
    }

    [Test]
    public Task DesignProbes_Snapshot()
    {
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 10));
        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 3)
            .Select(p => new { p.Sequence, p.Tm, p.GcContent, p.Score, p.Start, p.End })
            .ToList();

        return Verify(new { Probes = probes });
    }

    [Test]
    public Task ValidateProbe_Snapshot()
    {
        string probe = "ACGTACGTACGTACGTACGTACGT";
        var refs = new[] { "ACGTACGTACGTACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTTTTTT" };
        var validation = ProbeDesigner.ValidateProbe(probe, refs);

        return Verify(new
        {
            validation.IsValid,
            validation.SpecificityScore,
            validation.OffTargetHits,
            validation.SelfComplementarity,
            validation.HasSecondaryStructure
        });
    }
}
