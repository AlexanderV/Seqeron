// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
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
    public Task DesignPrimers_Snapshot()
    {
        string template = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 20));
        var dna = new DnaSequence(template);
        var result = ThermodynamicPrimerDesigner.DesignPrimers(dna, targetStart: 80, targetEnd: 160);

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
        var probes = ThermodynamicProbeDesigner.DesignProbes(target, maxProbes: 3)
            .Select(p => new { p.Sequence, p.Tm, p.GcContent, p.Score, p.Start, p.End })
            .ToList();

        return Verify(new { Probes = probes });
    }

    [Test]
    public Task ValidateProbe_Snapshot()
    {
        string probe = "ACGTACGTACGTACGTACGTACGT";
        var refs = new[] { "ACGTACGTACGTACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTTTTTT" };
        var validation = ThermodynamicProbeDesigner.ValidateProbe(probe, refs);

        return Verify(new
        {
            validation.IsValid,
            validation.SpecificityScore,
            validation.OffTargetHits,
            validation.SelfComplementarity,
            validation.HasSecondaryStructure,
            // Primer3 ntthal values at the probe conditions (primer3-py 2.3.1 calc_homodimer / calc_end_stability /
            // calc_hairpin, mv 50, dv 0, dntp 0, dna 50): 59.1857717189107 / 59.1857717189107 / 67.29188756961071.
            validation.ThermodynamicScreen,
            SelfDimerTm = Math.Round(validation.SelfDimerTm!.Value, 6),
            SelfEndDimerTm = Math.Round(validation.SelfEndDimerTm!.Value, 6),
            HairpinTm = Math.Round(validation.HairpinTm!.Value, 6),
            validation.Issues
        });
    }
}
