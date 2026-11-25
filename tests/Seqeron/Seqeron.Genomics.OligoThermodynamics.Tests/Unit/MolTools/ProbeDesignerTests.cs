// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Supplementary tests for ProbeDesigner utilities (oligo analysis, MW, extinction, concentration)
/// and smoke tests for validation (detailed in PROBE-VALID-001).
/// Probe design tests are in ProbeDesigner_ProbeDesign_Tests.cs (PROBE-DESIGN-001).
/// </summary>
[TestFixture]
public class ProbeDesignerTests
{

    [Test]
    [Category("Smoke")]
    public void ValidateProbe_BasicFunctionality_Smoke()
    {
        // Smoke test: Verify ValidateProbe returns valid result
        // Detailed tests in ProbeDesigner_ProbeValidation_Tests.cs (PROBE-VALID-001)
        string probe = "ACGTACGTACGTACGTACGT";
        var references = new[] { "NNNNACGTACGTACGTACGTACGTNNNN" };

        var validation = ThermodynamicProbeDesigner.ValidateProbe(probe, references);

        Assert.Multiple(() =>
        {
            Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0));
            Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0));
            Assert.That(validation.Issues, Is.Not.Null);
        });
    }




    [Test]
    public void AnalyzeOligo_CalculatesAllProperties()
    {
        string oligo = "ACGTACGTACGTACGTACGT";

        var (tm, gc, mw, extinction) = ThermodynamicProbeDesigner.AnalyzeOligo(oligo);

        Assert.That(tm, Is.GreaterThan(40));
        Assert.That(gc, Is.EqualTo(0.5).Within(0.01));
        Assert.That(mw, Is.GreaterThan(5000));
        Assert.That(extinction, Is.GreaterThan(100000));
    }
}
