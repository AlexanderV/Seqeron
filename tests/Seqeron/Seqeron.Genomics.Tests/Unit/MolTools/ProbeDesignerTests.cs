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
    public void CheckSpecificity_BasicFunctionality_Smoke()
    {
        // Smoke test: Verify CheckSpecificity works with suffix tree
        // Detailed tests in ProbeDesigner_ProbeValidation_Tests.cs (PROBE-VALID-001)
        string genome = "ACGTACGTACGTACGTACGTACGTACGT";
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        double specificity = ProbeDesigner.CheckSpecificity("ACGTACGT", genomeIndex);

        Assert.That(specificity, Is.InRange(0.0, 1.0));
    }

    [Test]
    public void CalculateMolecularWeight_20mer_ReasonableWeight()
    {
        string oligo = "ACGTACGTACGTACGTACGT";

        double mw = ProbeDesigner.CalculateMolecularWeight(oligo);

        // 20-mer should be around 6000-7000 Da
        Assert.That(mw, Is.InRange(5500, 7500));
    }

    [Test]
    public void CalculateExtinctionCoefficient_ReturnsPositive()
    {
        double extinction = ProbeDesigner.CalculateExtinctionCoefficient("ACGTACGT");

        Assert.That(extinction, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateConcentration_FromAbsorbance()
    {
        double extinction = 200000;
        double absorbance = 0.5;

        double concentration = ProbeDesigner.CalculateConcentration(absorbance, extinction);

        Assert.That(concentration, Is.GreaterThan(0));
        Assert.That(concentration, Is.EqualTo(2.5).Within(0.1)); // µM
    }

}
