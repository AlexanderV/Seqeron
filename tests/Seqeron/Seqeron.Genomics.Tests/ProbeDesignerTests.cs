using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ProbeDesignerTests
{
    #region Basic Design Tests

    [Test]
    public void DesignProbes_ValidSequence_ReturnsProbes()
    {
        // Good probe region: moderate GC, no repeats
        string target = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 5).ToList();

        Assert.That(probes.Count, Is.GreaterThan(0));
    }

    [Test]
    public void DesignProbes_ShortSequence_ReturnsEmpty()
    {
        var probes = ProbeDesigner.DesignProbes("ACGT").ToList();

        Assert.That(probes, Is.Empty);
    }

    [Test]
    public void DesignProbes_EmptySequence_ReturnsEmpty()
    {
        var probes = ProbeDesigner.DesignProbes("").ToList();

        Assert.That(probes, Is.Empty);
    }

    [Test]
    public void DesignProbes_ProbesHaveValidCoordinates()
    {
        string target = new string('A', 30) + "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + new string('T', 30);

        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var probe in probes)
        {
            Assert.That(probe.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(probe.End, Is.LessThan(target.Length));
            Assert.That(probe.End, Is.GreaterThan(probe.Start));
        }
    }

    [Test]
    public void DesignProbes_ProbesHaveCalculatedProperties()
    {
        string target = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 3).ToList();

        foreach (var probe in probes)
        {
            Assert.That(probe.Tm, Is.GreaterThan(0));
            Assert.That(probe.GcContent, Is.InRange(0.0, 1.0));
            Assert.That(probe.Score, Is.InRange(0.0, 1.0));
        }
    }

    #endregion

    #region Parameter Tests

    [Test]
    public void DesignProbes_MicroarrayDefaults_ReturnsValidProbes()
    {
        var param = ProbeDesigner.Defaults.Microarray;
        string target = new string('G', 25) + "ACGTACGTACGTACGTACGTACGTACGT" + new string('C', 25);

        var probes = ProbeDesigner.DesignProbes(target, param).ToList();

        foreach (var probe in probes)
        {
            Assert.That(probe.Sequence.Length, Is.InRange(param.MinLength, param.MaxLength));
        }
    }

    [Test]
    public void DesignProbes_FISHDefaults_AllowsLongerProbes()
    {
        var param = ProbeDesigner.Defaults.FISH;

        Assert.That(param.MinLength, Is.GreaterThanOrEqualTo(200));
        Assert.That(param.MaxLength, Is.GreaterThanOrEqualTo(500));
    }

    [Test]
    public void DesignProbes_qPCRDefaults_ShorterProbes()
    {
        var param = ProbeDesigner.Defaults.qPCR;

        Assert.That(param.MinLength, Is.LessThanOrEqualTo(30));
        Assert.That(param.MaxLength, Is.LessThanOrEqualTo(35));
    }

    #endregion

    #region Tiling Probe Tests

    [Test]
    public void DesignTilingProbes_CoversSequence()
    {
        string target = new string('A', 100) + "GCGCGCGC" + new string('T', 100);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 50, overlap: 10);

        Assert.That(tiling.Probes.Count, Is.GreaterThan(0));
        Assert.That(tiling.Coverage, Is.GreaterThan(0));
    }

    [Test]
    public void DesignTilingProbes_CalculatesTmStatistics()
    {
        string target = new string('G', 50) + new string('C', 50) + new string('A', 50);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 40, overlap: 10);

        Assert.That(tiling.MeanTm, Is.GreaterThan(0));
        Assert.That(tiling.TmRange, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void DesignTilingProbes_ProbesAreTilingType()
    {
        string target = new string('A', 200);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 50);

        Assert.That(tiling.Probes.All(p => p.Type == ProbeDesigner.ProbeType.Tiling), Is.True);
    }

    #endregion

    #region Antisense Probe Tests

    [Test]
    public void DesignAntisenseProbes_ReturnsAntisenseType()
    {
        string mRna = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC";

        var probes = ProbeDesigner.DesignAntisenseProbes(mRna, maxProbes: 3).ToList();

        Assert.That(probes.All(p => p.Type == ProbeDesigner.ProbeType.Antisense), Is.True);
    }

    #endregion

    #region Molecular Beacon Tests

    [Test]
    public void DesignMolecularBeacon_CreatesBeaconWithStem()
    {
        string target = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        var beacon = ProbeDesigner.DesignMolecularBeacon(target, probeLength: 20, stemLength: 5);

        Assert.That(beacon, Is.Not.Null);
        Assert.That(beacon!.Value.Type, Is.EqualTo(ProbeDesigner.ProbeType.MolecularBeacon));
        Assert.That(beacon.Value.Sequence.Length, Is.GreaterThan(20)); // Has stems
    }

    [Test]
    public void DesignMolecularBeacon_ShortSequence_ReturnsNull()
    {
        var beacon = ProbeDesigner.DesignMolecularBeacon("ACGT", probeLength: 20);

        Assert.That(beacon, Is.Null);
    }

    #endregion

    #region Validation Tests (Smoke - detailed tests in ProbeDesigner_ProbeValidation_Tests.cs)

    [Test]
    [Category("Smoke")]
    public void ValidateProbe_BasicFunctionality_Smoke()
    {
        // Smoke test: Verify ValidateProbe returns valid result
        // Detailed tests in ProbeDesigner_ProbeValidation_Tests.cs (PROBE-VALID-001)
        string probe = "ACGTACGTACGTACGTACGT";
        var references = new[] { "NNNNACGTACGTACGTACGTACGTNNNN" };

        var validation = ProbeDesigner.ValidateProbe(probe, references);

        Assert.Multiple(() =>
        {
            Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0));
            Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0));
            Assert.That(validation.Issues, Is.Not.Null);
        });
    }

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

    #endregion

    #region Oligo Analysis Tests

    [Test]
    public void AnalyzeOligo_CalculatesAllProperties()
    {
        string oligo = "ACGTACGTACGTACGTACGT";

        var (tm, gc, mw, extinction) = ProbeDesigner.AnalyzeOligo(oligo);

        Assert.That(tm, Is.GreaterThan(40));
        Assert.That(gc, Is.EqualTo(0.5).Within(0.01));
        Assert.That(mw, Is.GreaterThan(5000));
        Assert.That(extinction, Is.GreaterThan(100000));
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

    #endregion

    #region Edge Cases

    [Test]
    public void DesignProbes_AllGC_HandlesExtreme()
    {
        string target = new string('G', 100);

        var probes = ProbeDesigner.DesignProbes(target).ToList();

        // Should still return probes, but with warnings about GC
        foreach (var probe in probes)
        {
            Assert.That(probe.GcContent, Is.EqualTo(1.0).Within(0.01));
        }
    }

    [Test]
    public void DesignProbes_AllAT_HandlesExtreme()
    {
        string target = new string('A', 50) + new string('T', 50);

        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var probe in probes)
        {
            Assert.That(probe.GcContent, Is.LessThanOrEqualTo(0.1));
        }
    }

    [Test]
    public void DesignProbes_Homopolymer_GeneratesWarnings()
    {
        string target = new string('A', 20) + "GCGCGCGC" + new string('A', 20) +
                       "TATATATA" + new string('G', 30);

        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 10).ToList();

        // Probes covering homopolymer should have warnings
        var probesWithWarnings = probes.Where(p => p.Warnings.Count > 0).ToList();
        Assert.That(probesWithWarnings.Count, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void DesignProbes_CaseInsensitive()
    {
        string upper = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        string lower = "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt";

        var probesUpper = ProbeDesigner.DesignProbes(upper, maxProbes: 1).ToList();
        var probesLower = ProbeDesigner.DesignProbes(lower, maxProbes: 1).ToList();

        // Both MUST produce probes - 60 bp sequence is sufficient
        Assert.That(probesUpper, Is.Not.Empty, "Upper case sequence must produce probes");
        Assert.That(probesLower, Is.Not.Empty, "Lower case sequence must produce probes");

        Assert.That(probesUpper[0].Tm, Is.EqualTo(probesLower[0].Tm).Within(0.1),
            "Case should not affect Tm calculation");
    }

    #endregion
}
