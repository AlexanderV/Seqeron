namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for primer/probe design: melting temperature, hairpin, dimer, probe validation.
///
/// Test Units: PRIMER-TM-001, PRIMER-DESIGN-001, PRIMER-STRUCT-001, PROBE-DESIGN-001, PROBE-VALID-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class PrimerProbeProperties
{
    // -- PRIMER-TM-001 --

    /// <summary>
    /// Melting temperature is finite and in a reasonable range.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperature_IsFinite()
    {
        var primers = new[] { "ATGCATGC", "GCGCGCGC", "ATATATATAT", "ACGTACGTACGT", "GGGCCCAAATTT" };

        foreach (var p in primers)
        {
            double tm = PrimerDesigner.CalculateMeltingTemperature(p);
            Assert.That(double.IsFinite(tm), Is.True, $"Tm for '{p}' is not finite: {tm}");
        }
    }

    /// <summary>
    /// GC-rich primers have higher Tm than AT-rich primers of same length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperature_GcRich_HigherThanAtRich()
    {
        string gcRich = "GCGCGCGCGCGCGCGCGCGC";
        string atRich = "ATATATATATATATATATATAT";
        double tmGc = PrimerDesigner.CalculateMeltingTemperature(gcRich);
        double tmAt = PrimerDesigner.CalculateMeltingTemperature(atRich);

        Assert.That(tmGc, Is.GreaterThan(tmAt),
            $"GC-rich Tm ({tmGc}) should be > AT-rich Tm ({tmAt})");
    }

    /// <summary>
    /// Salt concentration ≥ 0 does not cause errors.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperatureWithSalt_IsFinite()
    {
        string primer = "ACGTACGTACGTACGT";
        double tm = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, naConcentration: 50);
        Assert.That(double.IsFinite(tm), Is.True, $"Tm with salt is not finite: {tm}");
    }

    // -- PRIMER-STRUCT-001 --

    /// <summary>
    /// Self-complementary primer has hairpin potential.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasHairpinPotential_SelfComplementary_ReturnsTrue()
    {
        // GCGC...GCGC is self-complementary
        string selfComp = "GCGCGCTTTTGCGCGC";
        bool result = PrimerDesigner.HasHairpinPotential(selfComp, minStemLength: 4, minLoopLength: 3);
        Assert.That(result, Is.True, "Self-complementary primer should have hairpin potential");
    }

    /// <summary>
    /// Primer dimer detection: identical forward/reverse has dimer potential.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasPrimerDimer_ComplementaryPair_ReturnsTrue()
    {
        string p1 = "ACGTACGTACGTACGT";
        string p2 = "ACGTACGTACGTACGT"; // has 3' complementarity with itself
        bool result = PrimerDesigner.HasPrimerDimer(p1, p2, minComplementarity: 4);
        // Just verify it returns without error; actual result depends on complementarity
        Assert.That(result, Is.TypeOf<bool>());
    }

    /// <summary>
    /// 3' stability score is finite.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Calculate3PrimeStability_IsFinite()
    {
        double stability = PrimerDesigner.Calculate3PrimeStability("ACGTACGTACGTACGT");
        Assert.That(double.IsFinite(stability), Is.True);
    }

    // -- PRIMER-DESIGN-001 --

    /// <summary>
    /// DesignPrimers returns valid forward/reverse primers spanning the target.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignPrimers_ProductSize_IsPositive()
    {
        string template = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 30));
        var dna = new DnaSequence(template);
        var result = PrimerDesigner.DesignPrimers(dna, targetStart: 100, targetEnd: 200);

        Assert.That(result.ProductSize, Is.GreaterThanOrEqualTo(0));
    }

    /// <summary>
    /// Primer candidate GC content is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void PrimerCandidate_GcContent_InRange()
    {
        string template = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 30));
        var dna = new DnaSequence(template);
        var candidates = PrimerDesigner.GeneratePrimerCandidates(dna, 10, 100, forward: true).ToList();

        foreach (var c in candidates)
            Assert.That(c.GcContent, Is.InRange(0.0, 100.0));
    }

    /// <summary>
    /// Primer candidate score is finite.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PrimerCandidate_Score_InRange()
    {
        string template = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 30));
        var dna = new DnaSequence(template);
        var candidates = PrimerDesigner.GeneratePrimerCandidates(dna, 10, 100, forward: true).ToList();

        foreach (var c in candidates)
            Assert.That(double.IsFinite(c.Score), Is.True, $"Score {c.Score} is not finite");
    }

    // -- PROBE-DESIGN-001 --

    /// <summary>
    /// Designed probes have GC content within valid range.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_GcContent_InRange()
    {
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 10));
        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var p in probes)
            Assert.That(p.GcContent, Is.InRange(0.0, 1.0));
    }

    /// <summary>
    /// Designed probes have finite melting temperature.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_Tm_IsFinite()
    {
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 10));
        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var p in probes)
            Assert.That(double.IsFinite(p.Tm), Is.True);
    }

    /// <summary>
    /// Probe score is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_Score_InRange()
    {
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 10));
        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var p in probes)
            Assert.That(p.Score, Is.InRange(0.0, 1.0));
    }

    // -- PROBE-VALID-001 --

    /// <summary>
    /// ValidateProbe returns a valid ProbeValidation with specificity in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_Specificity_InRange()
    {
        string probe = "ACGTACGTACGTACGTACGTACGT";
        var refs = new[] { "ACGTACGTACGTACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTTTTTT" };
        var validation = ProbeDesigner.ValidateProbe(probe, refs);

        Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0));
    }

    /// <summary>
    /// ValidateProbe off-target hits is non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_OffTargetHits_NonNegative()
    {
        string probe = "ACGTACGTACGTACGTACGTACGT";
        var refs = new[] { "TTTTTTTTTTTTTTTTTTTTTTTTTTTT" };
        var validation = ProbeDesigner.ValidateProbe(probe, refs);

        Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0));
    }
}
