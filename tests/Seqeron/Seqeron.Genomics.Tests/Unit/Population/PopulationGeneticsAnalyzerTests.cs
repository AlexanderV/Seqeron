namespace Seqeron.Genomics.Tests.Unit.Population;

[TestFixture]
public class PopulationGeneticsAnalyzerTests
{
    // Note: Allele Frequency tests (POP-FREQ-001) have been moved to
    // PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs

    // Note: Diversity Statistics tests (POP-DIV-001) have been moved to
    // PopulationGeneticsAnalyzer_Diversity_Tests.cs

    // Note: Hardy-Weinberg tests (POP-HW-001) have been moved to
    // PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs

    // Note: F-Statistics tests (POP-FST-001) have been moved to
    // PopulationGeneticsAnalyzer_FStatistics_Tests.cs

    // Note: Linkage Disequilibrium tests (POP-LD-001) have been moved to
    // PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs

    #region Selection Tests

    [Test]
    public void CalculateIHS_BalancedEHH_ReturnsNearZero()
    {
        var ehh0 = new List<double> { 1.0, 0.8, 0.5, 0.2, 0.1 };
        var ehh1 = new List<double> { 1.0, 0.8, 0.5, 0.2, 0.1 };
        var positions = new List<int> { 0, 1000, 2000, 3000, 4000 };

        double ihs = PopulationGeneticsAnalyzer.CalculateIHS(ehh0, ehh1, positions);

        Assert.That(Math.Abs(ihs), Is.LessThan(0.1));
    }

    [Test]
    public void CalculateIHS_ExtendedDerived_ReturnsPositive()
    {
        var ehh0 = new List<double> { 1.0, 0.5, 0.1, 0.05, 0.01 }; // Rapid decay
        var ehh1 = new List<double> { 1.0, 0.9, 0.8, 0.7, 0.6 };   // Extended
        var positions = new List<int> { 0, 1000, 2000, 3000, 4000 };

        double ihs = PopulationGeneticsAnalyzer.CalculateIHS(ehh0, ehh1, positions);

        Assert.That(ihs, Is.GreaterThan(0));
    }

    [Test]
    public void ScanForSelection_NegativeTajimaD_DetectsSignal()
    {
        var regions = new List<(string, int, int, double, double, double)>
        {
            ("Region1", 0, 10000, -2.5, 0.1, 0.5) // Significant Tajima's D
        };

        var signals = PopulationGeneticsAnalyzer.ScanForSelection(
            regions, tajimaDThreshold: -2.0).ToList();

        Assert.That(signals, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(signals.Any(s => s.TestType == "TajimasD"), Is.True);
    }

    [Test]
    public void ScanForSelection_HighFst_DetectsSignal()
    {
        var regions = new List<(string, int, int, double, double, double)>
        {
            ("Region1", 0, 10000, 0, 0.5, 0) // High Fst
        };

        var signals = PopulationGeneticsAnalyzer.ScanForSelection(
            regions, fstThreshold: 0.25).ToList();

        Assert.That(signals.Any(s => s.TestType == "Fst"), Is.True);
    }

    [Test]
    public void ScanForSelection_NoSignificantSignals_ReturnsEmpty()
    {
        var regions = new List<(string, int, int, double, double, double)>
        {
            ("Region1", 0, 10000, 0, 0.1, 0.5) // All neutral
        };

        var signals = PopulationGeneticsAnalyzer.ScanForSelection(regions).ToList();

        Assert.That(signals, Is.Empty);
    }

    #endregion

    // Ancestry Analysis tests moved to the canonical unit file
    // PopulationGeneticsAnalyzer_EstimateAncestry_Tests.cs (POP-ANCESTRY-001).

    // Inbreeding (ROH / F_ROH) tests moved to the canonical unit file
    // PopulationGeneticsAnalyzer_FindROH_Tests.cs (POP-ROH-001).

    #region Edge Cases Tests

    // Note: Diversity edge cases (CalculateNucleotideDiversity, CalculateWattersonTheta,
    // CalculateTajimasD) have been moved to PopulationGeneticsAnalyzer_Diversity_Tests.cs

    // Note: CalculateFst edge case tests (POP-FST-001) have been moved to
    // PopulationGeneticsAnalyzer_FStatistics_Tests.cs

    // Note: CalculateLD and FindHaplotypeBlocks edge case tests (POP-LD-001) have been moved to
    // PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs

    [Test]
    public void CalculateIHS_InsufficientData_ReturnsZero()
    {
        double ihs = PopulationGeneticsAnalyzer.CalculateIHS(
            new List<double> { 1.0 },
            new List<double> { 1.0 },
            new List<int> { 0 });

        Assert.That(ihs, Is.EqualTo(0));
    }

    #endregion
}
