namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for population genetics.
/// Verifies Hardy-Weinberg, allele frequency, and diversity output stability.
///
/// Test Units: POP-HW-001, POP-FREQ-001, POP-DIV-001, POP-LD-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("PopulationGenetics")]
public class PopulationSnapshotTests
{
    [Test]
    public Task HardyWeinberg_KnownGenotypes_MatchesSnapshot()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("SNP_rs123", 40, 45, 15);
        return Verify(new
        {
            result.VariantId,
            result.ObservedAA,
            result.ObservedAa,
            result.Observedaa,
            ExpectedAA = Math.Round(result.ExpectedAA, 2),
            ExpectedAa = Math.Round(result.ExpectedAa, 2),
            Expectedaa = Math.Round(result.Expectedaa, 2),
            ChiSquare = Math.Round(result.ChiSquare, 4),
            PValue = Math.Round(result.PValue, 4),
            result.InEquilibrium
        });
    }

    [Test]
    public Task AlleleFrequencies_KnownCounts_MatchesSnapshot()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(45, 40, 15);
        return Verify(new
        {
            MajorFreq = Math.Round(major, 4),
            MinorFreq = Math.Round(minor, 4)
        });
    }

    [Test]
    public Task DiversityStatistics_KnownSequences_MatchesSnapshot()
    {
        var sequences = new IReadOnlyList<char>[]
        {
            "ACGTACGTAC".ToList(),
            "ACGTACGTAT".ToList(),
            "ACGTACGTAG".ToList(),
            "ACGTACGTAC".ToList(),
            "ACGTACGTTC".ToList()
        };
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);
        return Verify(new
        {
            Pi = Math.Round(stats.NucleotideDiversity, 6),
            Theta = Math.Round(stats.WattersonTheta, 6),
            TajimasD = Math.Round(stats.TajimasD, 4),
            stats.SegregratingSites,
            stats.SampleSize,
            HetObs = Math.Round(stats.HeterozygosityObserved, 4),
            HetExp = Math.Round(stats.HeterozygosityExpected, 4)
        });
    }

    [Test]
    public Task LinkageDisequilibrium_KnownGenotypes_MatchesSnapshot()
    {
        // Simulate genotype pairs: 0 = AA, 1 = Aa, 2 = aa
        var genotypes = new[]
        {
            (Geno1: 0, Geno2: 0), (Geno1: 0, Geno2: 0),
            (Geno1: 1, Geno2: 1), (Geno1: 1, Geno2: 1),
            (Geno1: 2, Geno2: 2), (Geno1: 2, Geno2: 2),
            (Geno1: 0, Geno2: 1), (Geno1: 1, Geno2: 0),
            (Geno1: 0, Geno2: 0), (Geno1: 2, Geno2: 2)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("rs123", "rs456", genotypes, distance: 5000);

        return Verify(new
        {
            ld.Variant1,
            ld.Variant2,
            DPrime = Math.Round(ld.DPrime, 4),
            RSquared = Math.Round(ld.RSquared, 4),
            ld.Distance
        });
    }
}
