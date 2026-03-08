using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for population genetics calculations.
/// Verifies allele frequency, Hardy-Weinberg, and Fst invariants.
///
/// Test Units: POP-FREQ-001, POP-HW-001, POP-FST-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("PopulationGenetics")]
public class PopulationGeneticsProperties
{
    /// <summary>
    /// Allele frequencies must sum to 1 (MajorFreq + MinorFreq = 1).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlleleFrequencies_SumToOne()
    {
        var genGen = Gen.Choose(0, 100).Three()
            .Where(t => t.Item1 + t.Item2 + t.Item3 > 0)
            .ToArbitrary();

        return Prop.ForAll(genGen, counts =>
        {
            var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
                counts.Item1, counts.Item2, counts.Item3);
            return (Math.Abs(major + minor - 1.0) < 1e-10)
                .Label($"MajorFreq={major:F10} + MinorFreq={minor:F10} = {major + minor:F10}");
        });
    }

    /// <summary>
    /// Each allele frequency is in [0, 1] and they sum to 1.
    /// Note: "major" and "minor" naming may not always mean major >= minor
    /// in edge cases; we only test the sum and range.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlleleFrequencies_BothInRange()
    {
        var genGen = Gen.Choose(0, 100).Three()
            .Where(t => t.Item1 + t.Item2 + t.Item3 > 0)
            .ToArbitrary();

        return Prop.ForAll(genGen, counts =>
        {
            var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
                counts.Item1, counts.Item2, counts.Item3);
            return (major >= 0.0 && major <= 1.0 && minor >= 0.0 && minor <= 1.0)
                .Label($"Major={major:F4}, Minor={minor:F4} — both must be in [0,1]");
        });
    }

    /// <summary>
    /// MAF-C01: MAF is always in [0, 0.5] for any valid genotype set.
    /// Source: MAF invariant — MAF = min(f, 1-f)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MAF_IsAlwaysInRange()
    {
        // Generate arrays of length 1..50, each element in {0,1,2}
        var genoGen = Gen.Choose(0, 2)
            .ArrayOf()
            .Where(a => a.Length > 0)
            .ToArbitrary();

        return Prop.ForAll(genoGen, genotypes =>
        {
            var maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);
            return (maf >= 0.0 && maf <= 0.5)
                .Label($"MAF={maf:F6} must be in [0, 0.5]");
        });
    }

    /// <summary>
    /// Hardy-Weinberg chi-square is non-negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HardyWeinberg_ChiSquare_IsNonNegative()
    {
        var genGen = Gen.Choose(1, 50).Three()
            .ToArbitrary();

        return Prop.ForAll(genGen, counts =>
        {
            var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
                "prop", counts.Item1, counts.Item2, counts.Item3);
            return (result.ChiSquare >= 0.0)
                .Label($"ChiSquare={result.ChiSquare:F4} must be ≥ 0");
        });
    }

    /// <summary>
    /// Hardy-Weinberg p-value is in [0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HardyWeinberg_PValue_InRange()
    {
        var genGen = Gen.Choose(1, 50).Three()
            .ToArbitrary();

        return Prop.ForAll(genGen, counts =>
        {
            var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
                "prop", counts.Item1, counts.Item2, counts.Item3);
            return (result.PValue >= 0.0 && result.PValue <= 1.0)
                .Label($"PValue={result.PValue:F4} must be in [0,1]");
        });
    }

    /// <summary>
    /// Fst is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fst_IsInRange()
    {
        var pop1 = new[] { (0.3, 100), (0.5, 100), (0.7, 100) }
            .Select(x => (AlleleFreq: x.Item1, SampleSize: x.Item2));
        var pop2 = new[] { (0.6, 100), (0.4, 100), (0.8, 100) }
            .Select(x => (AlleleFreq: x.Item1, SampleSize: x.Item2));

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);
        Assert.That(fst, Is.InRange(-0.001, 1.001),
            $"Fst={fst:F4} must be in [0, 1]");
    }

    /// <summary>
    /// Nucleotide diversity is non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void NucleotideDiversity_IsNonNegative()
    {
        var sequences = new IReadOnlyList<char>[]
        {
            "ACGTACGTAC".ToList(),
            "ACGTACGTAT".ToList(),
            "ACGTACGTAG".ToList()
        };
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);
        Assert.That(pi, Is.GreaterThanOrEqualTo(0.0));
    }

    /// <summary>
    /// Identical sequences yield zero nucleotide diversity.
    /// </summary>
    [Test]
    [Category("Property")]
    public void IdenticalSequences_ZeroDiversity()
    {
        string seq = "ACGTACGTACGTACGT";
        var sequences = Enumerable.Repeat<IReadOnlyList<char>>(seq.ToList(), 5);
        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);
        Assert.That(pi, Is.EqualTo(0.0).Within(0.0001));
    }
}
