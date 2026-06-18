using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for population genetics calculations.
/// Verifies allele frequency, Hardy-Weinberg, and Fst invariants.
///
/// Test Units: POP-FREQ-001, POP-DIV-001, POP-HW-001, POP-FST-001, POP-LD-001, POP-ANCESTRY-001, POP-ROH-001, POP-SELECT-001
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

    #region POP-FREQ-001 (extended): D: deterministic

    /// <summary>
    /// INV-F3: Allele frequency calculation is deterministic (pure function).
    /// Evidence: CalculateAlleleFrequencies performs a simple arithmetic computation
    /// with no random state.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlleleFrequencies_AreDeterministic()
    {
        var genGen = Gen.Choose(0, 100).Three()
            .Where(t => t.Item1 + t.Item2 + t.Item3 > 0)
            .ToArbitrary();

        return Prop.ForAll(genGen, counts =>
        {
            var r1 = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
                counts.Item1, counts.Item2, counts.Item3);
            var r2 = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
                counts.Item1, counts.Item2, counts.Item3);
            return (r1.MajorFreq == r2.MajorFreq && r1.MinorFreq == r2.MinorFreq)
                .Label("CalculateAlleleFrequencies must be deterministic");
        });
    }

    #endregion

    #region POP-DIV-001: R: π ≥ 0; R: θ ≥ 0; M: more diverse → higher π; D: deterministic

    /// <summary>
    /// INV-D1: Watterson's theta (θ_W) is non-negative for any valid parameters.
    /// Evidence: θ = S / (a₁ · L) where S ≥ 0, a₁ > 0, L > 0.
    /// Source: Watterson (1975) "On the number of segregating sites".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property WattersonTheta_IsNonNegative()
    {
        var gen = Gen.Choose(0, 50).Two()
            .Select(t => (segregatingSites: t.Item1, sampleSize: t.Item2 + 2))
            .ToArbitrary();

        return Prop.ForAll(gen, data =>
        {
            double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(
                data.segregatingSites, data.sampleSize, 100);
            return (theta >= 0.0)
                .Label($"θ={theta:F6} must be ≥ 0 (S={data.segregatingSites}, n={data.sampleSize})");
        });
    }

    /// <summary>
    /// INV-D2: More divergent sequences produce higher nucleotide diversity.
    /// Evidence: π measures average proportion of pairwise differences;
    /// sequences with more substitutions yield higher π.
    /// Source: Nei &amp; Li (1979).
    /// </summary>
    [Test]
    [Category("Property")]
    public void NucleotideDiversity_MoreDiverse_HigherPi()
    {
        var lowDiv = new IReadOnlyList<char>[]
        {
            "ACGTACGTAC".ToList(),
            "ACGTACGTAC".ToList(),
            "ACGTACGTAT".ToList()
        };

        var highDiv = new IReadOnlyList<char>[]
        {
            "ACGTACGTAC".ToList(),
            "TGCATGCATG".ToList(),
            "GATCGATCGA".ToList()
        };

        double piLow = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(lowDiv);
        double piHigh = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(highDiv);

        Assert.That(piHigh, Is.GreaterThan(piLow),
            $"π(diverse)={piHigh:F4} should exceed π(similar)={piLow:F4}");
    }

    /// <summary>
    /// INV-D3: Diversity statistics calculation is deterministic.
    /// Evidence: CalculateDiversityStatistics is a pure function.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DiversityStatistics_AreDeterministic()
    {
        var seqs = new IReadOnlyList<char>[]
        {
            "ACGTACGTAC".ToList(),
            "ACGTACGTAT".ToList(),
            "ACGTACGTAG".ToList()
        };

        var stats1 = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(seqs);
        var stats2 = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(seqs);

        Assert.That(stats1.NucleotideDiversity, Is.EqualTo(stats2.NucleotideDiversity));
        Assert.That(stats1.WattersonTheta, Is.EqualTo(stats2.WattersonTheta));
        Assert.That(stats1.TajimasD, Is.EqualTo(stats2.TajimasD));
    }

    /// <summary>
    /// INV-D4: Nucleotide diversity is non-negative for random sequence sets.
    /// Evidence: π = average(differences) / length ≥ 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideDiversity_IsNonNegative_Property()
    {
        var seqGen = Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= 10)
            .Select(a => new string(a, 0, 10))
            .ToArbitrary();

        return Prop.ForAll(seqGen, seqGen, seqGen, (s1, s2, s3) =>
        {
            var sequences = new IReadOnlyList<char>[] { s1.ToList(), s2.ToList(), s3.ToList() };
            double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);
            return (pi >= 0.0)
                .Label($"π={pi:F6} must be ≥ 0");
        });
    }

    #endregion

    #region POP-HW-001 (extended): R: p² + 2pq + q² = 1.0; D: deterministic

    /// <summary>
    /// INV-HW3: Expected genotype frequencies under HWE sum to total sample size.
    /// Evidence: Under Hardy-Weinberg, p² + 2pq + q² = (p+q)² = 1,
    /// so ExpectedAA + ExpectedAa + Expectedaa = N.
    /// Source: Hardy (1908), Weinberg (1908).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HardyWeinberg_ExpectedFrequencies_SumToN()
    {
        var genGen = Gen.Choose(1, 50).Three().ToArbitrary();

        return Prop.ForAll(genGen, counts =>
        {
            int n = counts.Item1 + counts.Item2 + counts.Item3;
            var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
                "prop", counts.Item1, counts.Item2, counts.Item3);
            double expectedSum = result.ExpectedAA + result.ExpectedAa + result.Expectedaa;
            return (Math.Abs(expectedSum - n) < 0.01)
                .Label($"Expected sum={expectedSum:F4}, N={n}");
        });
    }

    /// <summary>
    /// INV-HW4: Hardy-Weinberg test is deterministic.
    /// Evidence: TestHardyWeinberg is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HardyWeinberg_IsDeterministic()
    {
        var genGen = Gen.Choose(1, 50).Three().ToArbitrary();

        return Prop.ForAll(genGen, counts =>
        {
            var r1 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
                "prop", counts.Item1, counts.Item2, counts.Item3);
            var r2 = PopulationGeneticsAnalyzer.TestHardyWeinberg(
                "prop", counts.Item1, counts.Item2, counts.Item3);
            return (r1.ChiSquare == r2.ChiSquare && r1.PValue == r2.PValue)
                .Label("TestHardyWeinberg must be deterministic");
        });
    }

    #endregion

    #region POP-FST-001: R: Fst ∈ [0,1]; S: Fst(A,B)=Fst(B,A); I: Fst(A,A)=0; D: deterministic

    /// <summary>
    /// INV-FST1: Fst ∈ [0, 1] for random allele frequencies.
    /// Evidence: Fst = σ²_S / p̄(1-p̄) is non-negative (numerator/denominator ≥ 0)
    /// and ≤ 1 when subpopulation frequencies lie in [0,1].
    /// Source: Wright (1965) Evolution 19:395-420.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Fst_IsInRange_Property()
    {
        var freqGen = Gen.Choose(1, 99).ArrayOf()
            .Where(a => a.Length >= 1 && a.Length <= 5)
            .ToArbitrary();

        return Prop.ForAll(freqGen, freqGen, (freqs1, freqs2) =>
        {
            int len = Math.Min(freqs1.Length, freqs2.Length);
            var pop1 = freqs1.Take(len).Select(f => (AlleleFreq: f / 100.0, SampleSize: 50));
            var pop2 = freqs2.Take(len).Select(f => (AlleleFreq: f / 100.0, SampleSize: 50));
            double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);
            return (fst >= -0.001 && fst <= 1.001)
                .Label($"Fst={fst:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-FST2: Fst(A,B) = Fst(B,A) — Fst is symmetric.
    /// Evidence: The variance formula uses pBar = weighted average of p1, p2
    /// and the variance term (p1-pBar)²n1 + (p2-pBar)²n2 is symmetric under swap.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fst_IsSymmetric()
    {
        var pop1 = new[] { (0.3, 100), (0.5, 80), (0.7, 120) }
            .Select(x => (AlleleFreq: x.Item1, SampleSize: x.Item2)).ToList();
        var pop2 = new[] { (0.6, 90), (0.4, 110), (0.8, 100) }
            .Select(x => (AlleleFreq: x.Item1, SampleSize: x.Item2)).ToList();

        double fstAB = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);
        double fstBA = PopulationGeneticsAnalyzer.CalculateFst(pop2, pop1);

        Assert.That(fstAB, Is.EqualTo(fstBA).Within(1e-10),
            $"Fst(A,B)={fstAB:F6} ≠ Fst(B,A)={fstBA:F6}");
    }

    /// <summary>
    /// INV-FST3: Fst(A,A) = 0 — identical populations have zero differentiation.
    /// Evidence: When p1 = p2 for all loci, variance = 0, so Fst = 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fst_SamePopulation_IsZero()
    {
        var pop = new[] { (0.3, 100), (0.5, 100), (0.7, 100) }
            .Select(x => (AlleleFreq: x.Item1, SampleSize: x.Item2)).ToList();

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop, pop);

        Assert.That(fst, Is.EqualTo(0.0).Within(1e-10),
            $"Fst(A,A)={fst:F6} should be 0");
    }

    /// <summary>
    /// INV-FST4: Fst calculation is deterministic.
    /// Evidence: CalculateFst is a pure function.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fst_IsDeterministic()
    {
        var pop1 = new[] { (0.3, 100), (0.5, 80) }
            .Select(x => (AlleleFreq: x.Item1, SampleSize: x.Item2)).ToList();
        var pop2 = new[] { (0.6, 90), (0.4, 110) }
            .Select(x => (AlleleFreq: x.Item1, SampleSize: x.Item2)).ToList();

        double fst1 = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);
        double fst2 = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst1, Is.EqualTo(fst2));
    }

    #endregion

    #region POP-LD-001: R: D' ∈ [0,1]; R: r² ∈ [0,1]; S: LD(a,b)=LD(b,a); D: deterministic

    /// <summary>
    /// INV-LD1: r² ∈ [0, 1] for random genotype data.
    /// Evidence: r² is the squared Pearson correlation of genotype values, bounded [0,1].
    /// Source: Hill &amp; Robertson (1968).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LD_RSquared_InRange()
    {
        var genoGen = Gen.Choose(0, 2).ArrayOf()
            .Where(a => a.Length >= 10 && a.Length <= 50)
            .ToArbitrary();

        return Prop.ForAll(genoGen, genoGen, (g1, g2) =>
        {
            int len = Math.Min(g1.Length, g2.Length);
            var pairs = g1.Take(len).Zip(g2.Take(len), (a, b) => (Geno1: a, Geno2: b));
            var ld = PopulationGeneticsAnalyzer.CalculateLD("v1", "v2", pairs, 100);
            return (ld.RSquared >= -0.001 && ld.RSquared <= 1.001)
                .Label($"r²={ld.RSquared:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-LD2: D' ∈ [0, 1] for random genotype data.
    /// Evidence: D' = |D|/D_max, clamped to [0,1] by construction.
    /// Source: Lewontin (1964).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LD_DPrime_InRange()
    {
        var genoGen = Gen.Choose(0, 2).ArrayOf()
            .Where(a => a.Length >= 10 && a.Length <= 50)
            .ToArbitrary();

        return Prop.ForAll(genoGen, genoGen, (g1, g2) =>
        {
            int len = Math.Min(g1.Length, g2.Length);
            var pairs = g1.Take(len).Zip(g2.Take(len), (a, b) => (Geno1: a, Geno2: b));
            var ld = PopulationGeneticsAnalyzer.CalculateLD("v1", "v2", pairs, 100);
            return (ld.DPrime >= -0.001 && ld.DPrime <= 1.001)
                .Label($"D'={ld.DPrime:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-LD3: LD is symmetric — LD(v1,v2) = LD(v2,v1) in r² and D'.
    /// Evidence: Pearson correlation and covariance are symmetric operations;
    /// D_max formula is also symmetric under variable swap.
    /// </summary>
    [Test]
    [Category("Property")]
    public void LD_IsSymmetric()
    {
        var genotypes = new (int Geno1, int Geno2)[]
        {
            (0, 0), (0, 1), (1, 0), (1, 1), (2, 2),
            (0, 2), (1, 1), (2, 0), (0, 1), (1, 2),
            (0, 0), (1, 0), (2, 1), (0, 2), (1, 1)
        };

        var ldForward = PopulationGeneticsAnalyzer.CalculateLD(
            "v1", "v2", genotypes, 100);
        var ldReverse = PopulationGeneticsAnalyzer.CalculateLD(
            "v2", "v1", genotypes.Select(g => (Geno1: g.Geno2, Geno2: g.Geno1)), 100);

        Assert.That(ldForward.RSquared, Is.EqualTo(ldReverse.RSquared).Within(1e-10),
            $"r² not symmetric: {ldForward.RSquared:F6} vs {ldReverse.RSquared:F6}");
        Assert.That(ldForward.DPrime, Is.EqualTo(ldReverse.DPrime).Within(1e-10),
            $"D' not symmetric: {ldForward.DPrime:F6} vs {ldReverse.DPrime:F6}");
    }

    /// <summary>
    /// INV-LD4: LD calculation is deterministic.
    /// Evidence: CalculateLD is a pure function.
    /// </summary>
    [Test]
    [Category("Property")]
    public void LD_IsDeterministic()
    {
        var genotypes = new (int Geno1, int Geno2)[]
        {
            (0, 0), (0, 1), (1, 0), (1, 1), (2, 2),
            (0, 2), (1, 1), (2, 0), (0, 1), (1, 2)
        };

        var r1 = PopulationGeneticsAnalyzer.CalculateLD("v1", "v2", genotypes, 100);
        var r2 = PopulationGeneticsAnalyzer.CalculateLD("v1", "v2", genotypes, 100);

        Assert.That(r1.RSquared, Is.EqualTo(r2.RSquared));
        Assert.That(r1.DPrime, Is.EqualTo(r2.DPrime));
    }

    #endregion

    #region POP-ANCESTRY-001: R: each proportion ∈ [0,1]; P: Σ proportions = 1.0; D: deterministic

    // EstimateAncestry runs the FRAPPE EM with fixed reference allele frequencies (Alexander et al.
    // 2009); each individual's ancestry fractions are a probability simplex summing to 1.

    private static Arbitrary<(List<(string, IReadOnlyList<int>)> inds, List<(string, IReadOnlyList<double>)> refs)>
        AncestryArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int k = 2 + rng.Next(2);   // 2..3 reference populations
            int m = 5 + rng.Next(6);   // 5..10 markers
            var refs = new List<(string, IReadOnlyList<double>)>();
            for (int p = 0; p < k; p++)
                refs.Add(($"pop{p}", (IReadOnlyList<double>)Enumerable.Range(0, m).Select(_ => rng.NextDouble()).ToList()));
            var inds = new List<(string, IReadOnlyList<int>)>();
            for (int i = 0; i < 1 + rng.Next(3); i++)
                inds.Add(($"ind{i}", (IReadOnlyList<int>)Enumerable.Range(0, m).Select(_ => rng.Next(0, 3)).ToList()));
            return (inds, refs);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): each individual's ancestry proportions are in [0,1], sum to 1, and are keyed by
    /// the reference population ids.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ancestry_ProportionsAreSimplex()
    {
        return Prop.ForAll(AncestryArbitrary(), input =>
        {
            var (inds, refs) = input;
            var popIds = refs.Select(r => r.Item1).ToHashSet();
            var results = PopulationGeneticsAnalyzer.EstimateAncestry(inds, refs).ToList();
            bool ok = results.All(a =>
                a.Proportions.Values.All(v => v is >= -1e-9 and <= 1.0 + 1e-9) &&
                Math.Abs(a.Proportions.Values.Sum() - 1.0) < 1e-6 &&
                a.Proportions.Keys.ToHashSet().SetEquals(popIds));
            return ok.Label("ancestry proportions were not a simplex over the reference populations");
        });
    }

    /// <summary>
    /// INV-2 (D): Ancestry estimation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ancestry_IsDeterministic()
    {
        return Prop.ForAll(AncestryArbitrary(), input =>
        {
            var (inds, refs) = input;
            var a = PopulationGeneticsAnalyzer.EstimateAncestry(inds, refs).SelectMany(x => x.Proportions.OrderBy(p => p.Key).Select(p => p.Value)).ToList();
            var b = PopulationGeneticsAnalyzer.EstimateAncestry(inds, refs).SelectMany(x => x.Proportions.OrderBy(p => p.Key).Select(p => p.Value)).ToList();
            return a.SequenceEqual(b).Label("EstimateAncestry must be deterministic");
        });
    }

    #endregion

    #region POP-ROH-001: R: ROH start < end; M: lower minLength → ≥ ROH; P: homozygous within run; D: deterministic

    // FindROH detects runs of homozygosity: maximal runs that begin/end on a homozygous SNP and
    // tolerate at most maxHeterozygotes interior heterozygous calls (PLINK ROH model).

    /// <summary>Generates 10..30 SNPs at increasing positions, mostly homozygous (0/2) with some het (1).</summary>
    private static Arbitrary<(int Position, int Genotype)[]> RohGenotypesArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 10 + rng.Next(21);
            var g = new (int, int)[n];
            for (int i = 0; i < n; i++)
            {
                int geno = rng.Next(10) < 7 ? (rng.Next(2) == 0 ? 0 : 2) : 1; // ~70% homozygous
                g[i] = (i * 1000, geno);
            }
            return g;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): every ROH starts and ends on a homozygous SNP at Start &lt; End, spans ≥ minSnps
    /// SNPs, and tolerates at most maxHeterozygotes interior heterozygotes.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Roh_RunsAreHomozygous()
    {
        return Prop.ForAll(RohGenotypesArbitrary(), genotypes =>
        {
            const int maxHet = 1;
            var byPos = genotypes.ToDictionary(g => g.Position, g => g.Genotype);
            var rohs = PopulationGeneticsAnalyzer.FindROH(genotypes, minSnps: 2, minLength: 0, maxHeterozygotes: maxHet, maxGap: 100_000).ToList();
            bool ok = rohs.All(r =>
            {
                if (!(r.Start < r.End && r.Start >= 0 && r.SnpCount >= 2)) return false;
                if (byPos[r.Start] == 1 || byPos[r.End] == 1) return false; // bounded by homozygous SNPs
                int het = genotypes.Count(g => g.Position >= r.Start && g.Position <= r.End && g.Genotype == 1);
                return het <= maxHet;
            });
            return ok.Label("a ROH violated the homozygosity/boundary/tolerance rules");
        });
    }

    /// <summary>
    /// INV-2 (M): a lower minimum length never reports fewer ROH.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Roh_LowerMinLength_MoreRuns()
    {
        return Prop.ForAll(RohGenotypesArbitrary(), genotypes =>
        {
            int loose = PopulationGeneticsAnalyzer.FindROH(genotypes, 2, 0, 1, 100_000).Count();
            int strict = PopulationGeneticsAnalyzer.FindROH(genotypes, 2, 15_000, 1, 100_000).Count();
            return (loose >= strict).Label($"lower minLength gave fewer ROH ({loose} < {strict})");
        });
    }

    /// <summary>
    /// INV-3 (D + positive control): an all-homozygous stretch is one ROH; detection is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Roh_AllHomozygous_IsOneRun_AndDeterministic()
    {
        var genotypes = Enumerable.Range(0, 10).Select(i => (i * 1000, 0)).ToArray();
        var a = PopulationGeneticsAnalyzer.FindROH(genotypes, 2, 0, 0, 100_000).ToList();
        var b = PopulationGeneticsAnalyzer.FindROH(genotypes, 2, 0, 0, 100_000).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(a, Has.Count.EqualTo(1), "a fully homozygous stretch is a single ROH");
            Assert.That(a[0].Start, Is.EqualTo(0));
            Assert.That(a[0].End, Is.EqualTo(9000));
            Assert.That(b.Select(r => (r.Start, r.End)), Is.EqualTo(a.Select(r => (r.Start, r.End))), "deterministic");
        });
    }

    #endregion

    #region POP-SELECT-001: R: statistic finite; M: stronger differentiation → higher signal; D: deterministic

    // ScanForSelection flags regions whose Tajima's D / Fst / iHS exceed selection thresholds; the
    // signal Score carries the test statistic. CalculateTajimasD is a finite summary statistic.

    /// <summary>
    /// INV-1 (R): Tajima's D is finite for valid samples (S ≥ 1, n ≥ 4, kHat ≥ 0).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TajimasD_IsFinite()
    {
        var gen = Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            return (kHat: rng.NextDouble() * 10, s: 1 + rng.Next(50), n: 4 + rng.Next(47));
        }).ToArbitrary();

        return Prop.ForAll(gen, d =>
            double.IsFinite(PopulationGeneticsAnalyzer.CalculateTajimasD(d.kHat, d.s, d.n))
                .Label("Tajima's D must be finite"));
    }

    /// <summary>
    /// INV-2 (R): every emitted selection signal has a finite score, ordered coordinates, and a
    /// non-empty test type.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Selection_SignalsAreWellFormed()
    {
        var gen = Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            return Enumerable.Range(0, 1 + rng.Next(4))
                .Select(i => ($"r{i}", i * 1000, i * 1000 + 999,
                    (rng.NextDouble() - 0.5) * 6, rng.NextDouble(), (rng.NextDouble() - 0.5) * 6))
                .ToList();
        }).ToArbitrary();

        return Prop.ForAll(gen, regions =>
        {
            var signals = PopulationGeneticsAnalyzer.ScanForSelection(regions).ToList();
            return signals.All(s => double.IsFinite(s.Score) && s.Start <= s.End && !string.IsNullOrEmpty(s.TestType))
                .Label("a selection signal was malformed");
        });
    }

    /// <summary>
    /// INV-3 (M): stronger differentiation (higher Fst) above the threshold is flagged, and a higher
    /// Fst gives a higher Fst-signal score; a below-threshold Fst is not flagged.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Selection_HigherFst_StrongerSignal()
    {
        (string, int, int, double, double, double) Region(string id, double fst) => (id, 0, 999, 0.0, fst, 0.0);

        var lowFst = PopulationGeneticsAnalyzer.ScanForSelection(new[] { Region("r", 0.10) }).Where(s => s.TestType == "Fst").ToList();
        var midFst = PopulationGeneticsAnalyzer.ScanForSelection(new[] { Region("r", 0.30) }).Where(s => s.TestType == "Fst").ToList();
        var highFst = PopulationGeneticsAnalyzer.ScanForSelection(new[] { Region("r", 0.60) }).Where(s => s.TestType == "Fst").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(lowFst, Is.Empty, "Fst below threshold is not flagged");
            Assert.That(midFst, Has.Count.EqualTo(1), "Fst above threshold is flagged");
            Assert.That(highFst.Single().Score, Is.GreaterThan(midFst.Single().Score), "higher Fst → higher signal score");
        });
    }

    /// <summary>
    /// INV-4 (D): Selection scanning is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Selection_IsDeterministic()
    {
        var regions = new[] { ("r0", 0, 999, -2.5, 0.4, 3.0), ("r1", 1000, 1999, 0.5, 0.1, 0.2) };
        var a = PopulationGeneticsAnalyzer.ScanForSelection(regions).Select(s => (s.TestType, s.Score)).ToList();
        var b = PopulationGeneticsAnalyzer.ScanForSelection(regions).Select(s => (s.TestType, s.Score)).ToList();
        Assert.That(b, Is.EqualTo(a));
    }

    #endregion
}
