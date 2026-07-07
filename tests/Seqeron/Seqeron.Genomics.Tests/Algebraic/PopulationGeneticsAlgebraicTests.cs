using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the PopGen area (allele frequencies, Hardy–Weinberg,
/// Fst, linkage disequilibrium).
///
/// Algebraic testing pins the conservation identities (frequencies and HWE
/// genotype proportions summing to one), the metric/identity behaviour of Fst,
/// and the locus-swap symmetry of pairwise LD.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 43, 45, 46, 47.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("PopGen")]
public class PopulationGeneticsAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-FREQ-001 — Allele frequencies (PopGen)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 43.
    //
    // Model: from diploid genotype counts (AA, Aa, aa) the major/minor allele
    //        frequencies are (2·AA+Aa)/2N and (2·aa+Aa)/2N. The two frequencies
    //        partition the allele pool.
    //   — docs/algorithms/Population_Genetics; PopulationGeneticsAnalyzer.CalculateAlleleFrequencies.
    //
    // Laws (row 43): ID — a monomorphic sample (all AA) → major frequency 1.0.
    //                DIST — majorFreq + minorFreq = 1 for any non-empty sample.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>ID: an all-homozygous-major sample has major allele frequency 1.0.</summary>
    [FsCheck.NUnit.Property]
    public Property AlleleFreq_Identity_MonomorphicIsOne()
    {
        return Prop.ForAll(Gen.Choose(1, 1000).ToArbitrary(), n =>
        {
            var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(n, 0, 0);
            return (major == 1.0 && minor == 0.0).Label($"n={n}: major={major}, minor={minor}");
        });
    }

    /// <summary>DIST: majorFreq + minorFreq = 1 for any non-empty genotype sample.</summary>
    [FsCheck.NUnit.Property]
    public Property AlleleFreq_Distributive_SumsToOne()
    {
        var counts = (from aa in Gen.Choose(0, 500)
                      from ab in Gen.Choose(0, 500)
                      from bb in Gen.Choose(0, 500)
                      where aa + ab + bb > 0
                      select (aa, ab, bb)).ToArbitrary();
        return Prop.ForAll(counts, t =>
        {
            var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(t.aa, t.ab, t.bb);
            return (Math.Abs(major + minor - 1.0) < 1e-12).Label($"major+minor={major + minor}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-HW-001 — Hardy–Weinberg equilibrium (PopGen)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 45.
    //
    // Model: under HWE the expected genotype proportions are p², 2pq, q² with
    //        p+q = 1, so p² + 2pq + q² = (p+q)² = 1 and the expected counts sum to N.
    //   — docs/algorithms/Population_Genetics; PopulationGeneticsAnalyzer.TestHardyWeinberg.
    //
    // Laws (row 45): DIST — p² + 2pq + q² = 1 and Σ expected = N.
    //                ID — a monomorphic sample (all AA) puts every individual in one
    //                genotype class and is trivially in equilibrium.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>DIST: the HWE genotype proportions sum to 1 and expected counts sum to N.</summary>
    [FsCheck.NUnit.Property]
    public Property HardyWeinberg_Distributive_ProportionsSumToOne()
    {
        var counts = (from aa in Gen.Choose(0, 500)
                      from ab in Gen.Choose(0, 500)
                      from bb in Gen.Choose(0, 500)
                      where aa + ab + bb > 0
                      select (aa, ab, bb)).ToArbitrary();
        return Prop.ForAll(counts, t =>
        {
            int n = t.aa + t.ab + t.bb;
            double p = (2.0 * t.aa + t.ab) / (2.0 * n);
            double q = 1 - p;
            var hw = PopulationGeneticsAnalyzer.TestHardyWeinberg("v", t.aa, t.ab, t.bb);
            double propSum = p * p + 2 * p * q + q * q;
            double expSum = hw.ExpectedAA + hw.ExpectedAa + hw.Expectedaa;
            return (Math.Abs(propSum - 1.0) < 1e-12 && Math.Abs(expSum - n) < 1e-9)
                .Label($"propSum={propSum}, expSum={expSum}, n={n}");
        });
    }

    /// <summary>ID: a monomorphic sample concentrates all expectation in one genotype.</summary>
    [Test]
    public void HardyWeinberg_Identity_MonomorphicIsSingleGenotype()
    {
        var hw = PopulationGeneticsAnalyzer.TestHardyWeinberg("v", observedAA: 200, observedAa: 0, observedaa: 0);
        hw.ExpectedAA.Should().BeApproximately(200, 1e-9);
        hw.ExpectedAa.Should().BeApproximately(0, 1e-9);
        hw.Expectedaa.Should().BeApproximately(0, 1e-9);
        hw.InEquilibrium.Should().BeTrue();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-FST-001 — Wright's Fst (PopGen)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 46.
    //
    // Model: Fst = weighted variance of allele frequencies between populations /
    //        mean heterozygosity. Identical populations have zero between-population
    //        variance, and the statistic is symmetric in the two populations.
    //   — docs/algorithms/Population_Genetics; PopulationGeneticsAnalyzer.CalculateFst.
    //
    // Laws (row 46): ID — Fst(identical pops) = 0.  COMM — Fst(a,b) = Fst(b,a).
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<List<(double AlleleFreq, int SampleSize)>> PopulationArbitrary() =>
        (from k in Gen.Choose(1, 8)
         from freqs in Gen.Choose(0, 100).Select(x => x / 100.0).ArrayOf(k)
         from sizes in Gen.Choose(1, 200).ArrayOf(k)
         select freqs.Zip(sizes, (f, s) => (f, s)).ToList())
        .ToArbitrary();

    /// <summary>ID: Fst of a population against an identical copy is 0.</summary>
    [FsCheck.NUnit.Property]
    public Property Fst_Identity_IdenticalPopulationsIsZero()
    {
        return Prop.ForAll(PopulationArbitrary(), pop =>
        {
            double fst = PopulationGeneticsAnalyzer.CalculateFst(pop, pop);
            return (Math.Abs(fst) < 1e-12).Label($"Fst(identical)={fst}");
        });
    }

    /// <summary>COMM: Fst(a, b) = Fst(b, a).</summary>
    [FsCheck.NUnit.Property]
    public Property Fst_Commutative_Symmetric()
    {
        return Prop.ForAll(PopulationArbitrary(), PopulationArbitrary(), (a, b) =>
        {
            if (a.Count != b.Count) return true.ToProperty(); // Fst requires matching loci counts
            double ab = PopulationGeneticsAnalyzer.CalculateFst(a, b);
            double ba = PopulationGeneticsAnalyzer.CalculateFst(b, a);
            return (Math.Abs(ab - ba) < 1e-12).Label($"Fst(a,b)={ab} != Fst(b,a)={ba}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-LD-001 — Linkage disequilibrium (PopGen)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 47.
    //
    // Model: pairwise LD between two loci from diploid 0/1/2 genotypes: r² is the
    //        squared genotype correlation and D' the normalized covariance. Both
    //        are symmetric functions of the two loci, and independent (orthogonal)
    //        loci have zero covariance hence D' = r² = 0.
    //   — docs/algorithms/Population_Genetics; PopulationGeneticsAnalyzer.CalculateLD.
    //
    // Laws (row 47): ID — independent loci → D' ≈ 0.  COMM — LD(a,b) = LD(b,a).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: a fully crossed orthogonal design (each (0/2)×(0/2) cell equally
    /// represented) has zero genotype covariance, so D' = 0 and r² = 0.
    /// </summary>
    [Test]
    public void Ld_Identity_IndependentLociIsZero()
    {
        var genotypes = new List<(int, int)> { (0, 0), (0, 2), (2, 0), (2, 2) };
        var ld = PopulationGeneticsAnalyzer.CalculateLD("v1", "v2", genotypes, distance: 100);
        ld.DPrime.Should().BeApproximately(0.0, 1e-12);
        ld.RSquared.Should().BeApproximately(0.0, 1e-12);
    }

    /// <summary>COMM: swapping the two loci leaves D' and r² unchanged.</summary>
    [FsCheck.NUnit.Property]
    public Property Ld_Commutative_LocusSwapInvariant()
    {
        var genoPairs = (from g1 in Gen.Choose(0, 2)
                         from g2 in Gen.Choose(0, 2)
                         select (g1, g2)).ArrayOf().Where(a => a.Length >= 2).ToArbitrary();
        return Prop.ForAll(genoPairs, pairs =>
        {
            var ab = PopulationGeneticsAnalyzer.CalculateLD("a", "b", pairs.Select(p => (p.Item1, p.Item2)), 1);
            var ba = PopulationGeneticsAnalyzer.CalculateLD("b", "a", pairs.Select(p => (p.Item2, p.Item1)), 1);
            return (Math.Abs(ab.DPrime - ba.DPrime) < 1e-12 && Math.Abs(ab.RSquared - ba.RSquared) < 1e-12)
                .Label($"D'={ab.DPrime}/{ba.DPrime}, r2={ab.RSquared}/{ba.RSquared}");
        });
    }
}
