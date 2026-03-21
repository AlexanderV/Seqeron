using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for metagenomics diversity analysis:
/// alpha diversity indices and beta diversity metrics.
///
/// Test Units: META-ALPHA-001, META-BETA-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Metagenomics")]
public class MetagenomicsProperties
{
    #region Generators

    /// <summary>
    /// Generates a random abundance dictionary with N species, values normalized to sum to 1.0.
    /// </summary>
    private static Arbitrary<IReadOnlyDictionary<string, double>> AbundanceArbitrary(int minSpecies = 2) =>
        Gen.Choose(minSpecies, 20)
            .SelectMany(n =>
                Gen.Choose(1, 100).ArrayOf(n)
                    .Select(counts =>
                    {
                        double total = counts.Sum();
                        var dict = new Dictionary<string, double>();
                        for (int i = 0; i < n; i++)
                            dict[$"species_{i + 1}"] = counts[i] / total;
                        return (IReadOnlyDictionary<string, double>)dict;
                    }))
            .ToArbitrary();

    /// <summary>
    /// Creates a uniform abundance for N species (equal proportions).
    /// </summary>
    private static IReadOnlyDictionary<string, double> UniformAbundance(int n)
    {
        var dict = new Dictionary<string, double>();
        for (int i = 0; i < n; i++)
            dict[$"species_{i + 1}"] = 1.0 / n;
        return dict;
    }

    /// <summary>
    /// Creates a single-species dominance abundance.
    /// </summary>
    private static IReadOnlyDictionary<string, double> SingleSpeciesAbundance()
    {
        return new Dictionary<string, double> { ["species_1"] = 1.0 };
    }

    #endregion

    #region META-ALPHA-001: R: Shannon ≥ 0; R: Simpson ∈ [0,1]; M: more species → higher diversity; D: deterministic

    /// <summary>
    /// INV-1: Shannon index H ≥ 0 for any abundance distribution.
    /// Evidence: H = −∑ pᵢ ln(pᵢ) where pᵢ ∈ (0,1] → each term is ≥ 0, so H ≥ 0.
    /// Source: Shannon (1948); Wikipedia "Diversity index".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_ShannonIndex_NonNegative()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.ShannonIndex >= -1e-9)
                .Label($"Shannon={result.ShannonIndex:F6} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2: Simpson index λ ∈ [0, 1] for any abundance distribution.
    /// Evidence: λ = ∑ pᵢ² where pᵢ ∈ [0,1] and ∑pᵢ = 1; 
    /// by Cauchy-Schwarz: 1/S ≤ λ ≤ 1.
    /// Source: Simpson (1949); Wikipedia "Diversity index".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_SimpsonIndex_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.SimpsonIndex >= -1e-9 && result.SimpsonIndex <= 1.0 + 1e-9)
                .Label($"Simpson={result.SimpsonIndex:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-3: Inverse Simpson ≥ 1.0 for any distribution with at least one species.
    /// Evidence: D = 1/λ where λ ≤ 1, so D ≥ 1.
    /// Source: Hill (1973); Wikipedia "Diversity index".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_InverseSimpson_AtLeast1()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.InverseSimpson >= 1.0 - 1e-9)
                .Label($"InverseSimpson={result.InverseSimpson:F4} must be ≥ 1");
        });
    }

    /// <summary>
    /// INV-4: Shannon entropy is maximized (H = ln S) when all species are equally abundant.
    /// Evidence: By Gibbs' inequality, uniform distribution maximizes entropy.
    /// Source: Shannon (1948); Wikipedia "Diversity index".
    /// </summary>
    [TestCase(2)]
    [TestCase(5)]
    [TestCase(10)]
    [Category("Property")]
    public void AlphaDiversity_UniformDistribution_MaxShannon(int n)
    {
        var abundances = UniformAbundance(n);
        var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        double expectedMaxH = Math.Log(n);
        Assert.That(result.ShannonIndex, Is.EqualTo(expectedMaxH).Within(0.01),
            $"Uniform distribution of {n} species should have H = ln({n}) = {expectedMaxH:F4}, got {result.ShannonIndex:F4}");
    }

    /// <summary>
    /// INV-5: Single species → Shannon H = 0 (no uncertainty).
    /// Evidence: H = −(1 × ln 1) = 0.
    /// Source: Shannon (1948).
    /// </summary>
    [Test]
    [Category("Property")]
    public void AlphaDiversity_SingleSpecies_ShannonIsZero()
    {
        var abundances = SingleSpeciesAbundance();
        var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.That(result.ShannonIndex, Is.EqualTo(0.0).Within(0.001),
            $"Single species should have Shannon H = 0, got {result.ShannonIndex:F4}");
    }

    /// <summary>
    /// INV-6: Pielou evenness J ∈ [0, 1] for distributions with &gt; 1 species.
    /// Evidence: J = H / ln(S) where H ≤ ln(S), so J ∈ [0, 1].
    /// Source: Pielou (1966); Wikipedia "Species evenness".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_PielouEvenness_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(2), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.PielouEvenness >= -1e-9 && result.PielouEvenness <= 1.0 + 1e-9)
                .Label($"Pielou J={result.PielouEvenness:F4} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-7: Alpha diversity is deterministic.
    /// Evidence: CalculateAlphaDiversity is a pure mathematical function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_IsDeterministic()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            var r2 = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            bool same = Math.Abs(r1.ShannonIndex - r2.ShannonIndex) < 1e-10 &&
                        Math.Abs(r1.SimpsonIndex - r2.SimpsonIndex) < 1e-10 &&
                        r1.ObservedSpecies == r2.ObservedSpecies;
            return same.Label("CalculateAlphaDiversity must be deterministic");
        });
    }

    /// <summary>
    /// INV-8: ObservedSpecies equals the number of taxa with non-zero abundance.
    /// Evidence: By definition, observed species count = |{i : pᵢ &gt; 0}|.
    /// Source: Wikipedia "Species richness".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_ObservedSpecies_EqualsNonZeroTaxa()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            int expected = abundances.Values.Count(v => v > 0);
            return (result.ObservedSpecies == expected)
                .Label($"ObservedSpecies={result.ObservedSpecies}, expected={expected}");
        });
    }

    #endregion

    #region META-BETA-001: S: dist(a,b)=dist(b,a); R: dist ∈ [0,1]; I: dist(x,x)=0; D: deterministic

    /// <summary>
    /// INV-1: Bray-Curtis ∈ [0, 1].
    /// Evidence: BC = 1 − 2C/(S₁+S₂) where 0 ≤ C ≤ min(S₁,S₂).
    /// Source: Bray &amp; Curtis (1957).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_BrayCurtis_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            return (result.BrayCurtis >= -1e-9 && result.BrayCurtis <= 1.0 + 1e-9)
                .Label($"BrayCurtis={result.BrayCurtis:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-2: Jaccard distance ∈ [0, 1].
    /// Evidence: J = 1 − |A∩B|/|A∪B| where |A∩B| ≤ |A∪B|.
    /// Source: Jaccard (1901).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_JaccardDistance_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            return (result.JaccardDistance >= -1e-9 && result.JaccardDistance <= 1.0 + 1e-9)
                .Label($"Jaccard={result.JaccardDistance:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-3: Symmetry — BrayCurtis(a,b) = BrayCurtis(b,a).
    /// Evidence: The formula is symmetric in its inputs.
    /// Source: Bray &amp; Curtis (1957).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_BrayCurtis_IsSymmetric()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            var r2 = MetagenomicsAnalyzer.CalculateBetaDiversity("s2", a2, "s1", a1);
            return (Math.Abs(r1.BrayCurtis - r2.BrayCurtis) < 1e-9)
                .Label($"BrayCurtis not symmetric: {r1.BrayCurtis:F6} vs {r2.BrayCurtis:F6}");
        });
    }

    /// <summary>
    /// INV-4: Symmetry — Jaccard(a,b) = Jaccard(b,a).
    /// Evidence: Jaccard distance is symmetric by definition.
    /// Source: Jaccard (1901).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_JaccardDistance_IsSymmetric()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            var r2 = MetagenomicsAnalyzer.CalculateBetaDiversity("s2", a2, "s1", a1);
            return (Math.Abs(r1.JaccardDistance - r2.JaccardDistance) < 1e-9)
                .Label($"Jaccard not symmetric: {r1.JaccardDistance:F6} vs {r2.JaccardDistance:F6}");
        });
    }

    /// <summary>
    /// INV-5: Identity — distance of a sample to itself is 0.
    /// Evidence: BC = 1 − 2C/(S+S) = 1 − 1 = 0 when both samples are identical.
    /// Jaccard = 0 when all species are shared.
    /// Source: Bray &amp; Curtis (1957); Jaccard (1901).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_SelfDistance_IsZero()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", abundances, "s2", abundances);
            return (Math.Abs(result.BrayCurtis) < 1e-9 && Math.Abs(result.JaccardDistance) < 1e-9)
                .Label($"Self-distance must be 0: BC={result.BrayCurtis:F6}, Jaccard={result.JaccardDistance:F6}");
        });
    }

    /// <summary>
    /// INV-6: Disjoint samples → maximum Jaccard distance = 1.0.
    /// Evidence: When A ∩ B = ∅, Jaccard = 1 − 0/|A∪B| = 1.
    /// Source: Jaccard (1901).
    /// </summary>
    [Test]
    [Category("Property")]
    public void BetaDiversity_DisjointSamples_JaccardIsOne()
    {
        var s1 = new Dictionary<string, double> { ["A"] = 0.5, ["B"] = 0.5 };
        var s2 = new Dictionary<string, double> { ["C"] = 0.5, ["D"] = 0.5 };

        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", s1, "s2", s2);
        Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(0.001),
            $"Disjoint samples should have Jaccard = 1.0, got {result.JaccardDistance:F4}");
    }

    /// <summary>
    /// INV-7: Beta diversity is deterministic.
    /// Evidence: CalculateBetaDiversity is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_IsDeterministic()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            var r2 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            bool same = Math.Abs(r1.BrayCurtis - r2.BrayCurtis) < 1e-10 &&
                        Math.Abs(r1.JaccardDistance - r2.JaccardDistance) < 1e-10 &&
                        r1.SharedSpecies == r2.SharedSpecies;
            return same.Label("CalculateBetaDiversity must be deterministic");
        });
    }

    /// <summary>
    /// INV-8: SharedSpecies + UniqueToSample1 + UniqueToSample2 = total unique species.
    /// Evidence: Every species across both samples is classified into exactly one category.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_SpeciesCounts_SumToTotal()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            int total = a1.Keys.Union(a2.Keys)
                .Count(k => (a1.GetValueOrDefault(k) > 0) || (a2.GetValueOrDefault(k) > 0));
            int sum = result.SharedSpecies + result.UniqueToSample1 + result.UniqueToSample2;
            return (sum == total)
                .Label($"Species count: shared+unique1+unique2={sum}, total={total}");
        });
    }

    #endregion
}
