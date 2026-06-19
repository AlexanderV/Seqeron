using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Population;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Population Genetics area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-FREQ-001 — allele frequency calculation (PopGen).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 43.
///
/// APIs under test:
///   CalculateAlleleFrequencies(homMajor, het, homMinor) → (MajorFreq, MinorFreq), where
///     MajorFreq = (2·homMajor + het) / 2(homMajor + het + homMinor), MinorFreq symmetric;
///     (0,0) when there are no individuals.
///   CalculateMAF(genotypes) where 0/1/2 = hom-ref / het / hom-alt: alt allele frequency
///     folded to the minor side, min(altFreq, 1 − altFreq).
///
/// Relations (derived from these ratio definitions, NOT from output):
///   • COMP (frequencies sum to 1): major and minor allele counts partition the 2N gene
///          copies, so MajorFreq + MinorFreq = 1 whenever the sample is non-empty.
///   • INV (count scaling): a frequency is a ratio, so multiplying every genotype count by
///          the same factor leaves both frequencies unchanged.
///   • INV (sample reordering): MAF depends only on the allele-count sum, so permuting the
///          genotype list does not change it; duplicating the list (a special scaling) also
///          leaves it unchanged.
///   • COMP (cross-API consistency): mapping hom-major/het/hom-minor to genotypes 0/1/2 makes
///          CalculateMAF equal min(MajorFreq, MinorFreq) from CalculateAlleleFrequencies.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class PopulationGeneticsMetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(20260619);

    /// <summary>A genotype-count triple (homMajor, het, homMinor) with at least one individual.</summary>
    private static IEnumerable<(int Major, int Het, int Minor)> CountTriples()
    {
        yield return (10, 5, 2);
        yield return (1, 0, 0);
        yield return (0, 7, 0);
        yield return (3, 3, 3);
        yield return (50, 0, 50);
        for (int t = 0; t < 6; t++)
            yield return (Rng.Next(0, 40), Rng.Next(0, 40), Rng.Next(1, 40));
    }

    private static List<int> GenotypesFromCounts(int major, int het, int minor) =>
        Enumerable.Repeat(0, major)
            .Concat(Enumerable.Repeat(1, het))
            .Concat(Enumerable.Repeat(2, minor))
            .ToList();

    private static List<int> Shuffle(IEnumerable<int> items)
    {
        var list = items.ToList();
        for (int i = list.Count - 1; i > 0; i--)
        {
            int j = Rng.Next(i + 1);
            (list[i], list[j]) = (list[j], list[i]);
        }
        return list;
    }

    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    private static char Flip(char c) => c == 'A' ? 'C' : 'A';

    /// <summary>Returns a copy of <paramref name="seq"/> with exactly <paramref name="count"/> spread-out positions mutated.</summary>
    private static string MutatePositions(string seq, int count)
    {
        var chars = seq.ToCharArray();
        for (int n = 0; n < count; n++)
        {
            int pos = (int)((long)n * seq.Length / count);
            chars[pos] = Flip(chars[pos]);
        }
        return new string(chars);
    }

    private static List<char[]> Panel(params string[] sequences) =>
        sequences.Select(s => s.ToCharArray()).ToList();

    #endregion

    #region COMP — major and minor frequencies sum to 1

    [Test]
    [Description("COMP: for a non-empty sample MajorFreq + MinorFreq = 1, since the two allele counts partition the 2N gene copies.")]
    public void AlleleFrequencies_SumToOne()
    {
        foreach (var (major, het, minor) in CountTriples())
        {
            var (maj, min) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(major, het, minor);

            (maj + min).Should().BeApproximately(1.0, 1e-12,
                because: "every gene copy is counted as either the major or the minor allele, so the two frequencies partition probability mass 1");
        }
    }

    #endregion

    #region INV — scaling all counts preserves the frequencies

    [Test]
    [Description("INV: multiplying every genotype count by the same factor leaves both allele frequencies unchanged, because each is a ratio.")]
    public void AlleleFrequencies_ScalingCounts_IsInvariant()
    {
        foreach (var (major, het, minor) in CountTriples())
        {
            var baseline = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(major, het, minor);

            foreach (int k in new[] { 2, 3, 5 })
            {
                var scaled = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(k * major, k * het, k * minor);

                scaled.MajorFreq.Should().BeApproximately(baseline.MajorFreq, 1e-12,
                    because: $"scaling counts by {k} multiplies numerator and denominator equally, leaving the major frequency unchanged");
                scaled.MinorFreq.Should().BeApproximately(baseline.MinorFreq, 1e-12,
                    because: $"scaling counts by {k} multiplies numerator and denominator equally, leaving the minor frequency unchanged");
            }
        }
    }

    #endregion

    #region INV — MAF is invariant to sample order and to duplication

    [Test]
    [Description("INV: permuting or duplicating the genotype list does not change the MAF, which depends only on the allele-count sum.")]
    public void Maf_SampleReorderingAndDuplication_IsInvariant()
    {
        foreach (var (major, het, minor) in CountTriples())
        {
            var genotypes = GenotypesFromCounts(major, het, minor);
            double baseline = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

            for (int trial = 0; trial < 3; trial++)
            {
                double shuffled = PopulationGeneticsAnalyzer.CalculateMAF(Shuffle(genotypes));
                shuffled.Should().BeApproximately(baseline, 1e-12,
                    because: "MAF sums alternate-allele dosages over a fixed denominator, so sample order is irrelevant");
            }

            double doubled = PopulationGeneticsAnalyzer.CalculateMAF(genotypes.Concat(genotypes));
            doubled.Should().BeApproximately(baseline, 1e-12,
                because: "duplicating every sample scales both the allele sum and the denominator equally, leaving the folded frequency unchanged");
        }
    }

    #endregion

    #region COMP — MAF equals min(major, minor) from the allele-frequency API

    [Test]
    [Description("COMP: mapping hom-major/het/hom-minor to genotypes 0/1/2 makes CalculateMAF agree with min(MajorFreq, MinorFreq).")]
    public void Maf_EqualsMinorAlleleFrequency()
    {
        foreach (var (major, het, minor) in CountTriples())
        {
            var (majFreq, minFreq) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(major, het, minor);
            double maf = PopulationGeneticsAnalyzer.CalculateMAF(GenotypesFromCounts(major, het, minor));

            maf.Should().BeApproximately(Math.Min(majFreq, minFreq), 1e-12,
                because: "the minor allele frequency is by definition the smaller of the two allele frequencies");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: POP-DIV-001 — diversity statistics: nucleotide diversity π and Watterson θ (PopGen).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 44.
    //
    // APIs under test:
    //   CalculateNucleotideDiversity(sequences) = mean over all C(n,2) sequence pairs of the
    //     per-site Hamming difference = totalPairwiseDiffs / (comparisons · length).
    //   CalculateWattersonTheta(S, n, L) = S / (a₁ · L), with a₁ = Σ_{i=1}^{n−1} 1/i.
    //
    // Relations (derived from these formulas, NOT from output):
    //   • MON (more diverse ⇒ higher π): adding substitutions raises the pairwise difference
    //          total over a fixed denominator, so π increases. For a two-sequence panel
    //          differing at m sites, π = m/L exactly and is strictly increasing in m; a
    //          homogeneous panel has π = 0 while any polymorphism gives π > 0.
    //   • MON (more segregating sites ⇒ higher θ): θ is linear in S with a positive slope
    //          1/(a₁·L), so it is strictly increasing in S and 0 at S = 0.
    //   • INV (sample reordering): π sums over UNORDERED pairs, so permuting the sequence panel
    //          leaves it unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region MON — more substitutions strictly increase nucleotide diversity π

    [Test]
    [Description("MON: for a two-sequence panel differing at m sites, π = m/L and strictly increases with m; a homogeneous panel has π = 0.")]
    public void NucleotideDiversity_MoreDivergence_IncreasesPi()
    {
        const int len = 50;
        string baseSeq = RandomDna(len);

        PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(Panel(baseSeq, baseSeq))
            .Should().Be(0.0, because: "two identical sequences differ nowhere, so π = 0");

        double previous = double.NegativeInfinity;
        for (int m = 0; m <= 8; m++)
        {
            string variant = MutatePositions(baseSeq, m);
            double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(Panel(baseSeq, variant));

            pi.Should().BeApproximately((double)m / len, 1e-12,
                because: $"the single sequence pair differs at exactly {m} of {len} sites, so π = {m}/{len}");
            pi.Should().BeGreaterThan(previous,
                because: "each extra substitution adds a pairwise difference over a fixed denominator, so π strictly increases");
            previous = pi;
        }
    }

    [Test]
    [Description("MON: a panel with a polymorphism has strictly greater π than the otherwise-identical homogeneous panel.")]
    public void NucleotideDiversity_AddingPolymorphism_IncreasesPi()
    {
        const int len = 40;
        string s = RandomDna(len);

        double homogeneous = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(Panel(s, s, s, s));
        double withMutant = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(Panel(s, s, s, MutatePositions(s, 3)));

        homogeneous.Should().Be(0.0, because: "four identical sequences are monomorphic");
        withMutant.Should().BeGreaterThan(homogeneous,
            because: "replacing one sequence with a 3-substitution variant introduces pairwise differences, raising π above 0");
    }

    #endregion

    #region MON — more segregating sites strictly increase Watterson's θ

    [Test]
    [Description("MON: θ = S/(a₁·L) is linear in S with positive slope, so it is 0 at S = 0 and strictly increases with the number of segregating sites.")]
    public void WattersonTheta_MoreSegregatingSites_IncreasesTheta()
    {
        foreach (int n in new[] { 2, 5, 10 })
        {
            foreach (int length in new[] { 100, 500 })
            {
                double a1 = Enumerable.Range(1, n - 1).Sum(i => 1.0 / i);
                double previous = double.NegativeInfinity;

                for (int s = 0; s <= 10; s++)
                {
                    double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(s, n, length);

                    theta.Should().BeApproximately(s / (a1 * length), 1e-12,
                        because: $"θ = S/(a₁·L) with S={s}, a₁={a1:F4}, L={length}");
                    theta.Should().BeGreaterThan(previous,
                        because: "θ is linear in S with the positive slope 1/(a₁·L), so it strictly increases with segregating sites");
                    previous = theta;
                }
            }
        }
    }

    #endregion

    #region INV — permuting the sequence panel preserves π

    [Test]
    [Description("INV: π sums over unordered sequence pairs, so permuting the panel leaves it unchanged.")]
    public void NucleotideDiversity_SampleReordering_IsInvariant()
    {
        const int len = 60;
        string s = RandomDna(len);
        var sequences = new[]
        {
            s,
            MutatePositions(s, 2),
            MutatePositions(s, 5),
            MutatePositions(s, 9),
            MutatePositions(s, 1),
        };

        double baseline = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(Panel(sequences));

        for (int trial = 0; trial < 5; trial++)
        {
            var order = Enumerable.Range(0, sequences.Length).ToList();
            for (int i = order.Count - 1; i > 0; i--)
            {
                int j = Rng.Next(i + 1);
                (order[i], order[j]) = (order[j], order[i]);
            }
            var permuted = Panel(order.Select(idx => sequences[idx]).ToArray());

            PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(permuted)
                .Should().BeApproximately(baseline, 1e-12,
                    because: "the mean over unordered pairs does not depend on the order of the sequences");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: POP-HW-001 — Hardy–Weinberg equilibrium test (PopGen).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 45.
    //
    // API under test (PopulationGeneticsAnalyzer.TestHardyWeinberg):
    //   From observed genotype counts (AA, Aa, aa) it derives p = (2·AA + Aa)/2n, q = 1−p,
    //   expected counts (p²n, 2pqn, q²n) and the 1-df χ² = Σ (obs−exp)²/exp.
    //
    // Relations (derived from these formulas, NOT from output):
    //   • COMP (genotype proportions sum to 1): p² + 2pq + q² = (p+q)² = 1, so the expected
    //          counts sum to exactly n.
    //   • INV (sample-size scaling): scaling every genotype count by k leaves p (hence the
    //          expected PROPORTIONS) unchanged and scales χ² by exactly k — both obs and exp
    //          scale by k, so each (obs−exp)²/exp term scales by k.
    //   • MON (larger deviation ⇒ larger χ²): at a FIXED allele frequency (p = ½, fixed n),
    //          moving the genotype configuration further from the HWE expectation strictly
    //          increases χ²; the perfectly-HWE configuration gives χ² = 0.
    // ───────────────────────────────────────────────────────────────────────────

    #region COMP — expected genotype counts sum to n (p² + 2pq + q² = 1)

    [Test]
    [Description("COMP: the HWE expected genotype counts sum to n, because p² + 2pq + q² = (p+q)² = 1.")]
    public void HardyWeinberg_ExpectedCounts_SumToN()
    {
        foreach (var (aa, het, bb) in CountTriples())
        {
            int n = aa + het + bb;
            var hw = PopulationGeneticsAnalyzer.TestHardyWeinberg("v", aa, het, bb);

            (hw.ExpectedAA + hw.ExpectedAa + hw.Expectedaa).Should().BeApproximately(n, 1e-9,
                because: "the expected proportions p², 2pq, q² sum to (p+q)² = 1, so the expected counts sum to the sample size n");
        }
    }

    #endregion

    #region INV — scaling all counts preserves proportions and scales χ² by the factor

    [Test]
    [Description("INV: scaling every genotype count by k leaves the expected proportions unchanged and scales χ² by exactly k.")]
    public void HardyWeinberg_ScalingCounts_PreservesProportionsAndScalesChiSquare()
    {
        foreach (var (aa, het, bb) in CountTriples())
        {
            int n = aa + het + bb;
            var baseHw = PopulationGeneticsAnalyzer.TestHardyWeinberg("v", aa, het, bb);

            foreach (int k in new[] { 2, 3, 5 })
            {
                int nk = k * n;
                var scaled = PopulationGeneticsAnalyzer.TestHardyWeinberg("v", k * aa, k * het, k * bb);

                (scaled.ExpectedAA / nk).Should().BeApproximately(baseHw.ExpectedAA / n, 1e-12,
                    because: $"scaling by {k} leaves the allele frequency p unchanged, so the expected AA proportion p² is unchanged");
                (scaled.Expectedaa / nk).Should().BeApproximately(baseHw.Expectedaa / n, 1e-12,
                    because: "the expected aa proportion q² depends only on the (scale-invariant) allele frequency");

                scaled.ChiSquare.Should().BeApproximately(k * baseHw.ChiSquare,
                    Math.Abs(k * baseHw.ChiSquare) * 1e-9 + 1e-9,
                    because: $"each (obs−exp)²/exp term scales by {k} when both obs and exp scale by {k}");
            }
        }
    }

    #endregion

    #region MON — larger deviation from HWE gives a larger χ²

    [Test]
    [Description("MON: at a fixed allele frequency (p = ½), moving the genotype configuration further from the HWE expectation strictly increases χ²; the HWE configuration gives χ² = 0.")]
    public void HardyWeinberg_LargerDeviation_IncreasesChiSquare()
    {
        const int n = 1000;   // p fixed at 1/2 by setting aa = AA and Aa = n − 2·AA

        double previous = double.NegativeInfinity;

        // AA descending from the HWE value 250 ⇒ deviation |AA−250| ascending.
        foreach (int homAA in new[] { 250, 200, 150, 100, 50, 0 })
        {
            int het = n - 2 * homAA;   // keeps 2·AA + Aa = n ⇒ p = 1/2
            int homaa = homAA;         // keeps the sample size at n

            var hw = PopulationGeneticsAnalyzer.TestHardyWeinberg("v", homAA, het, homaa);

            // p must stay exactly 1/2 ⇒ expected counts are the HWE-balanced 250/500/250.
            hw.ExpectedAa.Should().BeApproximately(0.5 * n, 1e-9,
                because: "the construction holds the allele frequency at p = 1/2, so the expected heterozygote count is 2pq·n = n/2");

            if (homAA == 250)
                hw.ChiSquare.Should().BeApproximately(0.0, 1e-9,
                    because: "the genotype counts equal the HWE expectation, so there is no deviation");

            hw.ChiSquare.Should().BeGreaterThan(previous,
                because: "with the allele frequency fixed, a genotype configuration further from the HWE expectation has a strictly larger χ²");
            previous = hw.ChiSquare;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: POP-FST-001 — Wright's Fst between populations (PopGen).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 46.
    //
    // API under test (PopulationGeneticsAnalyzer.CalculateFst):
    //   Per locus, p̄ = (n1·p1 + n2·p2)/(n1+n2), the among-population variance σ²_S and the
    //   total heterozygosity p̄(1−p̄); Fst = Σ σ²_S / Σ p̄(1−p̄).
    //
    // Relations (derived from the formula, NOT from output):
    //   • SYM (symmetry): p̄, the variance and the heterozygosity are all symmetric in the two
    //          (freq, size) populations, so Fst(A,B) = Fst(B,A).
    //   • COMP (identical ⇒ 0): if both populations have identical per-locus frequencies then
    //          p1 = p2 = p̄ at every locus, the variance is 0, and Fst = 0.
    //   • MON (more differentiated ⇒ higher Fst): for a single locus with symmetric frequencies
    //          ½∓δ/2 and equal sample sizes, p̄ = ½ and p̄(1−p̄) = ¼ are fixed while σ²_S = δ²/4,
    //          so Fst = δ² — strictly increasing in the differentiation δ.
    // ───────────────────────────────────────────────────────────────────────────

    #region Helpers (Fst)

    private static List<(double AlleleFreq, int SampleSize)> RandomPopulation(int loci)
    {
        var pop = new List<(double, int)>(loci);
        for (int i = 0; i < loci; i++)
            pop.Add((0.05 + 0.9 * Rng.NextDouble(), Rng.Next(10, 100)));   // freq kept in (0,1) ⇒ het > 0
        return pop;
    }

    #endregion

    #region SYM — Fst(A,B) = Fst(B,A)

    [Test]
    [Description("SYM: Fst is symmetric, since p̄, the among-population variance and the heterozygosity are all symmetric in the two populations.")]
    public void Fst_IsSymmetric()
    {
        for (int trial = 0; trial < 20; trial++)
        {
            int loci = Rng.Next(1, 6);
            var a = RandomPopulation(loci);
            var b = RandomPopulation(loci);

            double ab = PopulationGeneticsAnalyzer.CalculateFst(a, b);
            double ba = PopulationGeneticsAnalyzer.CalculateFst(b, a);

            ba.Should().BeApproximately(ab, 1e-12,
                because: "swapping the two populations leaves every per-locus term of the Fst ratio unchanged");
        }
    }

    #endregion

    #region COMP — identical populations give Fst = 0

    [Test]
    [Description("COMP: two populations with identical per-locus allele frequencies have zero among-population variance, so Fst = 0.")]
    public void Fst_IdenticalPopulations_IsZero()
    {
        for (int trial = 0; trial < 20; trial++)
        {
            var pop = RandomPopulation(Rng.Next(1, 6));

            PopulationGeneticsAnalyzer.CalculateFst(pop, pop).Should().BeApproximately(0.0, 1e-12,
                because: "when p1 = p2 = p̄ at every locus the among-population variance vanishes, so Fst = 0");
        }
    }

    #endregion

    #region MON — greater differentiation gives a strictly larger Fst

    [Test]
    [Description("MON: for symmetric single-locus frequencies ½∓δ/2 with equal sample sizes, Fst = δ², strictly increasing in the differentiation δ.")]
    public void Fst_MoreDifferentiation_IncreasesFst()
    {
        const int n = 50;
        double previous = double.NegativeInfinity;

        foreach (double delta in new[] { 0.0, 0.1, 0.2, 0.4, 0.6, 0.8 })
        {
            var pop1 = new List<(double, int)> { (0.5 - delta / 2, n) };
            var pop2 = new List<(double, int)> { (0.5 + delta / 2, n) };

            double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

            fst.Should().BeApproximately(delta * delta, 1e-12,
                because: $"with p̄ = ½ and p̄(1−p̄) = ¼ held fixed, σ²_S = δ²/4 gives Fst = δ² = {delta * delta}");
            fst.Should().BeGreaterThan(previous,
                because: "a larger allele-frequency difference between the populations gives a strictly larger Fst");
            previous = fst;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: POP-LD-001 — linkage disequilibrium D' and r² (PopGen).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 47.
    //
    // API under test (PopulationGeneticsAnalyzer.CalculateLD):
    //   From diploid genotype pairs (0/1/2 dosage at each of two loci): r² = squared Pearson
    //   correlation of the dosages; D = Cov(X₁,X₂)/2, normalised to D' = |D|/D_max ∈ [0,1].
    //
    // Relations (derived from the formula, NOT from output):
    //   • SYM (locus symmetry): covariance and correlation are symmetric, and D_max is a min
    //          over symmetric products, so swapping the two loci leaves D' and r² unchanged.
    //   • COMP (independent loci ⇒ D' ≈ 0): a full 3×3 factorial of genotype pairs makes the
    //          two dosages statistically independent, so Cov = 0 exactly ⇒ D' = 0 and r² = 0.
    //   • MON/boundary (perfect linkage ⇒ D' = 1): HWE-proportioned, perfectly co-inherited
    //          genotypes give the maximal D' = 1 (and r² = 1) — the upper bound of the measure.
    // ───────────────────────────────────────────────────────────────────────────

    #region SYM — LD is symmetric in the two loci

    [Test]
    [Description("SYM: swapping the two loci leaves D' and r² unchanged, as covariance, correlation and D_max are all symmetric.")]
    public void Ld_IsSymmetric()
    {
        for (int trial = 0; trial < 20; trial++)
        {
            int n = Rng.Next(5, 30);
            var pairs = Enumerable.Range(0, n).Select(_ => (Rng.Next(0, 3), Rng.Next(0, 3))).ToList();
            var swapped = pairs.Select(p => (p.Item2, p.Item1)).ToList();

            var ab = PopulationGeneticsAnalyzer.CalculateLD("a", "b", pairs, 100);
            var ba = PopulationGeneticsAnalyzer.CalculateLD("b", "a", swapped, 100);

            ba.DPrime.Should().BeApproximately(ab.DPrime, 1e-12,
                because: "D' is symmetric: D_max is a minimum over symmetric allele-frequency products");
            ba.RSquared.Should().BeApproximately(ab.RSquared, 1e-12,
                because: "r² is the squared correlation, which is symmetric in its two variables");
        }
    }

    #endregion

    #region COMP — independent loci give D' = 0 and r² = 0

    [Test]
    [Description("COMP: a full 3×3 factorial of genotype pairs makes the two dosages independent, so D' = 0 and r² = 0 exactly.")]
    public void Ld_IndependentLoci_AreZero()
    {
        var factorial = new List<(int, int)>();
        for (int g1 = 0; g1 <= 2; g1++)
            for (int g2 = 0; g2 <= 2; g2++)
                factorial.Add((g1, g2));

        var ld = PopulationGeneticsAnalyzer.CalculateLD("a", "b", factorial, 100);

        ld.DPrime.Should().BeApproximately(0.0, 1e-12,
            because: "in a full factorial the two dosages are independent, so their covariance — and hence D — is exactly 0");
        ld.RSquared.Should().BeApproximately(0.0, 1e-12,
            because: "zero covariance gives zero correlation, so r² = 0");
    }

    #endregion

    #region MON/boundary — perfectly linked loci reach D' = 1

    [Test]
    [Description("MON/boundary: HWE-proportioned, perfectly co-inherited genotypes attain the maximal D' = 1 and r² = 1.")]
    public void Ld_PerfectLinkage_ReachesOne()
    {
        // HWE proportions at p = 1/2 (0.25/0.5/0.25) ⇒ dosage variance 2pq = 1/2 = D_max scaling,
        // with the second locus identical to the first (perfect co-inheritance).
        var perfect = new List<(int, int)> { (0, 0), (1, 1), (1, 1), (2, 2) };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("a", "b", perfect, 100);

        ld.DPrime.Should().BeApproximately(1.0, 1e-9,
            because: "perfectly co-inherited HWE genotypes give D = D_max, so D' attains its maximum of 1");
        ld.RSquared.Should().BeApproximately(1.0, 1e-9,
            because: "identical dosages are perfectly correlated, so r² = 1");
    }

    #endregion
}
