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
}
