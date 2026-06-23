using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Chromosome;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 4): exact-value pins for the internal arithmetic of the
/// comparison / completeness / repeat / suspicious-region methods (shared-kmer fractions, fragmented vs
/// complete classification, contig length boundaries, k-mer size/peak edges, imperfect tandem purity,
/// low-complexity score formula).
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers4_Tests
{
    private const double Tol = 1e-9;

    private static string RandomDna(int length, int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";
        return string.Concat(Enumerable.Range(0, length).Select(_ => bases[rng.Next(4)]));
    }

    [Test]
    public void CompareAssemblies_PartialOverlap_ExactFractions()
    {
        // a1 4-mers {ACGT,CGTA,GTAC,TACG}; a2 4-mers {ACGT,CGTA,GTAA,TAAA,AAAA}; shared = 2.
        var a1 = new[] { ("s", "ACGTACGT") };
        var a2 = new[] { ("s", "ACGTAAAA") };
        var cmp = CompareAssemblies(a1, a2, kmerSize: 4);
        Assert.That(cmp.AlignedFraction1, Is.EqualTo(2.0 / 4).Within(Tol)); // 0.5
        Assert.That(cmp.AlignedFraction2, Is.EqualTo(2.0 / 5).Within(Tol)); // 0.4
        Assert.That(cmp.SequenceIdentity, Is.EqualTo((0.5 + 0.4) / 2).Within(Tol)); // 0.45
    }

    [Test]
    public void AssessCompleteness_PartialCoverage_IsFragmented()
    {
        // 100-bp unique marker; assembly carries only its first 80 bp ⇒ ~71% k-mer coverage:
        // above the relaxed hit threshold but below the 0.9 completeness cutoff ⇒ Fragmented.
        string marker = RandomDna(100, 11);
        var r = AssessCompleteness(new[] { ("scaf", marker.Substring(0, 80)) }, new[] { ("g", marker) });
        Assert.That(r.Fragmented, Is.EqualTo(1));
        Assert.That(r.Complete, Is.EqualTo(0));
        Assert.That(r.Missing, Is.EqualTo(0));
    }

    [Test]
    public void ExtractContigs_LengthBoundaryIsInclusive()
    {
        // A 3-bp contig is kept at minContigLength 3 (inclusive '>=') but dropped at 4.
        Assert.That(ExtractContigs(new[] { ("s", "AAAN") }, minContigLength: 3).Single(),
            Is.EqualTo(("s_contig1", "AAA")));
        Assert.That(ExtractContigs(new[] { ("s", "AAAN") }, minContigLength: 4).ToList(), Is.Empty);
        // Trailing contig after a gap also honours the inclusive bound.
        Assert.That(ExtractContigs(new[] { ("s", "NAAA") }, minContigLength: 3).Single(),
            Is.EqualTo(("s_contig1", "AAA")));
    }

    [Test]
    public void EstimateCompletenessFromKmers_AllSingletons_AndExpectedCoverage()
    {
        // No non-singleton k-mers ⇒ (0, 1.0, 0).
        var allSingleton = EstimateCompletenessFromKmers(new[] { ("A", 1), ("B", 1) });
        Assert.That(allSingleton.Completeness, Is.EqualTo(0.0).Within(Tol));
        Assert.That(allSingleton.ErrorRate, Is.EqualTo(1.0).Within(Tol));
        Assert.That(allSingleton.EstimatedGenomeSize, Is.EqualTo(0L));

        // expectedCoverage overrides the computed peak: size = total 11 / 2 = 5.
        var forced = EstimateCompletenessFromKmers(
            new[] { ("AAA", 4), ("CCC", 4), ("GGG", 2), ("TTT", 1) }, expectedCoverage: 2);
        Assert.That(forced.EstimatedGenomeSize, Is.EqualTo(5L));
    }

    [Test]
    public void FindTandemRepeats_ImperfectCopy_ExactPurity()
    {
        // AAAAA·AAAAA·AAAAC (then a trailing G so length 16 > minUnit·minCopies = 15): the 3rd
        // copy matches 4/5 = 80% ⇒ counted; matches = 5 + 4 = 9.
        // purity = (unitLen + matches)/(copies·unitLen) = (5 + 9)/15 = 14/15.
        var r = FindTandemRepeats(new[] { ("s", "AAAAAAAAAAAAAACG") }, minUnitLength: 5, minCopies: 3).First();
        Assert.That(r.Unit, Is.EqualTo("AAAAA"));
        Assert.That(r.Copies, Is.EqualTo(3));
        Assert.That(r.Purity, Is.EqualTo(14.0 / 15.0).Within(Tol));
    }

    [Test]
    public void FindSuspiciousRegions_PolyA_ExactLowComplexityScore()
    {
        // 600-bp poly(A): a 500-bp window has linguistic complexity 1/256, giving the peak
        // low-complexity score 1 − (1/256)/0.3.
        var r = FindSuspiciousRegions(new[] { ("s", new string('A', 600)) }).Single();
        Assert.That(r.Score, Is.EqualTo(1.0 - (1.0 / 256.0) / 0.3).Within(1e-9));
    }
}
