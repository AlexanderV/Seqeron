using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Chromosome;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 3): exact-value tests for assembly comparison,
/// syntenic-block clustering, suspicious-region boundaries, duplicated-marker completeness and
/// imperfect tandem repeats — the remaining uncovered branches.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers3_Tests
{
    private const double Tol = 1e-9;

    // Deterministic non-repetitive 120-bp DNA (fixed seed ⇒ unique 21-mers).
    private static string RandomDna(int length, int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";
        return string.Concat(Enumerable.Range(0, length).Select(_ => bases[rng.Next(4)]));
    }

    #region CompareAssemblies

    [Test]
    public void CompareAssemblies_IdenticalAndDisjoint()
    {
        var a = new[] { ("s", "ACGTACGTACGT") };
        // Identical ⇒ every k-mer shared ⇒ aligned fractions and identity all 1.0.
        var same = CompareAssemblies(a, a, kmerSize: 4);
        Assert.That(same.AlignedFraction1, Is.EqualTo(1.0).Within(Tol));
        Assert.That(same.AlignedFraction2, Is.EqualTo(1.0).Within(Tol));
        Assert.That(same.SequenceIdentity, Is.EqualTo(1.0).Within(Tol));

        // Disjoint alphabets ⇒ no shared k-mers ⇒ identity 0.
        var disjoint = CompareAssemblies(a, new[] { ("s", "GGGGGGGGGGGG") }, kmerSize: 4);
        Assert.That(disjoint.SequenceIdentity, Is.EqualTo(0.0).Within(Tol));
    }

    #endregion

    #region FindSyntenicBlocks

    [Test]
    public void FindSyntenicBlocks_CollinearBlockBetweenIdenticalSequences()
    {
        // Two identical 100-bp sequences with unique 21-mers: anchors at positions 0,21,42,63
        // cluster into a single collinear (non-inverted) block spanning [0, 63] on both.
        string seq = RandomDna(100, 7);
        var blocks = FindSyntenicBlocks(
            new[] { ("s1", seq) }, new[] { ("s2", seq) }, minBlockSize: 50, kmerSize: 21).ToList();

        Assert.That(blocks, Has.Count.EqualTo(1));
        var b = blocks[0];
        Assert.That(b.Seq1, Is.EqualTo("s1"));
        Assert.That(b.Start1, Is.EqualTo(0));
        Assert.That(b.End1, Is.EqualTo(63));   // last anchor at i = 63 (step 21)
        Assert.That(b.Seq2, Is.EqualTo("s2"));
        Assert.That(b.Start2, Is.EqualTo(0));
        Assert.That(b.End2, Is.EqualTo(63));
        Assert.That(b.IsInverted, Is.False);
    }

    #endregion

    #region FindSuspiciousRegions exact boundary

    [Test]
    public void FindSuspiciousRegions_PolyA_ExactRegionAndReason()
    {
        // 600-bp poly(A): every window is low-complexity ⇒ a single region spanning the whole
        // sequence labelled "Low complexity".
        var r = FindSuspiciousRegions(new[] { ("s", new string('A', 600)) }).Single();
        Assert.That(r.SequenceId, Is.EqualTo("s"));
        Assert.That(r.Start, Is.EqualTo(0));
        Assert.That(r.End, Is.EqualTo(599));      // final-region branch ⇒ sequence.Length - 1
        Assert.That(r.Reason, Is.EqualTo("Low complexity"));
    }

    #endregion

    #region AssessCompleteness duplicated marker

    [Test]
    public void AssessCompleteness_DuplicatedMarker_CountsBothCopies()
    {
        // A marker that appears in TWO assembly sequences ⇒ complete & duplicated (not single-copy).
        string marker = new string('A', 40) + new string('C', 40);
        var r = AssessCompleteness(
            new[] { ("scafA", marker), ("scafB", marker) },
            new[] { ("g1", marker) });

        Assert.That(r.Complete, Is.EqualTo(1));
        Assert.That(r.CompleteDuplicated, Is.EqualTo(1));
        Assert.That(r.CompleteSingleCopy, Is.EqualTo(0));
        Assert.That(r.DuplicationPercent, Is.EqualTo(100.0).Within(Tol));
    }

    #endregion

    #region FindTandemRepeats imperfect

    [Test]
    public void FindTandemRepeats_ImperfectCopyStillCountedAt80Percent()
    {
        // (ATG)×3 then a 4th copy "ATC" (2/3 = 66% < 80% ⇒ not extended): exactly 3 copies.
        var r = FindTandemRepeats(new[] { ("s", "ATGATGATGATC") }, minUnitLength: 3, minCopies: 3).First();
        Assert.That(r.Unit, Is.EqualTo("ATG"));
        Assert.That(r.Copies, Is.EqualTo(3));
        Assert.That(r.Start, Is.EqualTo(0));
        Assert.That(r.End, Is.EqualTo(8)); // 3 copies × 3 = 9 bp ⇒ [0,8]
        Assert.That(r.Purity, Is.EqualTo(1.0).Within(Tol)); // the 3 counted copies are perfect
    }

    #endregion
}
