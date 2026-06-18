// KMER-GENERATE-001 — K-mer Generation (all possible k-mers over an alphabet)
// Evidence: docs/Evidence/KMER-GENERATE-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-GENERATE-001.md
// Source: Wikipedia — K-mer (https://en.wikipedia.org/wiki/K-mer);
//         Clavijo BJ (2018), BioInfoLogics — k-mer counting, part I
//         (https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/);
//         Python Std Library — itertools.product (https://docs.python.org/3/library/itertools.html)

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for KMER-GENERATE-001: KmerAnalyzer.GenerateAllKmers.
///
/// Expected values are derived from the sources, not from the implementation:
/// the k-mer universe over an n-letter alphabet has n^k members (Wikipedia/BioInfoLogics),
/// enumerated as the k-fold Cartesian product (itertools.product); for a sorted alphabet
/// the output is lexicographic with the rightmost position advancing fastest.
/// </summary>
[TestFixture]
public class KmerAnalyzer_GenerateAllKmers_Tests
{
    // Default DNA alphabet {A,C,G,T} is already in sorted order -> lexicographic output.
    private const string Dna = "ACGT";

    // 20 standard amino acids, sorted; |alphabet| = 20 -> universe size 20^k.
    private const string Protein = "ACDEFGHIKLMNPQRSTVWY";

    #region GenerateAllKmers — MUST

    // M1 — n^k with n=4, k=1: the four DNA monomers A,C,G,T in lexicographic order.
    [Test]
    public void GenerateAllKmers_K1Dna_ReturnsFourMonomersInOrder()
    {
        var result = KmerAnalyzer.GenerateAllKmers(1).ToList();

        var expected = new[] { "A", "C", "G", "T" };
        Assert.That(result, Is.EqualTo(expected),
            "k=1 over {A,C,G,T} must be exactly the 4 monomers (4^1=4) in lexicographic order.");
    }

    // M2 — n^k with n=4, k=2: all 16 (4^2) 2-mers, lexicographic odometer order AA..TT.
    [Test]
    public void GenerateAllKmers_K2Dna_ReturnsSixteenTwoMersLexicographic()
    {
        var result = KmerAnalyzer.GenerateAllKmers(2).ToList();

        // Cartesian product {A,C,G,T} x {A,C,G,T}, rightmost position fastest (itertools.product).
        var expected = new[]
        {
            "AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT"
        };
        Assert.That(result, Is.EqualTo(expected),
            "k=2 must be all 16 (4^2) 2-mers in lexicographic order per the Cartesian-product odometer ordering.");
    }

    // M3 — k=3: count is 4^3=64; boundaries first=AAA, second=AAC, last=TTT (odometer order).
    [Test]
    public void GenerateAllKmers_K3Dna_HasSixtyFourWithCorrectBoundaries()
    {
        var result = KmerAnalyzer.GenerateAllKmers(3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(64),
                "k=3 universe size is 4^3 = 64 (Wikipedia n^k).");
            Assert.That(result[0], Is.EqualTo("AAA"),
                "First k-mer in lexicographic order is AAA.");
            Assert.That(result[1], Is.EqualTo("AAC"),
                "Second is AAC: rightmost position advances first (odometer ordering).");
            Assert.That(result[^1], Is.EqualTo("TTT"),
                "Last k-mer in lexicographic order is TTT.");
        });
    }

    // M4 — universe size equals 4^k for k = 1..6 (default DNA alphabet).
    [Test]
    public void GenerateAllKmers_DnaVariousK_CountEqualsFourToTheK()
    {
        Assert.Multiple(() =>
        {
            for (int k = 1; k <= 6; k++)
            {
                int expected = (int)Math.Pow(4, k); // 4,16,64,256,1024,4096
                Assert.That(KmerAnalyzer.GenerateAllKmers(k).Count(), Is.EqualTo(expected),
                    $"k={k}: number of all possible DNA k-mers must be 4^{k} = {expected}.");
            }
        });
    }

    // M5 — n^k generalises to any alphabet: 20 amino acids, k=2 -> 20^2 = 400.
    [Test]
    public void GenerateAllKmers_ProteinAlphabetK2_ReturnsFourHundred()
    {
        var result = KmerAnalyzer.GenerateAllKmers(2, Protein).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(400),
                "Universe size for a 20-letter alphabet at k=2 is 20^2 = 400 (n^k).");
            Assert.That(result.Distinct().Count(), Is.EqualTo(400),
                "All 400 protein 2-mers must be distinct.");
        });
    }

    // M6 — output is exactly the Cartesian-product set: no duplicates, distinct count = 4^k.
    [Test]
    public void GenerateAllKmers_K4Dna_NoDuplicatesEqualsCartesianSet()
    {
        var result = KmerAnalyzer.GenerateAllKmers(4).ToList();

        // Independent reference: build the 4-fold Cartesian product of {A,C,G,T}.
        var reference =
            from a in Dna
            from b in Dna
            from c in Dna
            from d in Dna
            select string.Concat(a, b, c, d);
        var referenceSet = reference.ToHashSet();

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(256),
                "k=4 universe size is 4^4 = 256.");
            Assert.That(result.Distinct().Count(), Is.EqualTo(256),
                "Output must contain no duplicate k-mers (each is a unique length-k tuple).");
            Assert.That(result.ToHashSet(), Is.EquivalentTo(referenceSet),
                "Output set must equal the independently built 4-fold Cartesian product of {A,C,G,T}.");
        });
    }

    #endregion

    #region GenerateAllKmers — SHOULD (edge cases)

    // S1 — single-letter alphabet: 1^k = 1, the homopolymer only.
    [Test]
    public void GenerateAllKmers_SingleLetterAlphabet_ReturnsOnlyHomopolymer()
    {
        var result = KmerAnalyzer.GenerateAllKmers(4, "A").ToList();

        Assert.That(result, Is.EqualTo(new[] { "AAAA" }),
            "A 1-letter alphabet yields exactly one k-mer (1^4 = 1): the homopolymer AAAA.");
    }

    // S2 — INV-03: every k-mer has length exactly k and uses only alphabet characters.
    [Test]
    public void GenerateAllKmers_K2Dna_EveryKmerHasLengthKFromAlphabet()
    {
        var alphabetSet = Dna.ToHashSet();

        var result = KmerAnalyzer.GenerateAllKmers(2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result.All(km => km.Length == 2), Is.True,
                "Every generated k-mer must have length exactly k=2.");
            Assert.That(result.All(km => km.All(alphabetSet.Contains)), Is.True,
                "Every character of every k-mer must come from the alphabet {A,C,G,T}.");
        });
    }

    // S3 — k must be positive: k=0 and k<0 throw ArgumentOutOfRangeException.
    [Test]
    public void GenerateAllKmers_NonPositiveK_ThrowsArgumentOutOfRange()
    {
        Assert.Multiple(() =>
        {
            // Enumerate to trigger the validation (deferred-execution method).
            Assert.Throws<ArgumentOutOfRangeException>(
                () => KmerAnalyzer.GenerateAllKmers(0).ToList(),
                "k=0 is not a valid k-mer length.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => KmerAnalyzer.GenerateAllKmers(-1).ToList(),
                "Negative k is not a valid k-mer length.");
        });
    }

    // S4 — null/empty alphabet: no symbols means no k-mers can be formed -> ArgumentException.
    [Test]
    public void GenerateAllKmers_NullOrEmptyAlphabet_ThrowsArgumentException()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => KmerAnalyzer.GenerateAllKmers(2, "").ToList(),
                "An empty alphabet cannot form any k-mer.");
            Assert.Throws<ArgumentException>(
                () => KmerAnalyzer.GenerateAllKmers(2, null!).ToList(),
                "A null alphabet cannot form any k-mer.");
        });
    }

    #endregion

    #region GenerateAllKmers — COULD

    // C1 — ordering follows the supplied alphabet, not a forced sort: "TGCA", k=1 -> T,G,C,A.
    [Test]
    public void GenerateAllKmers_UnsortedAlphabetK1_FollowsAlphabetOrder()
    {
        var result = KmerAnalyzer.GenerateAllKmers(1, "TGCA").ToList();

        Assert.That(result, Is.EqualTo(new[] { "T", "G", "C", "A" }),
            "Output order follows the alphabet's own order; lexicographic order holds only for a sorted alphabet.");
    }

    #endregion
}
