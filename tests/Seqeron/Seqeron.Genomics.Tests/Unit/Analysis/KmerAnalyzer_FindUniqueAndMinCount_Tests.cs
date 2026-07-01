// KMER-UNIQUE-001 — Unique K-mers / K-mers with Minimum Count
// Evidence: docs/Evidence/KMER-UNIQUE-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-UNIQUE-001.md
// Source: BioInfoLogics — k-mer counting, part I (2018)
//         (https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/);
//         Wikipedia — K-mer (https://en.wikipedia.org/wiki/K-mer);
//         Compeau P, Pevzner P (2015). Bioinformatics Algorithms (2nd ed.).

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Tests for KMER-UNIQUE-001: Unique K-mers and K-mers with Minimum Count.
///
/// Canonical methods: KmerAnalyzer.FindUniqueKmers (occurrence count = 1) and
/// KmerAnalyzer.FindKmersWithMinCount (occurrence count >= minCount, ordered desc).
/// Expected values are derived from the published worked tables in the Evidence
/// (ATCGATCAC k=3 -> 5 unique; AGAT k=2 all distinct; ACGTACGT k=4 occurrences),
/// not from the implementation output.
/// </summary>
[TestFixture]
public class KmerAnalyzer_FindUniqueAndMinCount_Tests
{
    // BioInfoLogics worked example: ATCGATCAC, k=3 -> 7 total, 6 distinct, 5 unique.
    private const string AtcgatcacSeq = "ATCGATCAC";

    // Wikipedia AGAT example: all three 2-mers are distinct (each count 1).
    private const string AgatSeq = "AGAT";

    // ACGTACGT, k=4: ACGT occurs twice (pos 0,4); CGTA, GTAC, TACG once each.
    private const string Acgtacgt = "ACGTACGT";

    #region FindUniqueKmers — MUST

    // M1 — BioInfoLogics worked table: unique 3-mers of ATCGATCAC are exactly
    // {TCG, CGA, GAT, TCA, CAC}; ATC appears twice (count 2) and is excluded.
    [Test]
    public void FindUniqueKmers_AtcgatcacK3_ReturnsFiveUniqueExcludingRepeated()
    {
        var result = KmerAnalyzer.FindUniqueKmers(AtcgatcacSeq, 3).ToHashSet();

        var expected = new HashSet<string> { "TCG", "CGA", "GAT", "TCA", "CAC" };
        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EquivalentTo(expected),
                "Unique 3-mers of ATCGATCAC must be the 5 count-1 k-mers per BioInfoLogics worked table.");
            Assert.That(result, Does.Not.Contain("ATC"),
                "ATC occurs twice (positions 0 and 4) so it is distinct but NOT unique.");
        });
    }

    // M2 — Wikipedia AGAT example: all 2-mers AG, GA, AT are distinct, each count 1.
    [Test]
    public void FindUniqueKmers_AgatK2_ReturnsAllThreeTwoMers()
    {
        var result = KmerAnalyzer.FindUniqueKmers(AgatSeq, 2).ToHashSet();

        Assert.That(result, Is.EquivalentTo(new[] { "AG", "GA", "AT" }),
            "AGAT has 2-mers AG, GA, AT each appearing once (Wikipedia AGAT example).");
    }

    // M3 — Homopolymer: AAAAA k=3 has a single distinct 3-mer AAA with count 3 (>1),
    // so there are zero unique k-mers (count>1 => not unique, BioInfoLogics).
    [Test]
    public void FindUniqueKmers_HomopolymerK3_ReturnsEmpty()
    {
        var result = KmerAnalyzer.FindUniqueKmers("AAAAA", 3).ToList();

        Assert.That(result, Is.Empty,
            "AAAAA k=3 has only AAA (count 3); a count>1 k-mer is not unique, so the set is empty.");
    }

    // M11 — INV-1 property: every k-mer returned by FindUniqueKmers has count exactly 1,
    // cross-checked independently against the CountKmers frequency map.
    [Test]
    public void FindUniqueKmers_EveryReturnedKmer_HasCountOne()
    {
        const string seq = "ATCGATCACGATCG";
        var counts = KmerAnalyzer.CountKmers(seq, 3);

        var unique = KmerAnalyzer.FindUniqueKmers(seq, 3).ToList();

        Assert.That(unique, Is.Not.Empty, "Sequence must contain at least one count-1 3-mer.");
        Assert.That(unique.All(k => counts[k] == 1), Is.True,
            "INV-1: every k-mer FindUniqueKmers returns must have occurrence count exactly 1.");
        // Independent check: the unique set equals {k : count==1} from the map.
        var expected = counts.Where(kvp => kvp.Value == 1).Select(kvp => kvp.Key).ToHashSet();
        Assert.That(unique.ToHashSet(), Is.EquivalentTo(expected),
            "Unique set must equal exactly the count-1 keys of the frequency map.");
    }

    #endregion

    #region FindUniqueKmers — edge / SHOULD / COULD

    // M7 — empty sequence and k > length both yield no k-mers (L-k+1<=0).
    [Test]
    public void FindUniqueKmers_EmptyOrKExceedsLength_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.FindUniqueKmers("", 3).ToList(), Is.Empty,
                "Empty sequence has no k-mers.");
            Assert.That(KmerAnalyzer.FindUniqueKmers("ACG", 5).ToList(), Is.Empty,
                "k=5 > length 3 gives L-k+1 < 0, so no k-mers.");
        });
    }

    // M9 — k <= 0 is invalid (k-mer length must be positive).
    [Test]
    public void FindUniqueKmers_NonPositiveK_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => KmerAnalyzer.FindUniqueKmers("ACGT", 0).ToList(),
            "k must be positive; k=0 is invalid.");
    }

    // S1 — case normalisation: lower-case input yields the same unique set as upper-case.
    [Test]
    public void FindUniqueKmers_LowerCaseInput_MatchesUpperCase()
    {
        var lower = KmerAnalyzer.FindUniqueKmers("atcgatcac", 3).ToHashSet();

        Assert.That(lower, Is.EquivalentTo(new[] { "TCG", "CGA", "GAT", "TCA", "CAC" }),
            "Input is upper-cased internally; lower-case ATCGATCAC has the same unique set.");
    }

    // C1 — monomers (k=1): a letter is unique iff it appears exactly once.
    [Test]
    public void FindUniqueKmers_MonomersK1_ReturnsLettersAppearingOnce()
    {
        // AGAT: A appears twice (not unique), G once, T once.
        var result = KmerAnalyzer.FindUniqueKmers(AgatSeq, 1).ToHashSet();

        Assert.That(result, Is.EquivalentTo(new[] { "G", "T" }),
            "k=1: G and T appear once; A appears twice so is excluded.");
    }

    #endregion

    #region FindKmersWithMinCount — MUST

    // M4 — recurrent filter: ACGTACGT k=4 minCount=2 -> only ACGT (count 2).
    [Test]
    public void FindKmersWithMinCount_Acgtacgt_K4_Min2_ReturnsOnlyRepeated()
    {
        var result = KmerAnalyzer.FindKmersWithMinCount(Acgtacgt, 4, 2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1),
                "Only ACGT occurs >= 2 times in ACGTACGT (k=4).");
            Assert.That(result[0].Kmer, Is.EqualTo("ACGT"),
                "ACGT is the only recurrent 4-mer.");
            Assert.That(result[0].Count, Is.EqualTo(2),
                "ACGT occurs exactly twice (positions 0 and 4).");
        });
    }

    // M5 + M6 — minCount=1 returns all distinct 4-mers with correct counts,
    // ordered by count descending (ACGT count 2 first).
    [Test]
    public void FindKmersWithMinCount_Acgtacgt_K4_Min1_ReturnsAllDistinctOrderedDesc()
    {
        var result = KmerAnalyzer.FindKmersWithMinCount(Acgtacgt, 4, 1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(4),
                "ACGTACGT k=4 has 4 distinct 4-mers: ACGT, CGTA, GTAC, TACG.");
            Assert.That(result[0].Kmer, Is.EqualTo("ACGT"),
                "Count-descending order puts ACGT (count 2) first.");
            Assert.That(result[0].Count, Is.EqualTo(2),
                "ACGT occurs twice.");
            var asMap = result.ToDictionary(p => p.Kmer, p => p.Count);
            Assert.That(asMap["CGTA"], Is.EqualTo(1), "CGTA occurs once.");
            Assert.That(asMap["GTAC"], Is.EqualTo(1), "GTAC occurs once.");
            Assert.That(asMap["TACG"], Is.EqualTo(1), "TACG occurs once.");
        });
    }

    // M12 — INV-5: with minCount=1 the returned keys equal the distinct k-mer set,
    // cross-checked against CountKmers.
    [Test]
    public void FindKmersWithMinCount_Min1Keys_EqualDistinctKmerSet()
    {
        const string seq = "ATCGATCACGATCG";
        var distinct = KmerAnalyzer.CountKmers(seq, 3).Keys.ToHashSet();

        var keys = KmerAnalyzer.FindKmersWithMinCount(seq, 3, 1).Select(p => p.Kmer).ToHashSet();

        Assert.That(keys, Is.EquivalentTo(distinct),
            "INV-5: minCount=1 selects every distinct k-mer (count>=1 for all).");
    }

    // S2 — ordering invariant INV-4: counts are non-increasing across the result.
    [Test]
    public void FindKmersWithMinCount_OutputCounts_AreNonIncreasing()
    {
        var counts = KmerAnalyzer.FindKmersWithMinCount("AAAACGTAAA", 2, 1)
            .Select(p => p.Count).ToList();

        for (int i = 1; i < counts.Count; i++)
        {
            Assert.That(counts[i], Is.LessThanOrEqualTo(counts[i - 1]),
                "INV-4: FindKmersWithMinCount must order k-mers by count descending.");
        }
    }

    // S3 — threshold above the maximum observed count yields empty.
    [Test]
    public void FindKmersWithMinCount_ThresholdAboveMaxCount_ReturnsEmpty()
    {
        // Max count of any 4-mer in ACGTACGT is 2; minCount=3 selects nothing.
        var result = KmerAnalyzer.FindKmersWithMinCount(Acgtacgt, 4, 3).ToList();

        Assert.That(result, Is.Empty,
            "No 4-mer occurs 3 times in ACGTACGT, so minCount=3 returns empty.");
    }

    #endregion

    #region FindKmersWithMinCount — edge

    // M8 — empty sequence and k > length yield empty.
    [Test]
    public void FindKmersWithMinCount_EmptyOrKExceedsLength_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.FindKmersWithMinCount("", 3, 1).ToList(), Is.Empty,
                "Empty sequence has no k-mers.");
            Assert.That(KmerAnalyzer.FindKmersWithMinCount("ACG", 5, 1).ToList(), Is.Empty,
                "k > length gives no k-mers.");
        });
    }

    // M10 — k <= 0 is invalid.
    [Test]
    public void FindKmersWithMinCount_NonPositiveK_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => KmerAnalyzer.FindKmersWithMinCount("ACGT", -1, 1).ToList(),
            "k must be positive; k=-1 is invalid.");
    }

    #endregion
}
