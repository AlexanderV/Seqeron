// GENOMIC-COMMON-001 — Common Region Detection (Longest Common Substring + all common regions)
// Evidence: docs/Evidence/GENOMIC-COMMON-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-COMMON-001.md
// Source: Gusfield D. (1997). Algorithms on Strings, Trees and Sequences. Cambridge Univ. Press, ISBN 0-521-58519-8;
//         Wikipedia "Longest common substring"; GeeksforGeeks "Suffix Tree Application 5" — all accessed 2026-06-13.

using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class GenomicAnalyzer_FindCommonRegion_Tests
{
    #region FindLongestCommonRegion

    // M1 — Definition (Wikipedia): longest CONTIGUOUS substring of both. ACGTACGT vs TTACGTGG -> TACGT (len 5).
    // Brute-force-verified (enumerate s[i:j] in t): maximal is TACGT at seq1 pos 3, seq2 pos 1.
    [Test]
    public void FindLongestCommonRegion_UniqueLcs_ReturnsTacgt()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("TACGT"),
                "Longest contiguous substring common to ACGTACGT and TTACGTGG is TACGT.");
            Assert.That(region.Length, Is.EqualTo(5), "TACGT has length 5; no length-6 common substring exists.");
            Assert.That(region.PositionInFirst, Is.EqualTo(3), "TACGT starts at 0-based index 3 in ACGTACGT.");
            Assert.That(region.PositionInSecond, Is.EqualTo(1), "TACGT starts at 0-based index 1 in TTACGTGG.");
        });
    }

    // M2 — Tie (Wikipedia: BADANAT/CANADAS share ADA & ANA). DNA analogue CACAGAG vs TACATAGAT shares
    // two distinct length-3 substrings ACA and AGA; documented tie-break returns the one first found in
    // sequence2 (TACATAGAT): ACA ends earlier than AGA -> ACA.
    [Test]
    public void FindLongestCommonRegion_Tie_ReturnsFirstInOther()
    {
        var seq1 = new DnaSequence("CACAGAG");
        var seq2 = new DnaSequence("TACATAGAT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("ACA"),
                "Both ACA and AGA are maximal (len 3); tie-break returns the one first found in sequence2 -> ACA.");
            Assert.That(region.Length, Is.EqualTo(3), "Maximal common substring length is 3.");
            Assert.That(region.PositionInFirst, Is.EqualTo(1), "ACA starts at index 1 in CACAGAG.");
            Assert.That(region.PositionInSecond, Is.EqualTo(1), "ACA starts at index 1 in TACATAGAT.");
        });
    }

    // M3 — Identical sequences: a string is a substring of itself, so the whole sequence is the LCS at 0/0.
    [Test]
    public void FindLongestCommonRegion_Identical_ReturnsWholeSequence()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("ACGT"), "Identical sequences share the whole sequence.");
            Assert.That(region.Length, Is.EqualTo(4), "Whole 4-base sequence is the LCS.");
            Assert.That(region.PositionInFirst, Is.EqualTo(0), "Starts at 0 in the first sequence.");
            Assert.That(region.PositionInSecond, Is.EqualTo(0), "Starts at 0 in the second sequence.");
        });
    }

    // M4 — No shared character: only the empty string qualifies -> length-0 LCS -> CommonRegion.None.
    [Test]
    public void FindLongestCommonRegion_NoCommon_ReturnsNone()
    {
        var seq1 = new DnaSequence("AAAA");
        var seq2 = new DnaSequence("GGGG");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.IsEmpty, Is.True, "No shared character -> empty LCS (length 0).");
            Assert.That(region.Sequence, Is.EqualTo(string.Empty), "Empty substring is returned.");
            Assert.That(region.Length, Is.EqualTo(0), "Length is 0.");
            Assert.That(region.PositionInFirst, Is.EqualTo(-1), "Sentinel position -1 when no region.");
            Assert.That(region.PositionInSecond, Is.EqualTo(-1), "Sentinel position -1 when no region.");
        });
    }

    // S1 — Empty first sequence -> CommonRegion.None (length-0 LCS).
    [Test]
    public void FindLongestCommonRegion_EmptyFirst_ReturnsNone()
    {
        var seq1 = new DnaSequence("");
        var seq2 = new DnaSequence("ACGT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.That(region.IsEmpty, Is.True, "Empty first sequence shares nothing -> CommonRegion.None.");
    }

    // S2 — Empty second sequence -> CommonRegion.None.
    [Test]
    public void FindLongestCommonRegion_EmptySecond_ReturnsNone()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.That(region.IsEmpty, Is.True, "Empty second sequence shares nothing -> CommonRegion.None.");
    }

    // S4 — Single-character overlap: minimal non-empty LCS. "A" vs "TTTAT" -> A at seq1 pos 0, seq2 pos 3.
    [Test]
    public void FindLongestCommonRegion_SingleCharOverlap_ReturnsA()
    {
        var seq1 = new DnaSequence("A");
        var seq2 = new DnaSequence("TTTAT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("A"), "Only the single base A is shared.");
            Assert.That(region.Length, Is.EqualTo(1), "Length-1 LCS.");
            Assert.That(region.PositionInFirst, Is.EqualTo(0), "A is at index 0 in the first sequence.");
            Assert.That(region.PositionInSecond, Is.EqualTo(3), "First A in TTTAT is at index 3.");
        });
    }

    // C1 — INV-01 property: the returned substring occurs in both sequences at the reported positions.
    [Test]
    public void FindLongestCommonRegion_ReturnedSubstring_OccursInBothAtPositions()
    {
        var seq1 = new DnaSequence("GATTACACGT");
        var seq2 = new DnaSequence("CCTTACAGG");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.IsEmpty, Is.False, "These sequences share TTACA.");
            Assert.That(seq1.Sequence.Substring(region.PositionInFirst, region.Length), Is.EqualTo(region.Sequence),
                "INV-01: substring must appear at PositionInFirst in sequence1.");
            Assert.That(seq2.Sequence.Substring(region.PositionInSecond, region.Length), Is.EqualTo(region.Sequence),
                "INV-01: substring must appear at PositionInSecond in sequence2.");
        });
    }

    #endregion

    #region FindCommonRegions

    // M6 — Right-maximal common regions of length >= 3 in ACGTACGT vs TTACGTGG.
    // FindCommonRegions returns, per start position in seq2, the SINGLE LONGEST substring also in seq1
    // (deduplicated), NOT every common substring. Hand-enumerating the longest match at each start of
    // TTACGTGG (positions 0..5) and deduplicating gives exactly {TACGT@(3,1), ACGT@(0,2), CGT@(1,3)}.
    // The FULL set of all common substrings >=3 would also include the prefixes TAC, TACG, ACG — those
    // are intentionally OMITTED (see FindCommonRegions_DoesNotReturnPrefixesOfLongerMatches below).
    [Test]
    public void FindCommonRegions_MinLengthThree_ReturnsThreeRegions()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 3).ToList();
        var tuples = regions.Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();

        Assert.That(tuples, Is.EquivalentTo(new[]
        {
            ("TACGT", 3, 1),
            ("ACGT", 0, 2),
            ("CGT", 1, 3),
        }), "minLength=3 yields exactly the right-maximal set TACGT@(3,1), ACGT@(0,2), CGT@(1,3).");
    }

    // M5 — Same pair, minLength=4 drops CGT (length 3); leaves the right-maximal TACGT and ACGT.
    // Hand-verified: at each start of TTACGTGG, the longest substring also in ACGTACGT with length >=4 is
    // TACGT (start 1) and ACGT (start 2); all other starts have no match >=4.
    [Test]
    public void FindCommonRegions_MinLengthFour_ReturnsTwoRegions()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 4).ToList();
        var tuples = regions.Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();

        Assert.That(tuples, Is.EquivalentTo(new[]
        {
            ("TACGT", 3, 1),
            ("ACGT", 0, 2),
        }), "minLength=4 yields exactly TACGT@(3,1) and ACGT@(0,2).");
    }

    // M7 — Contract lock: FindCommonRegions reports ONLY the longest match per start position in seq2
    // (right-maximal), NOT every common substring. For ACGTACGT vs TTACGTGG (minLength 3) the prefixes
    // TAC, TACG (prefixes of TACGT) and ACG (prefix of ACGT) are common substrings of length >= 3 but are
    // intentionally NOT returned. Hand-verified by full enumeration of all common substrings >= 3:
    // {TAC, TACG, TACGT, ACG, ACGT, CGT}; the returned set is the right-maximal subset {TACGT, ACGT, CGT}.
    [Test]
    public void FindCommonRegions_DoesNotReturnPrefixesOfLongerMatches()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        var returned = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 3)
            .Select(r => r.Sequence)
            .ToList();

        Assert.Multiple(() =>
        {
            // These ARE common substrings of length >= 3, yet are deliberately omitted (not right-maximal).
            Assert.That(seq1.Sequence.Contains("TAC") && seq2.Sequence.Contains("TAC"), Is.True,
                "Sanity: TAC is a genuine common substring.");
            Assert.That(returned, Does.Not.Contain("TAC"), "Prefix TAC of TACGT must be omitted.");
            Assert.That(returned, Does.Not.Contain("TACG"), "Prefix TACG of TACGT must be omitted.");
            Assert.That(returned, Does.Not.Contain("ACG"), "Prefix ACG of ACGT must be omitted.");
            // The returned set is exactly the right-maximal one.
            Assert.That(returned, Is.EquivalentTo(new[] { "TACGT", "ACGT", "CGT" }),
                "Only the longest match per start position in seq2 is reported.");
        });
    }

    // S5 — minLength < 1 is treated as 1 (a region is non-empty). For ACGT vs ACGT the right-maximal
    // match at each start of seq2 is ACGT@(0,0), CGT@(1,1), GT@(2,2), T@(3,3). Hand-verified; minLength=0
    // and minLength=-5 must give the same result as minLength=1 (no empty/length-0 regions).
    [Test]
    public void FindCommonRegions_MinLengthBelowOne_TreatedAsOne()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGT");

        var expected = new[] { ("ACGT", 0, 0), ("CGT", 1, 1), ("GT", 2, 2), ("T", 3, 3) };

        var atZero = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 0)
            .Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();
        var atNeg = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: -5)
            .Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(atZero, Is.EquivalentTo(expected), "minLength=0 is treated as 1; no empty regions.");
            Assert.That(atNeg, Is.EquivalentTo(expected), "minLength<0 is treated as 1; no empty regions.");
            Assert.That(atZero.All(t => t.Sequence.Length >= 1), Is.True,
                "No length-0 / empty region is ever returned.");
        });
    }

    // S3 — No shared substring -> empty enumeration.
    [Test]
    public void FindCommonRegions_NoMatch_ReturnsEmpty()
    {
        var seq1 = new DnaSequence("AAAA");
        var seq2 = new DnaSequence("GGGG");

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 2).ToList();

        Assert.That(regions, Is.Empty, "No shared base -> no common regions.");
    }

    // C2 — INV-04 property (O(n*m) algorithm): every region has length >= minLength and is contained in sequence1.
    [Test]
    public void FindCommonRegions_AllRegions_SatisfyMinLengthAndContainment()
    {
        var seq1 = new DnaSequence("GATTACAGATTACA");
        var seq2 = new DnaSequence("TTACAGGGATTAC");
        const int MinLength = 3;

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, MinLength).ToList();

        Assert.That(regions, Is.Not.Empty, "Sequences share several substrings of length >= 3.");
        Assert.Multiple(() =>
        {
            foreach (var r in regions)
            {
                Assert.That(r.Length, Is.GreaterThanOrEqualTo(MinLength),
                    $"INV-04: region '{r.Sequence}' must be at least minLength={MinLength}.");
                Assert.That(seq1.Sequence.Contains(r.Sequence), Is.True,
                    $"INV-04: region '{r.Sequence}' must be a substring of sequence1.");
                Assert.That(seq2.Sequence.Substring(r.PositionInSecond, r.Length), Is.EqualTo(r.Sequence),
                    $"INV-01: region '{r.Sequence}' must appear at PositionInSecond in sequence2.");
            }
        });
    }

    #endregion
}
