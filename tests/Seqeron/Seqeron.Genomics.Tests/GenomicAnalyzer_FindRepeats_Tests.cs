// GENOMIC-REPEAT-001 — Repeat Detection (Longest Repeated Substring + all repeats)
// Evidence: docs/Evidence/GENOMIC-REPEAT-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-REPEAT-001.md
// Source: CMU 15-451 Lecture #10 §2.1 (2017); Wikipedia "Longest repeated substring problem" (2026);
//         GeeksforGeeks "Suffix Tree Application 3" — all accessed 2026-06-13.

using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class GenomicAnalyzer_FindRepeats_Tests
{
    #region FindLongestRepeat

    // M1 — Wikipedia worked example: ATCGATCGA$ -> ATCGA (deepest internal node), occurrences at 0 and 4.
    [Test]
    public void FindLongestRepeat_WikipediaExample_ReturnsAtcga()
    {
        var dna = new DnaSequence("ATCGATCGA");

        RepeatInfo lrs = GenomicAnalyzer.FindLongestRepeat(dna);

        Assert.Multiple(() =>
        {
            Assert.That(lrs.Sequence, Is.EqualTo("ATCGA"),
                "Wikipedia: longest repeated substring of ATCGATCGA is ATCGA.");
            Assert.That(lrs.Length, Is.EqualTo(5),
                "ATCGA has length 5 (root-to-deepest-internal-node string depth).");
            Assert.That(lrs.Count, Is.EqualTo(2),
                "ATCGA occurs exactly twice (a repeat is an internal node with >=2 leaves).");
            Assert.That(lrs.Positions, Is.EqualTo(new[] { 0, 4 }),
                "ATCGA starts at 0-based positions 0 and 4.");
        });
    }

    // M2 — GeeksforGeeks: AAAAAAAAAA (len 10) -> AAAAAAAAA (len 9); overlapping occurrences at 0 and 1.
    [Test]
    public void FindLongestRepeat_OverlappingRun_ReturnsLengthNine()
    {
        var dna = new DnaSequence("AAAAAAAAAA"); // 10 A's

        RepeatInfo lrs = GenomicAnalyzer.FindLongestRepeat(dna);

        Assert.Multiple(() =>
        {
            Assert.That(lrs.Sequence, Is.EqualTo("AAAAAAAAA"),
                "GeeksforGeeks: LRS of 10 A's is 9 A's (overlap permitted by the definition).");
            Assert.That(lrs.Length, Is.EqualTo(9), "9 A's.");
            Assert.That(lrs.Count, Is.EqualTo(2), "Occurs twice (overlapping).");
            Assert.That(lrs.Positions, Is.EqualTo(new[] { 0, 1 }),
                "The two overlapping occurrences start at 0 and 1.");
        });
    }

    // M3 — GeeksforGeeks ABABABA -> ABABA analog mapped to DNA: ATATATA -> ATATA; overlap at 0 and 2.
    [Test]
    public void FindLongestRepeat_OverlappingPeriodTwo_ReturnsAtata()
    {
        var dna = new DnaSequence("ATATATA"); // analog of ABABABA

        RepeatInfo lrs = GenomicAnalyzer.FindLongestRepeat(dna);

        Assert.Multiple(() =>
        {
            Assert.That(lrs.Sequence, Is.EqualTo("ATATA"),
                "ABABABA->ABABA analog: ATATATA's LRS is ATATA.");
            Assert.That(lrs.Length, Is.EqualTo(5), "ATATA has length 5.");
            Assert.That(lrs.Count, Is.EqualTo(2), "Occurs twice (overlapping).");
            Assert.That(lrs.Positions, Is.EqualTo(new[] { 0, 2 }),
                "The two overlapping occurrences start at 0 and 2.");
        });
    }

    // M4 — GeeksforGeeks ABCDEFG (no repeat) analog: ACGT has no repeated substring.
    [Test]
    public void FindLongestRepeat_NoRepeat_ReturnsNone()
    {
        var dna = new DnaSequence("ACGT"); // all distinct chars -> no internal node with >=2 leaves

        RepeatInfo lrs = GenomicAnalyzer.FindLongestRepeat(dna);

        Assert.Multiple(() =>
        {
            Assert.That(lrs.IsEmpty, Is.True,
                "No substring repeats in ACGT, so the result is RepeatInfo.None (ABCDEFG analog).");
            Assert.That(lrs.Sequence, Is.Empty, "None carries an empty sequence.");
            Assert.That(lrs.Count, Is.EqualTo(0), "None carries zero positions.");
        });
    }

    // M5 — Empty input: the empty string has no substring occurring twice.
    [Test]
    public void FindLongestRepeat_EmptySequence_ReturnsNone()
    {
        var dna = new DnaSequence("");

        RepeatInfo lrs = GenomicAnalyzer.FindLongestRepeat(dna);

        Assert.That(lrs.IsEmpty, Is.True,
            "No substring occurs twice in the empty sequence, so the result is RepeatInfo.None.");
    }

    #endregion

    #region FindRepeats

    private const string EnumInput = "ACGTACGTTTTTACGT";

    private static Dictionary<string, int[]> ToMap(IEnumerable<RepeatInfo> repeats) =>
        repeats.ToDictionary(r => r.Sequence, r => r.Positions.ToArray());

    // M6 — Full enumeration at minLength 3: exact set of repeated substrings with exact positions.
    // Each listed substring occurs >=2 times in ACGTACGTTTTTACGT (positions 0-based).
    [Test]
    public void FindRepeats_MinLengthThree_ReturnsExactSetWithPositions()
    {
        var dna = new DnaSequence(EnumInput);

        var map = ToMap(GenomicAnalyzer.FindRepeats(dna, minLength: 3));

        Assert.Multiple(() =>
        {
            Assert.That(map.Keys, Is.EquivalentTo(new[] { "ACGT", "CGT", "TACGT", "TTT", "TTTT" }),
                "Exactly these substrings occur >=2 times with length >=3 in " + EnumInput + ".");
            Assert.That(map["ACGT"], Is.EqualTo(new[] { 0, 4, 12 }), "ACGT occurs at 0,4,12.");
            Assert.That(map["CGT"], Is.EqualTo(new[] { 1, 5, 13 }), "CGT occurs at 1,5,13.");
            Assert.That(map["TACGT"], Is.EqualTo(new[] { 3, 11 }), "TACGT occurs at 3,11.");
            Assert.That(map["TTT"], Is.EqualTo(new[] { 7, 8, 9 }), "TTT occurs at 7,8,9 (overlapping).");
            Assert.That(map["TTTT"], Is.EqualTo(new[] { 7, 8 }), "TTTT occurs at 7,8 (overlapping).");
        });
    }

    // M7 — Every returned repeat occurs at least twice and meets minLength (INV-1, INV-5).
    [Test]
    public void FindRepeats_AllResults_OccurAtLeastTwiceAndMeetMinLength()
    {
        var dna = new DnaSequence(EnumInput);

        var repeats = GenomicAnalyzer.FindRepeats(dna, minLength: 3).ToList();

        Assert.That(repeats, Is.Not.Empty, "The input contains repeats of length >=3.");
        Assert.Multiple(() =>
        {
            foreach (var r in repeats)
            {
                Assert.That(r.Count, Is.GreaterThanOrEqualTo(2),
                    $"INV-1: '{r.Sequence}' must occur at least twice.");
                Assert.That(r.Length, Is.GreaterThanOrEqualTo(3),
                    $"INV-5: '{r.Sequence}' must have length >= minLength (3).");
            }
        });
    }

    // S1 — minLength above any repeat length yields no results.
    [Test]
    public void FindRepeats_MinLengthAboveAllRepeats_ReturnsEmpty()
    {
        var dna = new DnaSequence("ACGTACGT"); // only repeat is ACGT (len 4)

        var repeats = GenomicAnalyzer.FindRepeats(dna, minLength: 5).ToList();

        Assert.That(repeats, Is.Empty,
            "No repeated substring has length >=5 in ACGTACGT, so the result is empty.");
    }

    // S2 — minLength is inclusive (>=, not >): a repeat exactly equal to minLength is returned.
    [Test]
    public void FindRepeats_MinLengthEqualsRepeatLength_IncludesRepeat()
    {
        var dna = new DnaSequence("ACGTACGT");

        var map = ToMap(GenomicAnalyzer.FindRepeats(dna, minLength: 4));

        Assert.Multiple(() =>
        {
            Assert.That(map.ContainsKey("ACGT"), Is.True,
                "ACGT (length 4) must be returned when minLength == 4 (filter is inclusive).");
            Assert.That(map["ACGT"], Is.EqualTo(new[] { 0, 4 }), "ACGT occurs at 0 and 4.");
        });
    }

    // S3 — Property test (INV-1..INV-3, INV-6) over the O(n^2) enumeration path.
    [Test]
    public void FindRepeats_AllResults_SatisfyInvariants()
    {
        var dna = new DnaSequence(EnumInput);

        var repeats = GenomicAnalyzer.FindRepeats(dna, minLength: 2).ToList();

        Assert.That(repeats, Is.Not.Empty, "Input has repeats of length >=2.");
        Assert.Multiple(() =>
        {
            foreach (var r in repeats)
            {
                Assert.That(r.Count, Is.GreaterThanOrEqualTo(2),
                    $"INV-1: '{r.Sequence}' occurs >=2 times.");
                Assert.That(r.Length, Is.EqualTo(r.Sequence.Length),
                    $"INV-2: Length equals Sequence.Length for '{r.Sequence}'.");
                Assert.That(r.Count, Is.EqualTo(r.Positions.Count),
                    $"INV-2/3: Count equals number of positions for '{r.Sequence}'.");
                foreach (var p in r.Positions)
                {
                    Assert.That(EnumInput.Substring(p, r.Length), Is.EqualTo(r.Sequence),
                        $"INV-3: position {p} is a true occurrence of '{r.Sequence}'.");
                }
                Assert.That(r.Positions, Is.Ordered.Ascending,
                    $"INV-6: positions of '{r.Sequence}' are sorted ascending.");
            }
        });
    }

    // C1 — Degenerate minLength <= 0 still yields only substrings occurring >=2 times; no zero-length results.
    [Test]
    public void FindRepeats_MinLengthZero_ReturnsNoZeroLengthRepeats()
    {
        var dna = new DnaSequence("ACGTACGT");

        var repeats = GenomicAnalyzer.FindRepeats(dna, minLength: 0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(repeats.All(r => r.Length >= 1), Is.True,
                "minLength<=0 must not produce zero-length repeats; only real substrings qualify.");
            Assert.That(repeats.All(r => r.Count >= 2), Is.True,
                "INV-1: even with minLength<=0, every result occurs at least twice.");
            Assert.That(repeats.Any(r => r.Sequence == "ACGT" && r.Positions.SequenceEqual(new[] { 0, 4 })),
                Is.True, "ACGT@{0,4} is among the returned repeats.");
        });
    }

    #endregion
}
