// KMER-POSITIONS-001 — K-mer Positions
// Evidence: docs/Evidence/KMER-POSITIONS-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-POSITIONS-001.md
// Source: Rosalind BA1D (Find All Occurrences of a Pattern in a String),
//         https://rosalind.info/problems/ba1d/ ; Wikipedia "k-mer".

using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_FindKmerPositions_Tests
{
    #region FindKmerPositions

    // M1 — Rosalind BA1D sample: Pattern ATAT in GATATATGCATATACTT → 1 3 9 (0-based, overlapping).
    [Test]
    public void FindKmerPositions_RosalindBA1D_ReturnsExpectedPositions()
    {
        const string sequence = "GATATATGCATATACTT";
        const string kmer = "ATAT";

        var result = KmerAnalyzer.FindKmerPositions(sequence, kmer).ToList();

        // Exact Rosalind BA1D output; a non-overlapping or 1-based impl would NOT produce this.
        Assert.That(result, Is.EqualTo(new[] { 1, 3, 9 }),
            "Rosalind BA1D: ATAT occurs (overlapping, 0-based) at 1, 3 and 9");
    }

    // M2 — Overlapping self-occurrence: AA in AAAA → 0 1 2 (L-k+1 = 3).
    [Test]
    public void FindKmerPositions_OverlappingSelfOccurrence_ReturnsAllStarts()
    {
        var result = KmerAnalyzer.FindKmerPositions("AAAA", "AA").ToList();

        // Overlapping occurrences are all reported (Rosalind BA1D); a non-overlapping
        // scanner would return [0, 2] only.
        Assert.That(result, Is.EqualTo(new[] { 0, 1, 2 }),
            "AA overlaps itself at every start 0,1,2 in AAAA");
    }

    // M3 — Classic overlapping pattern: ana in banana → 1 3.
    [Test]
    public void FindKmerPositions_AnaInBanana_ReturnsOverlappingStarts()
    {
        var result = KmerAnalyzer.FindKmerPositions("banana", "ana").ToList();

        Assert.That(result, Is.EqualTo(new[] { 1, 3 }),
            "ana occurs at 0-based starts 1 and 3 in banana (overlapping)");
    }

    // M4 — Wikipedia AGAT 2-mers: AG@0, GA@1, AT@2.
    [Test]
    public void FindKmerPositions_AgatTwoMers_ReturnsEachStart()
    {
        const string sequence = "AGAT";

        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.FindKmerPositions(sequence, "AG").ToList(),
                Is.EqualTo(new[] { 0 }), "2-mer AG starts at index 0 (Wikipedia k-mer)");
            Assert.That(KmerAnalyzer.FindKmerPositions(sequence, "GA").ToList(),
                Is.EqualTo(new[] { 1 }), "2-mer GA starts at index 1");
            Assert.That(KmerAnalyzer.FindKmerPositions(sequence, "AT").ToList(),
                Is.EqualTo(new[] { 2 }), "2-mer AT starts at index 2");
        });
    }

    // M5 — Pattern absent: GG in ATATAT → empty.
    [Test]
    public void FindKmerPositions_PatternAbsent_ReturnsEmpty()
    {
        var result = KmerAnalyzer.FindKmerPositions("ATATAT", "GG").ToList();

        Assert.That(result, Is.Empty,
            "GG does not occur in ATATAT; only matching starts are reported");
    }

    // M6 — Ascending order invariant (INV-2) on an overlapping pattern.
    [Test]
    public void FindKmerPositions_OverlappingMatches_ReturnedInAscendingOrder()
    {
        var result = KmerAnalyzer.FindKmerPositions("ATATATAT", "ATAT").ToList();

        // ATAT overlaps at 0,2,4 in ATATATAT; verify exact values AND ascending order.
        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EqualTo(new[] { 0, 2, 4 }),
                "ATAT occurs (overlapping) at 0,2,4 in ATATATAT");
            Assert.That(result, Is.Ordered.Ascending,
                "positions must be returned in strictly ascending order (INV-2)");
        });
    }

    // S1 — Pattern longer than text: ACGT in AC → empty (L-k+1 ≤ 0).
    [Test]
    public void FindKmerPositions_PatternLongerThanText_ReturnsEmpty()
    {
        var result = KmerAnalyzer.FindKmerPositions("AC", "ACGT").ToList();

        Assert.That(result, Is.Empty,
            "No candidate start positions when |kmer| > |sequence| (L-k+1 ≤ 0)");
    }

    // S2 — Pattern equals whole sequence: ACGT in ACGT → [0].
    [Test]
    public void FindKmerPositions_PatternEqualsSequence_ReturnsZeroOnly()
    {
        var result = KmerAnalyzer.FindKmerPositions("ACGT", "ACGT").ToList();

        Assert.That(result, Is.EqualTo(new[] { 0 }),
            "A pattern equal to the whole sequence occurs exactly once, at index 0");
    }

    // S3 — Case-insensitive matching (repository convention): atat ≡ ATAT.
    [Test]
    public void FindKmerPositions_LowercaseKmer_MatchesCaseInsensitively()
    {
        var result = KmerAnalyzer.FindKmerPositions("GATATATGCATATACTT", "atat").ToList();

        Assert.That(result, Is.EqualTo(new[] { 1, 3, 9 }),
            "Matching is case-insensitive: lowercase atat matches ATAT (same as M1)");
    }

    // C1 — null/empty sequence → empty.
    [Test]
    public void FindKmerPositions_NullOrEmptySequence_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.FindKmerPositions(null!, "AT").ToList(), Is.Empty,
                "null sequence yields no positions");
            Assert.That(KmerAnalyzer.FindKmerPositions("", "AT").ToList(), Is.Empty,
                "empty sequence yields no positions");
        });
    }

    // C2 — null/empty kmer → empty.
    [Test]
    public void FindKmerPositions_NullOrEmptyKmer_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.FindKmerPositions("ACGT", null!).ToList(), Is.Empty,
                "null kmer yields no positions");
            Assert.That(KmerAnalyzer.FindKmerPositions("ACGT", "").ToList(), Is.Empty,
                "empty kmer yields no positions");
        });
    }

    #endregion
}
