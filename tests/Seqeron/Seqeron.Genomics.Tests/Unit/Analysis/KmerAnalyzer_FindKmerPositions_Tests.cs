// KMER-POSITIONS-001 — K-mer Positions
// Evidence: docs/Evidence/KMER-POSITIONS-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-POSITIONS-001.md
// Source: Rosalind BA1D (Find All Occurrences of a Pattern in a String),
//         https://rosalind.info/problems/ba1d/ ; Wikipedia "k-mer".

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

    #region Reference cross-checks (B06 review, KMER-POSITIONS-001)

    // R1 — Rosalind SUBS "Finding a Motif in DNA" (current sample, quoted in the
    // breezedu/rosalind FindingaMotifinDNA.java header): GATATATGCATATACTT / ATAT → 2 4 10
    // in SUBS's 1-based convention, i.e. 0-based 1 3 9 (= Rosalind BA1D; Python re lookahead
    // finditer and Biopython nt_search both give [1, 3, 9]).
    [Test]
    public void FindKmerPositions_RosalindSubsSample_OneBasedIs_2_4_10()
    {
        var zeroBased = KmerAnalyzer.FindKmerPositions("GATATATGCATATACTT", "ATAT").ToList();

        Assert.That(zeroBased.Select(p => p + 1), Is.EqualTo(new[] { 2, 4, 10 }),
            "Rosalind SUBS sample output (1-based)");
    }

    // R2 — Rosalind SUBS older sample (mtarbit/Rosalind-Problems e009-subs.py):
    // ACGTACGTACGTACGT / GTA → 3 7 11 (1-based) → 0-based 2 6 10 (re lookahead: [2, 6, 10]).
    [Test]
    public void FindKmerPositions_RosalindSubsLegacySample_ZeroBased_2_6_10()
    {
        Assert.That(KmerAnalyzer.FindKmerPositions("ACGTACGTACGTACGT", "GTA").ToList(),
            Is.EqualTo(new[] { 2, 6, 10 }));
    }

    // R3 — values produced by Python re.finditer("(?=P)") and Biopython 1.88
    // Bio.SeqUtils.nt_search (forward strand, 0-based); identical for every row.
    [TestCase("AAAAAAAAAA", "AAA", new[] { 0, 1, 2, 3, 4, 5, 6, 7 })]
    [TestCase("ACGTACGTACGT", "CGTA", new[] { 1, 5 })]
    [TestCase("GATATATGCATATACTT", "GCAT", new[] { 7 })]
    [TestCase("ATATATAT", "ATAT", new[] { 0, 2, 4 })]
    public void FindKmerPositions_MatchesBiopythonNtSearchAndRegexLookahead(
        string sequence, string kmer, int[] expected)
    {
        Assert.That(KmerAnalyzer.FindKmerPositions(sequence, kmer).ToList(), Is.EqualTo(expected));
    }

    // R4 — patterns with non-trivial borders exercise the KMP prefix-function fallback
    // (expected values from Python re lookahead finditer; ABABX case from the
    // TheAlgorithms/Python knuth_morris_pratt.py self-test text).
    [TestCase("ABABZABABYABABX", "ABABX", new[] { 10 })]
    [TestCase("AABAABAAB", "AABAAB", new[] { 0, 3 })]
    [TestCase("AABAACAADAABAABA", "AABA", new[] { 0, 9, 12 })]
    public void FindKmerPositions_BorderedPatterns_AllOverlappingStartsReported(
        string sequence, string kmer, int[] expected)
    {
        Assert.That(KmerAnalyzer.FindKmerPositions(sequence, kmer).ToList(), Is.EqualTo(expected));
    }

    // R5 — case folding applies to both arguments.
    [Test]
    public void FindKmerPositions_LowercaseSequenceMixedCaseKmer_MatchesCaseInsensitively()
    {
        Assert.That(KmerAnalyzer.FindKmerPositions("gatatatgcatatactt", "aTaT").ToList(),
            Is.EqualTo(new[] { 1, 3, 9 }));
    }

    // R6 — worst case for a window-by-window comparison: A^2,000,000 searched for
    // A^199,999·C costs ~4·10^11 character comparisons naively (Θ((L−k+1)·k)), while KMP
    // reads the text once (O(L + k); Knuth, Morris & Pratt 1977). Must finish quickly and
    // return the exact answer; the all-A pattern returns L − k + 1 overlapping starts.
    [Test]
    public void FindKmerPositions_AdversarialHomopolymer_IsLinearTime()
    {
        string text = new('A', 2_000_000);
        string miss = new string('A', 199_999) + "C";
        string hit = new('A', 200_000);

        var sw = System.Diagnostics.Stopwatch.StartNew();
        int missCount = KmerAnalyzer.FindKmerPositions(text, miss).Count();
        var hits = KmerAnalyzer.FindKmerPositions(text, hit).ToList();
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(missCount, Is.Zero);
            Assert.That(hits, Has.Count.EqualTo(2_000_000 - 200_000 + 1));
            Assert.That(hits[0], Is.Zero);
            Assert.That(hits[^1], Is.EqualTo(2_000_000 - 200_000));
            Assert.That(sw.Elapsed, Is.LessThan(TimeSpan.FromSeconds(5)),
                "O(L + k) matcher; a naive scan needs ~4e11 comparisons here");
        });
    }

    #endregion
}
