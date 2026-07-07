// COMPGEN-DOTPLOT-001 — Dot Plot Generation (word-match / k-tuple dot matrix)
// Evidence: docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-DOTPLOT-001.md
// Source: Gibbs AJ, McIntyre GA (1970). Eur J Biochem 16(1):1-11. DOI:10.1111/j.1432-1033.1970.tb01046.x
//         EMBOSS dottup (Rice, Longden & Bleasby 2000) — exact word-match dot plot, default wordsize 10.
//
// A dot at (x, y) marks an EXACT match of a length-wordSize word starting at sequence1[x] and
// sequence2[y] (dottup). Expected coordinate sets below are derived by hand from that match rule
// and the published worked example (AGCGT vs AT, k=1 => {(0,0),(4,1)}), NOT read back from the code.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ComparativeGenomics_GenerateDotPlot_Tests
{
    private static HashSet<(int x, int y)> Dots(string s1, string s2, int wordSize, int stepSize = 1)
        => ComparativeGenomics.GenerateDotPlot(s1, s2, wordSize, stepSize).ToHashSet();

    #region GenerateDotPlot

    // M1 — Huttley worked example: X=AGCGT (x-axis), Y=AT (y-axis), wordSize=1. The only equal
    // characters are the A's at (x=0,y=0) and the T's at (x=4,y=1). A wrong impl that matched on
    // index instead of residue, or transposed axes, would not produce exactly this set.
    [Test]
    public void GenerateDotPlot_HuttleyExample_ReturnsTwoDots()
    {
        var dots = Dots("AGCGT", "AT", wordSize: 1);

        Assert.That(dots, Is.EquivalentTo(new[] { (0, 0), (4, 1) }),
            "AGCGT vs AT (k=1) must give dots only where residues match: A@(0,0) and T@(4,1) (Huttley/dottup).");
    }

    // M2 — Exact word match, wordSize=4 on ACGTACGT vs itself. Re-derived by hand from the match
    // rule over EVERY length-4 word start x=0..4 (n-w = 4):
    //   x=0 "ACGT" -> y=0,4 ; x=1 "CGTA" -> y=1 ; x=2 "GTAC" -> y=2 ; x=3 "TACG" -> y=3 ;
    //   x=4 "ACGT" -> y=0,4. The repeated word "ACGT" gives the off-diagonal dots (0,4) and (4,0)
    // (all overlapping occurrences reported, dottup), while the unique internal words give the
    // diagonal (1,1),(2,2),(3,3). A first-occurrence-only or single-word impl would miss members.
    [Test]
    public void GenerateDotPlot_ExactWordMatch_ReturnsAllOverlappingOccurrences()
    {
        var dots = Dots("ACGTACGT", "ACGTACGT", wordSize: 4);

        Assert.That(dots, Is.EquivalentTo(new[] { (0, 0), (0, 4), (1, 1), (2, 2), (3, 3), (4, 0), (4, 4) }),
            "All length-4 words matched at all positions: ACGT->{0,4} at x in {0,4} plus the unique internal words on the diagonal (dottup, INV-01).");
    }

    // M3 — Self-comparison main diagonal: ACGT vs itself, wordSize=1. "The main diagonal represents
    // the sequence's alignment with itself" (Wikipedia; INV-02). ACGT has four DISTINCT residues, so
    // the only equal-character pairs are (i,i): the exact match set is the full main diagonal and
    // nothing else. Assert the exact set (not a superset) — a wrong impl emitting extra off-diagonal
    // dots, or missing a diagonal dot, would fail.
    [Test]
    public void GenerateDotPlot_SelfComparison_ContainsMainDiagonal()
    {
        var dots = Dots("ACGT", "ACGT", wordSize: 1);

        Assert.That(dots, Is.EquivalentTo(new[] { (0, 0), (1, 1), (2, 2), (3, 3) }),
            "Self-comparison of distinct-residue ACGT is exactly the full main diagonal (i,i) (Wikipedia main diagonal; INV-02).");
    }

    // M4 — No shared residues: AAAA vs CCCC, wordSize=1. A dot is placed only on a match, so with
    // disjoint alphabets the result is empty. A buggy impl that emitted a dot per position pair
    // regardless of equality would fail here.
    [Test]
    public void GenerateDotPlot_NoSharedResidues_ReturnsEmpty()
    {
        var dots = Dots("AAAA", "CCCC", wordSize: 1);

        Assert.That(dots, Is.Empty,
            "Disjoint alphabets share no residue, so no dots may be plotted (dot-on-match rule).");
    }

    // M5 — Word longer than a sequence: no length-4 word can start in a length-3 sequence, so no
    // dots can be formed (EMBOSS dottup / manpage; INV-03).
    [Test]
    public void GenerateDotPlot_WordLongerThanSequence_ReturnsEmpty()
    {
        var dots = Dots("ACG", "ACGT", wordSize: 4);

        Assert.That(dots, Is.Empty,
            "A word longer than sequence1 cannot start, so no dots are produced (dottup; INV-03).");
    }

    // M6 — Null/empty inputs yield no dots (INV-03).
    [Test]
    public void GenerateDotPlot_NullOrEmptyInput_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Dots(null!, "ACGT", 1), Is.Empty, "Null sequence1 must yield no dots (INV-03).");
            Assert.That(Dots("ACGT", null!, 1), Is.Empty, "Null sequence2 must yield no dots (INV-03).");
            Assert.That(Dots("", "ACGT", 1), Is.Empty, "Empty sequence1 must yield no dots (INV-03).");
            Assert.That(Dots("ACGT", "", 1), Is.Empty, "Empty sequence2 must yield no dots (INV-03).");
        });
    }

    // M7 — Non-positive wordSize is undefined for a word-match dot plot and must throw (INV-05).
    [Test]
    public void GenerateDotPlot_NonPositiveWordSize_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: 0).ToList(),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(), "wordSize=0 must throw (INV-05).");
            Assert.That(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: -1).ToList(),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(), "Negative wordSize must throw (INV-05).");
        });
    }

    // M8 — Non-positive stepSize is undefined and must throw (INV-05).
    [Test]
    public void GenerateDotPlot_NonPositiveStepSize_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: 2, stepSize: 0).ToList(),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(), "stepSize=0 must throw (INV-05).");
            Assert.That(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: 2, stepSize: -1).ToList(),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(), "Negative stepSize must throw (INV-05).");
        });
    }

    // S1 — stepSize=4 samples sequence1 starts at x in {0,4} only. With the M2 sequences and word
    // ACGT, the reachable starts are 0 and 4 anyway, so the result equals the M2 set; the point is
    // that x is always a multiple of stepSize (INV-04) and no x outside {0,4} appears.
    [Test]
    public void GenerateDotPlot_StepSizeFour_SamplesEveryFourthPosition()
    {
        var dots = Dots("ACGTACGT", "ACGTACGT", wordSize: 4, stepSize: 4);

        Assert.Multiple(() =>
        {
            Assert.That(dots, Is.EquivalentTo(new[] { (0, 0), (0, 4), (4, 0), (4, 4) }),
                "stepSize=4 keeps x in {0,4}; full overlapping match set at those starts (INV-04).");
            Assert.That(dots.Select(d => d.x), Is.All.Matches<int>(x => x % 4 == 0),
                "Every reported x must be a multiple of stepSize=4 (INV-04).");
        });
    }

    // S3 — Default wordSize locks the documented EMBOSS dottup default of 10 (dottup manpage:
    // "Default value: 10"). A length-10 self-comparison of "ACGTACGTAC" has exactly one length-10
    // word start (x=0), and that word occurs once in seq2 (x=0), so the default-parameter call must
    // yield exactly {(0,0)}. Guards against the default silently drifting from 10.
    [Test]
    public void GenerateDotPlot_DefaultWordSize_IsTenAndMatchesFullWord()
    {
        var dots = ComparativeGenomics.GenerateDotPlot("ACGTACGTAC", "ACGTACGTAC").ToHashSet();

        Assert.That(dots, Is.EquivalentTo(new[] { (0, 0) }),
            "Default wordSize=10 (EMBOSS dottup default): one length-10 word, one occurrence => exactly {(0,0)}.");
    }

    // S2 — Case-insensitivity: acgtacgt vs ACGTACGT, wordSize=4 must give the same match set as M2
    // because the implementation upper-cases both sequences before matching.
    [Test]
    public void GenerateDotPlot_CaseInsensitive_MatchesMixedCase()
    {
        var dots = Dots("acgtacgt", "ACGTACGT", wordSize: 4);

        Assert.That(dots, Is.EquivalentTo(new[] { (0, 0), (0, 4), (1, 1), (2, 2), (3, 3), (4, 0), (4, 4) }),
            "Lower-case query must match upper-case subject and yield the full M2 set (case-insensitive comparison).");
    }

    // C1 — Property (INV-04): on a constructed repeat, every reported x is a multiple of stepSize
    // and the dot count never exceeds (#sampled word-starts) * (#seq2 word-starts). Verifies the
    // O(n x m) bound and the sampling rule without hard-coding the full set.
    [Test]
    public void GenerateDotPlot_AllXCoordinatesAreMultiplesOfStep_Property()
    {
        const string seq = "ACGTACGTACGT";
        const int wordSize = 3;
        const int stepSize = 2;

        var dots = ComparativeGenomics.GenerateDotPlot(seq, seq, wordSize, stepSize).ToList();

        int sampledStarts = (seq.Length - wordSize) / stepSize + 1;
        int seq2Starts = seq.Length - wordSize + 1;

        Assert.Multiple(() =>
        {
            Assert.That(dots, Is.Not.Empty, "Self-comparison of a repeat must yield matches.");
            Assert.That(dots.Select(d => d.x), Is.All.Matches<int>(x => x % stepSize == 0),
                "Every x must be a multiple of stepSize (INV-04).");
            Assert.That(dots.Count, Is.LessThanOrEqualTo(sampledStarts * seq2Starts),
                "Dot count is bounded by sampled-starts x seq2-starts (O(n x m) bound, INV-04).");
            Assert.That(dots, Has.Member((0, 0)),
                "The main diagonal start (0,0) must be present in a self-comparison (INV-02).");
        });
    }

    #endregion
}
