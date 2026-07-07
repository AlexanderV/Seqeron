// RNA-INVERT-001 — RNA Inverted Repeats (potential stem regions)
// Evidence: docs/Evidence/RNA-INVERT-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-INVERT-001.md
// Source: Alamro et al. (2021) IUPACpal, BMC Bioinformatics 22:51, https://doi.org/10.1186/s12859-021-03983-2;
//         Ussery, Wassenaar & Borini (2008) via Wikipedia "Inverted repeat"; EMBOSS einverted (Rice et al. 2000).

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_FindInvertedRepeats_Tests
{
    private static List<(int Start1, int End1, int Start2, int End2, int Length)> Find(
        string seq, int minLength = 4, int minSpacing = 3, int maxSpacing = 100) =>
        FindInvertedRepeats(seq, minLength, minSpacing, maxSpacing).ToList();

    // Strict Watson-Crick / IUPAC complement used by the implementation (A<->U, C<->G).
    private static char Comp(char c) => Seqeron.Genomics.Core.SequenceExtensions.GetRnaComplementBase(c);

    #region FindInvertedRepeats

    // M1 — Known IR from Wikipedia/Ussery example TTACG..nnnnnn..CGTAA (RNA: UUACG..AAAAAA..CGUAA).
    // Right arm CGUAA is the reverse complement of left arm UUACG. Evidence: WGW^R (IUPACpal); INV-01,02,05.
    [Test]
    public void FindInvertedRepeats_WikipediaUssueryExample_ReturnsExactArmPositions()
    {
        var result = Find("UUACGAAAAAACGUAA");

        Assert.That(result, Is.EquivalentTo(new[] { (0, 4, 11, 15, 5) }),
            "UUACG (0-4) + loop AAAAAA (5-10) + CGUAA (11-15): right arm is the reverse complement of the left; exactly one perfect IR of arm length 5");
    }

    // M2 — Palindromic IR with minimal 3-nt loop, derived directly from WGW^R (IUPACpal).
    // GGCC + AAA + GGCC: right arm GGCC = revcomp(GGCC). Evidence: IUPACpal definition; INV-01..05.
    [Test]
    public void FindInvertedRepeats_PalindromicStemThreeNtLoop_ReturnsExactArms()
    {
        var result = Find("GGCCAAAGGCC");

        Assert.That(result, Is.EquivalentTo(new[] { (0, 3, 7, 10, 4) }),
            "GGCC (0-3) + loop AAA (4-6) + GGCC (7-10): arm length 4, right arm = reverse complement of left");
    }

    // M3 — Antiparallel REVERSE complement is required; a parallel direct repeat must be rejected.
    // AAGG + loop + AAGG: AAGG == AAGG read 5'->3' (parallel) but AAGG is NOT the reverse complement of AAGG
    // (revcomp(AAGG)=CCUU). Evidence: IUPACpal/Wikipedia "reverse complement".
    [Test]
    public void FindInvertedRepeats_ParallelDirectRepeat_IsNotReported()
    {
        var result = Find("AAGGAAAAAGG");

        Assert.That(result, Is.Empty,
            "AAGG/AAGG is a parallel direct repeat, not an inverted repeat (revcomp(AAGG)=CCUU); must not be reported");
    }

    // M4 — No inverted repeat: a homopolymer has no complementary arm (A pairs only with U). Evidence: definition.
    [Test]
    public void FindInvertedRepeats_NoComplementaryArm_ReturnsEmpty()
    {
        var result = Find("AAAAAAAAAAAA");

        Assert.That(result, Is.Empty, "poly-A has no reverse-complement arm; no WGW^R decomposition exists");
    }

    // M5 — Reverse-complement invariant (INV-01): for every reported repeat and every k,
    // complement(seq[Start2+Length-1-k]) == seq[Start1+k] (antiparallel). Evidence: IUPACpal; INV-01.
    [Test]
    public void FindInvertedRepeats_ReportedArms_SatisfyAntiparallelReverseComplement()
    {
        const string seq = "UUACGAAAAAACGUAA";
        var result = Find(seq);

        Assert.That(result, Is.Not.Empty, "the example contains a perfect IR");
        foreach (var (s1, e1, s2, e2, len) in result)
        {
            Assert.Multiple(() =>
            {
                Assert.That(e1 - s1 + 1, Is.EqualTo(len), "left arm length equals reported Length (INV-02)");
                Assert.That(e2 - s2 + 1, Is.EqualTo(len), "right arm length equals reported Length (INV-02)");
                Assert.That(s1 <= e1 && e1 < s2 && s2 <= e2, Is.True, "arms disjoint, left precedes right (INV-05)");
                for (int k = 0; k < len; k++)
                    Assert.That(Comp(seq[s2 + len - 1 - k]), Is.EqualTo(seq[s1 + k]),
                        $"position {s1 + k} must pair antiparallel with {s2 + len - 1 - k} (INV-01)");
            });
        }
    }

    // S1 — Arm shorter than minLength is not reported. GGCC arm is length 4; with minLength 5 -> empty. INV-04.
    [Test]
    public void FindInvertedRepeats_ArmShorterThanMinLength_ReturnsEmpty()
    {
        var result = Find("GGCCAAAGGCC", minLength: 5);

        Assert.That(result, Is.Empty, "the only stem has arm length 4; minLength 5 excludes it (INV-04)");
    }

    // S2 — Loop below minSpacing is excluded. GGCC+A+GGCC has loop length 1 < default minSpacing 3. INV-03.
    [Test]
    public void FindInvertedRepeats_LoopBelowMinSpacing_ReturnsEmpty()
    {
        var result = Find("GGCCAGGCC");

        Assert.That(result, Is.Empty, "loop length 1 is below minSpacing 3; the stem is not reported (INV-03)");
    }

    // S3 — Loop above maxSpacing is excluded. UUACG..AAAAAA(6)..CGUAA with maxSpacing 5 -> loop 6 > 5. INV-03.
    [Test]
    public void FindInvertedRepeats_LoopAboveMaxSpacing_ReturnsEmpty()
    {
        var result = Find("UUACGAAAAAACGUAA", maxSpacing: 5);

        Assert.That(result, Is.Empty, "loop length 6 exceeds maxSpacing 5; the stem is not reported (INV-03)");
    }

    // S4 — Null and empty inputs yield empty results without throwing (sibling-method convention).
    [Test]
    public void FindInvertedRepeats_NullOrEmpty_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Find(null!), Is.Empty, "null sequence -> empty, no throw");
            Assert.That(Find(""), Is.Empty, "empty sequence -> empty, no throw");
        });
    }

    // S5 — Sequence too short to hold two minLength arms plus the minimum loop -> empty (length guard).
    [Test]
    public void FindInvertedRepeats_TooShortSequence_ReturnsEmpty()
    {
        var result = Find("GGCC"); // needs >= 2*4 + 3 = 11

        Assert.That(result, Is.Empty, "length 4 < 2*minLength+minSpacing (11); cannot hold two arms and a loop");
    }

    // S6 — Out-of-range parameters yield an empty result with no throw, per the documented contract
    // (doc/TestSpec §3.1: minLength < 1, minSpacing < 0, maxSpacing < minSpacing -> empty). The input
    // UUACGAAAAAACGUAA contains a real IR under default params, so an empty result isolates the guard.
    [Test]
    public void FindInvertedRepeats_OutOfRangeParameters_ReturnEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Find("UUACGAAAAAACGUAA", minLength: 0), Is.Empty, "minLength < 1 -> empty (§3.1)");
            Assert.That(Find("UUACGAAAAAACGUAA", minSpacing: -1), Is.Empty, "minSpacing < 0 -> empty (§3.1)");
            Assert.That(Find("UUACGAAAAAACGUAA", minSpacing: 5, maxSpacing: 4), Is.Empty,
                "maxSpacing < minSpacing -> empty (§3.1)");
        });
    }

    // C1 — Maximal-arm extension: a perfect 5-nt arm is reported at full length, not truncated to minLength 4.
    // GGGGC + AAA + GCCCC: right arm GCCCC = revcomp(GGGGC). Evidence: einverted maximal local alignment.
    [Test]
    public void FindInvertedRepeats_PerfectArmLongerThanMinLength_ReportsFullLength()
    {
        var result = Find("GGGGCAAAGCCCC");

        Assert.That(result, Is.EquivalentTo(new[] { (0, 4, 8, 12, 5) }),
            "GGGGC (0-4) + loop AAA (5-7) + GCCCC (8-12): arm extends to length 5, not truncated to minLength 4");
    }

    // C2 — Property test (O(n^2) invariant): on several fixed sequences, every reported repeat must
    // satisfy INV-01..INV-05. Evidence: IUPACpal WGW^R; invariants in TestSpec §3.
    [Test]
    public void FindInvertedRepeats_AllReportedRepeats_SatisfyInvariants()
    {
        string[] sequences =
        {
            "UUACGAAAAAACGUAA",
            "GGCCAAAGGCC",
            "GGGGCAAAGCCCC",
            "AUGCAUGCGCAUUUGCAUGC",
            "CCCCGGGGAUACCCCGGGG",
        };

        foreach (var seq in sequences)
        {
            foreach (var (s1, e1, s2, e2, len) in Find(seq))
            {
                int loop = s2 - e1 - 1;
                Assert.Multiple(() =>
                {
                    Assert.That(len, Is.GreaterThanOrEqualTo(4), $"[{seq}] Length >= minLength (INV-04)");
                    Assert.That(e1 - s1 + 1, Is.EqualTo(len), $"[{seq}] left arm length (INV-02)");
                    Assert.That(e2 - s2 + 1, Is.EqualTo(len), $"[{seq}] right arm length (INV-02)");
                    Assert.That(loop, Is.InRange(3, 100), $"[{seq}] loop within [minSpacing,maxSpacing] (INV-03)");
                    Assert.That(s1 <= e1 && e1 < s2 && s2 <= e2, Is.True, $"[{seq}] arm ordering (INV-05)");
                    for (int k = 0; k < len; k++)
                        Assert.That(Comp(seq[s2 + len - 1 - k]), Is.EqualTo(seq[s1 + k]),
                            $"[{seq}] antiparallel reverse-complement pairing at k={k} (INV-01)");
                });
            }
        }
    }

    #endregion
}
