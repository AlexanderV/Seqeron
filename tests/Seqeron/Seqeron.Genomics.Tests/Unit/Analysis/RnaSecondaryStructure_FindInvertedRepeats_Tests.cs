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
        string seq, int minLength = 4, int minSpacing = 3, int maxSpacing = 100, int maxMismatches = 0,
        bool allowWobble = true) =>
        FindInvertedRepeats(seq, minLength, minSpacing, maxSpacing, maxMismatches, allowWobble).ToList();

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
                    Assert.That(CanPair(seq[s1 + k], seq[s2 + len - 1 - k]), Is.True,
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
                        Assert.That(CanPair(seq[s1 + k], seq[s2 + len - 1 - k]), Is.True,
                            $"[{seq}] antiparallel pairing at k={k} (INV-01)");
                });
            }
        }
    }

    #endregion

    #region Review 2026-09 (B12 F22): delegation to the canonical RepeatFinder engine, RNA pairing rule, oracle

    // Independent brute force of the stem model: S = every exact stem (L, R, len) whose pairs all pair, with
    // len >= minArm and loop in [minLoop, maxLoop]; a stem is reported iff no other member of S contains it in
    // BOTH arms (EMBOSS palindrome nesting filter). Written without reference to the implementation.
    private static List<(int, int, int, int, int)> BruteForce(
        string seq, int minArm, int minLoop, int maxLoop, bool wobble)
    {
        string s = seq.ToUpperInvariant().Replace('T', 'U');
        int n = s.Length;
        bool Pair(char a, char b) => (a, b) is ('A', 'U') or ('U', 'A') or ('G', 'C') or ('C', 'G')
            || (wobble && (a, b) is ('G', 'U') or ('U', 'G'));
        var all = new HashSet<(int L, int R, int Len)>();
        for (int l = 0; l < n; l++)
            for (int e = l + 1; e < n; e++)
                for (int len = 1; e - len + 1 - (l + len) >= 0; len++)
                {
                    int r = e - len + 1;
                    int loop = r - (l + len);
                    if (!Pair(s[l + len - 1], s[r])) break; // stems anchored at outer pair (l, e) grow inward
                    if (len >= minArm && loop >= minLoop && loop <= maxLoop) all.Add((l, r, len));
                }
        return all
            .Where(a => !all.Any(b => b != a && b.L <= a.L && a.L + a.Len <= b.L + b.Len
                                      && b.R <= a.R && a.R + a.Len <= b.R + b.Len))
            .OrderBy(a => a.L).ThenBy(a => a.R)
            .Select(a => (a.L, a.L + a.Len - 1, a.R, a.R + a.Len - 1, a.Len))
            .ToList();
    }

    // Completeness / no false positives / no duplicates / ordering vs the brute force, many seeds, all alphabets
    // (U/T/lower case/N/IUPAC), both pairing rules, minLength 1..4. Python oracle cross-check of the same 12 000
    // cases: 0 mismatches (Evidence §RNA-INVERT-001 F22).
    [Test]
    public void FindInvertedRepeats_RandomStrings_EqualBruteForceOracle(
        [Values(0, 1)] int wobbleFlag, [Values(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)] int seed)
    {
        string[] alphabets = { "ACGU", "ACGT", "ACGUN", "acguNR", "ACG", "GU", "ACGTUacgtu" };
        var rng = new Random(seed * 7919 + wobbleFlag);
        for (int iter = 0; iter < 300; iter++)
        {
            string a = alphabets[rng.Next(alphabets.Length)];
            string seq = new string(Enumerable.Range(0, rng.Next(0, 36)).Select(_ => a[rng.Next(a.Length)]).ToArray());
            int minArm = rng.Next(1, 5), minLoop = new[] { 0, 1, 3 }[rng.Next(3)], maxLoop = new[] { 3, 5, 100 }[rng.Next(3)];
            if (maxLoop < minLoop) continue;
            var got = Find(seq, minArm, minLoop, maxLoop, 0, wobbleFlag == 1);
            var expected = BruteForce(seq, minArm, minLoop, maxLoop, wobbleFlag == 1);
            Assert.That(got, Is.EqualTo(expected),
                $"seq={seq} minArm={minArm} loop=[{minLoop},{maxLoop}] wobble={wobbleFlag}");
        }
    }

    // G-U wobble (ViennaRNA CanPair) is a stem pair by default; allowWobble:false restricts to A-U / G-C.
    // GGGGAAAAUUUU: G-U x4 across a 4-nt loop (and the one-shifted diagonal (1..4)/(8..11)).
    [Test]
    public void FindInvertedRepeats_WobbleStem_ReportedOnlyWithWobble()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Find("GGGGAAAAUUUU"), Is.EqualTo(new[] { (0, 3, 8, 11, 4), (1, 4, 8, 11, 4) }));
            Assert.That(Find("GGGGAAAATTTT"), Is.EqualTo(new[] { (0, 3, 8, 11, 4), (1, 4, 8, 11, 4) }), "T is read as U");
            Assert.That(Find("ggggaaaauuuu"), Is.EqualTo(new[] { (0, 3, 8, 11, 4), (1, 4, 8, 11, 4) }), "case-insensitive");
            Assert.That(Find("GGGGAAAAUUUU", allowWobble: false), Is.Empty, "Watson-Crick only: G-U does not pair");
        });
    }

    // Alternative pairings are not suppressed (no unsourced greedy "non-overlapping" selection): the two stems
    // GGGG/CCCC (0-3,7-10) and CCCC/GGGG (7-10,14-17) share the CCCC block; the former greedy rule kept only one.
    [Test]
    public void FindInvertedRepeats_StemsSharingBases_AreBothReported()
    {
        Assert.That(Find("GGGGAAACCCCAAAGGGG"), Is.EqualTo(new[] { (0, 3, 7, 10, 4), (7, 10, 14, 17, 4) }));
    }

    // Ambiguity codes never pair (ViennaRNA CanPair); the former strict-IUPAC complement paired N-N, S-S, W-W.
    [Test]
    public void FindInvertedRepeats_AmbiguityCodes_NeverPair()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Find("NNNNAAANNNN"), Is.Empty);
            Assert.That(Find("SSSSAAASSSS"), Is.Empty);
            Assert.That(Find("WWWWGGGWWWW"), Is.Empty);
            Assert.That(Find("GGNNAAAAAANNCC", minLength: 2), Is.EqualTo(new[] { (0, 1, 12, 13, 2) }),
                "N·N never extends the GG/CC stem inward");
        });
    }

    // maxMismatches = k: IUPACpal k-mismatch / EMBOSS palindrome -nummismatches (interior mismatches only).
    // GGAGG AAA CCUCC with one interior mismatch: arms GGAGG vs revcomp(CCACC) -> A-A mismatch in the middle.
    [Test]
    public void FindInvertedRepeats_MaxMismatches_AllowsInteriorMismatchOnly()
    {
        // arms GGAGG / CCACC: pair (A,A) at the stem centre mismatches
        Assert.Multiple(() =>
        {
            Assert.That(Find("GGAGGAAACCACC", minLength: 5), Is.Empty, "k = 0: the central A-A mismatch breaks the stem");
            Assert.That(Find("GGAGGAAACCACC", minLength: 5, maxMismatches: 1), Is.EqualTo(new[] { (0, 4, 8, 12, 5) }));
            Assert.Throws<ArgumentOutOfRangeException>(() => FindInvertedRepeats("GGCC", maxMismatches: -1));
        });
    }

    // With k > 0 every reported stem still has pairing outer and inner pairs, <= k interior mismatches, equal arms.
    [Test]
    public void FindInvertedRepeats_MismatchResults_SatisfyStemInvariants()
    {
        var rng = new Random(20260925);
        for (int iter = 0; iter < 400; iter++)
        {
            string seq = new string(Enumerable.Range(0, rng.Next(10, 40)).Select(_ => "ACGU"[rng.Next(4)]).ToArray());
            int k = rng.Next(1, 3);
            var res = Find(seq, 4, 3, 30, k);
            Assert.That(res, Is.Ordered.Ascending.Using<(int, int, int, int, int)>((x, y) =>
                x.Item1 != y.Item1 ? x.Item1.CompareTo(y.Item1) : x.Item3.CompareTo(y.Item3)));
            Assert.That(res.Distinct().Count(), Is.EqualTo(res.Count), "no duplicates");
            foreach (var (s1, e1, s2, e2, len) in res)
            {
                int mism = 0;
                for (int j = 0; j < len; j++)
                    if (!CanPair(seq[s1 + j], seq[e2 - j])) mism++;
                Assert.That(mism, Is.LessThanOrEqualTo(k), seq);
                Assert.That(CanPair(seq[s1], seq[e2]) && CanPair(seq[e1], seq[s2]), Is.True, "mismatches only interior");
                Assert.That(len, Is.GreaterThanOrEqualTo(4));
                Assert.That(s2 - e1 - 1, Is.InRange(3, 30));
            }
        }
    }

    // Canonical-engine equivalence: for ACGT input the RNA wrapper reports exactly RepeatFinder's stems (WC mode),
    // so there is one implementation of the stem model (DUP_MAP RNA entry).
    [Test]
    public void FindInvertedRepeats_EqualsRepeatFinderCanonicalEngine()
    {
        var rng = new Random(77);
        for (int iter = 0; iter < 300; iter++)
        {
            string seq = new string(Enumerable.Range(0, rng.Next(8, 60)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
            int k = rng.Next(0, 2);
            foreach (bool wobble in new[] { false, true })
            {
                var rna = Find(seq, 4, 3, 50, k, wobble);
                var dna = Seqeron.Genomics.Analysis.RepeatFinder.FindInvertedRepeats(seq, 4, 50, 3, k, int.MaxValue, wobble)
                    .Select(r => (r.LeftArmStart, r.LeftArmStart + r.ArmLength - 1, r.RightArmStart,
                                  r.RightArmStart + r.ArmLength - 1, r.ArmLength)).ToList();
                Assert.That(rna, Is.EqualTo(dna), seq);
            }
        }
    }

    #endregion
}
