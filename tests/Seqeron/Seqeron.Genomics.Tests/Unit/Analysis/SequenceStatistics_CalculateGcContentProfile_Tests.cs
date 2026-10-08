// SEQ-GC-PROFILE-001 — GC Content Profile (sliding-window GC content)
// Evidence: docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-GC-PROFILE-001.md
// Source: Wikipedia, GC-content (citing primary literature),
//         https://en.wikipedia.org/wiki/GC-content (accessed 2026-06-14).
//         Biopython Bio.SeqUtils.gc_fraction, Cock P.J.A. et al. (2009)
//         Bioinformatics 25(11):1422-1423, doi:10.1093/bioinformatics/btp163.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateGcContentProfile_Tests
{
    // Expected GC% values are derived by hand from GC% = (G+C)/(A+T+G+C)×100,
    // computed independently of the implementation, per the Evidence datasets.
    private const double Tolerance = 1e-10;

    #region CalculateGcContentProfile (sliding window)

    // M1 — all-GC window: 10/10×100 = 100.0 (not the fraction 1.0).
    // Evidence: Wikipedia (G+C)/(A+T+G+C)×100.
    [Test]
    public void CalculateGcContentProfile_AllGcWindow_Returns100Percent()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGGGGGGGGG", 10).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1),
                "n=10,w=10,step=1 → exactly one window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(100.0).Within(Tolerance),
                "GGGGGGGGGG → 10/10×100 = 100.0 GC% (INV-01)");
        });
    }

    // M2 — half-GC windows: every ATGC window is 2/4×100 = 50.0.
    // Evidence: Wikipedia formula; Biopython gc_fraction("ACTG")=0.50 ×100.
    [Test]
    public void CalculateGcContentProfile_AtgcRepeats_Returns50PercentPerWindow()
    {
        // "ATGCATGCATGC", window 4, step 4 → three disjoint ATGC windows.
        var profile = SequenceStatistics.CalculateGcContentProfile("ATGCATGCATGC", 4, 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(3),
                "n=12,w=4,step=4 → offsets 0,4,8 = 3 windows (INV-03)");
            Assert.That(profile.All(p => System.Math.Abs(p - 50.0) <= Tolerance), Is.True,
                "each ATGC window → 2 GC / 4 bases × 100 = 50.0 GC%");
        });
    }

    // M3 — exact mixed profile: GGGAAATGCC, w=4, step=3 → GGGA, AAAT, TGCC.
    // Evidence: Wikipedia formula per window. Values 75/0/75 are not what a
    // fraction (0.75/0/0.75) or an off-by-one window would produce.
    [Test]
    public void CalculateGcContentProfile_MixedSequence_ReturnsExactProfile()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 3).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(3),
                "n=10,w=4,step=3 → offsets 0,3,6 = 3 windows (INV-03)");
            Assert.That(profile[0], Is.EqualTo(75.0).Within(Tolerance),
                "GGGA → 3 GC / 4 × 100 = 75.0 GC%");
            Assert.That(profile[1], Is.EqualTo(0.0).Within(Tolerance),
                "AAAT → 0 GC / 4 × 100 = 0.0 GC%");
            Assert.That(profile[2], Is.EqualTo(75.0).Within(Tolerance),
                "TGCC → 3 GC / 4 × 100 = 75.0 GC%");
        });
    }

    // M4 — N excluded from the denominator: GGAN → 2 GC / 3 standard bases × 100.
    // Evidence: Biopython gc_fraction("ACTGN")=0.50 under default "remove" (N removed
    // from length). Counting N in the denominator would give 2/4×100 = 50.0 (wrong).
    [Test]
    public void CalculateGcContentProfile_AmbiguousN_ExcludedFromDenominator()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGAN", 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1), "single window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(200.0 / 3.0).Within(Tolerance),
                "GGAN → 2 GC / 3 standard bases (N excluded) × 100 = 66.66… GC% (INV-02)");
        });
    }

    // M5 — RNA U is a non-GC base equivalent to T: GGAU → 2 GC / 4 × 100 = 50.0.
    // Evidence: Biopython gc_fraction("GGAUCUUCGGAUCU")=0.50 (U non-GC).
    [Test]
    public void CalculateGcContentProfile_RnaUracil_TreatedAsNonGc()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGAU", 4).ToArray();

        Assert.That(profile[0], Is.EqualTo(50.0).Within(Tolerance),
            "GGAU → G,G are GC, A and U are non-GC → 2/4×100 = 50.0 GC% (INV-04)");
    }

    // M6 — window count obeys INV-03 for several step sizes.
    // Evidence: sliding-window enumeration count = ⌊(n-w)/step⌋+1.
    [Test]
    public void CalculateGcContentProfile_WindowCount_MatchesInvariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 1).Count(),
                Is.EqualTo(7), "n=10,w=4,step=1 → 7 windows (INV-03)");
            Assert.That(SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 3).Count(),
                Is.EqualTo(3), "n=10,w=4,step=3 → offsets 0,3,6 = 3 windows (INV-03)");
            Assert.That(SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 2).Count(),
                Is.EqualTo(4), "n=10,w=4,step=2 → offsets 0,2,4,6 = 4 windows (INV-03)");
        });
    }

    #endregion

    #region Edge cases and invariants

    // S1 — windowSize greater than sequence length yields an empty profile.
    // Evidence: INV-03 (no full window when W > n).
    [Test]
    public void CalculateGcContentProfile_WindowLargerThanSequence_ReturnsEmpty()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("ATGC", 100).ToArray();

        Assert.That(profile, Is.Empty,
            "W=100 > n=4 → no full window → empty profile (INV-03)");
    }

    // S2 — windowSize equal to length yields exactly one whole-sequence window.
    // Evidence: INV-03; GGCC → 4/4×100 = 100.0.
    [Test]
    public void CalculateGcContentProfile_WindowEqualsLength_ReturnsSingleValue()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGCC", 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1), "W == n → one window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(100.0).Within(Tolerance),
                "GGCC → 4 GC / 4 × 100 = 100.0 GC%");
        });
    }

    // S3 — null and empty input yield empty profiles.
    // Evidence: guarded input (§3.3).
    [Test]
    public void CalculateGcContentProfile_NullOrEmpty_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateGcContentProfile(null!, 4), Is.Empty,
                "null sequence → empty profile");
            Assert.That(SequenceStatistics.CalculateGcContentProfile("", 4), Is.Empty,
                "empty sequence → empty profile");
        });
    }

    // S4 — window with no standard base (all-N) → 0 (zero-division convention).
    // Evidence: Assumption A1 (repository convention, matches SEQ-GC-ANALYSIS-001).
    [Test]
    public void CalculateGcContentProfile_AllAmbiguousWindow_ReturnsZero()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("NNNN", 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1), "single window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(0.0).Within(Tolerance),
                "no standard base → 0 (zero-division convention, INV-05)");
        });
    }

    // C1 — case-insensitivity: lowercase input produces the same profile as uppercase.
    // Evidence: implementation case-folds before counting (§3.3).
    [Test]
    public void CalculateGcContentProfile_LowercaseInput_MatchesUppercase()
    {
        var upper = SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 1).ToArray();
        var lower = SequenceStatistics.CalculateGcContentProfile("gggaaatgcc", 4, 1).ToArray();

        Assert.That(lower, Is.EqualTo(upper).Within(Tolerance),
            "case-folded counting → lowercase profile equals uppercase profile");
    }

    // C2 — INV-01: every window value is bounded in [0, 100].
    // Evidence: numerator ≤ denominator, then ×100 (Wikipedia formula).
    [Test]
    public void CalculateGcContentProfile_AllWindows_BoundedZeroToHundred()
    {
        var profile = SequenceStatistics
            .CalculateGcContentProfile("ATGCGGGGAAAACCCCATGCATGC", 6, 1).ToArray();

        Assert.That(profile, Is.Not.Empty, "profile must contain windows");
        Assert.That(profile.All(p => p >= -Tolerance && p <= 100.0 + Tolerance), Is.True,
            "every GC% value in [0,100] (INV-01)");
    }

    #endregion

    #region Reference cross-checks (B03 review 2026-09, F17/D6)

    // EMBOSS 6.6.0 isochore -window 4 (executed; fraction 0.750 0.500 0.250 0.000 0.250 0.500 0.750)
    // and Biopython 1.88 gc_fraction(window, "remove")·100 on every window: identical.
    [Test]
    public void CalculateGcContentProfile_MatchesEmbossIsochoreAndBiopythonGcFraction()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 1).ToArray();
        Assert.That(profile, Is.EqualTo(new[] { 75.0, 50.0, 25.0, 0.0, 25.0, 50.0, 75.0 }).Within(Tolerance));
    }

    // Biopython 1.88 [gc_fraction(s[i:i+4], "remove")*100 for i in range(0, len(s)-3, 2)]
    // on lowercase RNA with N runs; all-N window → 0 (Biopython returns 0 for a zero denominator).
    [Test]
    public void CalculateGcContentProfile_RnaWithN_Step2_MatchesBiopythonGcFractionRemove()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("ggnnAAUGCCnnnn", 4, 2).ToArray();
        Assert.That(profile, Is.EqualTo(new[] { 100.0, 0.0, 25.0, 75.0, 100.0, 0.0 }).Within(Tolerance));
    }

    // Biopython 1.88 gc_fraction "remove" on a polluted 60-mer (mixed case, N, gaps; random.seed(7)),
    // window 10, step 7 (partial trailing window not reported), percent and fraction forms.
    [Test]
    public void CalculateGcContentProfile_PollutedSequence_MatchesBiopythonGcFractionRemove()
    {
        const string seq = "cGgACNCc-ANTACggCTCNgA-CT-A--gATANGagGNC-aNGC--TcCNC-A-TtNgc";
        double[] expected = { 75.0, 62.5, 75.0, 33.33333333333333, 50.0, 71.42857142857143, 83.33333333333334, 50.0 };

        var percent = SequenceStatistics.CalculateGcContentProfile(seq, 10, 7).ToArray();
        var fraction = SequenceStatistics.CalculateGcContentProfile(seq, 10, 7, fraction: true).ToArray();

        Assert.That(percent, Is.EqualTo(expected).Within(1e-12));
        Assert.That(fraction, Is.EqualTo(expected.Select(v => v / 100.0).ToArray()).Within(1e-14));
    }

    // D6 — the profile delegates to the canonical SequenceExtensions.CalculateGcFraction; locks
    // behaviour preservation against the former inline kernel (ToUpperInvariant, G/C over A/C/G/T/U)
    // over every UTF-16 code unit and 2000 random polluted inputs.
    [Test]
    public void CalculateGcContentProfile_EqualsFormerInlineKernelAndCanonicalGcFraction()
    {
        static double Former(string w)
        {
            int gc = 0, total = 0;
            foreach (char ch in w.ToUpperInvariant())
            {
                if (ch == 'G' || ch == 'C') { gc++; total++; }
                else if (ch == 'A' || ch == 'T' || ch == 'U') total++;
            }
            return total > 0 ? (double)gc / total * 100.0 : 0;
        }

        var all = new string(Enumerable.Range(0, 0x10000).Select(c => (char)c).ToArray());
        var single = SequenceStatistics.CalculateGcContentProfile(all, 1, 1).ToArray();
        for (int c = 0; c < 0x10000; c++)
            Assert.That(single[c], Is.EqualTo(Former(((char)c).ToString())), $"U+{c:X4}");

        var rng = new Random(20260928);
        const string alphabet = "ACGTUacgtuNnRYSW-* .";
        for (int iter = 0; iter < 2000; iter++)
        {
            var s = new string(Enumerable.Range(0, rng.Next(1, 80)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
            int w = rng.Next(1, s.Length + 1), step = rng.Next(1, 6);
            var got = SequenceStatistics.CalculateGcContentProfile(s, w, step).ToArray();
            var exp = new List<double>();
            for (int i = 0; i + w <= s.Length; i += step)
            {
                exp.Add(Former(s.Substring(i, w)));
                Assert.That(got[exp.Count - 1], Is.EqualTo(s.Substring(i, w).AsSpan().CalculateGcContent()));
            }
            Assert.That(got, Is.EqualTo(exp), s);
        }
    }

    // F19 — Biopython 1.88 [gc_fraction(s[i:i+6], mode)*100 for i in range(0, len(s)-5, 3)] with
    // s = "ACGTSSWWNNRYacgusw" (executed) for ambiguous = remove / ignore / weighted.
    [TestCase(SequenceExtensions.GcAmbiguityMode.Remove, new[] { 66.66666666666666, 40.0, 0.0, 66.66666666666666, 50.0 })]
    [TestCase(SequenceExtensions.GcAmbiguityMode.Ignore, new[] { 66.66666666666666, 33.33333333333333, 0.0, 33.33333333333333, 50.0 })]
    [TestCase(SequenceExtensions.GcAmbiguityMode.Weighted, new[] { 66.66666666666666, 41.66666666666667, 33.33333333333333, 58.333333333333336, 50.0 })]
    public void CalculateGcContentProfile_AmbiguityMode_MatchesBiopythonGcFraction(
        SequenceExtensions.GcAmbiguityMode mode, double[] expected)
    {
        var percent = SequenceStatistics.CalculateGcContentProfile("ACGTSSWWNNRYacgusw", 6, 3, false, mode).ToArray();
        var fraction = SequenceStatistics.CalculateGcContentProfile("ACGTSSWWNNRYacgusw", 6, 3, true, mode).ToArray();

        Assert.That(percent, Is.EqualTo(expected).Within(1e-12));
        Assert.That(fraction, Is.EqualTo(expected.Select(v => v / 100.0).ToArray()).Within(1e-14));
    }

    [Test]
    public void CalculateGcContentProfile_AmbiguityMode_InvalidWindowThrows_EmptyAndShortAreEmpty()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.CalculateGcContentProfile("ACGT", 0, 1, false, SequenceExtensions.GcAmbiguityMode.Remove));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.CalculateGcContentProfile("ACGT", 2, 0, false, SequenceExtensions.GcAmbiguityMode.Weighted));
        Assert.That(SequenceStatistics.CalculateGcContentProfile("", 2, 1, false, SequenceExtensions.GcAmbiguityMode.Remove), Is.Empty);
        Assert.That(SequenceStatistics.CalculateGcContentProfile("ACG", 4, 1, false, SequenceExtensions.GcAmbiguityMode.Remove), Is.Empty);
    }

    // F17 — window/step below 1 rejected eagerly (step 0 previously never terminated; W 0 returned n+1 zeros).
    [TestCase(0, 1)]
    [TestCase(-5, 1)]
    [TestCase(4, 0)]
    [TestCase(4, -1)]
    public void CalculateGcContentProfile_WindowOrStepBelowOne_Throws(int windowSize, int stepSize)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.CalculateGcContentProfile("ACGTACGT", windowSize, stepSize));
    }

    #endregion
}
